#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This defines the methods that work with DSSP files.

Fábio Madeira, 2017+

"""


import os
import json
import pandas as pd

from string import digits
from string import ascii_uppercase

from prointvar.mmcif import MMCIFreader
from prointvar.mmcif import MMCIFwriter

from prointvar.utils import compute_rsa
from prointvar.utils import flash
from prointvar.utils import get_rsa_class
from prointvar.utils import row_selector
# TODO FIX
from prointvar.utils import logging_out

from prointvar.library import dssp_types

from prointvar.config import config

"""
Useful info from http://web.expasy.org/docs/userman.html#FT_line

Secondary structure (HELIX, STRAND, TURN) - The feature table of sequence entries of proteins whose tertiary
structure is known experimentally contains the secondary structure information corresponding to that protein.
The secondary structure assignment is made according to DSSP (see Kabsch W., Sander C.; Biopolymers,
22:2577-2637(1983)) and the information is extracted from the coordinate data sets of the Protein Data Bank (PDB).

In the feature table only three types of secondary structure are specified: helices (key HELIX), beta-strands
(key STRAND) and turns (key TURN). Residues not specified in one of these classes are in a 'loop' or
'random-coil' structure. Because the DSSP assignment has more than the three common secondary structure classes,
we have converted the different DSSP assignments to HELIX, STRAND and TURN as shown in the table below.

DSSP    code	                                            DSSP definition	Swiss-Prot assignment
H	    Alpha-helix	                                        HELIX
G	    3(10) helix	                                        HELIX
I	    Pi-helix	                                        HELIX
E	    Hydrogen-bonded beta-strand (extended strand)	    STRAND
B	    Residue in an isolated beta-bridge	                STRAND
T	    H-bonded turn (3-turn, 4-turn or 5-turn)	        TURN
S	    Bend (five-residue bend centered at residue i)	    Not specified

One should be aware of the following facts:
Segment length. For helices (alpha and 3-10), the residue just before and just after the helix as given by
DSSP participates in the helical hydrogen-bonding pattern with a single H-bond. For practical purposes,
one can extend the HELIX range by one residue on each side, e.g. HELIX 25-35 instead of HELIX 26-34.
Also, the ends of secondary structure segments are less well defined for lower-resolution structures.
A fluctuation of one residue is common.
Missing segments. In low-resolution structures, badly formed helices or strands may be omitted in the
DSSP definition.
Special helices and strands. Helices of length three are 3-10 helices, those of length four and longer
are either alpha-helices or 3-10 helices (pi helices are extremely rare). A strand of one residue
corresponds to a residue in an isolated beta-bridge. Such bridges can be structurally important.
Missing secondary structure. No secondary structure is currently given in the feature table in the following cases:
No sequence data in the PDB entry;
Structure for which only C-alpha coordinates are in PDB;
NMR structure with more than one coordinate data set;
Model (i.e. theoretical) structure.
"""


def parse_dssp_from_file(inputfile, excluded=(), add_full_chain=True, add_ss_reduced=False,
                         add_rsa=True, method="Sander", add_rsa_class=False, verbose=False):
    """
    Parse lines of the DSSP file to get entries for every Residue
    in each CHAIN. The hierachy is maintained. CHAIN->RESIDUE->[...].

    :param inputfile: path to the DSSP file
    :param excluded: option to exclude DSSP columns
    :param add_full_chain: boolean
    :param add_ss_reduced: boolean
    :param add_rsa: boolean
    :param add_rsa_class: boolean
    :param method: name of the method
    :param verbose: boolean
    :return: returns a pandas DataFrame
    """

    if verbose:
        flash("Parsing DSSP from lines...")

    # example lines with some problems
    """
      #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
        1    1 A M              0   0  127      0, 0.0   345,-0.1     0, 0.0     3,-0.1   0.000 360.0 360.0 360.0 162.0  -18.7   21.6  -55.4
        2    2 A R        +     0   0  117      1,-0.1    28,-0.4   343,-0.1     2,-0.3   0.455 360.0  81.5-136.8 -28.7  -17.0   22.3  -52.1

      381  394 A K              0   0  125     -2,-0.4   -21,-0.1   -21,-0.2    -2,-0.0  -0.421 360.0 360.0  64.1 360.0  -22.5   44.2  -25.4
      382        !*             0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0
      383    1 A M              0   0  127      0, 0.0   345,-0.1     0, 0.0     3,-0.1   0.000 360.0 360.0 360.0 162.0  -10.0   71.4  -55.4

    10278  103 H H  E     -XZ1023010269W  69     -9,-2.3    -9,-2.2    -2,-0.3     2,-1.0  -0.884  22.6-128.4-108.1 141.6  -97.0   28.7  112.2
    10279  104 H I  E     +XZ1022910268W   0    -50,-2.2   -50,-0.6    -2,-0.4   -11,-0.3  -0.801  30.6 175.4 -90.4  95.6  -98.5   32.0  111.3
    10280  105 H L  E     +     0   0   21    -13,-1.7   -55,-2.5    -2,-1.0     2,-0.3   0.812  62.6   4.9 -70.5 -35.5  -96.3   34.5  113.1

    # missing segment break
      145        !              0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0

    # chain break
      382        !*             0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0
    """

    if not os.path.isfile(inputfile):
        raise IOError("{} not available or could not be read...".format(inputfile))

    # column width descriptors
    header = ("LINE", "RES", "CHAIN", "AA", "SS", "STRUCTURE",
              "BP1", "BP2", "BP2_CHAIN", "ACC",
              "NH_O_1", "NH_O_1_nrg", "O_HN_1", "O_HN_1_nrg",
              "NH_O_2", "NH_O_2_nrg", "O_HN_2", "O_HN_2_nrg",
              "TCO", "KAPPA", "ALPHA", "PHI", "PSI",
              "X-CA", "Y-CA", "Z-CA")

    widths = ((0, 5), (5, 11), (11, 12), (12, 15), (16, 17), (17, 25),
              (25, 29), (29, 33), (33, 34), (34, 38),
              (38, 45), (46, 50), (50, 56), (57, 61),
              (61, 67), (68, 72), (72, 78), (79, 84),
              (85, 91), (91, 97), (97, 103), (103, 109), (109, 115),
              (115, 123), (123, 130), (130, 137))

    all_str = {key: str for key in header}
    table = pd.read_fwf(inputfile, skiprows=28, names=header, colspecs=widths,
                        compression=None, converters=all_str, keep_default_na=False)

    # table modular extensions
    if add_full_chain:
        table = add_dssp_full_chain(table)

    table['SS'] = table.SS.fillna('-')
    if add_ss_reduced:
        table = add_dssp_ss_reduced(table)

    if add_rsa:
        table = add_dssp_rsa(table, method=method)

    if add_rsa_class:
        table = add_dssp_rsa_class(table)

    # drop missing residues ("!")  and chain breaks ("!*")
    table = table[table['AA'] != '!']
    table = table[table['AA'] != '!*']

    if excluded is not None:
        assert type(excluded) is tuple
        try:
            table = table.drop(list(excluded), axis=1)
        except ValueError:
            # most likely theses are not in there
            pass

    # enforce some specific column types
    for col in table:
        if col in dssp_types:
            try:
                table[col] = table[col].astype(dssp_types[col])
            except ValueError:
                # there are some NaNs in there
                pass

    if table.empty:
        raise ValueError('{} resulted in an empty DataFrame...'.format(inputfile))

    return table


def get_dssp_selected_from_table(data, chain=None, chain_full=None, res=None):
    """
    Utility that filters a pandas DataFrame by the input tuples.

    :param data: pandas DataFrame object
    :param chain: (tuple) chain IDs or None
    :param chain_full: (tuple) alternative chain IDs or None
    :param res: (tuple) res IDs or None
    :return: returns a modified pandas DataFrame
    """

    # excluding rows
    table = data
    if chain is not None:
        table = row_selector(table, 'CHAIN', chain, method="isin")

    if chain_full is not None:
        table = row_selector(table, 'CHAIN_FULL', chain_full, method="isin")

    if res is not None:
        table = row_selector(table, 'RES', res, method="isin")

    return table


def add_dssp_full_chain(data):
    """
    Utility that adds a new column to the table.
    Specific to DSSP outputs that are generated from mmCIF files containing
    multiple char chain IDs (e.g. 'AA' and 'BA'). These are found in the
    Biological Unit mmCIF structures.

    :param data: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    # BioUnits chain naming seems to follow the pattern:
    # chain A becomes AA then
    # AA->AZ then A0->A9 [A-Z then 0-9] and then AAA->AAZ and AA0->AA9
    # then ABA->ABZ and AB0->AB9
    alpha1 = [k for k in ascii_uppercase + digits]
    alpha2 = ['A' + k for k in alpha1]
    alpha3 = ['B' + k for k in alpha1]
    new_alphabet = alpha1 + alpha2 + alpha3

    table = data
    chains_full = []
    c = -1
    for ix in table.index:
        chain_id = table.loc[ix, "CHAIN"]
        aa_id = table.loc[ix, "AA"]
        if aa_id == "!*":
            if table.loc[ix - 1, "CHAIN"] == table.loc[ix + 1, "CHAIN"]:
                c += 1
            else:
                c = -1
        if c != -1 and aa_id != "!*" and aa_id != "!":
            if c >= len(new_alphabet):
                raise IndexError('Alphabet needs update to accommodate '
                                 'such high number of chains...')
            chain_id += new_alphabet[c]
        chains_full.append(chain_id)
    if not chains_full:
        table["CHAIN_FULL"] = table["CHAIN"]
    else:
        table["CHAIN_FULL"] = chains_full
    return table


def add_dssp_rsa(data, method="Sander"):
    """
    Utility that adds a new column to the table.
    Adds a new column with Relative Solvent Accessibility (RSA).

    :param data: pandas DataFrame object
    :param method: name of the method
    :return: returns a modified pandas DataFrame
    """

    table = data
    rsas = []
    for i in table.index:
        rsas.append(compute_rsa(table.loc[i, "ACC"], table.loc[i, "AA"],
                                method=method))
    table["RSA"] = rsas
    return table


def add_dssp_rsa_class(data, rsa_col='RSA'):
    """
    Utility that adds a new column to the table.
    Adds a new column with Relative Solvent Accessibility (RSA) classes.

    :param data: pandas DataFrame object
    :param rsa_col: column name
    :return: returns a modified pandas DataFrame
    """

    table = data
    rsas_class = []
    for i in table.index:
        rsas_class.append(get_rsa_class(table.loc[i, "{}".format(rsa_col)]))
    table["{}_CLASS".format(rsa_col)] = rsas_class
    return table


def add_dssp_ss_reduced(data):
    """
    Utility that adds a new column to the table.
    Adds a reduced-stated Secondary Structure (SS).

    :param data: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    table = data
    alphas = ['H']
    betas = ['E']
    coils = ['G', 'I', 'B', 'C', 'T', 'S', '', ' ']

    # replace some NaN with custom strings
    # table['SS'] = table.SS.fillnan('-')
    sss = []
    for ix in table.index:
        ss = table.loc[ix, "SS"]
        if ss in alphas:
            ss = 'H'
        elif ss in betas:
            ss = 'E'
        elif ss in coils:
            ss = 'C'
        else:
            ss = '-'
        sss.append(ss)

    table["SS_CLASS"] = sss

    return table


class DSSPreader(object):
    def __init__(self, inputfile, verbose=False):
        """
        :param inputfile: Needs to point to a valid DSSP file.
        :param verbose: boolean
        """
        self.inputfile = inputfile
        self.verbose = verbose
        self.data = None
        self.excluded = ("LINE", "STRUCTURE", "BP1", "BP2", "BP2_CHAIN",
                         "NH_O_1", "NH_O_1_nrg", "O_HN_1", "O_HN_1_nrg",
                         "NH_O_2", "NH_O_2_nrg", "O_HN_2", "O_HN_2_nrg",
                         "X-CA", "Y-CA", "Z-CA")

        if not os.path.isfile(inputfile):
            raise IOError("{} not available or could not be read...".format(inputfile))

    def read(self, **kwargs):
        return self.residues(**kwargs)

    def residues(self, excluded=None, add_full_chain=True, add_ss_reduced=False,
                 add_rsa=True, method="Sander", add_rsa_class=False):
        if excluded is None:
            excluded = self.excluded
        self.data = parse_dssp_from_file(self.inputfile, excluded=excluded,
                                         add_full_chain=add_full_chain,
                                         add_ss_reduced=add_ss_reduced,
                                         add_rsa=add_rsa, method=method,
                                         add_rsa_class=add_rsa_class,
                                         verbose=self.verbose)
        return self.data

    def to_json(self, pretty=True):
        if self.data is not None:
            if type(self.data) is pd.core.frame.DataFrame:
                data = self.data.to_dict(orient='records')
            else:
                data = self.data
            if pretty:
                return json.dumps(data, sort_keys=False, indent=4)
            else:
                return json.dumps(data)
        else:
            flash('No DSSP data parsed...')


class DSSPgenerator(object):
    def __init__(self, inputfile, outputfile=None, verbose=False):
        """
        :param inputfile: Needs to point to a valid PDB or mmCIF file.
        :param outputfile: if not provided will use the same file name and <.dssp> extension
        :param verbose: boolean
        """
        self.inputfile = inputfile
        self.outputfile = outputfile
        self.verbose = verbose
        self.data = None

        # generate outputfile if missing
        self._generate_output()

        if not os.path.isfile(inputfile):
            raise IOError("{} not available or could not be read...".format(inputfile))

        # inputfile needs to be in PDB or mmCIF format
        filename, extension = os.path.splitext(inputfile)
        if extension not in ['.pdb', '.ent', '.cif']:
            raise ValueError("{} is expected to be in mmCIF or PDB format..."
                             "".format(inputfile))

    def _generate_output(self):
        if not self.outputfile:
            filename, extension = os.path.splitext(self.inputfile)
            self.outputfile = filename + ".dssp"

    def _run(self, dssp_bin):
        cmd = "{} -i {} -o {}".format(dssp_bin, self.inputfile, self.outputfile)
        os.system(cmd)
        if not os.path.isfile(self.outputfile):
            raise IOError("DSSP output not generated for {}".format(self.outputfile))

    def run(self, override=False):
        if not os.path.exists(self.outputfile) or override:
            if os.path.isfile(config.dssp_bin):
                dssp_bin = config.dssp_bin
            elif os.path.isfile(config.dssp_bin_local):
                dssp_bin = config.dssp_bin_local
            else:
                raise IOError('DSSP executable is not available...')

            # run dssp and generate output
            self._run(dssp_bin)

        else:
            flash('DSSP for {} already available...'.format(self.outputfile))
        return


def dssp_runner_split_chains(values, bio=False, override=False, logger=None, verbose=False):
    """
    Runs DSSP for mmCIF files (Asymmetric Units or BioUnits). The difference here
    is that this splits the mmCIF in its chains and runs dssp for each chain independently.

    :param values: List of PDB IDs
    :param bio: (boolean) use standard mmCIF files or the BioUnits
    :param override: boolean
    :param logger: standard logger or None
    :param verbose: boolean
    :return: (side effects) computes DSSP and generates new files <pdbid>_<chainid>.dssp
    """

    for pdbid in values:
        if bio:
            inputcif = "{}{}{}.cif".format(config.db_root, config.db_cif_biounit, pdbid)
        else:
            inputcif = "{}{}{}.cif".format(config.db_root, config.db_cif, pdbid)

        # read cif file and get available chains
        r = MMCIFreader(inputcif)
        data = r.read(add_res_full=False, add_contacts=False)
        chains = [k for k in data.loc[:, 'label_asym_id'].unique()]
        for chain in chains:
            # since len(chain) > 1 (e.g. 'AA' or 'BA') are repetitions these are skipped here
            if len(chain) == 1:
                # write out the new cif file with the current chain
                filename, extension = os.path.splitext(inputcif)
                outputcif = filename + '_{}.cif'.format(chain)
                if not os.path.isfile(outputcif) or override:
                    message = "Generating mmCIF for {}_{}...".format(pdbid, chain)
                    logging_out(message, 'info', logger, verbose)
                    w = MMCIFwriter(inputcif, outputcif)
                    w.run(data=data, chain=(chain,), res=None, atom=None, lines=None,
                          override=override)
                    message = "Generated mmCIF for {}_{}...".format(pdbid, chain)
                    logging_out(message, 'info', logger, verbose)
                else:
                    message = "mmCIF for {}_{} already available...".format(pdbid, chain)
                    logging_out(message, 'info', logger, verbose)
                # generating the dssp output for the current chain
                outputdssp = "{}{}{}_{}.dssp".format(config.db_root, config.db_dssp_generated, pdbid, chain)
                if not os.path.isfile(outputdssp) or override:
                    message = "Generating DSSP for {}_{}...".format(pdbid, chain)
                    logging_out(message, 'info', logger, verbose)
                    DSSPgenerator(outputcif, outputdssp, verbose=verbose).run(override=override)
                    message = "Generated DSSP for {}_{}...".format(pdbid, chain)
                    logging_out(message, 'info', logger, verbose)
                else:
                    message = "DSSP for {}_{} already available...".format(pdbid, chain)
                    logging_out(message, 'info', logger, verbose)
    return

if __name__ == '__main__':
    pdbid = "2pah"
    # pdbid = "1cg2"
    # pdbid = "3kic"

    # new error cases
    # pdbid = "2rie"
    # pdbid = "3f1i"
    pdbid = "3j6l"

    from prointvar.config import config as c

    # inputcif = "{}{}{}.cif".format(c.db_root, c.db_cif, pdbid)
    inputdssp = "{}{}{}_bio.dssp".format(c.db_root, c.db_dssp_generated, pdbid)
    # DSSPgenerator(inputcif, inputdssp).run()

    d = DSSPreader(inputdssp)
    d.read()
    nd = d.data
    # nd = d.to_json(pretty=True)
    # print(nd)
    # nd = json.loads(nd)
    # print(nd[0])
    # print([k for k in nd[0]])
    # print(nd.loc[0, "CHAIN"])
    # print(nd.loc[:, "CHAIN"].unique())
    # print(nd.loc[:, "RES"].unique())
    # print(nd.loc[nd['AA'] == '!*'])
    print(nd.tail())
    # print(get_dssp_selected_from_table(nd, chain_full=('BA',)).tail())
    pass
