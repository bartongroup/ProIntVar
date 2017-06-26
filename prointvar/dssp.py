#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This defines the methods that work with DSSP files.

Fábio Madeira, 2017+

"""


import os
import json
import logging
import pandas as pd
from io import StringIO
from string import digits
from string import ascii_uppercase

from prointvar.pdbx import PDBXreader
from prointvar.pdbx import PDBXwriter

from prointvar.utils import compute_rsa
from prointvar.utils import get_rsa_class
from prointvar.utils import row_selector
from prointvar.utils import lazy_file_remover
from prointvar.library import dssp_types

from prointvar.config import config

logger = logging.getLogger("prointvar")


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
                         add_rsa=True, method="Sander", add_rsa_class=False,
                         reset_res_id=False):
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
    :param reset_res_id: boolean
    :return: returns a pandas DataFrame
    """

    logger.info("Parsing DSSP from lines...")

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

    lines = []
    parse = False
    with open(inputfile) as inlines:
        for line in inlines:
            line = line.rstrip()
            if parse:
                lines.append(line)
            if line.startswith("  #"):
                parse = True
    lines = "\n".join(lines)

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
    table = pd.read_fwf(StringIO(lines), names=header, colspecs=widths, # skiprows=28
                        compression=None, converters=all_str, keep_default_na=False)

    # table modular extensions
    if add_full_chain:
        table = add_dssp_full_chain(table)
        logger.info("DSSP added full chain...")

    table['SS'] = table.SS.fillna('-')
    if add_ss_reduced:
        table = add_dssp_ss_reduced(table)
        logger.info("DSSP added reduced SS...")

    if add_rsa:
        table = add_dssp_rsa(table, method=method)
        logger.info("DSSP added RSA...")

    if add_rsa_class:
        table = add_dssp_rsa_class(table)
        logger.info("DSSP added RSA class...")

    # drop missing residues ("!")  and chain breaks ("!*")
    table = table[table['AA'] != '!']
    table = table[table['AA'] != '!*']

    if reset_res_id:
        table.reset_index(inplace=True)
        table = table.drop(['index'], axis=1)
        table['LINE'] = table.index + 1
        logger.info("DSSP reset residue number...")

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
        logger.info("DSSP table filtered by CHAIN...")

    if chain_full is not None:
        table = row_selector(table, 'CHAIN_FULL', chain_full, method="isin")
        logger.info("DSSP table filtered by CHAIN_FULL...")

    if res is not None:
        table = row_selector(table, 'RES', res, method="isin")
        logger.info("DSSP table filtered by RES...")

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
    def __init__(self, inputfile):
        """
        :param inputfile: Needs to point to a valid DSSP file.
        """
        self.inputfile = inputfile
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
                 add_rsa=True, method="Sander", add_rsa_class=False,
                 reset_res_id=False):
        if excluded is None:
            excluded = self.excluded
        self.data = parse_dssp_from_file(self.inputfile, excluded=excluded,
                                         add_full_chain=add_full_chain,
                                         add_ss_reduced=add_ss_reduced,
                                         add_rsa=add_rsa, method=method,
                                         add_rsa_class=add_rsa_class,
                                         reset_res_id=reset_res_id)
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
            logger.info('No DSSP data parsed...')


class DSSPrunner(object):
    def __init__(self, inputfile, outputfile=None):
        """
        :param inputfile: Needs to point to a valid PDB or mmCIF file.
        :param outputfile: if not provided will use the same file name and
          <.dssp> extension
        """
        self.inputfile = inputfile
        self.outputfile = outputfile
        self.data = None

        if not os.path.isfile(self.inputfile):
            raise IOError("{} not available or could not be read..."
                          "".format(self.inputfile))

        # inputfile needs to be in PDB or mmCIF format
        filename, extension = os.path.splitext(self.inputfile)
        if extension not in ['.pdb', '.ent', '.cif']:
            raise ValueError("{} is expected to be in mmCIF or PDB format..."
                             "".format(self.inputfile))

    def _generate_output(self, run_unbound=False):
        filename, extension = os.path.splitext(self.inputfile)
        if run_unbound:
            self.outputfile = filename + "_unbound.dssp"
        else:
            self.outputfile = filename + ".dssp"

    def _run(self, dssp_bin, run_unbound=False, override=False, save_new_input=False,
             clean_output=True, category='label'):
        """
        If run_unbound=True, the structure is split into its chains and dssp run for
        each independently.
        """

        if not run_unbound:
            cmd = "{} -i {} -o {}".format(dssp_bin, self.inputfile, self.outputfile)
            os.system(cmd)
            if not os.path.isfile(self.outputfile):
                raise IOError("DSSP output not generated for {}".format(self.outputfile))
        else:
            # read the file and get available chains
            r = PDBXreader(inputfile=self.inputfile)
            data = r.atoms(add_res_full=False, add_contacts=False, format_type=None)
            chains = [k for k in data.loc[:, '{}_asym_id'.format(category)].unique()]
            new_inputs = []
            new_outputs = []
            for chain in chains:
                # since len(chain) > 1 (e.g. 'AA' or 'BA') are repetitions and are skipped
                if len(chain) == 1:
                    # write out the new cif file with the current chain
                    filename, extension = os.path.splitext(self.inputfile)
                    outputpdb = filename + '_{}.pdb'.format(chain)
                    new_inputs.append(outputpdb)
                    if not os.path.isfile(outputpdb) or override:
                        w = PDBXwriter(inputfile=None, outputfile=outputpdb)
                        try:
                            w.run(data=data, chain=(chain,), res=None, atom=None,
                                  lines=('ATOM', ), override=override, format_type="pdb")
                        except ValueError:
                            # skipping only HETATM chains or (generally) empty tables
                            continue
                    else:
                        logger.info("PDB for {} already available...".format(outputpdb))
                    # generating the dssp output for the current chain
                    filename, extension = os.path.splitext(self.outputfile)
                    outputdssp = filename + '_{}.dssp'.format(chain)
                    new_outputs.append(outputdssp)
                    if not os.path.isfile(outputdssp) or override:
                        d = DSSPrunner(outputpdb, outputdssp)
                        d.run(override=override, run_unbound=False, save_new_input=False)
                    else:
                        logger.info("DSSP for {} already available...".format(outputdssp))

            # concat the DSSP output to a single file
            lines = ["  # DSSP generated by ProIntVar\n"]
            for outputdssp in new_outputs:
                parse = False
                with open(outputdssp, 'r') as inlines:
                    for line in inlines:
                        if parse:
                            lines.append(line)
                        if line.startswith("  #"):
                            parse = True

                if clean_output:
                    lazy_file_remover(outputdssp)

            with open(self.outputfile, 'w') as outlines:
                outlines.write(''.join(lines))

            if not save_new_input:
                for outputpdb in new_inputs:
                    lazy_file_remover(outputpdb)

    def run(self, run_unbound=False, override=False, save_new_input=False,
            clean_output=True):

        # generate outputfile if missing
        if not self.outputfile:
            self._generate_output(run_unbound=run_unbound)

        if not os.path.exists(self.outputfile) or override or run_unbound:
            if os.path.isfile(config.dssp_bin):
                dssp_bin = config.dssp_bin
            else:
                raise IOError('DSSP executable is not available...')

            # run dssp and generate output
            self._run(dssp_bin, run_unbound=run_unbound, override=override,
                      save_new_input=save_new_input, clean_output=clean_output)
        else:
            logger.info('DSSP for {} already available...'.format(self.outputfile))
        return


if __name__ == '__main__':
    pass
