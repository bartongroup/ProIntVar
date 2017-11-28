# -*- coding: utf-8 -*-

"""
This defines the methods that work with DSSP files.

Fábio Madeira, 2017+
"""

import os
import logging

from proteofav.structures import mmCIF, PDB, filter_structures
from proteofav.utils import GenericInputs, InputFileHandler

from prointvar.utils import lazy_file_remover
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


def dssp_generate_output_filename(filename, run_unbound=False):
    """
    Little helper function to generate the output filename,
    if it was missing.

    :param filename: path to input file
    :param run_unbound: boolean
    :return: (str)
    """

    filename, extension = os.path.splitext(filename)
    if run_unbound:
        filename_output = filename + "_unbound.dssp"
    else:
        filename_output = filename + ".dssp"
    return filename_output


def run_dssp(filename_input,  filename_output=None, dssp_bin=None,
             run_unbound=False, overwrite=False, save_new_input=False,
             clean_output=True, category='label'):
    """
    If run_unbound=True, the structure is split into its chains and dssp run for
    each independently.

    :param filename_input: Needs to point to a valid PDB or mmCIF file.
    :param filename_output: if not provided will use the same file name and
      <.dssp> extension

    :param filename_input: path to input file
    :param filename_output: path to output file
    :param dssp_bin: path to the DSSP binary/executable (overwrites the config)
    :param run_unbound: boolean
    :param overwrite: boolean
    :param save_new_input: boolean
    :param clean_output: boolean
    :param category: data category to be used as precedence in _atom_site.*_*
        asym_id, seq_id and atom_id
    :return: (side-effects) generates a DSSP file
    """

    InputFileHandler(filename_input)

    # input file needs to be in PDB or mmCIF format
    filename, extension = os.path.splitext(filename_input)
    if extension not in ['.pdb', '.ent', '.cif', '.mmcif']:
        raise ValueError("{} is expected to be in mmCIF or PDB format..."
                         "".format(filename_input))

    # generate output file if missing
    if not filename_output:
        dssp_generate_output_filename(filename=filename_input, run_unbound=run_unbound)

    if not os.path.exists(filename_output) or overwrite or run_unbound:
        if dssp_bin and os.path.isfile(dssp_bin):
            pass
        elif os.path.isfile(config.dssp_bin):
            dssp_bin = config.dssp_bin
        else:
            raise IOError('DSSP executable is not available...')

        if not run_unbound:
            cmd = "{} -i {} -o {}".format(dssp_bin, filename_input, filename_output)
            os.system(cmd)
            if not os.path.isfile(filename_output):
                raise IOError("DSSP output not generated for {}".format(filename_output))
        else:
            # read the file and get available chains
            r = mmCIF.read(filename=filename_input)
            table = filter_structures(r, add_res_full=False, add_contacts=False)
            chains = [k for k in table.loc[:, '{}_asym_id'.format(category)].unique()]
            new_inputs = []
            new_outputs = []
            for chain in chains:
                # since len(chain) > 1 (e.g. 'AA' or 'BA') are repetitions and are skipped
                if len(chain) == 1:
                    # write out the new cif file with the current chain
                    filename, extension = os.path.splitext(filename_input)
                    outputpdb = filename + '_{}.pdb'.format(chain)
                    new_inputs.append(outputpdb)
                    if not os.path.isfile(outputpdb) or overwrite:
                        try:
                            table = filter_structures(table, chains=(chain,), res=None, atoms=None,
                                                      lines='ATOM', category=category)
                            PDB.write(table=table, filename=outputpdb,
                                      output_format="pdb", overwrite=overwrite)
                        except ValueError:
                            # skipping only HETATM chains or (generally) empty tables
                            continue
                    else:
                        logger.info("PDB for %s already available...", outputpdb)
                    # generating the dssp output for the current chain
                    filename, extension = os.path.splitext(filename_output)
                    outputdssp = filename + '_{}.dssp'.format(chain)
                    new_outputs.append(outputdssp)
                    if not os.path.isfile(outputdssp) or overwrite:
                        DSSP.generate(outputpdb, outputdssp, run_unbound=False,
                                      overwrite=overwrite, save_new_input=save_new_input,
                                      clean_output=clean_output, category=category)
                    else:
                        logger.info("DSSP for %s already available...", outputdssp)

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

            with open(filename_output, 'w') as outlines:
                outlines.write(''.join(lines))

            if not save_new_input:
                for outputpdb in new_inputs:
                    lazy_file_remover(outputpdb)
    else:
        logger.info("DSSP for %s already available...", filename_output)
    return


class DSSP(GenericInputs):
    def generate(self, filename_input=None, filename_output=None, **kwargs):
        self.table = run_dssp(filename_input, filename_output, **kwargs)
        return self.table


DSSP = DSSP()

if __name__ == '__main__':
    pass
