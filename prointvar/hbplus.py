# -*- coding: utf-8 -*-

"""

This defines the methods that work with HBPLUS files.

Fábio Madeira, 2017+

"""

import os
import shutil
import logging
import pandas as pd

from proteofav.structures import PDB, mmCIF, filter_structures
from proteofav.utils import InputFileHandler, GenericInputs

from proteofav.utils import row_selector
from proteofav.utils import constrain_column_types
from proteofav.utils import exclude_columns

from prointvar.utils import lazy_file_remover
from prointvar.library import hbplus_types

from prointvar.config import config

logger = logging.getLogger("prointvar")


def parse_hb2_from_file(filename, excluded_cols=()):
    """
    Parse lines of the HBPLUS file to get entries for ...

    :param filename: path to the HBPLUS file
    :param excluded_cols: option to exclude HBPLUS columns
    :return: returns a pandas DataFrame
    """

    logger.info("Parsing HBPLUS from lines...")

    # example lines with some problems
    """
    HBPLUS Hydrogen Bond Calculator v 3.2            Jun 10 11:22:09 BST 2017
    (c) I McDonald, D Naylor, D Jones and J Thornton 1993 All Rights Reserved.
    Citing HBPLUS in publications that use these results is condition of use.
    2PAH <- Brookhaven Code "2pah.new" <- PDB file
    <---DONOR---> <-ACCEPTOR-->    atom                        ^               
    c    i                          cat <-CA-CA->   ^        H-A-AA   ^      H- 
    h    n   atom  resd res      DA  || num        DHA   H-A  angle D-A-AA Bond
    n    s   type  num  typ     dist DA aas  dist angle  dist       angle   num
    A0123-ARG N   A0127-GLU OE1 2.96 MS   4  6.71 166.9  1.98 114.1 116.0     1
    A0421-ILE N   A0123-ARG O   2.61 MM  -1  4.80 167.2  1.63 163.3 164.6     2
    A0127-GLU N   A0124-THR O   3.07 MM   3  5.10 150.6  2.15 100.2 109.4     3
    A0129-ASP N   A0126-GLN O   3.34 MM   3  6.00 154.7  2.41 113.3 116.4     4
    A0130-ARG N   A0127-GLU O   3.43 MM   3  6.08 158.5  2.48 110.6 112.9     5
    A0131-PHE N   A0128-LEU O   3.09 MM   3  5.48 153.3  2.17 106.5 113.6     6
    """

    if not os.path.isfile(filename):
        raise IOError("{} not available or could not be read...".format(filename))

    # column width descriptors
    header = ("CHAIN_D", "RES_D", "INSCODE_D", "COMP_D", "ATOM_D",
              "CHAIN_A", "RES_A", "INSCODE_A", "COMP_A", "ATOM_A",
              "DIST_DA", "CATEG_DA", "NUM_AAS", "DIST_CA-CA", "ANGLE_D-H-A",
              "DIST_H-A", "ANGLE_H-A-AA", "ANGLE_D-A-AA", "ID")

    widths = ((0, 1), (1, 5), (5, 6), (6, 9), (9, 14), (14, 15), (15, 19),
              (19, 20), (20, 23), (23, 27), (27, 33), (33, 36), (36, 40),
              (40, 45), (45, 52), (52, 58), (58, 64), (64, 70), (70, 78))

    all_str = {key: str for key in header}
    table = pd.read_fwf(filename, skiprows=8, names=header, colspecs=widths,
                        compression=None, converters=all_str, keep_default_na=False)

    # fix resid from 0123 to 123
    def fix_res_id(data, key):
        resid = data[key]
        if resid[0] == "0":
            resid = resid[1:]
        return resid

    table.is_copy = False
    table['RES_D'] = table.apply(fix_res_id, axis=1, args=('RES_D',))
    table['RES_A'] = table.apply(fix_res_id, axis=1, args=('RES_A',))
    logger.info("HBPLUS fixed residue ID...")

    # fix insertion codes
    table.INSCODE_D[table.INSCODE_D == "-"] = "?"
    table.INSCODE_A[table.INSCODE_A == "-"] = "?"

    # excluding columns
    table = exclude_columns(table, excluded=excluded_cols)

    # enforce some specific column types
    table = constrain_column_types(table, hbplus_types)

    if table.empty:
        raise ValueError('{} resulted in an empty DataFrame...'.format(filename))

    return table


def filter_hbplus(table, chain_A=None, chain_D=None,
                  res_A=None, res_D=None):
    """
    Utility that filters a pandas DataFrame by the input tuples.

    :param table: pandas DataFrame object
    :param chain_D: (tuple) chain IDs or None (donor)
    :param chain_A: (tuple) chain IDs or None (acceptor)
    :param res_D: (tuple) res IDs or None (donor)
    :param res_A: (tuple) res IDs or None (acceptor)
    :return: returns a modified pandas DataFrame
    """

    # excluding rows
    if chain_D is not None:
        table = row_selector(table, 'CHAIN_D', chain_D)
        logger.info("HBPLUS table filtered by CHAIN_D...")

    if chain_A is not None:
        table = row_selector(table, 'CHAIN_A', chain_A)
        logger.info("HBPLUS table filtered by CHAIN_A...")

    if res_D is not None:
        table = row_selector(table, 'RES_D', res_D)
        logger.info("HBPLUS table filtered by RES_D...")

    if res_A is not None:
        table = row_selector(table, 'RES_A', res_A)
        logger.info("HBPLUS table filtered by RES_A...")

    return table


def hbplus_generate_output_filename(filename, hydro_pdb_out=False):
    """
    Little helper function to generate the output filename,
    if it was missing.

    :param filename: path to input file
    :param hydro_pdb_out: boolean
    :return: (str)
    """

    filename, extension = os.path.splitext(filename)
    if hydro_pdb_out:
        filename_output = filename + ".h.pdb"
    else:
        filename_output = filename + ".hb2"
    return filename_output


def _run_hbplus(filename_input, filename_output,
                hbplus_bin, clean_bin=None, run_clean=True,
                clean_output=True, hydro_pdb_out=False):
    # clean and hbplus clip absolute paths to 78 chars
    root, filename = os.path.split(filename_input)
    if filename_input != filename:
        shutil.copyfile(filename_input, filename)
    basename, extension = os.path.splitext(filename)
    basename = os.path.join(".", basename)
    # output files (some optional)
    output_clean = basename + ".new"
    output_clean_log = basename + ".clean.log"
    output_hbplus = basename + ".hb2"
    output_hbplus_log = basename + ".hbplus.log"
    output_hbplus_h = basename + ".h"  # '-o' option
    # output_hbplus_hhb = basename + ".hhb" # '-L' option

    # it is recommended that the clean tool is run before calling hbplus
    if clean_bin and run_clean:
        cmd = "{} <<< {} > {}".format(clean_bin, filename, output_clean_log)
        os.system(cmd)
        if not os.path.isfile(output_clean):
            raise IOError("Clean (HBPLUS) output not generated for {}"
                          "".format(output_clean))
    else:
        output_clean = filename

    # run hbplus
    cmd = "{} -R {} -o > {}".format(hbplus_bin, output_clean, output_hbplus_log)
    os.system(cmd)
    if not os.path.isfile(output_hbplus):
        raise IOError("HBPLUS output not generated for {}".format(output_hbplus))

    # mv the automatically generated file -> to the provided outputfile
    if hydro_pdb_out:
        # this means all we care about is the '*.h.pdb' file
        if output_hbplus_h != filename_output:
            shutil.copyfile(output_hbplus_h, filename_output)
    else:
        if output_hbplus != filename_output:
            shutil.copyfile(output_hbplus, filename_output)

    if clean_output:
        # remove output files
        if filename_input != filename:
            lazy_file_remover(filename)
        if output_hbplus != filename_output:
            lazy_file_remover(output_hbplus)
        if output_hbplus_h != filename_output:
            lazy_file_remover(output_hbplus_h)
        lazy_file_remover(output_clean)
        lazy_file_remover(output_clean_log)
        lazy_file_remover(output_hbplus_log)
        # other files
        lazy_file_remover("./hbdebug.dat")
        lazy_file_remover("./fort.15")


def run_hbplus(filename_input, filename_output=None,
               run_clean=False, hydro_pdb_out=False,
               clean_output=True, save_new_input=False, overwrite=False):
    """
    Runs HBPLUS to add explicit Hydrogen atoms to a PDB
    structure.

    :param filename_input: path to input file.
      Needs to point to a valid PDB or mmCIF file.
    :param filename_output: path to output file
      if not provided will use the same file name and <*.hb2> extension
    :param run_clean: boolean
    :param hydro_pdb_out: boolean
    :param clean_output: boolean
    :param save_new_input: boolean
    :param overwrite: boolean
    :return: Runs HBPLUS on the provided PDB structure
    """

    InputFileHandler(filename_input)

    # inputfile needs to be in PDB or mmCIF format
    filename, extension = os.path.splitext(filename_input)
    if extension not in ['.pdb', '.ent', '.cif']:
        raise ValueError("{} is expected to be in mmCIF or PDB format..."
                         "".format(filename_input))

    # generate outputfile if missing
    if not filename_output:
        filename_output = hbplus_generate_output_filename(filename=filename_input,
                                                          hydro_pdb_out=hydro_pdb_out)

    if not os.path.exists(filename_output) or overwrite:
        if os.path.isfile(config.hbplus_bin):
            hbplus_bin = config.hbplus_bin
            clean_bin = None
            if run_clean and os.path.isfile(config.clean_bin):
                clean_bin = config.clean_bin
        else:
            raise IOError('HBPLUS executables are not available...')

        # inputfile needs to be in PDB format
        filename_input_back = filename_input
        filename, extension = os.path.splitext(filename_input)
        if extension == '.cif':
            filename, extension = os.path.splitext(filename_input)
            filename_input = filename + "_new.pdb"
            r = mmCIF.read(filename=filename_input_back)
            table = filter_structures(r, add_res_full=False, add_contacts=False,
                                      lines='ATOM', category='auth')
            PDB.write(table=table, filename=filename_input,
                      output_format="pdb", overwrite=overwrite)

        # run hbplus and generate output - also clean unnecessary output
        _run_hbplus(filename_input, filename_output,
                    hbplus_bin, clean_bin, run_clean=run_clean,
                    clean_output=clean_output, hydro_pdb_out=hydro_pdb_out)

        # clean the new PDB input file generated
        if not save_new_input:
            if filename_input != filename_input_back:
                lazy_file_remover(filename_input)
    else:
        logger.info("HBPLUS for %s already available...", filename_output)


class _HBPLUS(GenericInputs):
    def read(self, filename=None, **kwargs):
        self.table = parse_hb2_from_file(filename=filename, **kwargs)
        return self.table

    def run(self, filename_input=None, filename_output=None, **kwargs):
        self.table = run_hbplus(filename_input=filename_input,
                                filename_output=filename_output, **kwargs)
        return self.table


HBPLUS = _HBPLUS()
