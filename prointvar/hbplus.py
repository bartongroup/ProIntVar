#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This defines the methods that work with HBPLUS files.

Fábio Madeira, 2017+

"""


import os
import json
import shutil
import pandas as pd

from prointvar.pdbx import PDBXwriter

from prointvar.utils import flash
from prointvar.utils import row_selector
from prointvar.utils import lazy_file_remover
from prointvar.library import hbplus_types

from prointvar.config import config


def parse_hb2_from_file(inputfile, excluded=(), verbose=False):
    """
    Parse lines of the HBPLUS file to get entries for ...

    :param inputfile: path to the HBPLUS file
    :param excluded: option to exclude HBPLUS columns
    :param verbose: boolean
    :return: returns a pandas DataFrame
    """

    if verbose:
        flash("Parsing HBPLUS from lines...")

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

    if not os.path.isfile(inputfile):
        raise IOError("{} not available or could not be read...".format(inputfile))

    # column width descriptors
    header = ("CHAIN_D", "RES_D", "INSCODE_D", "COMP_D", "ATOM_D",
              "CHAIN_A", "RES_A", "INSCODE_A", "COMP_A", "ATOM_A",
              "DIST_DA", "CATEG_DA", "NUM_AAS", "DIST_CA-CA", "ANGLE_D-H-A",
              "DIST_H-A", "ANGLE_H-A-AA", "ANGLE_D-A-AA", "ID")

    widths = ((0, 1), (1, 5), (5, 6), (6, 9), (9, 14), (14, 15), (15, 19),
              (19, 20), (20, 23), (23, 27), (27, 33), (33, 36), (36, 40),
              (40, 45), (45, 52), (52, 58), (58, 64), (64, 70), (70, 78))

    all_str = {key: str for key in header}
    table = pd.read_fwf(inputfile, skiprows=8, names=header, colspecs=widths,
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

    # fix insertion codes
    table.INSCODE_D[table.INSCODE_D == "-"] = "?"
    table.INSCODE_A[table.INSCODE_A == "-"] = "?"

    if excluded is not None:
        assert type(excluded) is tuple
        try:
            table = table.drop(list(excluded), axis=1)
        except ValueError:
            # most likely theses are not in there
            pass

    # enforce some specific column types
    for col in table:
        if col in hbplus_types:
            try:
                table[col] = table[col].astype(hbplus_types[col])
            except ValueError:
                # there are some NaNs in there
                pass

    if table.empty:
        raise ValueError('{} resulted in an empty DataFrame...'.format(inputfile))

    return table


def get_hbplus_selected_from_table(data, chain_A=None, chain_D=None,
                                   res_A=None, res_D=None):
    """
    Utility that filters a pandas DataFrame by the input tuples.

    :param data: pandas DataFrame object
    :param chain_D: (tuple) chain IDs or None (donor)
    :param chain_A: (tuple) chain IDs or None (acceptor)
    :param res_D: (tuple) res IDs or None (donor)
    :param res_A: (tuple) res IDs or None (acceptor)
    :return: returns a modified pandas DataFrame
    """

    # excluding rows
    table = data
    if chain_D is not None:
        table = row_selector(table, 'CHAIN_D', chain_D, method="isin")

    if chain_A is not None:
        table = row_selector(table, 'CHAIN_A', chain_A, method="isin")

    if res_D is not None:
        table = row_selector(table, 'RES_D', res_D, method="isin")

    if res_A is not None:
        table = row_selector(table, 'RES_A', res_A, method="isin")

    return table


class HBPLUSreader(object):
    def __init__(self, inputfile, verbose=False):
        """
        :param inputfile: Needs to point to a valid HBPLUS file.
        :param verbose: boolean
        """
        self.inputfile = inputfile
        self.verbose = verbose
        self.data = None
        self.excluded = ("NUM_AAS", "DIST_CA-CA", "ANGLE_D-H-A",
                         "DIST_H-A", "ANGLE_H-A-AA", "ANGLE_D-A-AA")

        if not os.path.isfile(inputfile):
            raise IOError("{} not available or could not be read...".format(inputfile))

    def read(self, **kwargs):
        return self.residues(**kwargs)

    def residues(self, excluded=None):
        if excluded is None:
            excluded = self.excluded
        self.data = parse_hb2_from_file(self.inputfile, excluded=excluded,
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
            flash('No HBPLUS data parsed...')


class HBPLUSgenerator(object):
    def __init__(self, inputfile, outputfile=None, verbose=False):
        """
        :param inputfile: Needs to point to a valid PDB or mmCIF file.
        :param outputfile: if not provided will use the same file name and
          <.hb2> extension
        :param verbose: boolean
        """
        self.inputfile = inputfile
        self.inputfile_back = inputfile
        self.outputfile = outputfile
        self.verbose = verbose
        self.data = None

        if not os.path.isfile(self.inputfile):
            raise IOError("{} not available or could not be read..."
                          "".format(self.inputfile))

        # inputfile needs to be in PDB or mmCIF format
        filename, extension = os.path.splitext(self.inputfile)
        if extension not in ['.pdb', '.ent', '.cif']:
            raise ValueError("{} is expected to be in mmCIF or PDB format..."
                             "".format(self.inputfile))

    def _generate_output(self, hydro_pdb_out=False):
        filename, extension = os.path.splitext(self.inputfile)
        if hydro_pdb_out:
            self.outputfile = filename + ".h.pdb"
        else:
            self.outputfile = filename + ".hb2"

    def _generate_pdb(self, override=False):
        filename, extension = os.path.splitext(self.inputfile)
        self.inputfile = filename + "_new.pdb"
        w = PDBXwriter(inputfile=self.inputfile_back, outputfile=self.inputfile)
        w.run(format_type="pdb", override=override)

    def _run(self, hbplus_bin, clean_bin=None, run_clean=True,
             clean_output=True, hydro_pdb_out=False):
        # clean and hbplus clip absolute paths to 78 chars
        root, filename = os.path.split(self.inputfile)
        if self.inputfile != filename:
            shutil.copyfile(self.inputfile, filename)
        basename, extension = os.path.splitext(filename)
        basename = "./" + basename
        # output files (some optional)
        output_clean = basename + ".new"
        output_clean_log = basename + ".clean.log"
        output_hbplus = basename + ".hb2"
        output_hbplus_log = basename + ".hbplus.log"
        output_hbplus_h = basename + ".h" # '-o' option
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
            if output_hbplus_h != self.outputfile:
                shutil.copyfile(output_hbplus_h, self.outputfile)
        else:
            if output_hbplus != self.outputfile:
                shutil.copyfile(output_hbplus, self.outputfile)

        if clean_output:
            # remove output files
            if self.inputfile != filename:
                lazy_file_remover(filename)
            if output_hbplus != self.outputfile:
                lazy_file_remover(output_hbplus)
            if output_hbplus_h != self.outputfile:
                lazy_file_remover(output_hbplus_h)
            lazy_file_remover(output_clean)
            lazy_file_remover(output_clean_log)
            lazy_file_remover(output_hbplus_log)
            # other files
            lazy_file_remover("./hbdebug.dat")
            lazy_file_remover("./fort.15")

    def run(self, run_clean=False, hydro_pdb_out=False, override=False,
            clean_output=True, save_new_input=False):

        # generate outputfile if missing
        if not self.outputfile:
            self._generate_output(hydro_pdb_out=hydro_pdb_out)

        if not os.path.exists(self.outputfile) or override:
            if os.path.isfile(config.hbplus_bin):
                hbplus_bin = config.hbplus_bin
                clean_bin = None
                if run_clean and os.path.isfile(config.clean_bin):
                    clean_bin = config.clean_bin
            else:
                raise IOError('HBPLUS executables are not available...')

            # inputfile needs to be in PDB format
            filename, extension = os.path.splitext(self.inputfile)
            if extension == '.cif':
                self._generate_pdb(override=override)

            # run hbplus and generate output - also clean unnecessary output
            self._run(hbplus_bin, clean_bin, run_clean=run_clean,
                      clean_output=clean_output, hydro_pdb_out=hydro_pdb_out)

            # clean the new PDB input file generated
            if not save_new_input:
                if self.inputfile != self.inputfile_back:
                    lazy_file_remover(self.inputfile)

        else:
            flash('HBPLUS for {} already available...'.format(self.outputfile))
        return


if __name__ == '__main__':
    pass
