#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This defines the methods that work with ARPEGGIO files.

Needs the Edited Arpeggio Fork from: https://github.com/biomadeira/arpeggio

Propositions for things to work:
- Generate new PDB file from mmCIF and (even) PDB so that some issues (below) are cleared
    - Remove alternative locations and reset the atom sequential id
      MMCIFreader(<...>).atoms(remove_altloc=True, reset_atom_id=True)

    - Chain IDs need to be a single character so it is safer to run with category='auth'
      MMCIFwriter(<...>).run(category='auth')

    - Remove hydrogens
      MMCIFreader(<...>).atoms(remove_hydrogens=True)

Fábio Madeira, 2017+

"""

import os
import json
import shutil
import pandas as pd
from contextlib import suppress

from prointvar.mmcif import MMCIFwriter
from prointvar.mmcif import MMCIFreader

from prointvar.utils import flash
from prointvar.utils import row_selector
from prointvar.utils import string_split
from prointvar.library import arpeggio_types

from prointvar.config import config


def parse_arpeggio_from_file(inputfile, excluded=(), verbose=False):
    """
    Parse lines of the ARPEGGIO *contacts* file to get entries from...

    :param inputfile: path to the ARPEGGIO file
    :param excluded: option to exclude ARPEGGIO columns
    :param verbose: boolean
    :return: returns a pandas DataFrame
    """

    if verbose:
        flash("Parsing ARPEGGIO from lines...")

    # example lines
    # format documentation at https://github.com/biomadeira/arpeggio
    """
    B/376/ND2	B/374/O	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	4.971	1.901	INTRA_SELECTION
    B/376/CB	B/374/O	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	4.538	1.318	INTRA_SELECTION
    B/376/CA	B/374/O	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	4.202	0.982	INTRA_SELECTION
    B/376/N	B/374/O	0	0	0	0	1	0	0	0	0	0	0	0	0	1	0	3.298	0.228	INTRA_SELECTION
    B/374/C	B/376/CA	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	4.928	1.528	INTRA_SELECTION
    B/374/C	B/376/N	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	3.819	0.569	INTRA_SELECTION
    B/375/NE2	B/377/N	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	4.944	1.844	INTRA_SELECTION
    """

    if not os.path.isfile(inputfile):
        raise IOError("{} not available or could not be read...".format(inputfile))

    # column width descriptors
    header = ("ENTRY_A", "ENTRY_B", "STERIC_CLASH", "COVALENT", "VDW_CLASH", "VDW_INTER",
              "PROXIMAL", "HYDROGEN", "WEAK_HYDROGEN", "HALOGEN", "IONIC",
              "METAL_COMPLEX", "AROMATIC", "HYDROPHOBIC", "CARBONYL", "POLAR", "WEAK_POLAR",
              "DIST", "VDW_DIST", "ENTITIES")

    all_str = {key: str for key in header}
    table = pd.read_csv(inputfile, delim_whitespace=True, low_memory=False,
                        names=header, compression=None, converters=all_str,
                        keep_default_na=False)

    # split ENTRIES into CHAIN, RES, and ATOM
    table = add_arpeggio_res_split(table)

    if excluded is not None:
        assert type(excluded) is tuple
        try:
            table = table.drop(list(excluded), axis=1)
        except ValueError:
            # most likely theses are not in there
            pass

    # enforce some specific column types
    for col in table:
        if col in arpeggio_types:
            try:
                table[col] = table[col].astype(arpeggio_types[col])
            except ValueError:
                # there are some NaNs in there
                pass

    if table.empty:
        raise ValueError('{} resulted in an empty DataFrame...'.format(inputfile))

    return table


def add_arpeggio_res_split(data):
    """
    Utility that adds new columns to the table.
    Adds new columns from the 'full atom description' (e.g B/377/GLU/N).

    adds: 'CHAIN', 'RES', 'COMP', 'ATOM', and 'INSCODE'

    :param data: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """
    table = data

    def get_chain_id(data, key):
        return data[key].split('/')[0]

    def get_res_id(data, key):
        values = string_split(data[key].split('/')[1])
        return values[0]

    def get_comp_id(data, key):
        values = string_split(data[key].split('/')[2])
        return values[0]

    def get_atom_id(data, key):
        return data[key].split('/')[3]

    def get_icode_id(data, key):
        values = string_split(data[key].split('/')[1])
        if len(values) == 2:
            return values[1]
        else:
            return '?'

    table.is_copy = False
    table['CHAIN_A'] = table.apply(get_chain_id, axis=1, args=('ENTRY_A', ))
    table['RES_A'] = table.apply(get_res_id, axis=1, args=('ENTRY_A', ))
    table['COMP_A'] = table.apply(get_comp_id, axis=1, args=('ENTRY_A', ))
    table['ATOM_A'] = table.apply(get_atom_id, axis=1, args=('ENTRY_A', ))
    table['INSCODE_A'] = table.apply(get_icode_id, axis=1, args=('ENTRY_A', ))

    table['CHAIN_B'] = table.apply(get_chain_id, axis=1, args=('ENTRY_B', ))
    table['RES_B'] = table.apply(get_res_id, axis=1, args=('ENTRY_B', ))
    table['COMP_B'] = table.apply(get_comp_id, axis=1, args=('ENTRY_B', ))
    table['ATOM_B'] = table.apply(get_atom_id, axis=1, args=('ENTRY_B', ))
    table['INSCODE_B'] = table.apply(get_icode_id, axis=1, args=('ENTRY_B', ))
    return table


def get_arpeggio_selected_from_table(data, chain_A=None, res_A=None,
                                     chain_B=None, res_B=None):
    """
    Utility that filters a pandas DataFrame by the input tuples.

    :param data: pandas DataFrame object
    :param chain_B: (tuple) chain IDs or None (donor)
    :param chain_A: (tuple) chain IDs or None (acceptor)
    :param res_B: (tuple) res IDs or None (donor)
    :param res_A: (tuple) res IDs or None (acceptor)
    :return: returns a modified pandas DataFrame
    """

    # excluding rows
    table = data

    if chain_A is not None:
        table = row_selector(table, 'CHAIN_A', chain_A, method="isin")

    if chain_B is not None:
        table = row_selector(table, 'CHAIN_B', chain_B, method="isin")

    if res_A is not None:
        table = row_selector(table, 'RES_A', res_A, method="isin")

    if res_B is not None:
        table = row_selector(table, 'RES_B', res_B, method="isin")

    return table


class ARPEGGIOreader(object):
    def __init__(self, inputfile, verbose=False):
        """
        :param inputfile: Needs to point to a valid ARPEGGIO file.
        :param verbose: boolean
        """
        self.inputfile = inputfile
        self.verbose = verbose
        self.data = None
        self.excluded = ("ENTRY_A", "ENTRY_B")

        if not os.path.isfile(inputfile):
            raise IOError("{} not available or could not be read...".format(inputfile))

    def read(self, **kwargs):
        return self.residues(**kwargs)

    def residues(self, excluded=None):
        if excluded is None:
            excluded = self.excluded
        self.data = parse_arpeggio_from_file(self.inputfile, excluded=excluded,
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
            flash('No ARPEGGIO data parsed...')


class ARPEGGIOgenerator(object):
    def __init__(self, inputfile, outputfile=None, verbose=False):
        """
        :param inputfile: Needs to point to a valid PDB or mmCIF file.
        :param outputfile: if not provided will use the same file name and <.contacts> extension
        :param verbose: boolean
        """

        self.inputfile_back = inputfile
        self.inputfile = inputfile
        self.outputfile = outputfile
        self.verbose = verbose
        self.data = None

        # generate outputfile if missing
        self._generate_output()

        if not os.path.isfile(inputfile):
            raise IOError("{} not available or could not be read...".format(inputfile))

        # inputfile needs to be in PDB format
        filename, extension = os.path.splitext(inputfile)
        if extension in ['.pdb', '.ent']:
            self._generate_pdb_from_pdb()
        elif extension in ['.cif']:
            self._generate_pdb_from_mmcif()
        else:
            raise ValueError("{} is expected to be in mmCIF or PDB format..."
                             "".format(inputfile))

    def _generate_output(self):
        if not self.outputfile:
            filename, extension = os.path.splitext(self.inputfile)
            self.outputfile = filename + ".contacts"

    def _generate_pdb_from_mmcif(self):
        filename, extension = os.path.splitext(self.inputfile)
        w = MMCIFwriter(inputfile=self.inputfile, outputfile=filename + ".pdb")
        r = MMCIFreader(inputfile=self.inputfile)
        data = r.atoms(remove_altloc=True, reset_atom_id=True, format_type='mmcif')
        w.run(data=data, format_type="pdb", category="auth")
        self.inputfile = filename + ".pdb"

    def _generate_pdb_from_pdb(self):
        filename, extension = os.path.splitext(self.inputfile)
        w = MMCIFwriter(inputfile=None, outputfile=filename + "_clean.pdb")
        r = MMCIFreader(inputfile=self.inputfile)
        data = r.atoms(remove_altloc=True, reset_atom_id=True, format_type='pdb')
        w.run(data=data, format_type="pdb", category="auth")
        self.inputfile = filename + "_clean.pdb"

    def _run(self, python_path, python_exe, arpeggio_bin, clean_output=True):

        filename, extension = os.path.splitext(self.inputfile)
        output_arpeggio = filename + ".contacts"
        output_bs_contacts = filename + ".bs_contacts"
        output_atomtypes = filename + ".atomtypes"
        
        # run arpeggio
        cmd = 'PYTHONPATH={} {} {} {}'.format(python_path, python_exe,
                                              arpeggio_bin, self.inputfile)
        os.system(cmd)
        if not os.path.isfile(output_arpeggio):
            raise IOError("ARPEGGIO output not generated for {}".format(self.inputfile))

        # mv the automatically generated file -> to the provided outputfile
        if output_arpeggio != self.outputfile:
            shutil.copyfile(output_arpeggio, self.outputfile)

        if clean_output:
            # remove output files
            def lazy_remove_files(filename):
                with suppress(FileNotFoundError):
                    os.remove(filename)

            if output_arpeggio != self.outputfile:
                lazy_remove_files(output_arpeggio)
            lazy_remove_files(output_bs_contacts)
            lazy_remove_files(output_atomtypes)
            # if self.inputfile_back != self.inputfile:
            #     lazy_remove_files(self.inputfile)

    def run(self, override=False, clean_output=True, save_clean_pdb=False):

        if not os.path.exists(self.outputfile) or override:
            if os.path.isfile(config.python_exe) and os.path.exists(config.arpeggio_bin):
                arpeggio_bin = config.arpeggio_bin
                python_exe = config.python_exe
                python_path = config.python_path
            elif (os.path.isfile(config.python_exe_local) and
                    os.path.exists(config.arpeggio_bin_local)):
                arpeggio_bin = config.arpeggio_bin_local
                python_exe = config.python_exe_local
                python_path = config.python_path_local
            else:
                raise IOError('ARPEGGIO executable is not available...')

            # run arpeggio and generate output
            self._run(python_path, python_exe, arpeggio_bin, clean_output=clean_output)

        else:
            flash('ARPEGGIO for {} already available...'.format(self.outputfile))

        if self.inputfile.endswith('_clean.pdb') and not save_clean_pdb:
            os.remove(self.inputfile)
        return


if __name__ == '__main__':
    pass
