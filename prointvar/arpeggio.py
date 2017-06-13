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
from operator import attrgetter
from collections import Counter
from collections import namedtuple

from prointvar.mmcif import MMCIFwriter
from prointvar.mmcif import MMCIFreader

from prointvar.utils import flash
from prointvar.utils import lazy_file_remover
from prointvar.utils import row_selector
from prointvar.utils import string_split
from prointvar.library import arpeggio_types
from prointvar.library import aa_codes_3to1_extended

from prointvar.config import config


def parse_arpeggio_from_file(inputfile, excluded=(), add_res_split=True,
                             add_group_pdb=True, verbose=False):
    """
    Parse lines of the ARPEGGIO *contacts* file to get entries from...

    :param inputfile: path to the ARPEGGIO file
    :param excluded: option to exclude ARPEGGIO columns
    :param add_res_split: (boolean) splits ENTRY_* into 'CHAIN', 'ATOM', 'COMP', 'INSCODE'
    :param add_group_pdb: (boolean) adds a column with 'ATOM' or 'HETATM'
    :param verbose: boolean
    :return: returns a pandas DataFrame
    """

    if verbose:
        flash("Parsing ARPEGGIO from lines...")

    # example lines
    # format documentation at https://github.com/biomadeira/arpeggio
    """
    B/376/ASN/ND2	B/374/ILE/O	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	4.971	1.901	INTRA_SELECTION
    B/376/ASN/CB	B/374/ILE/O	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	4.538	1.318	INTRA_SELECTION
    B/376/ASN/CA	B/374/ILE/O	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	4.202	0.982	INTRA_SELECTION
    B/376/ASN/N	    B/374/ILE/O	0	0	0	0	1	0	0	0	0	0	0	0	0	1	0	3.298	0.228	INTRA_SELECTION
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
    if add_res_split:
        table = add_arpeggio_res_split(table)

    if add_group_pdb:
        table = add_arpeggio_group_pdb(table)

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
    Also flips the order or the contacts (i.e. from B->A to A->B)
    Atom-atom pairs are only provided A->B.

    adds: 'CHAIN', 'RES', 'RES_FULL', 'COMP', 'ATOM', and 'INSCODE'

    :param data: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """
    table = data

    # get most frequent chain in ENTRY_A and use it to define the inter. direction
    # for multiple chains use decreasing frequency
    chains_a = [v.split('/')[0] for v in table['ENTRY_A'].tolist()]
    Chains = namedtuple('Chains', 'key numb')
    freqs = [Chains(key=k, numb=n) for k, n in zip(Counter(chains_a).keys(),
                                                   Counter(chains_a).values())]
    freqs_dict = {ent.key: ent.numb for ent in sorted(freqs, key=attrgetter('numb'),
                                                      reverse=True)}

    def get_chain_id(entry):
        return entry.split('/')[0]

    def get_res_id(entry):
        values = string_split(entry.split('/')[1])
        return values[0]

    def get_icode_id(entry):
        values = string_split(entry.split('/')[1])
        if len(values) == 2:
            return values[1]
        else:
            return '?'

    def get_res_full_id(entry):
        return entry.split('/')[1]

    def get_comp_id(entry):
        values = string_split(entry.split('/')[2])
        return values[0]

    def get_atom_id(entry):
        return entry.split('/')[3]

    chain_a = []
    res_a = []
    inscode_a = []
    res_full_a = []
    comp_a = []
    atom_a = []
    chain_b = []
    res_b = []
    inscode_b = []
    res_full_b = []
    comp_b = []
    atom_b = []
    for ix in table.index:
        chain1 = get_chain_id(table.loc[ix, 'ENTRY_A'])
        chain2 = get_chain_id(table.loc[ix, 'ENTRY_B'])
        # sort the A->B order based on the frequency of each chain ID
        if freqs_dict[chain1] >= freqs_dict[chain2]:
            entry1 = 'ENTRY_A'
            entry2 = 'ENTRY_B'
        else:
            entry1 = 'ENTRY_B'
            entry2 = 'ENTRY_A'
        chain_a.append(get_chain_id(table.loc[ix, entry1]))
        res_a.append(get_res_id(table.loc[ix, entry1]))
        inscode_a.append(get_icode_id(table.loc[ix, entry1]))
        res_full_a.append(get_res_full_id(table.loc[ix, entry1]))
        comp_a.append(get_comp_id(table.loc[ix, entry1]))
        atom_a.append(get_atom_id(table.loc[ix, entry1]))
        chain_b.append(get_chain_id(table.loc[ix, entry2]))
        res_b.append(get_res_id(table.loc[ix, entry2]))
        inscode_b.append(get_icode_id(table.loc[ix, entry2]))
        res_full_b.append(get_res_full_id(table.loc[ix, entry2]))
        comp_b.append(get_comp_id(table.loc[ix, entry2]))
        atom_b.append(get_atom_id(table.loc[ix, entry2]))

    assert len(chain_a) == len(table)
    table['CHAIN_A'] = chain_a
    table['RES_A'] = res_a
    table['INSCODE_A'] = inscode_a
    table['RES_FULL_A'] = res_a
    table['COMP_A'] = comp_a
    table['ATOM_A'] = atom_a
    table['CHAIN_B'] = chain_b
    table['RES_B'] = res_b
    table['INSCODE_B'] = inscode_b
    table['RES_FULL_B'] = res_b
    table['COMP_B'] = comp_b
    table['ATOM_B'] = atom_b
    return table


def add_arpeggio_group_pdb(data):
    """
    Utility that adds new columns to the table.
    Adds new columns: 'GROUP_A' and 'GROUP_B'

    :param data: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    table = data

    def get_group_pdb(data, key):
        # this maybe a bit crude way of assigning this value
        if data[key] in aa_codes_3to1_extended:
            return 'ATOM'
        else:
            return 'HETATM'

    table.is_copy = False
    table['GROUP_A'] = table.apply(get_group_pdb, axis=1, args=('COMP_A', ))
    table['GROUP_B'] = table.apply(get_group_pdb, axis=1, args=('COMP_B', ))
    return table


def get_arpeggio_selected_from_table(data, chain_A=None, chain_B=None,
                                     res_A=None, res_B=None,
                                     res_full_A=None, res_full_B=None,
                                     group_A=None, group_B=None):
    """
    Utility that filters a pandas DataFrame by the input tuples.

    :param data: pandas DataFrame object
    :param chain_A: (tuple) chain IDs or None (donor)
    :param chain_B: (tuple) chain IDs or None (acceptor)
    :param res_A: (tuple) res IDs or None (donor)
    :param res_B: (tuple) res IDs or None (acceptor)
    :param res_full_A: (tuple) res IDs + inscode or None (donor)
    :param res_full_B: (tuple) res IDs + inscode or None (acceptor)
    :param group_A: (tuple) group_PDB or None (donor)
    :param group_B: (tuple) group_PDB or None (acceptor)
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

    if res_full_A is not None:
        table = row_selector(table, 'RES_FULL_A', res_full_A, method="isin")

    if res_full_B is not None:
        table = row_selector(table, 'RES_FULL_B', res_full_B, method="isin")

    if group_A is not None:
        table = row_selector(table, 'GROUP_A', res_A, method="equals")

    if group_B is not None:
        table = row_selector(table, 'GROUP_B', res_B, method="equals")

    return table


def residues_aggregation(data, agg_method='unique'):
    """
    Gets the contacts res-by-res, instead of atom-atom.

    :param data: pandas DataFrame object
    :param agg_method: current values: 'first', 'unique', and 'minimum'
    :return: returns a modified pandas DataFrame
    """
    table = data
    agg_generic = agg_method
    agg_method_origin = agg_method
    agg_cols = ['CHAIN_A', 'RES_FULL_A', 'COMP_A', 'CHAIN_B', 'RES_FULL_B', 'COMP_B']
    if agg_method not in ['first', 'unique', 'minimum', 'maximum']:
        raise ValueError('Method {} is not currently implemented...'
                         ''.format(agg_method))

    if agg_method != 'minimum' and agg_method != 'maximum':
        columns_to_agg = {col: agg_generic if table[col].dtype == 'object' else agg_method
                          for col in table.columns if col not in agg_cols}
    else:
        if agg_method_origin == 'minimum':
            # need the table sort by distance first: ascending
            table = table.sort_values(["DIST", "VDW_DIST"], ascending=[True, True])
            table.reset_index(inplace=True)
            table = table.drop(['index'], axis=1)
            agg_generic = 'first'
            agg_method = 'max'
            columns_to_agg = {col: agg_generic if table[col].dtype == 'object' else agg_method
                              for col in table.columns if col not in agg_cols}
            columns_to_agg['DIST'] = 'min'
            columns_to_agg['VDW_DIST'] = 'min'
        elif agg_method_origin == 'maximum':
            # need the table sort by distance first: descending
            table = table.sort_values(["DIST", "VDW_DIST"], ascending=[False, False])
            table.reset_index(inplace=True)
            table = table.drop(['index'], axis=1)
            agg_generic = 'first'
            agg_method = 'max'
            columns_to_agg = {col: agg_generic if table[col].dtype == 'object' else agg_method
                              for col in table.columns if col not in agg_cols}
    table = table.groupby(by=agg_cols, as_index=False).agg(columns_to_agg)
    return table


def interaction_modes(data, int_mode='inter-chain'):
    """
    Gets the contacts filtered base on the entities that are interacting.
    Interaction modes possible: inter-chain, intra-chain and hetatm.

    :param data: pandas DataFrame object
    :param int_mode: current values: 'inter-chain', 'intra-chain' and 'hetatm'
    :return: returns a modified pandas DataFrame
    """

    table = data
    if int_mode == 'inter-chain':
        table = table.loc[table['CHAIN_A'] != table['CHAIN_B']]
    elif int_mode == 'intra-chain':
        table = table.loc[table['CHAIN_A'] == table['CHAIN_B']]
    elif int_mode == 'hetatm':
        table = table.loc[(table['GROUP_A'] == 'HETATM') | (table['GROUP_B'] == 'HETATM')]
    else:
        raise ValueError('Interaction mode {} is not currently implemented...'
                         ''.format(int_mode))
    # FIXME optional?
    table.reset_index(inplace=True)
    table = table.drop(['index'], axis=1)
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
        self.excluded = ("ENTRY_A", "ENTRY_B", "ENTITIES")

        if not os.path.isfile(inputfile):
            raise IOError("{} not available or could not be read...".format(inputfile))

    def read(self, **kwargs):
        return self.contacts(**kwargs)

    def contacts(self, excluded=None, add_res_split=True, add_group_pdb=True,
                 residue_agg=False, agg_method='minimum',
                 int_filter=False, int_mode='inter-chain'):

        if excluded is None:
            excluded = self.excluded
        self.data = parse_arpeggio_from_file(self.inputfile, excluded=excluded,
                                             add_res_split=add_res_split,
                                             add_group_pdb=add_group_pdb,
                                             verbose=self.verbose)

        if residue_agg:
            self.data = residues_aggregation(self.data, agg_method=agg_method)

        if int_filter:
            self.data = interaction_modes(self.data, int_mode=int_mode)

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
            if output_arpeggio != self.outputfile:
                lazy_file_remover(output_arpeggio)
            lazy_file_remover(output_bs_contacts)
            lazy_file_remover(output_atomtypes)

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
            lazy_file_remover(self.inputfile)
        return


if __name__ == '__main__':
    pass
