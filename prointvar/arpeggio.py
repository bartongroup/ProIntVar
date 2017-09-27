#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This defines the methods that work with ARPEGGIO files.

Needs the Edited Arpeggio Fork from: https://bitbucket.org/biomadeira/arpeggio

Propositions for things to work:
- Generate new PDB file from mmCIF and (even) PDB so that some issues (below) are cleared
    - Remove alternative locations and reset the atom sequential id
      PDBXreader(<...>).atoms(remove_altloc=True, reset_atom_id=True)

    - Chain IDs need to be a single character so it is safer to run with category='auth'
      PDBXwriter(<...>).run(category='auth')

    - Remove hydrogens
      PDBXreader(<...>).atoms(remove_hydrogens=True)

Fábio Madeira, 2017+

"""

import os
import json
import shutil
import logging
import pandas as pd
from operator import attrgetter
from collections import Counter
from collections import namedtuple

from prointvar.pdbx import PDBXwriter
from prointvar.pdbx import PDBXreader
from prointvar.reduce import REDUCErunner
from prointvar.hbplus import HBPLUSrunner

from prointvar.utils import lazy_file_remover
from prointvar.utils import row_selector
from prointvar.utils import string_split
from prointvar.utils import constrain_column_types
from prointvar.utils import exclude_columns
from prointvar.library import arpeggio_types
from prointvar.library import arpeggio_col_renames

from prointvar.config import config

logger = logging.getLogger("prointvar")


def parse_arpeggio_from_file(inputfile, excluded=(), add_res_split=True,
                             parse_special=False):
    """
    Parse lines of the ARPEGGIO *contacts* file to get entries from...

    :param inputfile: path to the ARPEGGIO file
    :param excluded: option to exclude ARPEGGIO columns
    :param add_res_split: (boolean) splits ENTRY_* into 'CHAIN', 'ATOM', 'RES', 'INSCODE'
    :param parse_special: (boolean) tries to parse special contact types
    :return: returns a pandas DataFrame
    """

    logger.info("Parsing ARPEGGIO from lines...")

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
    if add_res_split:
        table = add_arpeggio_res_split(table)
        logger.info("Added Arpeggio full chain, res and atom information...")

    if parse_special:
        # tries to parse *.amam, *.amri, *.ari and *.ri
        table = add_special_cont_types(inputfile, table)
        logger.info("Parsed special contact-types...")

    # excluding columns
    table = exclude_columns(table, excluded=excluded)

    # enforce some specific column types
    table = constrain_column_types(table, arpeggio_types)

    if table.empty:
        raise ValueError('{} resulted in an empty DataFrame...'.format(inputfile))

    return table


def parse_arpeggio_spec_from_file(inputfile, excluded=(), add_res_split=True,
                                  int_type="res-res"):
    """
    Parse lines of the ARPEGGIO *.amam*, *.amri*, *.ri*  and *.ari* files..

    :param inputfile: path to the ARPEGGIO file
    :param excluded: option to exclude ARPEGGIO columns
    :param add_res_split: (boolean) splits ENTRY_* into 'CHAIN', 'ATOM', 'RES', 'INSCODE'
    :param int_type: (str) 'res-res' or 'atom-res'
    :return: returns a pandas DataFrame
    """

    # example lines
    # format documentation at https://github.com/biomadeira/arpeggio
    """
    # amam (res-res)
    37	A/159/	[23.1595, 16.125999, 35.171997]	40	A/161/	[19.171001, 15.512501, 36.269501]	AMIDEAMIDE	INTER_RESIDUE	INTRA_SELECTION
    44	A/165/	[14.835501, 12.307, 32.723]	46	A/167/	[18.681, 8.6389999, 31.563]	AMIDEAMIDE	INTER_RESIDUE	INTRA_SELECTION
    46	A/167/	[18.681, 8.6389999, 31.563]	44	A/165/	[14.835501, 12.307, 32.723]	AMIDEAMIDE	INTER_RESIDUE	INTRA_SELECTION
    94	A/212/	[6.1964998, 12.783, 58.588501]	96	A/214/	[9.8845005, 12.103001, 60.452499]	AMIDEAMIDE	INTER_RESIDUE	INTRA_SELECTION
    # amri (res-res)
    3	A/121/	[-2.691, 14.5735, 27.1605]	64	A/120/	[0.17149999999999999, 14.428833333333333, 24.741499999999995]	AMIDERING	INTER_RESIDUE	INTRA_SELECTION
    3	A/121/	[-2.691, 14.5735, 27.1605]	151	A/120/	[-0.33888888888888885, 14.43488888888889, 23.960888888888885]	AMIDERING	INTER_RESIDUE	INTRA_SELECTION
    6	A/124/	[-1.8935, 8.1655006, 31.7075]	81	A/260/	[-3.0968333333333331, 11.735499999999998, 33.67583333333333]	AMIDERING	INTER_RESIDUE	INTRA_SELECTION
    21	A/143/	[23.228001, 32.715, 33.210999]	80	A/149/	[26.669499999999999, 30.518000000000004, 32.972000000000001]	AMIDERING	INTER_RESIDUE	INTRA_SELECTION
    # ri (res-res)
    5	B/208/	[-39.389400000000002, 39.8444, 24.696200000000005]	67	B/204/	[-33.653833333333331, 40.87116666666666, 24.8565]	OT	INTER_RESIDUE	INTRA_SELECTION
    10	B/187/	[-49.272400000000005, 57.387, 31.215600000000002]	71	B/191/	[-46.081333333333333, 55.690666666666665, 27.043833333333332]	OT	INTER_RESIDUE	INTRA_SELECTION
    15	A/264/	[5.4510000000000005, 15.938800000000001, 35.380200000000002]	76	A/131/	[3.509666666666666, 16.270499999999998, 30.280166666666663]	OT	INTER_RESIDUE	INTRA_SELECTION
    25	A/326/	[2.8020000000000005, 33.325000000000003, 41.7986]	93	A/377/	[2.335, 35.671666666666667, 38.610999999999997]	FF	INTER_RESIDUE	INTRA_SELECTION
    25	A/326/	[2.8020000000000005, 33.325000000000003, 41.7986]	108	A/327/	[-1.7431666666666665, 36.167333333333332, 42.487833333333334]	EF	INTER_RESIDUE	INTRA_SELECTION
    # ari (atom-res)
    B/222/N	3	B/220/	[-51.032000000000011, 56.128399999999999, 22.770000000000003]	['DONORPI']	INTER_RESIDUE	INTRA_SELECTION
    B/200/CB	8	B/201/	[-29.119000000000007, 56.571000000000005, 27.233799999999999]	['CARBONPI']	INTER_RESIDUE	INTRA_SELECTION
    B/200/CG2	8	B/201/	[-29.119000000000007, 56.571000000000005, 27.233799999999999]	['CARBONPI']	INTER_RESIDUE	INTRA_SELECTION
    B/224/N	10	B/187/	[-49.272400000000005, 57.387, 31.215600000000002]	['DONORPI']	INTER_RESIDUE	INTRA_SELECTION
    B/224/CA	10	B/187/	[-49.272400000000005, 57.387, 31.215600000000002]	['CARBONPI']	INTER_RESIDUE	INTRA_SELECTION
    B/191/CB	10	B/187/	[-49.272400000000005, 57.387, 31.215600000000002]	['CARBONPI']	INTER_RESIDUE	INTRA_SELECTION
    B/224/CG1	10	B/187/	[-49.272400000000005, 57.387, 31.215600000000002]	['CARBONPI']	INTER_RESIDUE	INTRA_SELECTION
    """

    if not os.path.isfile(inputfile):
        raise IOError("{} not available or could not be read...".format(inputfile))

    # column width descriptors
    if int_type == "res-res":
        header = ("ID_A", "ENTRY_A", "COORDS_A", "ID_B", "ENTRY_B", "COORDS_B",
                  "CONT_TYPE", "INT_TYPE", "SELECT")
    elif int_type == "atom-res":
        header = ("ENTRY_A", "ID_B", "ENTRY_B", "COORDS_B",
                  "CONT_TYPE", "INT_TYPE", "SELECT")
    else:
        raise ValueError('Input Type {} is not currently implemented...'
                         ''.format(int_type))

    all_str = {key: str for key in header}
    table = pd.read_csv(inputfile, sep='\t', low_memory=False,
                        names=header, compression=None, converters=all_str,
                        keep_default_na=False)

    # split ENTRIES into CHAIN, RES, and ATOM
    if add_res_split:
        table = add_arpeggio_res_split(table)

    # excluding columns
    table = exclude_columns(table, excluded=excluded)

    if table.empty:
        raise ValueError('{} resulted in an empty DataFrame...'.format(inputfile))

    return table


def add_arpeggio_res_split(data):
    """
    Utility that adds new columns to the table.
    Adds new columns from the 'full atom description' (e.g B/377/GLU/N).
    Also flips the order or the contacts (i.e. from B->A to A->B)
    Atom-atom pairs are only provided A->B.

    adds: 'CHAIN', 'RES', 'RES_FULL', 'ATOM', and 'INSCODE'

    :param data: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """
    table = data

    #FIXME: This is no longer used... I don't think it was working as intended in the first place
    # # get most frequent chain in ENTRY_A and use it to define the inter. direction
    # # for multiple chains use decreasing frequency
    # chains_a = [v.split('/')[0] for v in table['ENTRY_A'].tolist()]
    # Chains = namedtuple('Chains', 'key numb')
    # freqs = [Chains(key=k, numb=n) for k, n in zip(Counter(chains_a).keys(),
    #                                                Counter(chains_a).values())]
    # freqs_dict = {ent.key: ent.numb for ent in sorted(freqs, key=attrgetter('numb'),
    #                                                   reverse=True)}

    def _parse_arpeggio_atom(s):
        # Note that this method does not reorder contact direction...
        table = s.str.split('/', expand=True)

        # Rename columns
        suffix = s.name[-1]
        column_name_dict = {k: v.format(suffix) for k, v in enumerate(['CHAIN_{}', 'RES_FULL_{}', 'ATOM_{}', 'X_{}'])}
        table = table.rename(columns=column_name_dict)

        # Parse RES_FULL_X to DataFrame
        res_full = table['RES_FULL_{}'.format(suffix)].str.split(r'(\d+)', expand=True)
        res_full = res_full.drop(0, axis=1)  # the split always returns 3 columns, we'll never need the first

        # Format residues with no insertion code
        res_full.loc[:, 2] = res_full.loc[:, 2].str.replace('', '?')

        # Rename columns
        column_name_dict = {k: v.format(suffix) for k, v in enumerate(['RES_{}', 'INSCODE_{}'], 1)}
        res_full = res_full.rename(columns=column_name_dict)

        table = table.join(res_full)

        return table

    if not table.empty:
        tbA = _parse_arpeggio_atom(table['ENTRY_A'])
        tbB = _parse_arpeggio_atom(table['ENTRY_B'])
        table = table.join(tbA.join(tbB))
    return table


def add_special_cont_types(inputfile, data):
    """
    Tries to parse *.amam, *.amri, *.ari and *.ri and add these data as new
     columns.

    :param inputfile: path to the ARPEGGIO file
    :param data: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """
    table = data

    excluded = ["ID_A", "ENTRY_A", "COORDS_A", "ID_B", "ENTRY_B", "COORDS_B",
                "CONT_TYPE", "INT_TYPE", "SELECT"]

    filename, extension = os.path.splitext(inputfile)
    input_amam = filename + ".amam"
    amam = parse_arpeggio_spec_from_file(input_amam, excluded=tuple(excluded),
                                         int_type="res-res")
    table = add_contact_info(table, amam, col_name="Amide-Amide", int_type="res-res")

    input_amri = filename + ".amri"
    amri = parse_arpeggio_spec_from_file(input_amri, excluded=tuple(excluded),
                                         int_type="res-res")
    table = add_contact_info(table, amri, col_name="Amide-Aromatic", int_type="res-res")

    input_ri = filename + ".ri"
    ri = parse_arpeggio_spec_from_file(input_ri, excluded=tuple(excluded),
                                       int_type="res-res")
    table = add_contact_info(table, ri, col_name="Aromatic-Aromatic", int_type="res-res")

    excluded.remove("ID_A")
    excluded.remove("COORDS_A")
    input_ari = filename + ".ari"
    ari = parse_arpeggio_spec_from_file(input_ari, excluded=tuple(excluded),
                                        int_type="atom-res")
    table = add_contact_info(table, ari, col_name="Atom-Ring", int_type="atom-res")

    return table


def add_contact_info(data, info, col_name="Amide-Amide", int_type="res-res"):
    """
    Identifies residues or atoms that are observed in both data and info,
     adding a new column with '1' or '0' as a proxy of observation.

    :param data: pandas DataFrame object
    :param info: pandas DataFrame object
    :param col_name: (str) name of the new column (contact type)
    :param int_type: (str) 'res-res' or 'atom-res'
    :return: returns a modified pandas DataFrame
    """

    table = data
    a2b = ["CHAIN_A", "RES_FULL_A", "ATOM_A", "CHAIN_B", "RES_FULL_B"]
    b2a = ["CHAIN_B", "RES_FULL_B", "ATOM_B", "CHAIN_A", "RES_FULL_A"]
    if int_type == 'res-res':
        a2b.remove("ATOM_A")
        b2a.remove("ATOM_B")

    def new_column(data, keys):
        return "_".join([data[k] for k in keys])

    table.is_copy = False
    table['a2b'] = table.apply(new_column, axis=1, args=(a2b,))
    table['b2a'] = table.apply(new_column, axis=1, args=(b2a,))
    info['a2b'] = info.apply(new_column, axis=1, args=(a2b,))

    conts = []
    for ix in info.index:
        # always A->B
        conts.append(info.loc[ix, 'a2b'])

    values = ['0'] * len(table.index)
    for ix in table.index:
        # if '_'.join([v for v in table.loc[ix, a2b]]) in conts:
        if table.loc[ix, 'a2b'] in conts:
            values[ix] = '1'
        elif table.loc[ix, 'b2a'] in conts:
            values[ix] = '1'
        else:
            values[ix] = '0'

    table[col_name] = values
    table = table.drop(['a2b', 'b2a'], axis=1)
    return table


def get_arpeggio_selected_from_table(data, chain_A=None, chain_B=None,
                                     res_A=None, res_B=None,
                                     res_full_A=None, res_full_B=None,
                                     atom_A=None, atom_B=None):
    """
    Utility that filters a pandas DataFrame by the input tuples.

    :param data: pandas DataFrame object
    :param chain_A: (tuple) chain IDs or None (donor)
    :param chain_B: (tuple) chain IDs or None (acceptor)
    :param res_A: (tuple) res IDs or None (donor)
    :param res_B: (tuple) res IDs or None (acceptor)
    :param res_full_A: (tuple) res IDs + inscode or None (donor)
    :param res_full_B: (tuple) res IDs + inscode or None (acceptor)
    :param atom_A: (tuple) atom IDs or None (donor)
    :param atom_B: (tuple) atom IDs or None (acceptor)
    :return: returns a modified pandas DataFrame
    """

    # excluding rows
    table = data

    if chain_A is not None:
        table = row_selector(table, 'CHAIN_A', chain_A)
        logger.info("Arpeggio table filtered by CHAIN_A...")

    if chain_B is not None:
        table = row_selector(table, 'CHAIN_B', chain_B)
        logger.info("Arpeggio table filtered by CHAIN_B...")

    if res_A is not None:
        table = row_selector(table, 'RES_A', res_A)
        logger.info("Arpeggio table filtered by RES_A...")

    if res_B is not None:
        table = row_selector(table, 'RES_B', res_B)
        logger.info("Arpeggio table filtered by RES_B...")

    if res_full_A is not None:
        table = row_selector(table, 'RES_FULL_A', res_full_A)
        logger.info("Arpeggio table filtered by RES_FULL_A...")

    if res_full_B is not None:
        table = row_selector(table, 'RES_FULL_B', res_full_B)
        logger.info("Arpeggio table filtered by RES_FULL_B...")

    if atom_A is not None:
        table = row_selector(table, 'ATOM_A', atom_A)
        logger.info("Arpeggio table filtered by ATOM_A...")

    if atom_B is not None:
        table = row_selector(table, 'ATOM_B', atom_B)
        logger.info("Arpeggio table filtered by ATOM_B...")

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
    agg_cols = ['CHAIN_A', 'RES_FULL_A', 'CHAIN_B', 'RES_FULL_B']
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
            # if contacts columns have been collapsed
            if 'Int_Types' in list(table):
                columns_to_agg['Int_Types'] = 'unique'
        elif agg_method_origin == 'maximum':
            # need the table sort by distance first: descending
            table = table.sort_values(["DIST", "VDW_DIST"], ascending=[False, False])
            table.reset_index(inplace=True)
            table = table.drop(['index'], axis=1)
            agg_generic = 'first'
            agg_method = 'max'
            columns_to_agg = {col: agg_generic if table[col].dtype == 'object' else agg_method
                              for col in table.columns if col not in agg_cols}
            # if contacts columns have been collapsed
            if 'Int_Types' in list(table):
                columns_to_agg['Int_Types'] = 'unique'
    table = table.groupby(by=agg_cols, as_index=False).agg(columns_to_agg)
    return table


def interaction_modes(data, int_mode='inter-chain'):
    """
    Gets the contacts filtered base on the entities that are interacting.
    Interaction modes possible: inter-chain and intra-chain.

    :param data: pandas DataFrame object
    :param int_mode: current values: 'inter-chain' and 'intra-chain'
    :return: returns a modified pandas DataFrame
    """

    table = data
    if int_mode == 'inter-chain':
        table = table.loc[table['CHAIN_A'] != table['CHAIN_B']]
    elif int_mode == 'intra-chain':
        table = table.loc[table['CHAIN_A'] == table['CHAIN_B']]
    else:
        raise ValueError('Interaction mode {} is not currently implemented...'
                         ''.format(int_mode))
    # FIXME optional?
    table.reset_index(inplace=True)
    table = table.drop(['index'], axis=1)
    return table


def collapsed_contacts(data, col_method='full'):
    """
    Collapses the various contact information columns into a single column.

    :param data: pandas DataFrame object
    :param col_method: collapse method: 'full' means all columns are used
        'minimal': only selected contact types are kept
    :return: returns a modified pandas DataFrame
    """

    table = data
    if col_method not in ['full', 'minimal']:
        raise ValueError('Method {} is not currently implemented...'
                         ''.format(col_method))
    # rename the columns
    table = table.rename(columns=arpeggio_col_renames)
    col_names = [k for k in list(arpeggio_col_renames.values()) if k in table]
    if col_method == 'minimal':
        col_min = (
            "Steric-Clash", "VDW-Bond", "Hydrogen-Bond", "Halogen-Bond", "Ionic-Bond",
            "Aromatic-Bond", "Hydrophobic-Bond", "Carbonyl-Bond", "Polar-Bond",
            "Amide-Amide", "Amide-Aromatic", "Aromatic-Aromatic", "Atom-Ring",
        )
        col_min = [k for k in col_min if k in table]
        excluded = [k for k in col_names if k not in col_min]
        col_names = list(col_min)
        table = table.drop(excluded, axis=1)
    # aggregate results
    melted = pd.melt(table.reset_index(), id_vars=['index'], value_vars=col_names, var_name='Int_Types')
    melted.query('value == 1', inplace=True)
    aggregated = melted.groupby(['index'])['Int_Types'].aggregate(lambda x: set(x))
    table = table.join(aggregated)
    # finally remove all the columns that are not needed anymore
    table = table.drop(col_names, axis=1)
    return table


def ignore_consecutive_residues(data, numb_res=3):
    """
    Drop atom-atom (or res-res) contacts that occur between n
    'numb_res' consecutive residues.

    If arpeggio has been generated from a 'pro_format' PDB generated file,
     which uses new_asym_id and new_seq_id, this method should work fine.

    :param data: pandas DataFrame object
    :param numb_res: (int) number of residues that are skipped
    :return: returns a modified pandas DataFrame
    """
    table = data
    # this only works if there are no ins_codes
    ins_codes_1 = [k for k in table.INSCODE_A.unique()]
    ins_codes_2 = [k for k in table.INSCODE_B.unique()]
    if (len(ins_codes_1) and len(ins_codes_2) and
            ins_codes_1[0] == '?' and ins_codes_2[0] == '?'):
        table = table.loc[((table["CHAIN_A"] != table["CHAIN_B"]) |
                           ((table["CHAIN_A"] == table["CHAIN_B"]) &
                            (abs(table["RES_A"].astype(int) -
                                 table["RES_B"].astype(int)) >= numb_res)))]
    else:
        message = ("Warning: Atoms in consecutive residues were not removed as there are "
                   "some with insertion codes. These are not handled at this time...")
        logger.debug(message)
    return table


class ARPEGGIOreader(object):
    def __init__(self, inputfile):
        """
        :param inputfile: Needs to point to a valid ARPEGGIO file.
        """
        self.inputfile = inputfile
        self.data = None
        self.excluded = ("ENTRY_A", "ENTRY_B", "ENTITIES")

        if not os.path.isfile(inputfile):
            raise IOError("{} not available or could not be read...".format(inputfile))

    def read(self, **kwargs):
        return self.contacts(**kwargs)

    def contacts(self, excluded=None, add_res_split=True,
                 residue_agg=False, agg_method='minimum',
                 int_filter=False, int_mode='inter-chain',
                 collapsed_cont=False, col_method='full',
                 ignore_consecutive=False, numb_res=3,
                 parse_special=False):

        if excluded is None:
            excluded = self.excluded
        self.data = parse_arpeggio_from_file(self.inputfile, excluded=excluded,
                                             add_res_split=add_res_split,
                                             parse_special=parse_special)
        if ignore_consecutive:
            self.data = ignore_consecutive_residues(self.data, numb_res=numb_res)

        if residue_agg:
            self.data = residues_aggregation(self.data, agg_method=agg_method)

        if collapsed_cont:
            self.data = collapsed_contacts(self.data, col_method=col_method)

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
            logger.info("No ARPEGGIO data parsed...")


class ARPEGGIOrunner(object):
    def __init__(self, inputfile, outputfile=None):
        """
        :param inputfile: Needs to point to a valid PDB or mmCIF file.
        :param outputfile: if not provided will use the same file name and
          <.contacts> extension
        """
        self.inputfile = inputfile
        self.inputfile_back = inputfile
        self.outputfile = outputfile
        self.data = None
        self.inputfile_h = None

        if not os.path.isfile(self.inputfile):
            raise IOError("{} not available or could not be read..."
                          "".format(self.inputfile))

        # inputfile needs to be in PDB or mmCIF format
        filename, extension = os.path.splitext(self.inputfile)
        if extension not in ['.pdb', '.ent', '.cif']:
            raise ValueError("{} is expected to be in mmCIF or PDB format..."
                             "".format(self.inputfile))

    def _generate_output(self):
        filename, extension = os.path.splitext(self.inputfile)
        self.outputfile = filename + ".contacts"

    def _generate_pdb(self, override=False, pro_format=False):
        filename, extension = os.path.splitext(self.inputfile)
        self.inputfile = filename + "_new.pdb"
        w = PDBXwriter(inputfile=None, outputfile=self.inputfile)
        r = PDBXreader(inputfile=self.inputfile_back)
        data = r.atoms(remove_altloc=True, reset_atom_id=True, add_new_pro_id=True,
                       remove_partial_res=True, format_type=None)
        w.run(data=data, format_type="pdb", category="auth",
              override=override, pro_format=pro_format)

    def _generate_pdb_with_hydrogens(self, hydro_method="hbplus", override=False):
        if hydro_method == "hbplus":
            w = HBPLUSrunner(inputfile=self.inputfile, outputfile=self.inputfile_h)
            w.run(hydro_pdb_out=True, override=override)
        elif hydro_method == "reduce":
            w = REDUCErunner(inputfile=self.inputfile, outputfile=self.inputfile_h)
            w.run(override=override)
        else:
            raise ValueError('Method {} is not currently implemented...'
                             ''.format(hydro_method))

    def _run(self, python_path, python_exe, arpeggio_bin, clean_output=True,
             hydro_method="arpeggio"):

        filename, extension = os.path.splitext(self.inputfile)
        input_arpeggio = filename + ".pdb"
        output_arpeggio = filename + ".contacts"  # atom-atom contact information
        output_bs_contacts = filename + ".bs_contacts"
        output_atomtypes = filename + ".atomtypes"  # atom types
        output_hydro = filename + "_hydrogenated.pdb"

        output_amam = filename + ".amam"  # res-res amide-amide
        output_amri = filename + ".amri"  # res-res amide-ring
        output_ari = filename + ".ari"  # atom-res atom-ring
        output_ri = filename + ".ri"  # res-res aromatic ring-aromatic ring
        output_rings = filename + ".rings"  # list of res rings
        output_residue_sifts = filename + ".residue_sifts"
        output_sift = filename + ".sift"
        output_siftmatch = filename + ".siftmatch"
        output_polarmatch = filename + ".polarmatch"
        output_specific_sift = filename + ".specific.sift"
        output_specific_siftmatch = filename + ".specific.siftmatch"
        output_specific_polarmatch = filename + ".specific.polarmatch"

        if hydro_method in ["hbplus", "reduce"]:
            input_arpeggio = filename + ".h.pdb"

        # run arpeggio
        cmd = 'PYTHONPATH={} {} {} -wh {}'.format(python_path, python_exe,
                                                  arpeggio_bin, input_arpeggio)
        os.system(cmd)
        if not os.path.isfile(output_arpeggio):
            raise IOError("ARPEGGIO output not generated for {}".format(input_arpeggio))

        # mv the automatically generated file -> to the provided outputfile
        if output_arpeggio != self.outputfile:
            shutil.copyfile(output_arpeggio, self.outputfile)
            nfilename, nextension = os.path.splitext(self.outputfile)
            new_output_amam = nfilename + ".amam"
            new_output_amri = nfilename + ".amri"
            new_output_ari = nfilename + ".ari"
            new_output_ri = nfilename + ".ri"
            shutil.copyfile(output_amam, new_output_amam)
            shutil.copyfile(output_amri, new_output_amri)
            shutil.copyfile(output_ari, new_output_ari)
            shutil.copyfile(output_ri, new_output_ri)

        if hydro_method == "arpeggio":
            shutil.copyfile(output_hydro, self.inputfile_h)

        if clean_output:
            # remove output files
            if output_arpeggio != self.outputfile:
                lazy_file_remover(output_arpeggio)
                lazy_file_remover(output_amam)
                lazy_file_remover(output_amri)
                lazy_file_remover(output_ari)
                lazy_file_remover(output_ri)

            lazy_file_remover(output_bs_contacts)
            lazy_file_remover(output_atomtypes)
            lazy_file_remover(output_hydro)
            lazy_file_remover(output_rings)
            lazy_file_remover(output_residue_sifts)
            lazy_file_remover(output_sift)
            lazy_file_remover(output_siftmatch)
            lazy_file_remover(output_polarmatch)
            lazy_file_remover(output_specific_sift)
            lazy_file_remover(output_specific_siftmatch)
            lazy_file_remover(output_specific_polarmatch)

    def write(self, **kwargs):
        return self.run(**kwargs)

    def run(self, hydro_method="arpeggio", override=False,
            clean_output=True, save_new_input=False, pro_format=False):

        # generate outputfile if missing
        if not self.outputfile:
            self._generate_output()

        if not os.path.exists(self.outputfile) or override:
            if os.path.isfile(config.python_exe) and os.path.exists(config.arpeggio_bin):
                arpeggio_bin = config.arpeggio_bin
                python_exe = config.python_exe
                python_path = config.python_path
            else:
                raise IOError('ARPEGGIO executable is not available...')

            # inputfile needs to be in PDB format
            filename, extension = os.path.splitext(self.inputfile)
            self._generate_pdb(override=override, pro_format=pro_format)

            # get PDB with explicit hydrogen atoms
            self.inputfile_h = filename + ".h.pdb"
            if hydro_method in ["hbplus", "reduce"]:
                self._generate_pdb_with_hydrogens(hydro_method=hydro_method,
                                                  override=override)

            # run arpeggio and generate output - also clean unnecessary output
            self._run(python_path, python_exe, arpeggio_bin, clean_output=clean_output,
                      hydro_method=hydro_method)

            # clean the new PDB input file generated
            if not save_new_input:
                if self.inputfile != self.inputfile_back:
                    lazy_file_remover(self.inputfile)
                    lazy_file_remover(self.inputfile_h)

        else:
            logger.info("ARPEGGIO for %s already available...", self.outputfile)
        return


if __name__ == '__main__':
    pass
