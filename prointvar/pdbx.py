#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""

This defines the methods that work with PDB/mmCIF files.

Fábio Madeira, 2017+

"""

import os
import json
import shlex
import logging
import pandas as pd
from io import StringIO
from scipy.spatial import cKDTree
from string import ascii_uppercase
from collections import OrderedDict

from prointvar.utils import row_selector
from prointvar.utils import string_split
from prointvar.utils import get_new_pro_ids
from prointvar.library import mmcif_types
from prointvar.library import aa_default_atoms

logger = logging.getLogger("prointvar")

_PDB_FORMAT = "%s%5i %-4s%c%3s %c%4s%c   %8.3f%8.3f%8.3f%s%6.2f      %4s%2s%2s\n"


def parse_mmcif_atoms_from_file(inputfile, excluded=(), add_res_full=True,
                                add_contacts=False, dist=5, first_model=True,
                                add_atom_altloc=False, remove_altloc=False,
                                remove_hydrogens=True, reset_atom_id=True,
                                add_new_pro_id=False, remove_partial_res=True):
    """
    Parse mmCIF ATOM and HETATM lines.

    :param inputfile: path to the mmCIF file
    :param excluded: option to exclude mmCIF columns
    :param add_res_full: option to extend the table with 'res_full'
        i.e. res_number + insertion_code (e.g. '12A')
    :param add_contacts: boolean
    :param dist: distance threshold in Angstrom
    :param first_model: boolean
    :param add_atom_altloc: boolean (new string join)
    :param remove_altloc: boolean
    :param remove_hydrogens: boolean
    :param reset_atom_id: boolean
    :param add_new_pro_id: (boolean) used for chain_id mapping
    :param remove_partial_res: (boolean) removes amino acids with missing atoms
    :return: returns a pandas DataFrame
    """

    logger.info("Parsing mmCIF atoms from lines...")

    # example lines with some problems
    """
    _atom_site.pdbx_PDB_model_num
    _atom_site.pdbe_label_seq_id
    _atom_site.orig_label_asym_id
    _atom_site.orig_auth_asym_id
    ATOM 1 N N . VAL A 1 1 ? -7.069 21.943 18.770 1.0 56.51 ? ? ? ? ? ? 118 VAL A N 1 1 A A
    ATOM 2 C CA . VAL A 1 1 ? -7.077 21.688 20.244 1.0 59.09 ? ? ? ? ? ? 118 VAL A CA 1 1 A A
    ATOM 3 C C . VAL A 1 1 ? -5.756 21.077 20.700 1.0 44.63 ? ? ? ? ? ? 118 VAL A C 1 1 A A
    ATOM 4 O O . VAL A 1 1 ? -5.346 20.029 20.204 1.0 59.84 ? ? ? ? ? ? 118 VAL A O 1 1 A A
    """

    if not os.path.isfile(inputfile):
        raise IOError("{} not available or could not be read...".format(inputfile))

    # parsing atom lines
    header = []
    lines = []
    with open(inputfile) as inlines:
        for line in inlines:
            if line.startswith("_atom_site."):
                header.append(line.split('.')[1].rstrip())
            elif line.startswith("ATOM"):
                lines.append(line)
            elif line.startswith("HETATM"):
                lines.append(line)
    lines = "".join(lines)

    all_str = {key: str for key in header}
    table = pd.read_table(StringIO(lines), delim_whitespace=True, low_memory=False,
                          names=header, compression=None, converters=all_str,
                          keep_default_na=False)

    # excluding columns
    if excluded is not None:
        assert type(excluded) is tuple
        try:
            table = table.drop(list(excluded), axis=1)
        except ValueError:
            # most likely theses are not in there
            pass

    # if only first model (>1 in NMR structures)
    if first_model:
        table = row_selector(table, key='pdbx_PDB_model_num', value=None,
                             method='first')

    # table modular extensions
    if add_contacts:
        table = add_mmcif_contacts(table, dist=dist)
        logger.info("PDBx added contact indexes...")

    if add_res_full:
        table = add_mmcif_res_full(table)
        logger.info("PDBx added full res (res + ins_code)...")

    if add_atom_altloc:
        table = add_mmcif_atom_altloc(table)
        logger.info("PDBx added full atom (atom + altloc)...")

    if add_new_pro_id:
        table = add_mmcif_new_pro_ids(table)
        logger.info("PDBx added 'new' chain and res ids...")

    if remove_altloc:
        table = remove_multiple_altlocs(table)
        reset_atom_id = True
        logger.info("PDBx removed altlocs...")

    if remove_hydrogens:
        table = row_selector(table, key='type_symbol', value='H', method='diffs')
        logger.info("PDBx removed existing hydrogens...")

    if remove_partial_res:
        table = remove_partial_residues(table)
        logger.info("PDBx removed incomplete residues...")

    if reset_atom_id:
        table.reset_index(inplace=True)
        table = table.drop(['index'], axis=1)
        table['id'] = table.index + 1
        logger.info("PDBx reset atom numbers...")

    # enforce some specific column types
    for col in table:
        if col in mmcif_types:
            try:
                table[col] = table[col].astype(mmcif_types[col])
            except ValueError:
                # there are some NaNs in there
                pass

    if table.empty:
        raise ValueError('{} resulted in an empty DataFrame...'.format(inputfile))

    return table


def parse_pdb_atoms_from_file(inputfile, excluded=(), add_contacts=False,
                              dist=5, first_model=True, add_atom_altloc=False,
                              remove_altloc=False, remove_hydrogens=True,
                              reset_atom_id=True, add_new_pro_id=False,
                              remove_partial_res=True):
    """
    Parse PDB ATOM and HETATM lines.

    :param inputfile: path to the PDB file
    :param excluded: option to exclude mmCIF columns
    :param add_contacts: boolean
    :param dist: distance threshold in Angstrom
    :param first_model: boolean 
    :param add_atom_altloc: boolean (new string join)
    :param remove_altloc: boolean
    :param remove_hydrogens: boolean
    :param reset_atom_id: boolean
    :param add_new_pro_id: (boolean) used for chain_id mapping
    :param remove_partial_res: (boolean) removes amino acids with missing atoms
    :return: returns a pandas DataFrame
    """

    logger.info("Parsing PDB atoms from lines...")

    # example lines
    """
    MODEL        1
    ATOM      0  N   SER A  -1     104.083  78.916  -1.349  1.00 61.47           N
    ATOM      1  N   VAL A 118      -7.069  21.943  18.770  1.00 56.51           N  
    ATOM      2  CA  VAL A 118      -7.077  21.688  20.244  1.00 59.09           C  
    ATOM      3  C   VAL A 118      -5.756  21.077  20.700  1.00 44.63           C  
    ATOM      4  O   VAL A 118      -5.346  20.029  20.204  1.00 59.84           O
    """

    if not os.path.isfile(inputfile):
        raise IOError("{} not available or could not be read...".format(inputfile))

    # parsing atom lines, converting it to mmcif-style headers
    lines = []
    modelnumb = '1'
    with open(inputfile) as inlines:
        for line in inlines:
            line = line.rstrip()
            line = line[0:78]
            if line.startswith("MODEL"):
                modelnumb = line.split()[1]
            elif line.startswith("ATOM"):
                lines.append(line + "%s" % modelnumb)
            elif line.startswith("HETATM"):
                lines.append(line + "%s" % modelnumb)
    lines = "\n".join(lines)

    header = ('group_PDB', 'id', 'label_atom_id', 'label_alt_id', 'label_comp_id',
              'label_asym_id', 'label_seq_id_full', 'label_seq_id',
              'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z',
              'occupancy', 'B_iso_or_equiv', 'type_symbol', 'auth_atom_id', 'auth_comp_id',
              'auth_asym_id', 'auth_seq_id_full', 'auth_seq_id', 'pdbx_PDB_model_num')

    # https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    widths = ((0, 6), (6, 11), (12, 16), (16, 17), (17, 20), (21, 22), (22, 27), (22, 26), (26, 27),
              (30, 38), (38, 46), (46, 54), (54, 60), (60, 66), (76, 78), # (72, 76), ('seg_id')
              (12, 16), (17, 20), (21, 22), (22, 27), (22, 26), (78, 79))

    all_str = {key: str for key in header}
    table = pd.read_fwf(StringIO(lines), names=header, colspecs=widths,
                        compression=None, converters=all_str, keep_default_na=False)

    # excluding columns
    if excluded is not None:
        assert type(excluded) is tuple
        try:
            table = table.drop(list(excluded), axis=1)
        except ValueError:
            # most likely theses are not in there
            pass

    # if only first model (>1 in NMR structures)
    if first_model:
        table = row_selector(table, key='pdbx_PDB_model_num', value=None,
                             method='first')

    # fixes the 'pdbx_PDB_ins_code'
    table = fix_pdb_ins_code(table)
    # fixes the 'label_alt_id
    table = fix_label_alt_id(table)
    # fixes 'type_symbol' if missing
    table = fix_type_symbol(table)

    # table modular extensions
    if add_contacts:
        table = add_mmcif_contacts(table, dist=dist)
        logger.info("PDBx added contact indexes...")

    if add_atom_altloc:
        table = add_mmcif_atom_altloc(table)
        logger.info("PDBx added full atom (atom + alt_loc)...")

    if add_new_pro_id:
        table = add_mmcif_new_pro_ids(table)
        logger.info("PDBx added 'new' chain and res ids...")

    if remove_altloc:
        table = remove_multiple_altlocs(table)
        reset_atom_id = True
        logger.info("PDBx removed altlocs...")

    if remove_hydrogens:
        table = row_selector(table, key='type_symbol', value='H', method='diffs')
        logger.info("PDBx removed existing hydrogens...")

    if remove_partial_res:
        table = remove_partial_residues(table)
        logger.info("PDBx removed incomplete residues...")

    if reset_atom_id:
        table.reset_index(inplace=True)
        table = table.drop(['index'], axis=1)
        table['id'] = table.index + 1
        logger.info("PDBx reset atom numbers...")

    # enforce some specific column types
    for col in table:
        if col in mmcif_types:
            try:
                table[col] = table[col].astype(mmcif_types[col])
            except ValueError:
                # there are some NaNs in there
                pass

    if table.empty:
        raise ValueError('{} resulted in an empty DataFrame...'.format(inputfile))

    return table


def _tokenize(handle):
    for line in handle:
        if line.startswith("#"):
            continue
        elif line.startswith(";"):
            token = line[1:].strip()
            for line in handle:
                line = line.strip()
                if line == ';':
                    break
                token += line
            yield token
        else:
            tokens = shlex.split(line)
            for token in tokens:
                yield token


def parse_mmcif_categories_from_file(inputfile, excluded=(), category=None):
    """
    Generic method that gets all categories and fields (except for the .atom_site*)

    :param inputfile: path to the mmCIF file
    :param excluded: option to exclude mmCIF categories
    :param category: data category to be parsed (otherwise all)
    :return: returns a nested dictionary
    """

    logger.info("Parsing mmCIF categories from lines...")

    if not os.path.isfile(inputfile):
        raise IOError("{} not available or could not be read...".format(inputfile))

    full_dict = OrderedDict()
    with open(inputfile) as handle:
        loop_flag = False
        key = None
        full_tokens = _tokenize(handle)
        full_token = next(full_tokens)
        full_dict[full_token[0:5]] = full_token[5:]
        i = 0
        n = 0
        for full_token in full_tokens:
            if full_token == "loop_":
                loop_flag = True
                keys = []
                i = 0
                n = 0
                continue
            elif loop_flag:
                if full_token.startswith("_"):
                    if i > 0:
                        loop_flag = False
                    else:
                        full_dict[full_token] = []
                        keys.append(full_token)
                        n += 1
                        continue
                else:
                    full_dict[keys[i % n]].append(full_token)
                    i += 1
                    continue
            if key is None:
                key = full_token
            else:
                full_dict[key] = [full_token]
                key = None

    mmcid_object = OrderedDict()
    for full_key in full_dict:
        if full_key != "data_":
            main_key = full_key.split(".")[0]
            lead_key = full_key.split(".")[1]
            if (category is not None and main_key == category) or category is None:
                if main_key not in excluded:
                    if main_key not in mmcid_object:
                        mmcid_object[main_key] = {}
                    mmcid_object[main_key][lead_key] = full_dict[full_key]

    return mmcid_object


def get_mmcif_selected_from_table(data, chain=None, res=None, res_full=None, comp=None,
                                  atom=None, lines=None, category='label'):
    """
    Utility that filters a pandas DataFrame by the input tuples. 
    There are various filtering options available:
        chains, residues, atoms, lines ('ATOM' or 'HETATM')

    :param data: pandas DataFrame object
    :param chain: (tuple) chain IDs or None
    :param res: (tuple) res IDs or None
    :param res_full: (tuple) full res IDs ('*_seq_id' + 'pdbx_PDB_ins_code') or None
    :param comp: (tuple) comp IDs or None
    :param atom: (tuple) atom IDs or None
    :param lines: (tuple) 'ATOM' or 'HETATM' or None (both)
    :param category: data category to be used as precedence in _atom_site.*_*
        asym_id, seq_id and atom_id
    :return: returns a modified pandas DataFrame
    """

    # excluding rows
    table = data
    if chain is not None:
        table = row_selector(table, '{}_asym_id'.format(category), chain, method="isin")
        logger.info("PDBx table filtered by %s_asym_id...", category)

    if res is not None:
        table = row_selector(table, '{}_seq_id'.format(category), res, method="isin")
        logger.info("PDBx table filtered by %s_seq_id...", category)

    if res_full is not None:
        table = row_selector(table, '{}_seq_id_full'.format(category), res_full, method="isin")
        logger.info("PDBx table filtered by %s_seq_id_full...", category)

    if comp is not None:
        table = row_selector(table, '{}_comp_id'.format(category), comp, method="isin")
        logger.info("PDBx table filtered by %s_comp_id...", category)

    if atom is not None:
        table = row_selector(table, '{}_atom_id'.format(category), atom, method="isin")
        logger.info("PDBx table filtered by %s_atom_id...", category)

    if lines is not None:
        table = row_selector(table, 'group_PDB', lines, method="isin")
        logger.info("PDBx table filtered by group_PDB...")

    return table


def residues_aggregation(data, agg_method='centroid', category='label'):
    """
    Gets the residues' atoms and their centroids (mean).

    :param data: pandas DataFrame object
    :param agg_method: current values: 'centroid', 'first', 'mean' and 'unique'
    :param category: data category to be used as precedence in _atom_site.*_*
        asym_id, seq_id and atom_id
    :return: returns a modified pandas DataFrame
    """
    table = data
    agg_generic = agg_method
    agg_cols = ['pdbx_PDB_model_num', '{}_asym_id'.format(category),
                '{}_seq_id'.format(category)]
    if agg_method not in ['centroid', 'first', 'unique', 'mean']:
        raise ValueError('Method {} is not currently implemented...'
                         ''.format(agg_method))
    if agg_method == 'centroid' or agg_method == 'mean':
        agg_generic = 'first'
        agg_method = 'mean'
    columns_to_agg = {col: agg_generic if table[col].dtype == 'object' else agg_method
                      for col in table.columns if col not in agg_cols}
    columns_to_agg['id'] = 'first'
    table = table.groupby(by=agg_cols, as_index=False).agg(columns_to_agg)
    table.sort_values(by='id').reset_index()
    return table


def write_mmcif_from_table(outputfile, data, override=False):
    """
    Generic method that writes 'atom' lines in mmCIF format.
    
    :param outputfile: path to the mmCIF file
    :param data: pandas DataFrame object
    :param override: boolean
    :return: (side effects) writes to file
    """

    table = data
    atom_lines = ['data_mmCIF_generated_by_ProIntVar', 'loop_']
    atom_lines += ["_atom_site.{}".format(v) for v in list(table)]
    for i in table.index:
        line = ' '.join([str(v) for v in list(table.loc[i, :])])
        atom_lines.append(line)

    # write the final output
    if not os.path.exists(outputfile) or override:
        with open(outputfile, 'w') as outlines:
            outlines.write("\n".join(atom_lines))
    else:
        logger.info("mmCIF for %s already available...", outputfile)
    return


def write_pdb_from_table(outputfile, data, override=False, pro_format=False):
    """
    Generic method that writes 'atom' lines in PDB format.
    
    :param outputfile: path to the PDB file
    :param data: pandas DataFrame object
    :param override: boolean
    :param pro_format: ProIntVar internal format where asym_id and seq_id
        are outputted from new_asyn_id and new_seq_id
    :return: (side effects) writes to file
    """

    table = data
    atom_lines = ['REMARK 100 PDB generated by ProIntVar\n']
    atom_number = 0
    for i in table.index:
        atom_number += 1
        atom_lines.append(get_atom_line(data=table, index=i,
                                        atom_number=atom_number,
                                        pro_format=pro_format))

    # write the final output
    if not os.path.exists(outputfile) or override:
        with open(outputfile, 'w') as outlines:
            outlines.write("".join(atom_lines))
    else:
        logger.info("PDB for %s already available...", outputfile)
    return


def get_contact_indexes_from_table(data, dist=5):
    """
    Gets the DataFrame indexes of the nearby atoms under a provided radius (dist).

    :param data: pandas DataFrame object
    :param dist: distance threshold in Angstrom
    :return: list of pandas DataFrame indexes
    """

    table = data
    tree = cKDTree(table[['Cartn_x', 'Cartn_y', 'Cartn_z']])
    query_point = table.loc[:, ['Cartn_x', 'Cartn_y', 'Cartn_z']]
    indexes = tree.query_ball_point(query_point, r=dist)

    return indexes


def add_mmcif_contacts(data, dist=5):
    """
    Utility that adds a new column to the table.
    Adds a new column with the contacting indexes (str)

    :param data: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    table = data
    indexes = get_contact_indexes_from_table(table, dist=dist)
    assert len(indexes) == len(table)
    nindexes = [','.join([str(z) for z in x]) for x in indexes]
    table['contact_indexes'] = nindexes

    return table


def add_mmcif_res_full(data):
    """
    Utility that adds a new column to the table.
    Adds a new column with the 'full res' (i.e. seq_id + ins_code).

    :param data: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    table = data
    # adds both 'label' and 'auth' entries
    if 'label_seq_id' in table:
        seqs_full = []
        for ix in table.index:
            seq = "{}{}".format(table.loc[ix, 'label_seq_id'],
                                table.loc[ix, 'pdbx_PDB_ins_code']).replace('?', '')
            seqs_full.append(seq)
        assert len(seqs_full) == len(table)
        table['label_seq_id_full'] = seqs_full
    if 'auth_seq_id' in table:
        seqs_full = []
        for ix in table.index:
            seq = "{}{}".format(table.loc[ix, 'auth_seq_id'],
                                table.loc[ix, 'pdbx_PDB_ins_code']).replace('?', '')
            seqs_full.append(seq)
        assert len(seqs_full) == len(table)
        table['auth_seq_id_full'] = seqs_full

    return table


def get_mmcif_res_split(data):
    """
    Utility that adds new columns to the table.
    Adds new columns from the 'full res' (i.e. seq_id + ins_code).
    
    adds: 'label_seq_id', 'auth_seq_id', and 'pdbx_PDB_ins_code'
    
    :param data: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    table = data
    seq_ids = []
    ins_codes = []
    for i in table.index:
        values = string_split(table.loc[i, "label_seq_id_full"])
        if len(values) == 1:
            values.append('?')
        elif len(values) == 2 and values[0] == '-':
            values = ["".join(values), '?']
        elif len(values) == 3 and values[0] == '-':
            values = ["".join(values[0] + values[1]), values[2]]
        seq_ids.append(values[0])
        ins_codes.append(values[1])
    table["label_seq_id"] = seq_ids
    table["auth_seq_id"] = seq_ids
    table["pdbx_PDB_ins_code"] = ins_codes

    return table


def fix_pdb_ins_code(data):
    """
    Utility that fixes the 'pdbx_PDB_ins_code' column to match is expected
    in the mmcif format.

    :param data: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    table = data
    ins_codes = []
    for i in table.index:
        value = table.loc[i, "pdbx_PDB_ins_code"]
        if value == '' or value == ' ' or value == '?':
            value = '?'
        ins_codes.append(value)
    table["pdbx_PDB_ins_code"] = ins_codes
    table['pdbx_PDB_ins_code'] = table['pdbx_PDB_ins_code'].fillna('?').astype(str)
    return table


def fix_label_alt_id(data):
    """
    Utility that fixes the 'label_alt_id' column to match what is
    expected in the mmCIF format.

    :param data: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    table = data
    alt_locs = []
    for i in table.index:
        value = table.loc[i, "label_alt_id"]
        if value == '' or value == ' ' or value == '?':
            value = '.'
        alt_locs.append(value)
    table["label_alt_id"] = alt_locs
    table['label_alt_id'] = table['label_alt_id'].fillna('.').astype(str)
    return table


def fix_type_symbol(data):
    """
    Utility that fixes the 'type_symbol' column to match what is
    expected in the mmCIF format - when missing in the Structure.

    :param data: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """
    table = data

    # print(table.head())
    def get_type_symbol(data, key, key_fix):
        # this maybe a bit crude way of assigning this value
        if data[key] != " " and data[key] != "" and len(data[key]):
            return data[key]
        else:
            return ''.join([x for x in data[key_fix] if x in ascii_uppercase])[0]
    table.is_copy = False
    table['type_symbol'] = table.apply(get_type_symbol, axis=1,
                                       args=('type_symbol', 'label_atom_id'))
    return table


def add_mmcif_atom_altloc(data):
    """
    Utility that adds new columns to the table.
    adds: 'label_atom_altloc_id', and 'auth_atom_altloc_id' which is a
        string join between '*_atom_id' + 'label_alt_id'
    
    :param data: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    def join_atom_altloc(data, category='label'):
        atom = data['{}_atom_id'.format(category)]
        altloc = data['label_alt_id']
        if altloc == "." or altloc == '' or altloc == ' ':
            return atom
        else:
            return atom + '.' + altloc

    data.is_copy = False
    data['label_atom_altloc_id'] = data.apply(join_atom_altloc,
                                              axis=1, args=('label', ))
    data['auth_atom_altloc_id'] = data.apply(join_atom_altloc,
                                             axis=1, args=('auth', ))
    return data


def add_mmcif_new_pro_ids(data, category='label'):
    """
    Adds a new column to the table with a new Entity/Chain ID to be used
    for mapping chains.

    :param data: pandas DataFrame object
    :param category: data category to be used as precedence in _atom_site.*_*
        asym_id, seq_id and atom_id
    :return: returns a modified pandas DataFrame
    """

    table = data
    try:
        new_pro_id = get_new_pro_ids()
        ochain_ids = table.loc[:, "{}_asym_id".format(category)].tolist()
        oseq_ids = table.loc[:, "{}_seq_id".format(category)].tolist()
        flat = [k + l for k, l in zip(ochain_ids, oseq_ids)]
        nchain_ids = []
        nseq_ids = []
        nchain = 'A'
        nseq = '1'
        prev_pro = "     "
        for pro in flat:
            if prev_pro != pro:
                prev_pro = pro
                nchain, nseq = next(new_pro_id)
            nchain_ids.append(nchain)
            nseq_ids.append(nseq)
        assert len(table.index) == len(nchain_ids)
        assert len(table.index) == len(nseq_ids)
        table['new_asym_id'] = nchain_ids
        table['new_seq_id'] = nseq_ids
    except StopIteration:
        message = ("This structure contains >9999 (seq_ids) * 62 (asym_ids), which causes "
                   "problems to work with DSSP and arpeggio. Please review the structure...")
        raise StopIteration(message)
    return table


def remove_multiple_altlocs(data):
    """
    Removes alternative locations (i.e. 'rows') leaving only the first.
    Needs to find rows with alt_id != '.' and then find following rows until
    '.' appears again (expects atoms to be consequent).

    :param data: pandas DataFrame object
    :return: returns a modified pandas DataFrame
    """

    table = data
    drop_ixs = []
    for ix in table.index:
        altloc = table.loc[ix, 'label_alt_id']
        if altloc != '.':
            # table.loc[ix, 'label_alt_id'] = '.'
            table.set_value(ix, 'label_alt_id', '.')
            atomid = table.loc[ix, 'label_atom_id']
            try:
                for nx in range(1, 100, 1):
                    altnx = table.loc[ix + nx, 'label_alt_id']
                    atomnx = table.loc[ix + nx, 'label_atom_id']
                    if altnx != '.' and atomnx == atomid:
                        # store indexes of the rows to be dropped
                        drop_ixs.append(ix + nx)
                    else:
                        break
            except KeyError:
                break
    return table.drop(table.index[drop_ixs])


def remove_partial_residues(data, category='label'):
    """
    Removes residues that contain missing atoms. Needs to check
    which atoms are available for each residue. Also removes residues with
    same '*_seq_id' as the previous residue.

    :param data: pandas DataFrame object
    :param category: data category to be used as precedence in _atom_site.*_*
        asym_id, seq_id and atom_id
    :return: returns a modified pandas DataFrame
    """
    table = data
    drop_ixs = []
    curr_ixs = []
    curr_atoms = []
    prev_res = ''
    prev_seq = ''
    next_res_for_rm = False
    table.reset_index(inplace=True)
    table = table.drop(['index'], axis=1)
    for ix in table.index:
        group = table.loc[ix, 'group_PDB']
        if group == 'ATOM':
            curr_res = table.loc[ix, '{}_comp_id'.format(category)]
            curr_seq = table.loc[ix, '{}_seq_id'.format(category)]
            if curr_res in aa_default_atoms:
                curr_atom = table.loc[ix, '{}_atom_id'.format(category)]
                if prev_res == curr_res and prev_seq == curr_seq:
                    curr_ixs.append(ix)
                    curr_atoms.append(curr_atom)
                else:
                    if curr_ixs:
                        # check available atoms
                        default_atoms = aa_default_atoms[prev_res]
                        intersection = list(set(default_atoms) - set(curr_atoms))
                        if intersection != [] or next_res_for_rm:
                            # missing atoms: means that there are atoms in the 'default' list
                            # that are not observed in the structure
                            drop_ixs += curr_ixs
                            next_res_for_rm = False
                        elif prev_seq == curr_seq:
                            # duplicated *_seq_id: means that the next residue has got the same
                            # *_seq_id as the previous res (could come from removing altlocs)
                            next_res_for_rm = True
                    # resetting variables
                    prev_res = curr_res
                    prev_seq = curr_seq
                    curr_ixs = [ix]
                    curr_atoms = [curr_atom]

    return table.drop(table.index[drop_ixs])


def get_atom_line(data, index, atom_number, pro_format=False,
                  coords=None, new_chain=None, category='auth'):
    """
    Returns an ATOM PDB-formatted string.
    (Based on code from the PDB module in Biopython.)

    :param data: pandas DataFrame object
    :param index: atom index
    :param atom_number: incremental number
    :param pro_format: ProIntVar internal format where asym_id and seq_id
        are outputted from new_asyn_id and new_seq_id
    :param coords: list of transformed coordinates
    :param new_chain: if true defaults to chain "A"
    :param category: data category to be used as precedence in _atom_site.*_*
        asym_id, seq_id and atom_id
    :return: returns a PDB-formatted ATOM/HETATM line
    """

    table = data
    ix = index

    record_type = table.loc[ix, 'group_PDB']
    if record_type == "ATOM":
        record_type = "ATOM  "

    """
    ATOM     16  CB  ASN A   2      22.780  31.612   8.039  1.00 97.98           C
    ATOM     17  CG  ASN A   2      23.735  31.870   9.167  1.00100.56           C
    ATOM     18  OD1 ASN A   2      23.345  32.366  10.218  1.00 97.84           O
    ATOM     19  ND2 ASN A   2      25.014  31.606   8.922  1.00106.62           N
    ATOM     20  H   ASN A   2      24.256  34.106   6.858  1.00  0.00           H
    ATOM     21 HD21 ASN A   2      25.654  31.751   9.644  1.00  0.00           H
    ATOM     22 HD22 ASN A   2      25.276  31.270   8.035  1.00  0.00           H
    """

    name = table.loc[ix, '{}_atom_id'.format(category)]
    if len(name) == 1:
        name = " {}  ".format(name.strip())
    elif len(name) == 2:
        name = " {} ".format(name.strip())
    elif len(name) == 3:
        name = " {}".format(name.strip())
    elif len(name) == 4:
        name = name.strip()

    altloc = table.loc[ix, 'label_alt_id']
    if altloc == ".":
        altloc = " "

    resname = table.loc[ix, '{}_comp_id'.format(category)]

    # FIXME is this needed?
    # overriding the original chain name
    if new_chain is not None:
        chain_id = new_chain
    else:
        if pro_format:
            chain_id = table.loc[ix, 'new_asym_id']
        else:
            chain_id = table.loc[ix, '{}_asym_id'.format(category)]
            chain_id = chain_id[0]

    if pro_format:
        resseq = str(table.loc[ix, 'new_seq_id'])
    else:
        resseq = str(table.loc[ix, '{}_seq_id'.format(category)])

    icode = table.loc[ix, 'pdbx_PDB_ins_code']
    if icode == "?":
        icode = " "

    # overriding the original coordinates
    if coords is not None:
        x = float(coords[0])
        y = float(coords[1])
        z = float(coords[2])
    else:
        x = float(table.loc[ix, 'Cartn_x'])
        y = float(table.loc[ix, 'Cartn_y'])
        z = float(table.loc[ix, 'Cartn_z'])

    occupancy_str = "%6.2f" % float(table.loc[ix, 'occupancy'])

    bfactor = float(table.loc[ix, 'B_iso_or_equiv'])

    segid = ""

    element = table.loc[ix, 'type_symbol']
    element = element.strip().upper()

    charge = "  "

    values = (record_type, atom_number, name, altloc, resname, chain_id,
              resseq, icode, x, y, z, occupancy_str, bfactor, segid,
              element, charge)

    return _PDB_FORMAT % values


class PDBXreader(object):
    def __init__(self, inputfile):
        """
        :param inputfile: Needs to point to a valid mmCIF file.
        """
        self.inputfile = inputfile
        self.data = None
        # defaults
        self.excluded = ("Cartn_x_esd", "Cartn_y_esd", "Cartn_z_esd",
                         "occupancy_esd", "B_iso_or_equiv_esd", "pdbx_formal_charge")

        if not os.path.isfile(inputfile):
            raise IOError("{} not available or could not be read...".format(inputfile))

    def read(self, **kwargs):
        return self.atoms(**kwargs)

    def atoms(self, excluded=None, add_res_full=True, add_contacts=False, dist=5,
              first_model=True, add_atom_altloc=False, remove_altloc=False,
              remove_hydrogens=True, reset_atom_id=True, format_type="mmcif",
              residue_agg=False, agg_method='centroid', category='label',
              add_new_pro_id=False, remove_partial_res=False):
        if excluded is None:
            excluded = self.excluded

        if not format_type:
            # try to guess the correct format
            if self.inputfile.endswith('.pdb') or self.inputfile.endswith('.ent'):
                format_type = "pdb"
            elif self.inputfile.endswith('.cif'):
                format_type = "mmcif"
            else:
                message = ("Could not guess the format of the input file... "
                           "Please define it by passing 'format_type'='<name>'")
                raise ValueError(message)

        if format_type == "mmcif":
            self.data = parse_mmcif_atoms_from_file(self.inputfile, excluded=excluded,
                                                    add_res_full=add_res_full,
                                                    add_contacts=add_contacts, dist=dist,
                                                    first_model=first_model,
                                                    add_atom_altloc=add_atom_altloc,
                                                    remove_altloc=remove_altloc,
                                                    remove_hydrogens=remove_hydrogens,
                                                    reset_atom_id=reset_atom_id,
                                                    add_new_pro_id=add_new_pro_id,
                                                    remove_partial_res=remove_partial_res)

        elif format_type == "pdb":
            self.data = parse_pdb_atoms_from_file(self.inputfile, excluded=excluded,
                                                  add_contacts=add_contacts, dist=dist,
                                                  first_model=first_model,
                                                  add_atom_altloc=add_atom_altloc,
                                                  remove_altloc=remove_altloc,
                                                  remove_hydrogens=remove_hydrogens,
                                                  reset_atom_id=reset_atom_id,
                                                  add_new_pro_id=add_new_pro_id,
                                                  remove_partial_res=remove_partial_res)
        else:
            message = 'The provided format {} is not implemented...'.format(format_type)
            raise ValueError(message)

        if residue_agg:
            self.data = residues_aggregation(self.data, agg_method=agg_method,
                                             category=category)

        return self.data

    def categories(self, excluded=None, category=None):
        if excluded is None:
            excluded = ('_atom_site', '_entity_poly_seq', '_pdbx_poly_seq_scheme')
        self.data = parse_mmcif_categories_from_file(self.inputfile, excluded=excluded,
                                                     category=category)
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
            logger.info("No mmCIF data parsed...")


class PDBXwriter(object):
    def __init__(self, inputfile=None, outputfile=None):
        """
        :param inputfile: Needs to point to a valid mmCIF file.
        :param outputfile: if not provided will use the same file name and
            <_filtered.cif> extension
        """
        self.inputfile = inputfile
        self.outputfile = outputfile
        self.data = None

        if inputfile is not None and not os.path.isfile(inputfile):
            raise IOError("{} not available or could not be read...".format(inputfile))

    def _generate_output(self, format_type="cif"):
        if not self.outputfile:
            filename, extension = os.path.splitext(self.inputfile)
            self.outputfile = filename + ".cif"
            if format_type == "pdb":
                self.outputfile = filename + ".cif"

    def run(self, data=None, chain=None, res=None, atom=None, lines=None, category='label',
            override=False, format_type="mmcif", pro_format=False):

        # generate outputfile if missing
        if self.outputfile is None:
            self._generate_output(format_type=format_type)

        # writes a new mmCIF with selected chains, residues, atoms, or lines
        if data is None:
            r = PDBXreader(inputfile=self.inputfile)
            # guess the input format as 'format_type' refers to the output format
            try:
                data = r.atoms(excluded=(), format_type="mmcif")
            except ValueError:
                data = r.atoms(excluded=(), format_type="pdb")
        self.data = get_mmcif_selected_from_table(data, chain=chain, res=res, atom=atom,
                                                  lines=lines, category=category)
        if format_type == "mmcif":
            write_mmcif_from_table(outputfile=self.outputfile, data=self.data,
                                   override=override)
        elif format_type == "pdb":
            write_pdb_from_table(outputfile=self.outputfile, data=self.data,
                                 override=override, pro_format=pro_format)
        else:
            message = 'The provided format {} is not implemented...'.format(format_type)
            raise ValueError(message)
        return


if __name__ == '__main__':
    pass
