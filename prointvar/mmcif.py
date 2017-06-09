#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""

This defines the methods that work with SIFTS files.

Fábio Madeira, 2017+

"""

import os
import json
import shlex
import pandas as pd
from io import StringIO
from collections import OrderedDict

from scipy.spatial import cKDTree

from prointvar.contacts import get_interaction_chains
from prointvar.contacts import get_interaction_molecules
from prointvar.contacts import get_interaction_topologies
from prointvar.contacts import get_interaction_properties
from prointvar.contacts import get_distance_between_atoms
from prointvar.contacts import get_distance_between_atoms_vdw

from prointvar.utils import flash
from prointvar.utils import row_selector
from prointvar.utils import string_split

from prointvar.library import mmcif_types

_PDB_FORMAT = "%s%5i %-4s%c%3s %c%4s%c   %8.3f%8.3f%8.3f%s%6.2f      %4s%2s%2s\n"


def parse_mmcif_atoms_from_file(inputfile, excluded=(), add_res_full=True,
                                add_contacts=False, dist=5, first_model=True,
                                add_atom_altloc=False, verbose=False):
    """
    Parse mmCIF ATOM and HETATM lines.

    :param inputfile: path to the mmCIF file
    :param excluded: option to exclude mmCIF columns
    :param add_res_full: option to extend the table with 'res_full'
        i.e. res_number + insertion_code (e.g. '12A')
    :param add_contacts: boolean
    :param dist: distance threshold in Angstrom
    :param first_model: boolean
    :param add_atom_altloc: boolean
    :param verbose: boolean
    :return: returns a pandas DataFrame
    """

    if verbose:
        flash("Parsing mmCIF atoms from lines...")

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

    if add_res_full:
        table = add_mmcif_res_full(table)

    if add_atom_altloc:
        table = add_mmcif_atom_altloc(table)

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


def parse_pdb_atoms_from_file(inputfile, excluded=(), add_contacts=False, dist=5,
                              first_model=True, add_atom_altloc=False, verbose=False):
    """
    Parse PDB ATOM and HETATM lines.

    :param inputfile: path to the PDB file
    :param excluded: option to exclude mmCIF columns
    :param add_contacts: boolean
    :param dist: distance threshold in Angstrom
    :param first_model: boolean 
    :param add_atom_altloc: boolean
    :param verbose: boolean
    :return: returns a pandas DataFrame
    """

    if verbose:
        flash("Parsing PDB atoms from lines...")

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
            if line.startswith("MODEL"):
                modelnumb = line.split()[1]
            elif line.startswith("ATOM"):
                lines.append(line + "%s" % modelnumb)
            elif line.startswith("HETATM"):
                lines.append(line + "%s" % modelnumb)
    lines = "\n".join(lines)

    header = ('group_PDB', 'id', 'label_atom_id', 'label_alt_id', 'label_comp_id',
              'label_asym_id', 'label_seq_id_full', 'Cartn_x', 'Cartn_y', 'Cartn_z',
              'occupancy', 'B_iso_or_equiv', 'type_symbol', 'auth_atom_id', 'auth_comp_id',
              'auth_asym_id', 'auth_seq_id_full', 'pdbx_PDB_model_num')

    widths = ((0, 6), (6, 11), (12, 16), (16, 17), (17, 20), (21, 22), (22, 26),
              (30, 38), (38, 46), (46, 54), (54, 60), (60, 66), (76, 78), # (72, 76), ('seg_id')
              (12, 16), (17, 20), (21, 22), (22, 26), (78, 79))

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

    # adding a 'label_seq_id' , 'auth_seq_id' and 'pdbx_PDB_ins_code'
    table = add_mmcif_res_split(table)
    # adding the 'label_alt_id
    table = add_label_alt_id(table)

    # table modular extensions
    if add_contacts:
        table = add_mmcif_contacts(table, dist=dist)

    if add_atom_altloc:
        table = add_mmcif_atom_altloc(table)

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


def parse_mmcif_categories_from_file(inputfile, excluded=(), category=None, verbose=False):
    """
    Generic method that gets all categories and fields (except for the .atom_site*)

    :param inputfile: path to the mmCIF file
    :param excluded: option to exclude mmCIF categories
    :param category: data category to be parsed (otherwise all)
    :param verbose: boolean
    :return: returns a nested dictionary
    """

    if verbose:
        flash("Parsing mmCIF categories from lines...")

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

    if res is not None:
        table = row_selector(table, '{}_seq_id'.format(category), res, method="isin")

    if res_full is not None:
        table = row_selector(table, '{}_seq_id_full'.format(category), res_full, method="isin")

    if comp is not None:
        table = row_selector(table, '{}_comp_id'.format(category), comp, method="isin")

    if atom is not None:
        table = row_selector(table, '{}_atom_id'.format(category), atom, method="isin")

    if lines is not None:
        table = row_selector(table, 'group_PDB', lines, method="isin")

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
    columns_to_agg['id'.format(category)] = 'first'
    table = table.groupby(by=agg_cols, as_index=False).agg(columns_to_agg)
    return table.sort_values(by='id').reset_index()


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


def get_mmcif_full_contacts_from_table(data, dist=5, add_contact_info=True,
                                       chain=None, res=None, atom=None, lines=None,
                                       category='label', ignore_same_chain=False,
                                       ignore_same_res=True, ignore_consecutive=3,
                                       left_sided=True):
    """
    Finding residues/atoms within a distance threshold. The new DataFrame only contains
    rows that participate in interactions.
    
    :param data: pandas DataFrame
    :param dist: distance threshold in Angstrom
    :param add_contact_info: boolean
    :param chain: (tuple) chain IDs or None
    :param res: (tuple) res IDs or None (needs chain)
    :param atom: (tuple) atom IDs or None (needs chain and res)
    :param lines: 'ATOM' or 'HETATM' or None (both)
    :param category: data category to be used as precedence in _atom_site.*_*
        asym_id, seq_id and atom_id
    :param dist_vdw: distance with Van der Waals threshold in Angstrom
    :param ignore_same_chain: ignores contacts within the same chain
    :param ignore_same_res: ignores contacts within the same res
    :param ignore_consecutive: number of surrounding residues (in sequence in both 
        directions) that will be ignored
    :param left_sided: (boolean) if True only A->B index are reported, otherwise
        both A->B and B->A pairs are reported
    :return: new pandas DataFrame
    """

    table = data

    if 'contact_indexes' not in table:
        table = add_mmcif_contacts(table, dist=dist)

    # filter based on chain, res, atom and lines
    ftable = get_mmcif_selected_from_table(table, chain=chain, res=res, atom=atom,
                                           lines=lines, category=category)
    ftable = ftable.copy()

    # get indexes from the table
    indexes = []
    for ix in ftable.index:
        indexes.append([int(i) for i in ftable.loc[ix, 'contact_indexes'].split(',')])

    flatten = []
    for i, ix in enumerate(ftable.index):
        pairs = [(ix, jx) for jx in indexes[i]]
        flatten += tuple(pairs)

    if left_sided:
        # if pairs are present in both directions (i.e. A->B and B->A)
        # remove one of the pairs from the bottom up
        for pair in flatten[::-1]:
            inv_pair = (pair[1], pair[1])
            if inv_pair in flatten:
                flatten.remove(pair)

    # get new DataFrame
    rows = []
    header1 = list(table)
    header2 = list(map(lambda x: '%s_2' % x, list(table)))
    header = header1 + header2
    for entry in flatten:
        i = entry[0]
        j = entry[1]
        values = list(table.loc[i, :]) + list(table.loc[j, :])
        assert len(header) == len(values)
        new_dict = {k: v for k, v in zip(header, values)}
        rows.append(new_dict)
    table = pd.DataFrame(rows)

    # filter based on the type of contacts
    if ignore_same_res:
        table = table.loc[((table['{}_asym_id'.format(category)] !=
                            table['{}_asym_id_2'.format(category)]) |
                           ((table['{}_asym_id'.format(category)] ==
                             table['{}_asym_id_2'.format(category)]) &
                            (table['{}_seq_id'.format(category)] !=
                             table['{}_seq_id_2'.format(category)])))]

    if ignore_same_chain:
        table = table.loc[(table['{}_asym_id'.format(category)] !=
                           table['{}_asym_id_2'.format(category)])]

    if ignore_consecutive:
        try:
            # FIXME use the query
            table = table.loc[((table['{}_asym_id'.format(category)] !=
                                table['{}_asym_id_2'.format(category)]) |
                               ((table['{}_asym_id'.format(category)] ==
                                 table['{}_asym_id_2'.format(category)]) &
                                (abs(table['auth_seq_id'].astype(int) -
                                     table['auth_seq_id_2'].astype(int)) >=
                                 ignore_consecutive)))]
        except ValueError as e:
            if 'pdbe_label_seq_id' in table and 'pdbe_label_seq_id_2' in table:
                table = table.loc[((table['{}_asym_id'.format(category)] !=
                                    table['{}_asym_id_2'.format(category)]) |
                                   ((table['{}_asym_id'.format(category)] ==
                                     table['{}_asym_id_2'.format(category)]) &
                                    (abs(table['auth_seq_id'].astype(int) -
                                         table['auth_seq_id_2'].astype(int)) >=
                                     ignore_consecutive)))]
            else:
                flash("Skipping 'ignore_consecutive' with {}...".format(e))

    # extend the table with information about the interaction
    table = table.apply(lambda x: x.fillna('-'))
    if add_contact_info:
        table = add_mmcif_contact_info(table)

    if table.empty:
        raise ValueError('Filtering resulted in an empty DataFrame...')

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
    atom_lines = ['data_Autogenerated_mmCIF', 'loop_']
    atom_lines += ["_atom_site.{}".format(v) for v in list(table)]
    for i in table.index:
        line = ' '.join([str(v) for v in list(table.loc[i, :])])
        atom_lines.append(line)

    # write the final output
    if not os.path.exists(outputfile) or override:
        with open(outputfile, 'w') as outlines:
            outlines.write("\n".join(atom_lines))
    else:
        flash('mmCIF for {} already available...'.format(outputfile))
    return


def write_pdb_from_table(outputfile, data, override=False):
    """
    Generic method that writes 'atom' lines in PDB format.
    
    :param outputfile: path to the PDB file
    :param data: pandas DataFrame object
    :param override: boolean
    :return: (side effects) writes to file
    """

    table = data
    atom_lines = []
    atom_number = 0
    for i in table.index:
        atom_number += 1
        atom_lines.append(get_atom_line(data=table, index=i,
                                        atom_number=atom_number))

    # write the final output
    if not os.path.exists(outputfile) or override:
        with open(outputfile, 'w') as outlines:
            outlines.write("".join(atom_lines))
    else:
        flash('PDB for {} already available...'.format(outputfile))
    return


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


def add_mmcif_res_split(data):
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


def add_label_alt_id(data):
    """
    Utility that adds fixes the label_alt_id column to match what is
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


def add_mmcif_atom_altloc(data):
    """
    Utility that adds new columns to the table.
    adds: 'label_atom_altloc_id', and 'auth_atom_altloc_id'
    
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


def add_mmcif_contact_info(data, category='label'):
    """
    Utility that adds new columns to the table. Columns added are:
    'distance', 'distance_vdw', 'int_types', 'int_atom' and 'int_res'.

    :param data: pandas DataFrame object
    :param category: data category to be used as precedence in _atom_site.*_*
        asym_id, seq_id and atom_id
    :return: returns a modified pandas DataFrame
    """

    table = data
    headers = [x for x in table.columns.values if x.endswith('_2')]
    if headers:
        table['distance'] = table.apply(get_distance_between_atoms, axis=1)
        table['distance_vdw'] = table.apply(get_distance_between_atoms_vdw,
                                            axis=1, args=(category, ))
        table['int_chains'] = table.apply(get_interaction_chains,
                                          axis=1, args=(category, ))
        table['int_molecules'] = table.apply(get_interaction_molecules,
                                             axis=1, args=(category, ))
        table['int_topologies'] = table.apply(get_interaction_topologies,
                                              axis=1, args=(category, ))
        table['int_properties'] = table.apply(get_interaction_properties,
                                              axis=1, args=(category, ))
        table['int_atom'] = 'A'
        table['int_res'] = '-'

        # temp table with used to find the closest atom-atom distances used for
        # reporting 'int_res'
        ntable = table.groupby(['{}_asym_id'.format(category),
                                '{}_asym_id_2'.format(category),
                                '{}_seq_id_full'.format(category),
                                '{}_seq_id_full_2'.format(category)],
                               as_index=False)['distance'].min()

        # add 'int_res' values for the closest atom-atom interactions
        for ix in ntable.index:
            # FIXME use query
            tmptable = table[((table['{}_asym_id'.format(category)] ==
                               ntable.loc[ix, '{}_asym_id'.format(category)]) &
                              (table['{}_asym_id_2'.format(category)] ==
                               ntable.loc[ix, '{}_asym_id_2'.format(category)]) &
                              (table['{}_seq_id_full'.format(category)] ==
                               ntable.loc[ix, '{}_seq_id_full'.format(category)]) &
                              (table['{}_seq_id_full_2'.format(category)] ==
                               ntable.loc[ix, '{}_seq_id_full_2'.format(category)]) &
                              (table['distance'] == ntable.loc[ix, 'distance']))]
            index = tmptable.index.tolist()[0]
            table.set_value(index, 'int_res', 'R')

    else:
        raise ValueError('This method expects a DataFrame with contact pairs...')
    return table


def get_atom_line(data, index, atom_number, coords=None, new_chain=None):
    """
    Returns an ATOM PDB-formatted string.
    (Based on code from the PDB module in Biopython.)

    :param data: pandas DataFrame object
    :param index: atom index
    :param atom_number: incremental number
    :param coords: list of transformed coordinates
    :param new_chain: if true defaults to chain "A"
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

    name = table.loc[ix, 'label_atom_id']
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

    resname = table.loc[ix, 'label_comp_id']

    # FIXME is this needed?
    # overriding the original chain name
    if new_chain is not None:
        chain_id = new_chain
    else:
        chain_id = table.loc[ix, 'label_asym_id']
        chain_id = chain_id[0]

    try:
        resseq = str(table.loc[ix, 'auth_seq_id'])
    except:
        resseq = str(table.loc[ix, 'label_seq_id'])

    icode = table.loc[ix, 'pdbx_PDB_ins_code']
    if icode == "?":
        icode = " "

    # FIXME is this needed?
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


class MMCIFreader(object):
    def __init__(self, inputfile, verbose=False):
        """
        :param inputfile: Needs to point to a valid mmCIF file.
        :param verbose: boolean
        """
        self.inputfile = inputfile
        self.verbose = verbose
        self.data = None
        # defaults
        self.excluded = ("Cartn_x_esd", "Cartn_y_esd", "Cartn_z_esd",
                         "occupancy_esd", "B_iso_or_equiv_esd", "pdbx_formal_charge")

        if not os.path.isfile(inputfile):
            raise IOError("{} not available or could not be read...".format(inputfile))

    def read(self, **kwargs):
        return self.atoms(**kwargs)

    def atoms(self, excluded=None, add_res_full=True, add_contacts=False, dist=5,
              first_model=True, add_atom_altloc=False, format_type="mmcif",
              residue_agg=False, agg_method='centroid', category='label'):
        if excluded is None:
            excluded = self.excluded
        if format_type == "mmcif":
            self.data = parse_mmcif_atoms_from_file(self.inputfile, excluded=excluded,
                                                    add_res_full=add_res_full,
                                                    add_contacts=add_contacts, dist=dist,
                                                    first_model=first_model,
                                                    add_atom_altloc=add_atom_altloc,
                                                    verbose=self.verbose)

        elif format_type == "pdb":
            self.data = parse_pdb_atoms_from_file(self.inputfile, excluded=excluded,
                                                  add_contacts=add_contacts, dist=dist,
                                                  first_model=first_model,
                                                  add_atom_altloc=add_atom_altloc,
                                                  verbose=self.verbose)
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
                                                     category=category, verbose=self.verbose)
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
            flash('No mmCIF data parsed...')


class MMCIFwriter(object):
    def __init__(self, inputfile=None, outputfile=None, verbose=False):
        """
        :param inputfile: Needs to point to a valid mmCIF file.
        :param outputfile: if not provided will use the same file name and
            <_filtered.cif> extension
        :param verbose: boolean
        """
        self.inputfile = inputfile
        self.outputfile = outputfile
        self.verbose = verbose
        self.data = None

        # generate outputfile if missing
        if self.outputfile is None:
            self._generate_output()

        if inputfile is not None and not os.path.isfile(inputfile):
            raise IOError("{} not available or could not be read...".format(inputfile))

    def _generate_output(self):
        if not self.outputfile:
            filename, extension = os.path.splitext(self.inputfile)
            self.outputfile = filename + "_filtered.cif"

    def run(self, data=None, chain=None, res=None, atom=None, lines=None, category='label',
            override=False, format_type="mmcif"):
        # writes a new mmCIF with selected chains, residues, atoms, or lines
        if data is None:
            r = MMCIFreader(inputfile=self.inputfile)
            # guess the input format as 'format_type' refers to the output format
            try:
                data = r.read(excluded=(), format_type="mmcif")
            except ValueError:
                data = r.read(excluded=(), format_type="pdb")
        self.data = get_mmcif_selected_from_table(data, chain=chain, res=res, atom=atom,
                                                  lines=lines, category=category)
        if format_type == "mmcif":
            write_mmcif_from_table(outputfile=self.outputfile, data=self.data,
                                   override=override)
        elif format_type == "pdb":
            write_pdb_from_table(outputfile=self.outputfile, data=self.data,
                                 override=override)
        else:
            message = 'The provided format {} is not implemented...'.format(format_type)
            raise ValueError(message)
        return


if __name__ == '__main__':
    pdbid = "2pah"
    # pdbid = "1cg2"
    # pdbid = "3kic"

    from prointvar.config import config as c

    inputcif = "{}{}{}.cif".format(c.db_root, c.db_cif, pdbid)

    # d = MMCIFreader(inputcif)
    # d.read()
    # nd = d.data
    # nd = d.to_json()
    # print(nd)
    # nd = json.loads(nd)
    # print([k for k in nd[0]])
    # print(nd.loc[0, "label_asym_id"])
    # print(nd.loc[:, "label_asym_id"].unique())
    # print(nd.loc[:, "auth_asym_id"].unique())
    # print(nd.tail())

    # d = MMCIFreader(inputcif)
    # d.categories()
    # d.categories(category="_citation")
    # d.categories(excluded=('_atom_site',))
    # d = d.to_json()
    # # print(d)
    # print(json.loads(d)["_citation"])

    # d = MMCIFreader(inputcif)
    # d.read()
    # d = d.to_json()
    # print(d)

    # outputcif = "{}{}{}.cif".format(c.db_root, c.tmp_dir_local, pdbid)
    # data = get_mmcif_selected_from_table(nd, chain=('A',), res=None, atom=None, lines=None)
    # data = get_mmcif_selected_from_table(nd, chain=None, res=('1', '2', '3'), atom=None, lines=None)
    # data = get_mmcif_selected_from_table(nd, chain=('A',), res=('1', '2', '3'), atom=None, lines=None)
    # data = get_mmcif_selected_from_table(nd, chain=None, res=None, atom=('CA',), lines=None)
    # data = get_mmcif_selected_from_table(nd, chain=None, res=None, atom=None, lines=('HETATM',))
    # data = get_mmcif_selected_from_table(nd, chain=('A',), res=None, atom=None, lines=('HETATM',))
    # write_mmcif_from_table(outputcif, data)

    # d = MMCIFreader(inputcif)
    # data = d.read(add_contacts=True)
    # w = MMCIFwriter(inputcif, outputcif)
    # w.run(data=data, chain=None, res=None, atom=None, lines=None)
    #
    # d = MMCIFreader(outputcif)
    # d.read(add_contacts=False)
    # nd = d.data
    # nd = d.to_json()
    # print(nd)

    inputcif = "{}{}{}.cif".format(c.db_root, c.db_cif_biounit, pdbid)
    d = MMCIFreader(inputcif)
    d = d.read()
    # get the auth_seq_id from the label_seq_id
    seqid = d[(d['auth_seq_id'] == '285') & (d['label_asym_id'] == 'A')]
    # seqid = mmcif[mmcif['label_asym_id'] == 'A']
    # seqid = mmcif[mmcif['auth_seq_id'] == 285 & ]
    seqid = seqid.loc[seqid.index[0], 'label_seq_id']
    seqids = (seqid,)

    cmmcif = get_mmcif_full_contacts_from_table(d, dist=5, add_contact_info=True,
                                                chain=('A',), res=seqids, atom=None,
                                                lines=None, category='label',
                                                ignore_same_chain=False,
                                                ignore_same_res=True,
                                                ignore_consecutive=3)

    table = cmmcif
    # print(table.groupby([k for k in table if k != 'distance'],
    #                     as_index=False)['distance'].min())
    # print(len(table))
    # print(table.loc[:, ('label_asym_id', 'label_seq_id', 'label_atom_id',
    #                     'label_asym_id_2', 'label_seq_id_2', 'label_atom_id_2')])
    print(table.loc[:, ('label_asym_id', 'label_seq_id', 'label_atom_id',
                        'label_asym_id_2', 'label_seq_id_2', 'label_atom_id_2',
                        'distance', 'distance_vdw', 'int_chains', 'int_molecules',
                        'int_topologies', 'int_properties', 'int_atom', 'int_res')])

    # inputcif = "{}/sgc_fragments/BRD1A/BRD1A-x038_event1.pdb".format(c.db_root)
    # d = MMCIFreader(inputcif)
    # d = d.atoms(format_type='pdb', add_atom_altloc=True)
    #
    # nd = row_selector(d, 'label_comp_id', 'LIG', 'equals')
    # print(nd.loc[:, 'label_atom_altloc_id'])

    # d = MMCIFreader(inputcif)
    # d = d.atoms()
    # nd = residues_aggregation(d, agg_method='unique').reset_index()
    # print(nd.loc[0, :])

    pass
