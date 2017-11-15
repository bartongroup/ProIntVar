# -*- coding: utf-8 -*-


"""
Extends the functionality of ProteoFAV.structures

Fábio Madeira, 2017+
"""

import os
import shlex
import logging
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from collections import OrderedDict

from prointvar.utils import get_new_pro_ids


logger = logging.getLogger("prointvar")


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


def add_mmcif_new_pro_ids(data, category='auth'):
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


def get_coordinates(table):
    """
    Get coordinates from the provided PDBx table
    and return a vector-set with all the coordinates.

    :param table: pandas DataFrame from PDBXreader
    :returns: coordinates - array (N,3) where N is number of atoms
    """

    # assumes it's a valid PDBx table
    if isinstance(table, pd.DataFrame):
        coords = [np.array([table.loc[i, "Cartn_x"],
                            table.loc[i, "Cartn_y"],
                            table.loc[i, "Cartn_z"]],
                           dtype=float) for i in table.index]
    else:
        return ValueError("Pandas DataFrame is not valid...")
    return coords

