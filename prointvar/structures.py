# -*- coding: utf-8 -*-


"""
Extends the functionality of ProteoFAV.structures

Fábio Madeira, 2017+
"""

import logging

from prointvar.utils import get_new_pro_ids


logger = logging.getLogger("prointvar")


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


