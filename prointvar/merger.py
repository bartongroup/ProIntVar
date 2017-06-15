#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""

This defines the methods that work with structures and sequences.

Fábio Madeira, 2017+

"""

from __future__ import print_function

import os
import pandas as pd

from prointvar.dssp import DSSPgenerator
from prointvar.dssp import DSSPreader
from prointvar.dssp import get_dssp_selected_from_table
from prointvar.sifts import SIFTSreader
from prointvar.sifts import get_sifts_selected_from_table
from prointvar.mmcif import MMCIFreader
from prointvar.mmcif import get_mmcif_selected_from_table
from prointvar.fetchers import fetch_best_structures_pdbe

from prointvar.utils import flash

from prointvar.config import config


class TableMergerError(Exception):
    pass


def mmcif_sifts_table_merger(mmcif_table, sifts_table):
    """
    Merge the mmCIF and SIFTS tables.
    
    :param mmcif_table: mmCIF pandas DataFrame
    :param sifts_table: SIFTS pandas DataFrame
    :return: merged pandas DataFrame
    """

    # bare minimal columns needed
    if ('auth_seq_id_full' in mmcif_table and 'label_asym_id' in mmcif_table and
            'PDB_dbResNum' in sifts_table and 'PDB_entityId' in sifts_table):

        # workaround for BioUnits
        if 'orig_label_asym_id' in mmcif_table:
            table = mmcif_table.merge(sifts_table, how='left',
                                      left_on=['auth_seq_id_full', 'orig_label_asym_id'],
                                      right_on=['PDB_dbResNum', 'PDB_entityId'])
        else:
            table = mmcif_table.merge(sifts_table, how='left',
                                      left_on=['auth_seq_id_full', 'label_asym_id'],
                                      right_on=['PDB_dbResNum', 'PDB_entityId'])

    else:
        raise TableMergerError('Not possible to merge mmCIF and SIFTS table! '
                               'Some of the necessary columns are missing...')
    return table


def mmcif_dssp_table_merger(mmcif_table, dssp_table):
    """
    Merge the mmCIF and DSSP tables.
    
    :param mmcif_table: mmCIF pandas DataFrame
    :param dssp_table: DSSP pandas DataFrame
    :return: merged pandas DataFrame
    """

    # bare minimal columns needed
    if ('new_seq_id' in mmcif_table and 'new_asym_id' in mmcif_table and
            'RES' in dssp_table and 'CHAIN' in dssp_table):

        table = mmcif_table.merge(dssp_table, how='left',
                                  left_on=['new_seq_id', 'new_asym_id'],
                                  right_on=['RES', 'CHAIN_FULL'])
    elif ('auth_seq_id_full' in mmcif_table and 'label_asym_id' in mmcif_table and
            'RES' in dssp_table and 'CHAIN_FULL' in dssp_table):

        table = mmcif_table.merge(dssp_table, how='left',
                                  left_on=['auth_seq_id_full', 'label_asym_id'],
                                  right_on=['RES', 'CHAIN_FULL'])
    else:
        raise TableMergerError('Not possible to merge mmCIF and DSSP table! '
                               'Some of the necessary columns are missing...')
    return table


def dssp_sifts_table_merger(dssp_table, sifts_table):
    """
    Merge the DSSP and SIFTS tables.
    
    :param dssp_table: mmCIF pandas DataFrame
    :param sifts_table: SIFTS pandas DataFrame
    :return: merged pandas DataFrame
    """

    # bare minimal columns needed
    if ('RES' in dssp_table and 'CHAIN' in dssp_table and
            'PDB_dbResNum' in sifts_table and 'PDB_entityId' in sifts_table):

        table = dssp_table.merge(sifts_table, how='left',
                                 left_on=['RES', 'CHAIN'],
                                 right_on=['PDB_dbResNum', 'PDB_entityId'])

    else:
        raise TableMergerError('Not possible to merge DSSP and SIFTS table! '
                               'Some of the necessary columns are missing...')
    return table


def dssp_dssp_table_merger(bound_table, unbound_table):
    """
    Merge the DSSP bound and unbound table.
    
    :param bound_table: DSSP bound pandas DataFrame
    :param unbound_table: DSSP unbound pandas DataFrame
    :return: merged pandas DataFrame
    """

    # assumes a pre-merge with mmcif was performed, and that bound dssp is observed
    if ('label_seq_id' in bound_table and 'label_asym_id' in bound_table and
            'RES' in unbound_table and 'CHAIN_FULL' in unbound_table and
            'RES' in bound_table and 'CHAIN_FULL' in bound_table):

        # rename the unbound columns
        cols = ('ACC', 'SS', 'SS_CLASS', 'RSA', 'RSA_CLASS', 'RES', 'CHAIN_FULL')
        excluded = tuple([k for k in unbound_table if k not in cols])
        if excluded is not None:
            unbound_table = unbound_table.drop(list(excluded), axis=1)
        unbound_table.columns = ['{}_UNB'.format(k) for k in unbound_table]

        table = bound_table.merge(unbound_table, how='left',
                                  left_on=['label_seq_id', 'label_asym_id'],
                                  right_on=['RES_UNB', 'CHAIN_FULL_UNB'])
    # bare minimal columns needed
    elif ('RES' in bound_table and 'CHAIN_FULL' in bound_table and
            'RES' in unbound_table and 'CHAIN_FULL' in unbound_table):

        # rename the unbound columns
        cols = ('ACC', 'SS', 'SS_CLASS', 'RSA', 'RSA_CLASS', 'RES', 'CHAIN_FULL')
        excluded = tuple([k for k in unbound_table if k not in cols])
        if excluded is not None:
            unbound_table = unbound_table.drop(list(excluded), axis=1)
        unbound_table.columns = ['{}_UNB'.format(k) for k in unbound_table]

        table = bound_table.merge(unbound_table, how='left',
                                  left_on=['RES', 'CHAIN_FULL'],
                                  right_on=['RES_UNB', 'CHAIN_FULL_UNB'])
    else:
        raise TableMergerError('Not possible to merge the DSSP tables! '
                               'Some of the necessary columns are missing...')

    return table


def table_merger(mmcif_table=None, dssp_table=None, sifts_table=None):
    """
    Merges the tables provided using appropriate columns.
    
    :param mmcif_table: (optional) mmCIF pandas DataFrame
    :param dssp_table: (optional) DSSP pandas DataFrame
    :param sifts_table: (optional) SIFTS pandas DataFrame
    :return: merged pandas DataFrame
    """

    if mmcif_table is not None and dssp_table is not None and sifts_table is not None:
        mmcif_table = mmcif_dssp_table_merger(mmcif_table, dssp_table)
        table = mmcif_sifts_table_merger(mmcif_table, sifts_table)

    elif mmcif_table is not None and dssp_table is not None:
        table = mmcif_dssp_table_merger(mmcif_table, dssp_table)

    elif mmcif_table is not None and sifts_table is not None:
        table = mmcif_sifts_table_merger(mmcif_table, sifts_table)

    elif dssp_table is not None and sifts_table is not None:
        table = dssp_sifts_table_merger(dssp_table, sifts_table)
    else:
        return TableMergerError('Nothing to merge...')

    return table


def table_generator(uniprot_id=None, pdb_id=None, chain=None, res=None,
                    site=None, atom=None, lines=None, bio=True, add_dssp=True):
    """
    Simplifies the process of generating tables and merging them.
    
    :param uniprot_id: (str) UniProt ID
    :param pdb_id: (str) PDB ID
    :param chain: (tuple) Chain ID ('*_asym_id')
    :param res: (tuple) PDB ResNum ('*_seq_id_full')
    :param site: (tuple) UniProt ResNum (positional index)
    :param atom: (tuple) atom IDs or None
    :param lines: 'ATOM' or 'HETATM' or None (both)
    :param bio: boolean for using AsymUnits or BioUnits
    :param add_dssp: boolean
    """

    # generates the tables for the uniprot_id, pdb_id
    # chain, res, site are only processed after
    if uniprot_id or pdb_id:
        # get pdb_id if only uniprot_id is provided
        if uniprot_id and not pdb_id:
            # fetch the best structures from the PDBe API
            data = fetch_best_structures_pdbe(uniprot_id)
            if data is not None:
                # uses the first structure
                pdb_id = data[0]['pdb_id']
                chain = (data[0]['chain_id'],)
            else:
                flash('Best structures not available from the PDBe API for '
                      '{}'.format(uniprot_id))
                raise TableMergerError('Nothing to merge...')

        # mmCIF table
        if bio:
            inputcif = "{}{}{}.cif".format(config.db_root, config.db_cif_biounit, pdb_id)
        else:
            inputcif = "{}{}{}.cif".format(config.db_root, config.db_cif, pdb_id)

        r = MMCIFreader(inputcif)
        mmcif_table = r.atoms(add_res_full=True, add_contacts=False)
        mmcif_table = get_mmcif_selected_from_table(mmcif_table, chain=chain, res_full=res,
                                                    atom=atom, lines=lines, category='auth')

        # SIFTS table
        inputsifts = "{}{}{}.xml".format(config.db_root, config.db_sifts_xml, pdb_id)

        r = SIFTSreader(inputsifts)
        sifts_table = r.residues(add_regions=True, add_dbs=False)
        sifts_table = get_sifts_selected_from_table(sifts_table, chain=chain, res=res,
                                                    uniprot=uniprot_id, site=site)

        # DSSP table
        if add_dssp:
            if bio:
                outputdssp = "{}{}{}_bio.dssp" \
                            "".format(config.db_root, config.db_dssp_generated, pdb_id)
            else:
                outputdssp = "{}{}{}.dssp" \
                            "".format(config.db_root, config.db_dssp_generated, pdb_id)
            w = DSSPgenerator(inputcif, outputdssp)
            w.run()
            if os.path.exists(outputdssp):
                r = DSSPreader(outputdssp)
                try:
                    dssp_table = r.residues(add_ss_reduced=True, add_rsa_class=True)
                except IndexError:
                    print(outputdssp)
                    raise IndexError
                dssp_table = get_dssp_selected_from_table(dssp_table, chain_full=chain, res=res)
            else:
                dssp_table = None
        else:
            dssp_table = None

        return mmcif_table, dssp_table, sifts_table

    else:
        raise ValueError('No UniProt ID or PDB ID provided...')


class TableMerger(object):
    def __init__(self, mmcif_table=None, dssp_table=None, sifts_table=None,
                 store=True, verbose=False):

        self.mmcif_table = mmcif_table
        self.dssp_table = dssp_table
        self.sifts_table = sifts_table
        self.store = store
        self.merged_table = None
        self.filename = None
        self.verbose = verbose

    # TODO move out these private methods
    def _dump_merged_table(self, outputfile, override=False):
        if not os.path.isfile(outputfile) or override:
            self.merged_table.to_pickle(outputfile)
            if self.verbose:
                flash('{} stored locally...'.format(outputfile))
        else:
            flash('{} already available...'.format(outputfile))

    def _load_merged_table(self, inputfile):
        if os.path.isfile(inputfile):
            if self.verbose:
                flash('{} is already available...'.format(inputfile))
            self.merged_table = pd.read_pickle(inputfile)
            return self.merged_table
        else:
            raise IOError("{} not available or could not be read..."
                          "".format(inputfile))

    @staticmethod
    def _get_filename(uniprot_id=None, pdb_id=None, chain=None, res=None,
                      site=None, atom=None, lines=None, bio=True, add_dssp=True):

        if uniprot_id is not None and pdb_id is None and chain is None:
            name = "{}".format(uniprot_id)
        elif uniprot_id is not None and pdb_id is not None and chain is None:
            name = "{}_{}".format(uniprot_id, pdb_id)
        elif uniprot_id is None and pdb_id is not None and chain is None:
            name = "{}".format(pdb_id)
        elif uniprot_id is None and pdb_id is not None and chain is not None:
            name = "{}_{}".format(pdb_id, "_".join(list(chain)))
        elif uniprot_id is None and pdb_id is not None and chain is not None and \
                res is not None:
            name = "{}_{}_{}".format(pdb_id, "_".join(list(chain)), "_".join(list(res)))
        elif uniprot_id is not None and pdb_id is not None and chain is not None and \
                site is not None:
            name = "{}_{}_{}_{}".format(uniprot_id, "_".join(list(site)),
                                        pdb_id, "_".join(list(chain)))
        elif uniprot_id is not None and pdb_id is not None and chain is not None and \
                site is None:
            name = "{}_{}_{}".format(uniprot_id, pdb_id, "_".join(list(chain)))
        elif uniprot_id is not None and pdb_id is None and chain is None and \
                site is not None:
            name = "{}_{}".format(uniprot_id, "_".join(list(site)))
        else:
            raise ValueError('No UniProt ID or PDB ID provided...')
        # TODO use pickled dir?
        # db_pickle
        if bio:
            filename = "{}{}{}_bio.pkl".format(config.db_root, config.tmp_dir_local, name)
        else:
            filename = "{}{}{}.pkl".format(config.db_root, config.tmp_dir_local, name)
        return filename

    def merge(self, outputfile=None, override=False):
        """Merges the provided tables and stores"""

        self.merged_table = table_merger(mmcif_table=self.mmcif_table,
                                         dssp_table=self.dssp_table,
                                         sifts_table=self.sifts_table)
        if self.store and outputfile is not None:
            self.filename = outputfile
            self._dump_merged_table(self.filename, override=override)
        return self.merged_table

    def run(self, uniprot_id=None, pdb_id=None, chain=None, res=None, site=None,
            atom=None, lines=None, bio=True, add_dssp=True, outputfile=None, override=False):
        """Generates the tables, merges and stores"""

        self.mmcif_table, self.dssp_table, self.sifts_table = \
            table_generator(uniprot_id=uniprot_id, pdb_id=pdb_id, chain=chain, res=res,
                            site=site, atom=atom, lines=lines, bio=bio, add_dssp=add_dssp)
        self.merged_table = self.merge(outputfile=outputfile, override=override)

        if self.store:
            if outputfile is not None:
                self.filename = outputfile
            else:
                self.filename = self._get_filename(uniprot_id=uniprot_id, pdb_id=pdb_id,
                                                   chain=chain, res=res, site=site,
                                                   atom=atom, lines=lines,
                                                   bio=bio, add_dssp=add_dssp)
            self._dump_merged_table(self.filename, override=override)

        return self.merged_table

    def load(self, uniprot_id=None, pdb_id=None, chain=None, res=None,
             site=None, atom=None, lines=None, bio=True, add_dssp=True, inputfile=None):
        """Loads a merged table"""

        if inputfile is not None:
            self.filename = inputfile
        else:
            self.filename = self._get_filename(uniprot_id=uniprot_id, pdb_id=pdb_id,
                                               chain=chain, res=res, site=site,
                                               atom=atom, lines=lines,
                                               bio=bio, add_dssp=add_dssp)
        self.merged_table = self._load_merged_table(self.filename)
        return self.merged_table


if __name__ == '__main__':
    pass
