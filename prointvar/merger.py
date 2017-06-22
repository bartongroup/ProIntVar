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


def table_merger(mmcif_table=None, dssp_table=None, sifts_table=None,
                 contacts_table=None):
    """
    Merges the tables provided using appropriate columns.
    
    :param mmcif_table: (optional) mmCIF pandas DataFrame
    :param dssp_table: (optional) DSSP pandas DataFrame
    :param sifts_table: (optional) SIFTS pandas DataFrame
    :param contacts_table: (optional) Arpeggio pandas DataFrame
    :return: merged pandas DataFrame
    """

    available = [mmcif_table, dssp_table, sifts_table, contacts_table]
    available = [k for k in available if k is not None]
    if len(available) < 2:
        return TableMergerError('Nothing to merge...')

    table = None
    if mmcif_table is not None:
        if dssp_table is not None:
            mmcif_table = mmcif_dssp_table_merger(mmcif_table, dssp_table)
        if sifts_table is not None:
            mmcif_table = mmcif_sifts_table_merger(mmcif_table, sifts_table)
        # if contacts_table is not None:
        #     mmcif_table = contacts_mmcif_table_merger(contacts_table, mmcif_table)
        table = mmcif_table

    elif dssp_table is not None:
        if sifts_table is not None:
            dssp_table = dssp_sifts_table_merger(dssp_table, sifts_table)
        table = dssp_table

    return table


def table_generator(uniprot_id=None, pdb_id=None, chain=None, res=None,
                    site=None, atom=None, lines=None, bio=True, sifts=True, dssp=True,
                    dssp_unbound=False, contacts=False, override=False):
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
    :param sifts: boolean
    :param dssp: boolean
    :param dssp_unbound: (boolean) if true runs both DSSP bound and unbound
    :param contacts: boolean
    :param override: boolean
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
        if sifts:
            inputsifts = "{}{}{}.xml".format(config.db_root, config.db_sifts_xml, pdb_id)

            r = SIFTSreader(inputsifts)
            sifts_table = r.residues(add_regions=True, add_dbs=False)
            sifts_table = get_sifts_selected_from_table(sifts_table, chain=chain, res=res,
                                                        uniprot=uniprot_id, site=site)
        else:
            sifts_table = None

        # DSSP table
        if dssp:
            # (default) bound DSSP
            if bio:
                outputdssp = "{}{}{}_bio.dssp" \
                             "".format(config.db_root, config.db_dssp_generated, pdb_id)
            else:
                outputdssp = "{}{}{}.dssp" \
                             "".format(config.db_root, config.db_dssp_generated, pdb_id)
            if not outputdssp or override:
                w = DSSPgenerator(inputcif, outputdssp)
                w.run(override=override)
            if os.path.exists(outputdssp):
                r = DSSPreader(outputdssp)
                dssp_table = r.residues(add_ss_reduced=True, add_rsa_class=True)
                dssp_table = get_dssp_selected_from_table(dssp_table,
                                                          chain_full=chain, res=res)
            else:
                dssp_table = None

            if dssp_unbound and dssp_table is not None:
                # unbound DSSP
                if bio:
                    outputdssp_unb = "{}{}{}_bio_unbound.dssp" \
                                     "".format(config.db_root, config.db_dssp_generated, pdb_id)
                else:
                    outputdssp_unb = "{}{}{}_unbound.dssp" \
                                     "".format(config.db_root, config.db_dssp_generated, pdb_id)

                if not outputdssp_unb or override:
                    w = DSSPgenerator(inputcif, outputdssp_unb)
                    w.run(override=override, run_unbound=dssp_unbound)
                if os.path.exists(outputdssp_unb):
                    r = DSSPreader(outputdssp_unb)
                    dssp_unb_table = r.residues(add_ss_reduced=True, add_rsa_class=True,
                                                reset_res_id=True)
                    dssp_unb_table = get_dssp_selected_from_table(dssp_unb_table,
                                                                  chain_full=chain, res=res)
                    # merge bound and unbound
                    dssp_table = dssp_dssp_table_merger(dssp_table, dssp_unb_table)
        else:
            dssp_table = None

        if contacts:
            contacts_table = None
        else:
            contacts_table = None

        return mmcif_table, dssp_table, sifts_table, contacts_table

    else:
        raise ValueError('No UniProt ID or PDB ID provided...')


class TableMerger(object):
    def __init__(self, mmcif_table=None, dssp_table=None, sifts_table=None,
                 contacts_table=None, store=True, verbose=False):

        self.mmcif_table = mmcif_table
        self.dssp_table = dssp_table
        self.sifts_table = sifts_table
        self.contacts_table = contacts_table
        self.store = store
        self.merged_table = None
        self.filename = None
        self.verbose = verbose

    @staticmethod
    def _get_filename(uniprot_id=None, pdb_id=None, chain=None, res=None,
                      site=None, atom=None, lines=None,
                      bio=False, dssp=False, sifts=True, contacts=False):

        if uniprot_id is None and pdb_id is None:
            raise ValueError('No UniProt ID or PDB ID provided...')

        name = [uniprot_id, pdb_id, chain, res, site, atom, lines]
        name = '_'.join(['-'.join([v for v in k]) if type(k) is tuple else k
                         for k in name if k is not None])
        if sifts:
            name += "_sifts"
        if dssp:
            name += "_dssp"
        if bio:
            name += "_bio"
        if contacts:
            name += "_cont"

        filename = "{}{}{}.pkl".format(config.db_root, config.db_pickle, name)
        return filename

    def merge(self, outputfile=None, override=False):
        """Merges the provided tables and stores"""

        self.merged_table = table_merger(mmcif_table=self.mmcif_table,
                                         dssp_table=self.dssp_table,
                                         sifts_table=self.sifts_table,
                                         contacts_table=self.contacts_table)
        if self.store and outputfile is not None:
            self.filename = outputfile
            dump_merged_table(self.merged_table, self.filename,
                              override=override, verbose=self.verbose)
        return self.merged_table

    def run(self, uniprot_id=None, pdb_id=None, chain=None, res=None, site=None,
            atom=None, lines=None, bio=False, sifts=True, dssp=False, dssp_unbound=False,
            contacts=False, outputfile=None, override=False):
        """Generates the tables, merges and stores"""

        self.mmcif_table, self.dssp_table, self.sifts_table, self.contacts_table = \
            table_generator(uniprot_id=uniprot_id, pdb_id=pdb_id, chain=chain, res=res,
                            site=site, atom=atom, lines=lines, bio=bio, sifts=sifts,
                            dssp=dssp, dssp_unbound=dssp_unbound, contacts=contacts,
                            override=override)

        self.merged_table = self.merge(outputfile=outputfile, override=override)

        if self.store:
            if outputfile is not None:
                self.filename = outputfile
            else:
                self.filename = self._get_filename(uniprot_id=uniprot_id, pdb_id=pdb_id,
                                                   chain=chain, res=res, site=site,
                                                   atom=atom, lines=lines, bio=bio,
                                                   sifts=sifts, dssp=dssp, contacts=contacts)
            dump_merged_table(self.merged_table, self.filename,
                              override=override, verbose=self.verbose)

        return self.merged_table

    def load(self, uniprot_id=None, pdb_id=None, chain=None, res=None,
             site=None, atom=None, lines=None, bio=False, sifts=True,
             dssp=False, contacts=False, inputfile=None):
        """Loads a merged table"""

        if inputfile is not None:
            self.filename = inputfile
        else:
            self.filename = self._get_filename(uniprot_id=uniprot_id, pdb_id=pdb_id,
                                               chain=chain, res=res, site=site,
                                               atom=atom, lines=lines, bio=bio,
                                               sifts=sifts, dssp=dssp, contacts=contacts)
        self.merged_table = load_merged_table(self.filename, verbose=self.verbose)
        return self.merged_table


def dump_merged_table(data, outputfile, override=False, verbose=False):
    """
    Dumps the table to a file.

    :param data: pandas DataFrame object - merged table
    :param outputfile: path to the pickled file
    :param override: boolean
    :param verbose: boolean
    :return: (side-effects) writes a picked object to a file
    """
    if not os.path.isfile(outputfile) or override:
        data.to_pickle(outputfile)
        if verbose:
            flash('{} stored locally...'.format(outputfile))
    else:
        if verbose:
            flash('{} already available...'.format(outputfile))


def load_merged_table(inputfile, verbose=False):
    """
    Loads a table from a file.

    :param inputfile: path to the pickled file
    :param verbose: boolean
    :return: returns a picked object to a file
    """
    if os.path.isfile(inputfile):
        if verbose:
            flash('{} is already available...'.format(inputfile))
        merged_table = pd.read_pickle(inputfile)
        return merged_table
    else:
        raise IOError("{} not available or could not be read..."
                      "".format(inputfile))


if __name__ == '__main__':
    pass
