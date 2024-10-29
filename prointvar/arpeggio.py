#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This defines the methods that work with ARPEGGIO files.

Propositions for things to work:
- Generate new PDB file from mmCIF and (even) PDB so that some issues (below) are cleared
    - Remove alternative locations and reset the atom sequential id
      PDBXreader(<...>).atoms(remove_altloc=True, reset_atom_id=True)

    - Chain IDs need to be a single character so it is safer to run with category='auth'
      PDBXwriter(<...>).run(category='auth')

    - Remove hydrogens
      PDBXreader(<...>).atoms(remove_hydrogens=True)

Fábio Madeira, 2017+
Stuart MacGowan, 2024+

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
from prointvar.library import arpeggio_types
from prointvar.library import arpeggio_col_renames

from prointvar.config import config

logger = logging.getLogger("prointvar")

def parse_arpeggio_from_file(inputfile, parse_special=False):
    """
    Parse lines of the ARPEGGIO *json* file to get entries from...

    :param inputfile: path to the ARPEGGIO file
    :param parse_special: (boolean) tries to parse special contact types
    :return: returns a pandas DataFrame
    """
    
    def swap_bgn_end_in_columns(df):
        new_columns = {}
        for col in df.columns:
            if col.startswith('bgn.'):
                new_columns[col] = col.replace('bgn.', 'end.', 1)
            elif col.startswith('end.'):
                new_columns[col] = col.replace('end.', 'bgn.', 1)
        return df.rename(columns=new_columns)

    logger.info("Parsing ARPEGGIO from lines...")

    if not os.path.isfile(inputfile):
        raise IOError(f"{inputfile} not available or could not be read...")

    with open(inputfile, 'r') as f:
        f.seek(0)
        data = json.load(f)
        table = pd.json_normalize(data)
        
    # Process special contact types 'atom-plane', 'plane-plane' and 'group-group'
    special_table = table[table['type'] != 'atom-atom']
    if not special_table.empty:
        if parse_special:
            logger.info("Retaining special contact-types...")
            # table = table[table['type'] != 'atom-atom']
            # table = pd.concat([table, special_table])
        else:
            logger.info("Dropping special contact-types...")
            table = table[table['type'] == 'atom-atom']
    
    # Reorder contact direction to enforce bgn < end
    table_A = table.loc[table['bgn.auth_seq_id'] <= table['end.auth_seq_id']]
    table_B = table.loc[table['bgn.auth_seq_id'] > table['end.auth_seq_id']]
    table_B = swap_bgn_end_in_columns(table_B)
    table = pd.concat([table_A, table_B])
    
    # Format residues with no insertion code
    table.loc[:, "bgn.pdbx_PDB_ins_code"] = table.loc[:, "bgn.pdbx_PDB_ins_code"].str.replace(' ', '?')
    table.loc[:, "end.pdbx_PDB_ins_code"] = table.loc[:, "end.pdbx_PDB_ins_code"].str.replace(' ', '?')

    # Enforce some specific column types
    table = constrain_column_types(table, arpeggio_types)

    if table.empty:
        raise ValueError(f'{inputfile} resulted in an empty DataFrame...')

    return table


def get_arpeggio_selected_from_table(data, chain_A=None, chain_B=None,
                                     res_A=None, res_B=None,
                                     atom_A=None, atom_B=None):
    """
    Utility that filters a pandas DataFrame by the input tuples.

    :param data: pandas DataFrame object
    :param chain_A: (tuple) chain IDs or None (donor)
    :param chain_B: (tuple) chain IDs or None (acceptor)
    :param res_A: (tuple) res IDs or None (donor)
    :param res_B: (tuple) res IDs or None (acceptor)
    :param atom_A: (tuple) atom IDs or None (donor)
    :param atom_B: (tuple) atom IDs or None (acceptor)
    :return: returns a modified pandas DataFrame
    """

    # excluding rows
    table = data

    if chain_A is not None:
        table = row_selector(table, 'bgn.auth_asym_id', chain_A)
        logger.info("Arpeggio table filtered by CHAIN_A...")

    if chain_B is not None:
        table = row_selector(table, 'end.auth_asym_id', chain_B)
        logger.info("Arpeggio table filtered by CHAIN_B...")

    if res_A is not None:
        table = row_selector(table, 'bgn.auth_seq_id', res_A)
        logger.info("Arpeggio table filtered by RES_A...")

    if res_B is not None:
        table = row_selector(table, 'end.auth_seq_id', res_B)
        logger.info("Arpeggio table filtered by RES_B...")

    if atom_A is not None:
        table = row_selector(table, 'bgn.auth_atom_id', atom_A)
        logger.info("Arpeggio table filtered by ATOM_A...")

    if atom_B is not None:
        table = row_selector(table, 'end.auth_atom_id', atom_B)
        logger.info("Arpeggio table filtered by ATOM_B...")

    return table


def residues_aggregation(data, agg_method='unique'):
    """
    Gets the contacts res-by-res, instead of atom-atom.

    :param data: pandas DataFrame object
    :param agg_method: current values: 'first', 'unique', and 'minimum'
    :return: returns a modified pandas DataFrame
    """
    
    def contact_aggregator():
        return lambda all_atom_contacts: set(contact for atom_contacts in all_atom_contacts for contact in atom_contacts)

    table = data
    agg_generic = agg_method
    agg_method_origin = agg_method
    agg_cols = ['bgn.auth_asym_id', 'bgn.auth_seq_id', 'end.auth_asym_id', 'end.auth_seq_id']
    if agg_method not in ['first', 'unique', 'minimum', 'maximum']:
        raise ValueError('Method {} is not currently implemented...'
                         ''.format(agg_method))

    if agg_method != 'minimum' and agg_method != 'maximum':
        columns_to_agg = {col: agg_generic if table[col].dtype == 'object' else agg_method
                          for col in table.columns if col not in agg_cols}
    else:
        if agg_method_origin == 'minimum':
            # need the table sort by distance first: ascending
            table = table.sort_values("distance", ascending=True)
            table.reset_index(inplace=True)
            table = table.drop(['index'], axis=1)
            agg_generic = 'first'
            agg_method = 'max'
            columns_to_agg = {col: agg_generic if table[col].dtype == 'object' else agg_method
                              for col in table.columns if col not in agg_cols}
            columns_to_agg['distance'] = 'min'
            columns_to_agg['contact'] = contact_aggregator()
        elif agg_method_origin == 'maximum':
            # need the table sort by distance first: descending
            table = table.sort_values("distance", ascending=False)
            table.reset_index(inplace=True)
            table = table.drop(['index'], axis=1)
            agg_generic = 'first'
            agg_method = 'max'
            columns_to_agg = {col: agg_generic if table[col].dtype == 'object' else agg_method
                              for col in table.columns if col not in agg_cols}
            # if contacts columns have been collapsed
            columns_to_agg['contact'] = contact_aggregator()
    table = table.groupby(by=agg_cols, as_index=False).agg(columns_to_agg)
    
    # TODO: Drop atom columns after residue aggregation? Or implement suitable aggregation? Approach has implications for merging.
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
        table = table.loc[table['bgn.auth_asym_id'] != table['end.auth_asym_id']]
    elif int_mode == 'intra-chain':
        table = table.loc[table['bgn.auth_asym_id'] == table['end.auth_asym_id']]
    else:
        raise ValueError('Interaction mode {} is not currently implemented...'
                         ''.format(int_mode))
    # FIXME optional?
    table.reset_index(inplace=True)
    table = table.drop(['index'], axis=1)
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
    ins_codes_1 = [k for k in table["bgn.pdbx_PDB_ins_code"].unique()]
    ins_codes_2 = [k for k in table["end.pdbx_PDB_ins_code"].unique()]
    if (len(ins_codes_1) and len(ins_codes_2) and
            ins_codes_1[0] == '?' and ins_codes_2[0] == '?'):
        table = table.loc[((table["bgn.auth_asym_id"] != table["end.auth_asym_id"]) |
                           ((table["bgn.auth_asym_id"] == table["end.auth_asym_id"]) &
                            (abs(table["bgn.auth_seq_id"].astype(int) -
                                 table["end.auth_seq_id"].astype(int)) >= numb_res)))]
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

        if not os.path.isfile(inputfile):
            raise IOError(f"{inputfile} not available or could not be read...")

    def read(self, **kwargs):
        return self.contacts(**kwargs)

    def contacts(self, residue_agg=False, agg_method='minimum',
                 int_filter=False, int_mode='inter-chain',
                 ignore_consecutive=False, numb_res=3,
                 parse_special=False):

        self.data = parse_arpeggio_from_file(self.inputfile, parse_special=parse_special)
        if ignore_consecutive:
            self.data = ignore_consecutive_residues(self.data, numb_res=numb_res)

        if residue_agg:
            self.data = residues_aggregation(self.data, agg_method=agg_method)

        if int_filter:
            self.data = interaction_modes(self.data, int_mode=int_mode)

        return self.data

    def to_json(self, pretty=True):
        if self.data is not None:
            if isinstance(self.data, pd.DataFrame):
                data = self.data.to_dict(orient='records')
            else:
                data = self.data
            return json.dumps(data, sort_keys=False, indent=4) if pretty else json.dumps(data)
        else:
            logger.info("No ARPEGGIO data parsed...")


class ARPEGGIOrunner(object):
    def __init__(self, inputfile, outputfile=None):
        """
        :param inputfile: Needs to point to a valid PDB or mmCIF file.
        :param outputfile: if not provided will use the same file name and <.contacts> extension
        """
        self.inputfile = inputfile
        self.inputfile_back = inputfile
        self.outputfile = outputfile
        self.data = None
        self.inputfile_h = None

        if not os.path.isfile(self.inputfile):
            raise IOError(f"{self.inputfile} not available or could not be read...")

        # inputfile needs to be in PDB or mmCIF format
        filename, extension = os.path.splitext(self.inputfile)
        if extension not in ['.pdb', '.ent', '.cif']:
            raise ValueError(f"{self.inputfile} is expected to be in mmCIF or PDB format...")

    def _generate_output(self):
        filename, extension = os.path.splitext(self.inputfile)
        self.outputfile = filename + ".json"

    def _generate_mmcif(self, override=False, pro_format=False):
        filename, extension = os.path.splitext(self.inputfile)
        self.inputfile = filename + "_new.cif"
        w = PDBXwriter(inputfile=None, outputfile=self.inputfile)
        r = PDBXreader(inputfile=self.inputfile_back)
        data = r.atoms(remove_altloc=True, reset_atom_id=True, add_new_pro_id=True,
                       remove_partial_res=True, format_type=None)
        w.run(data=data, format_type="mmcif", category="auth",
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

    def _run(self, arpeggio_bin, clean_output=True,
             hydro_method="arpeggio"):

        filename, extension = os.path.splitext(self.inputfile)
        input_arpeggio = filename + ".cif"
        output_path = os.path.dirname(input_arpeggio)
        output_arpeggio = filename + ".json"  # atom-atom contact information
        output_hydro = filename + "_hydrogenated.mmcif"

        if hydro_method in ["hbplus", "reduce"]:
            input_arpeggio = filename + ".h.cif"

        # run arpeggio
        cmd = f'{arpeggio_bin} -wh {input_arpeggio}'
        if output_path:
            cmd = f'{cmd} --output {output_path}'
        os.system(cmd)
        if not os.path.isfile(output_arpeggio):
            raise IOError(f"ARPEGGIO output not generated for {input_arpeggio}")

        # mv the automatically generated file -> to the provided outputfile
        if output_arpeggio != self.outputfile:
            shutil.copyfile(output_arpeggio, self.outputfile)

        if hydro_method == "arpeggio":
            shutil.copyfile(output_hydro, self.inputfile_h)

        if clean_output:
            # remove output files
            if output_arpeggio != self.outputfile:
                lazy_file_remover(output_arpeggio)

    def write(self, **kwargs):
        return self.run(**kwargs)

    def run(self, hydro_method="arpeggio", override=False,
            clean_output=True, save_new_input=False, pro_format=False):

        # generate outputfile if missing
        if not self.outputfile:
            self._generate_output()

        if not os.path.exists(self.outputfile) or override:
            if os.path.exists(config.arpeggio_bin):
                arpeggio_bin = config.arpeggio_bin
            else:
                raise IOError('ARPEGGIO executable is not available...')

            # get PDB with explicit hydrogen atoms
            filename, extension = os.path.splitext(self.inputfile)
            self.inputfile_h = filename + ".h.pdb"
            if hydro_method in ["hbplus", "reduce"]:
                self._generate_pdb_with_hydrogens(hydro_method=hydro_method,
                                                  override=override)
                
            # input file needs to be in mmCIF format
            if extension == '.pdb' or extension == '.ent':
                # FIXME: causes ValueError: Missing _chem_comp. category in mmcif in ARPEGGIO
                self._generate_mmcif(override=override, pro_format=pro_format)

            # run arpeggio and generate output - also clean unnecessary output
            self._run(arpeggio_bin, clean_output=clean_output,
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
