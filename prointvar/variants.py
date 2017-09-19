#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""

This is where methods that handle genetic variants.

Fábio Madeira, 2017+

"""

import logging
import pandas as pd

from prointvar.fetchers import fetch_uniprot_variants_ebi
from prointvar.fetchers import fetch_uniprot_species_from_id
from prointvar.fetchers import fetch_ensembl_uniprot_ensembl_mapping
from prointvar.fetchers import fetch_ensembl_ensembl_uniprot_mapping
from prointvar.fetchers import fetch_ensembl_transcript_variants
from prointvar.fetchers import fetch_ensembl_somatic_variants
from prointvar.fetchers import InvalidEnsemblSpecies

from prointvar.merger import uniprot_vars_ensembl_vars_merger

from prointvar.utils import row_selector
from prointvar.utils import splitting_up_by_key
from prointvar.utils import merging_down_by_key
from prointvar.utils import flatten_nested_structure
from prointvar.utils import refactor_key_val_singletons
from prointvar.utils import constrain_column_types
from prointvar.utils import exclude_columns
from prointvar.library import uni_ens_var_types
from prointvar.library import update_ensembl_to_uniprot

logger = logging.getLogger("prointvar")


def flatten_uniprot_variants_ebi(data, excluded=()):
    """
    Flattens the json output obtained from the Proteins API variants
     endpoint.

    :param data: original response (json output)
    :param excluded: option to exclude VAR columns
    :return: returns a pandas DataFrame
    """

    try:
        data = data.json()
    except AttributeError:
        assert type(data) is dict

    var_rows = []
    for entry in data["features"]:
        entries = {key: val for key, val in data.items() if key != 'features'}

        flatten_nested_structure(entry, entries)
        var_rows.append(refactor_key_val_singletons(entries))

    table = pd.DataFrame(var_rows)

    # excluding columns
    table = exclude_columns(table, excluded=excluded)

    # enforce some specific column types
    table = constrain_column_types(table, uni_ens_var_types)

    # split multi id rows
    table = splitting_up_by_key(table, key='xrefs_id')

    # merge down multi rows with same id
    table = merging_down_by_key(table, key='xrefs_id')

    if table.empty:
        raise ValueError('Variants collapsing resulted in an empty DataFrame...')

    return table


def flatten_ensembl_variants(data, excluded=(), synonymous=True):
    """
    Flattens the json output obtained from the Proteins API variants
     endpoint.

    :param data: original response (json output)
    :param excluded: option to exclude VAR columns
    :param synonymous: (boolean)
    :return: returns a pandas DataFrame
    """

    try:
        data = data.json()
    except AttributeError:
        assert type(data) is dict

    table = pd.DataFrame(data)
    # rename columns
    table.rename(columns=update_ensembl_to_uniprot, inplace=True)

    # excluding columns
    table = exclude_columns(table, excluded=excluded)

    # enforce some specific column types
    table = constrain_column_types(table, uni_ens_var_types)

    # split multi id rows
    table = splitting_up_by_key(table, key='xrefs_id')

    # merge down multi rows with same id
    table = merging_down_by_key(table, key='xrefs_id')

    # filter synonymous
    if not synonymous:
        table = row_selector(table, key='consequenceType',
                             value='synonymous_variant', method="diffs")
    return table


def get_ensembl_protein_id_from_mapping(data):
    """
    Gets a list of Ensembl IDs from a 'xrefs/symbol/' mapping.

    :param data: Requests object from the Ensembl-UniProt Mapping
    :return: list of Ensembl Protein IDs
    """
    ensps = []
    for entry in data:
        if 'type' in entry and 'id' in entry:
            if entry['type'] == 'translation':
                if entry['id'] not in ensps:
                    ensps.append(entry['id'])
    return ensps


def get_uniprot_id_from_mapping(data, full_entry=False, uniprot_id=None):
    """
    Gets a list of UniProt IDs from a '"xrefs/id/"' mapping.

    :param data: Requests object from the Ensembl-UniProt Mapping
    :param full_entry: (boolean) if True gets dictionary instead of just
        the UniProt IDs
    :param uniprot_id: if not None means that we want the data for a specific
        UniProt ID
    :return: list of UniProt IDs
    """
    uniprots = []
    for entry in data:
        if 'dbname' in entry and 'primary_id' in entry:
            if uniprot_id is not None and entry['primary_id'] == uniprot_id:
                if full_entry:
                    uniprots = [entry]
                else:
                    uniprots = [entry['primary_id']]
                break
            elif entry['dbname'] == 'Uniprot/SWISSPROT':
                if entry['primary_id'] not in uniprots:
                    if full_entry:
                        uniprots.append(entry)
                    else:
                        uniprots.append(entry['primary_id'])
            elif entry['dbname'] == 'Uniprot/SPTREMBL':
                if entry['primary_id'] not in uniprots:
                    if full_entry:
                        uniprots.append(entry)
                    else:
                        uniprots.append(entry['primary_id'])
    return uniprots


def get_preferred_uniprot_id_from_mapping(data):
    """
    Takes a list of Ensembl xrefs/ids mapped from a UniProt ID
    and gets the preferred entry (Many-to-one), based on seq
    identity and coverage.

    :param data: list of dictionaries
    :return: (str) preferred UniProt ID
    """

    best_match = None
    curr_ix = -1
    prev_identity = 0
    prev_coverage = 0
    prev_id = "-" * 100
    for ix, entry in enumerate(data):
        if ('ensembl_identity' in entry and 'xref_identity' in entry and
                'xref_start' in entry and 'xref_end' in entry):
            identity = entry['ensembl_identity'] + entry['xref_identity']
            coverage = entry['xref_end'] - entry['xref_start']
            if identity + coverage >= prev_identity + prev_coverage:
                prev_identity = identity
                prev_coverage = coverage
                # preferring the smallest UniProt ID (for getting variants)
                if len(entry['primary_id']) < len(prev_id):
                    prev_id = entry['primary_id']
                    curr_ix = ix
    if curr_ix != -1 and 'primary_id' in data[curr_ix]:
        best_match = data[curr_ix]['primary_id']
    return best_match


def get_preferred_ensembl_id_from_mapping(identifiers, cached=False,
                                          uniprot_id=None):
    """
    Takes a list of Ensembl xrefs/ids mapped from a UniProt ID
    and gets the preferred entry (Many-to-one), based on seq
    identity and coverage.

    :param identifiers: list of Ensembl IDs
    :param cached: (boolean) if True, stores a pickle file locally
    :param uniprot_id: if not None means that we want the data for a specific
        UniProt ID
    :return: (str) preferred Ensembl ID
    """

    best_match = None
    curr_ix = -1
    prev_identity = 0
    prev_coverage = 0
    for ix, ensp in enumerate(identifiers):
        info = fetch_ensembl_ensembl_uniprot_mapping(ensp,
                                                     cached=cached).json()
        # gets the mapping for a specific uniprot
        data = get_uniprot_id_from_mapping(info, full_entry=True,
                                           uniprot_id=uniprot_id)
        for entry in data:
            if ('ensembl_identity' in entry and 'xref_identity' in entry and
                    'xref_start' in entry and 'xref_end' in entry):

                identity = entry['ensembl_identity'] + entry['xref_identity']
                coverage = entry['xref_end'] - entry['xref_start']
                if identity + coverage > prev_identity + prev_coverage:
                    prev_identity = identity
                    prev_coverage = coverage
                    curr_ix = ix
    if curr_ix != -1:
        best_match = identifiers[curr_ix]
    return best_match


def get_ensembl_species_from_uniprot(data):
    """
    Gets a Species Name from a UniProt organism lookup.

    :param data: Request object from the UniProt Query endpoint.
    :return: (str) formatted species name
    """
    organism = str(data.content, encoding='utf-8').split('\n')[1]
    species = '_'.join(organism.split()[0:2]).lower()
    return species


class VariantsAgreggator(object):
    def __init__(self, identifier, uniprot=True, cached=False):
        """
        Aggregates Variants from UniProt Proteins API and Ensembl REST API.

        :param identifier: UniProt or Ensembl Protein ID
        :param uniprot: (boolean) if True assumes the inputted ID is from UniProt
        :param cached: (boolean) passed to the fetcher methods
        """

        self.cached = cached
        self.data = None
        if uniprot:
            self.id_source = 'UniProt'
            self.uniprot_id = identifier
            self.ensembl_id = self._ensembl_id_from_uniprot()
        else:
            self.id_source = 'Ensembl'
            self.ensembl_id = identifier
            self.uniprot_id = self._uniprot_id_from_ensembl()

    def _get_uniprot_species(self):
        info = fetch_uniprot_species_from_id(self.uniprot_id,
                                             cached=self.cached)
        self.species = get_ensembl_species_from_uniprot(info)

    def _ensembl_id_from_uniprot(self):

        self._get_uniprot_species()
        try:
            info = fetch_ensembl_uniprot_ensembl_mapping(self.uniprot_id,
                                                         cached=self.cached,
                                                         species=self.species).json()
        except InvalidEnsemblSpecies:
            logger.info('Provided species {} is not valid'.format(self.species))
            return None
        ensps = get_ensembl_protein_id_from_mapping(info)
        best_match = get_preferred_ensembl_id_from_mapping(ensps, cached=self.cached,
                                                           uniprot_id=self.uniprot_id)
        return best_match

    def _uniprot_id_from_ensembl(self):

        info = fetch_ensembl_ensembl_uniprot_mapping(self.ensembl_id,
                                                     cached=self.cached).json()
        data = get_uniprot_id_from_mapping(info, full_entry=True)
        best_match = get_preferred_uniprot_id_from_mapping(data)
        return best_match

    def run(self, synonymous=True, uniprot_vars=True,
            ensembl_transcript_vars=True, ensembl_somatic_vars=True):

        uni_vars = None
        trans_vars = None
        som_vars = None
        ens_vars = None

        if self.uniprot_id is not None and uniprot_vars:
            r = fetch_uniprot_variants_ebi(self.uniprot_id, cached=self.cached)
            if r is not None:
                uni_vars = flatten_uniprot_variants_ebi(r)

        if (self.ensembl_id is not None and
                (ensembl_transcript_vars or ensembl_somatic_vars)):

            if ensembl_transcript_vars:
                r = fetch_ensembl_transcript_variants(self.ensembl_id,
                                                      cached=self.cached)
                if r is not None:
                    trans_vars = flatten_ensembl_variants(r, synonymous=synonymous)

            if ensembl_somatic_vars:
                r = fetch_ensembl_somatic_variants(self.ensembl_id,
                                                   cached=self.cached)
                if r is not None:
                    som_vars = flatten_ensembl_variants(r, synonymous=synonymous)

            if isinstance(trans_vars, pd.DataFrame) and isinstance(som_vars, pd.DataFrame):
                ens_vars = pd.concat([trans_vars, som_vars]).reset_index(drop=True)
            elif isinstance(trans_vars, pd.DataFrame):
                ens_vars = trans_vars
            elif isinstance(som_vars, pd.DataFrame):
                ens_vars = som_vars

        if isinstance(uni_vars, pd.DataFrame) and isinstance(ens_vars, pd.DataFrame):
            self.data = uniprot_vars_ensembl_vars_merger(uni_vars, ens_vars)
        elif isinstance(uni_vars, pd.DataFrame):
            self.data = uni_vars
        elif isinstance(ens_vars, pd.DataFrame):
            self.data = ens_vars
        else:
            logger.info('No variants found...')
            self.data = pd.DataFrame()

        return self.data


if __name__ == '__main__':
    pass
