#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""

This is where methods that handle genetic variants.

Fábio Madeira, 2017+

"""

import logging
import pandas as pd

from prointvar.fetchers import fetch_uniprot_species_from_id
from prointvar.fetchers import fetch_ensembl_uniprot_ensembl_mapping
from prointvar.fetchers import fetch_ensembl_ensembl_uniprot_mapping
from prointvar.fetchers import InvalidEnsemblSpecies

from prointvar.library import uni_var_types

logger = logging.getLogger("prointvar")


def flatten_json(data, entries, ix=1, prefix=None):
    """
    Flattens a bunch of nested dictionaries

    :param data: "features" Dict
    :param entries: Dict
    :param ix: Count index
    :param prefix: name prefix
    :return: updates the 'entries' dictionary
    """
    for key, val in data.items():
        if type(val) not in (dict, list):
            if prefix is None:
                entries["%s_%s" % (ix, key)] = val
            else:
                entries["%s_%s_%s" % (ix, prefix, key)] = val
        else:
            if type(val) is dict:
                for j, k in enumerate(val.keys()):
                    nkey = "%s_%s_%s_%s" % (ix, prefix, key, k)
                    if type(val[k]) is dict:
                        flatten_json(val[k], entries, ix + j, prefix=nkey)
                    else:
                        entries[nkey] = val[k]

            elif type(val) is list:
                for l, v in enumerate(val):
                    flatten_json(v, entries, ix + l, prefix=key)


def collapse_unique_values(entries):
    """
    Collapses the Dictionary generated by flatten_json to
     unique columns, where columns with multiple times are
     aggregated to the same column, and the items grouped by order
     in a list

    :param entries: Dict
    :return: returns a column aggregated dictionary
        (similar to the 'unique' group_by aggregation in pandas)
    """
    unique_names = []
    for key in entries.keys():
        name = '_'.join(key.split('_')[1:])
        if name not in unique_names:
            unique_names.append(name)

    col_counts = {}
    for key in entries.keys():
        name = '_'.join(key.split('_')[1:])
        if name in unique_names:
            if name not in col_counts:
                col_counts[name] = 1
            else:
                col_counts[name] += 1

    new_entries = {}
    for key, val in entries.items():
        col = '_'.join(key.split('_')[1:])

        if col not in new_entries:
            if col_counts[col] > 1:
                new_entries[col] = [val]
            else:
                new_entries[col] = val
        else:
            new_entries[col].append(val)
    return new_entries


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
        entries = {"1_%s" % k: v for k, v in data.items() if k != "features"}

        flatten_json(entry, entries)
        var_rows.append(collapse_unique_values(entries))

    table = pd.DataFrame(var_rows)

    if excluded is not None:
        assert type(excluded) is tuple
        try:
            table = table.drop(list(excluded), axis=1)
        except ValueError:
            # most likely theses are not in there
            pass

    # enforce some specific column types
    for col in table:
        if col in uni_var_types:
            try:
                table[col] = table[col].astype(uni_var_types[col])
            except ValueError:
                # there are some NaNs in there
                pass
        else:
            print('missing', col)

    if table.empty:
        raise ValueError('Variants collapsing resulted in an empty DataFrame...')

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
                    uniprots.append(entry)
                else:
                    uniprots.append(entry['primary_id'])
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
    for ix, entry in enumerate(data):
        if ('ensembl_identity' in entry and 'xref_identity' in entry and
                'xref_start' in entry and 'xref_end' in entry):
            identity = entry['ensembl_identity'] + entry['xref_identity']
            coverage = entry['xref_end'] - entry['xref_start']
            if identity + coverage > prev_identity + prev_coverage:
                prev_identity = identity
                prev_coverage = coverage
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
        if uniprot:
            self.id_source = 'UniProt'
            self.uniprot_id = identifier
            self.ensembl_id = self._ensembl_id_from_uniprot()
        else:
            self.id_source = 'Ensembl'
            self.uniprot_id = self._uniprot_id_from_ensembl()
            self.ensembl_id = identifier

    def _ensembl_id_from_uniprot(self):

        info = fetch_uniprot_species_from_id(self.uniprot_id,
                                             cached=self.cached)
        species = get_ensembl_species_from_uniprot(info)
        try:
            info = fetch_ensembl_uniprot_ensembl_mapping(self.uniprot_id,
                                                         cached=self.cached,
                                                         species=species).json()
        except InvalidEnsemblSpecies:
            logger.info('Provided species {} is not valid'.format(species))
            return None
        ensps = get_ensembl_protein_id_from_mapping(info)
        best_match = get_preferred_ensembl_id_from_mapping(ensps, cached=self.cached,
                                                           uniprot_id=self.uniprot_id)
        return best_match

    def _uniprot_id_from_ensembl(self):

        info = fetch_ensembl_ensembl_uniprot_mapping(self.ensembl_id,
                                                     cached=self.cached)
        data = get_uniprot_id_from_mapping(info, full_entry=True)
        best_match = get_preferred_uniprot_id_from_mapping(data)
        return best_match


if __name__ == '__main__':
    pass
