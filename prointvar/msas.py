#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""

This defines the methods that work with Multiple sequence alignments.

Fábio Madeira, 2017+

"""

import os
import re
import copy
import json
import logging
import pandas as pd

from Bio import AlignIO

from prointvar.fetchers import fetch_uniprot_id_from_name
from prointvar.utils import constrain_column_types

logger = logging.getLogger("prointvar")


def read_alignment(inputfile, aln_format=None):
    """
    Reads an input multiple sequence alignment (MSA).
    Formats recognised by Biopython:
        "clustal", "emboss", "nexus", "fasta", "phylip" and "stockholm"

    :param inputfile: Input MSA (read by biopython)
    :param aln_format: (str) or None
    :return: returns the Biopython alignment object
    """

    if not os.path.isfile(inputfile):
        raise IOError("{} not available or could not be read...".format(inputfile))

    if aln_format is not None:
        align_format = aln_format
    elif inputfile.endswith('.fasta') or inputfile.endswith('.fa'):
        align_format = 'fasta'
    elif inputfile.endswith('.sto') or inputfile.endswith('.sth'):
        align_format = 'stockholm'
    elif inputfile.endswith('.aln') or inputfile.endswith('.clw'):
        align_format = 'clustal'
    else:
        raise ValueError("Alignment format unrecognised...")

    alignment = AlignIO.read(inputfile, align_format)
    return alignment


def parse_msa_sequences_from_file(inputfile, excluded=(), get_uniprot_id=False,
                                  cached=False):
    """
    Reads a Pfam/CATH MSA and returns a pandas table with a
    collection of protein IDs and sequences

    :param inputfile: path to the MSA file
    :param excluded: option to exclude some columns
    :param get_uniprot_id: (boolean)
    :param cached: (boolean)
    :return: returns a pandas DataFrame
    """

    rows = []
    alignment = read_alignment(inputfile)
    for record in alignment:
        seq = str(record.seq)
        desc = str(record.description)
        # get cross-reference information for each entry
        entry = {}
        entry['Sequence'] = seq
        entry = parse_sequence_info_from_description(desc, entry, cached=cached,
                                                     get_uniprot_id=get_uniprot_id)

        rows.append(entry)

    table = pd.DataFrame(rows)

    # excluding columns
    if excluded is not None:
        assert type(excluded) is tuple
        try:
            table = table.drop(list(excluded), axis=1)
        except ValueError:
            # most likely theses are not in there
            pass

    # enforce some specific column types
    msa_types = {key: str for key in list(table) if key != 'Start' and key != 'End'}
    table = constrain_column_types(table, msa_types)

    if table.empty:
        raise ValueError('{} resulted in an empty DataFrame...'.format(inputfile))

    return table


def parse_sequence_info_from_description(desc, entry, get_uniprot_id=False,
                                         cached=False):
    """
    Parses the Biopython alignment sequence description and tries to guess
    the content. (only works for known formats).

    :param desc: Biopython's 'description' field
    :param entry: Mutable dictionary
    :param get_uniprot_id: (boolean)
    :param cached: (boolean)
    :return: (updated) Dictionary
    """

    prev_entry = copy.deepcopy(entry)
    # trying the UniProt fasta seq description
    parse_uniprot_fasta_seq_description(desc, entry)
    if entry != prev_entry:
        return entry

    # trying the Pfam Stockholm seq description
    parse_pfam_sth_seq_description(desc, entry, cached=cached,
                                   get_uniprot_id=get_uniprot_id)
    if entry != prev_entry:
        return entry

    # trying the CATH fasta seq description
    parse_cath_fasta_seq_description(desc, entry, cached=cached,
                                     get_uniprot_id=get_uniprot_id)

    if entry != prev_entry:
        return entry

    # trying a generic sequence description
    parse_generic_seq_description(desc, entry, cached=cached,
                                  get_uniprot_id=get_uniprot_id)

    if entry != prev_entry:
        return entry

    logger.debug("Nothing parsed from the MSA sequence description...")
    return entry


def parse_uniprot_fasta_seq_description(desc, entry):
    """
    Pattern: <source>|<Accession_ID>|<Accession_Name> ++
    Example: sp|P00439|PH4H_HUMAN Phenylalanine-4-hydroxylase (...)
        OS=Homo sapiens GN=PAH PE=1 SV=1

    :param desc: Biopython's 'description' field
    :param entry: Mutable dictionary
    :return: (updated) Dictionary
    """

    # trying the UniProt fasta seq description
    pattern = re.compile("([a-zA-Z])+\|([A-Z0-9])+\|([A-Z0-9])+_([A-Z0-9])+")
    match = re.search(pattern, desc, flags=0)
    if match:
        match = match.group()
        # pattern matching
        pattern = re.compile("([a-zA-Z])+\|")
        source = re.search(pattern, match, flags=0)
        if source:
            source = source.group().rstrip('|')
            entry['Collection'] = source

        pattern = re.compile("\|([a-zA-Z0-9])+\|")
        identifier = re.search(pattern, match, flags=0)
        if identifier:
            identifier = identifier.group().lstrip('|').rstrip('|')
            entry['Accession'] = identifier

        pattern = re.compile("\|([A-Z0-9])+_([A-Z0-9])+")
        name = re.search(pattern, match, flags=0)
        if name:
            name = name.group().lstrip('|')
            entry['Name'] = name

        entry['Source'] = 'UniProt'

        # remaining description
        if desc != match:
            desc = desc.replace(match, "")
            entry['Description'] = desc.strip()

    return entry


def parse_pfam_sth_seq_description(desc, entry, get_uniprot_id=False,
                                   cached=False):
    """
    Pattern: <Accession_Name>/<Start>-<End>
    Example: C7P4T5_HALMD/44-372

    :param desc: Biopython's 'description' field
    :param entry: Mutable dictionary
    :param get_uniprot_id: (boolean)
    :param cached: (boolean)
    :return: (updated) Dictionary
    """

    pattern = re.compile("([A-Z0-9])+_([A-Z0-9])+\/[\-]?([0-9])+-[\-]?([0-9])+")
    match = re.search(pattern, desc, flags=0)
    if match:
        match = match.group()
        entry = parse_generic_seq_description(match, entry, cached=cached,
                                              get_uniprot_id=get_uniprot_id)

        entry['Source'] = 'Pfam'

        # remaining description
        if desc != match:
            desc = desc.replace(match, "")
            entry['Description'] = desc.strip()

    return entry


def parse_cath_fasta_seq_description(desc, entry, get_uniprot_id=False,
                                     cached=False):
    """
    Pattern: cath|<version>|<Accession_ID>/<Start>-<End>
    Example:
    # sequence from CATH alignments - with structure domain
    # cath|4.1.0|1rwcA01/4-372 (...)
        CATH_S35=1.50.10.100.1;GENE=P84141_ARTAU;GO=GO:0005576,GO:0005975;
        MDA=1.50.10.100;ORG=Arthrobacter_aurescens;
        TAXON=43663;UNIPROT=P84141

    # sequence from CATH alignments - without structure domain
    # biomap|4.1.0|b7f2808bde19a485bdd0a90545c40d99/29-337 (...)
        EC=4.2.2.5;GENE=CSLA_PEDHD;GO=GO:0005576,GO:0005975;
        MDA=1.50.10.100_2.70.98.10_2.60.220.10;
        ORG=Pedobacter_heparinus;SWISSPROT=1;
        TAXON=485917;UNIPROT=Q59288

    :param desc: Biopython's 'description' field
    :param entry: Mutable dictionary
    :param get_uniprot_id: (boolean)
    :param cached: (boolean)
    :return: (updated) Dictionary
    """
    # trying the CATH fasta seq description
    pattern = re.compile("([a-zA-Z])+\|([0-9])(.|-)([0-9])(.|-)([0-9])\|"
                         "([a-zA-Z0-9])+\/[\-]?([0-9])+-[\-]?([0-9])+")
    match = re.search(pattern, desc, flags=0)
    if match:
        match = match.group()
        # pattern matching
        pattern = re.compile("([a-zA-Z])+\|")
        source = re.search(pattern, match, flags=0)
        if source:
            source = source.group().rstrip('|')
            entry['Collection'] = source

        pattern = re.compile("\|([0-9])(.|-)([0-9])(.|-)([0-9])\|")
        version = re.search(pattern, match, flags=0)
        if version:
            version = version.group().lstrip('|').rstrip('|')
            entry['Version'] = version

        pattern = re.compile("([A-Z0-9])+[\_]?([a-zA-Z0-9])+\/[\-]?([0-9])+-[\-]?([0-9])+")
        # same as the generic seq description
        nmatch = re.search(pattern, match, flags=0)
        if nmatch:
            entry = parse_generic_seq_description(nmatch.group(),
                                                  entry, cached=cached,
                                                  get_uniprot_id=get_uniprot_id)

        entry['Source'] = 'CATH'

        # remaining description
        if desc != match:
            desc = desc.replace(match, "")
            entry['Description'] = desc.strip()

    return entry


def parse_generic_seq_description(desc, entry, get_uniprot_id=False,
                                  cached=False):
    """
    Pattern: <Accession_ID>/<Start>-<End> or <Accession_Name>/<Start>-<End>
    Example: P00439/24-145

    :param desc: Biopython's 'description' field
    :param entry: Mutable dictionary
    :param get_uniprot_id: (boolean)
    :param cached: (boolean)
    :return: (updated) Dictionary
    """

    # trying a generic sequence description
    pattern = re.compile("([A-Z0-9])+[\_]?([a-zA-Z0-9])+\/[\-]?([0-9])+-[\-]?([0-9])+")
    match = re.search(pattern, desc, flags=0)
    if match:
        match = match.group()
        # pattern matching
        pattern = re.compile("([A-Z0-9])+_([a-zA-Z0-9])+\/")
        name = re.search(pattern, match, flags=0)
        if name:
            name = name.group().rstrip('/')
            entry['Name'] = name

        pattern = re.compile("([a-zA-Z0-9])+\/")
        identifier = re.search(pattern, match, flags=0)
        if identifier:
            identifier = identifier.group().rstrip('/')
            entry['Accession'] = identifier

        pattern = re.compile("\/[\-]?([0-9])+")
        start = re.search(pattern, match, flags=0)
        if start:
            start = int(start.group().lstrip('/'))
            entry['Start'] = start

        pattern = re.compile("-[\-]?([0-9])+")
        end = re.search(pattern, match, flags=0)
        if end:
            end = int((end.group())[1:])
            entry['End'] = end

        # UniProt ID missing
        if name and get_uniprot_id:
            r = fetch_uniprot_id_from_name(name, cached=cached)
            try:
                identifier = r.json()[0]["id"]
            except Exception:
                # JSONDecodeError
                # new - to match the update in the UniProt API endpoint
                identifier = str(r.content, encoding='utf-8').strip()
            entry['Accession'] = identifier

        entry['Source'] = 'GenericParser'

        # remaining description
        if desc != match:
            desc = desc.replace(match, "")
            entry['Description'] = desc.strip()

    return entry


class MSAreader(object):
    def __init__(self, inputfile):
        """
        :param inputfile: Needs to point to a valid MSA file.
        """
        self.inputfile = inputfile
        self.data = None
        self.excluded = ()

        if not os.path.isfile(inputfile):
            raise IOError("{} not available or could not be read...".format(inputfile))

    def read(self, **kwargs):
        return self.msas(**kwargs)

    def msas(self, excluded=None, get_uniprot_id=False, cached=False):
        if excluded is None:
            excluded = self.excluded
        self.data = parse_msa_sequences_from_file(self.inputfile, excluded=excluded,
                                                  get_uniprot_id=get_uniprot_id,
                                                  cached=cached)
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
            logger.info("No MSA data parsed...")


if __name__ == '__main__':
    pass
