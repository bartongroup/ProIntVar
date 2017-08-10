#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""

This defines the methods that work with Multiple sequence alignments.

Fábio Madeira, 2017+

"""

import re
import logging
import pandas as pd

from Bio import AlignIO

from prointvar.fetchers import fetch_uniprot_id_from_name

logger = logging.getLogger("prointvar")


def read_alignment(inputfile, aln_format=None):
    """
    Reads an input multiple sequence alignment (MSA).
    Formats recognised by Biopython:
        "clustal", "emboss", "nexus", "fasta", "phylip" and "stockholm"

    :param inputfile: Input MSA (read by biopython)
    :param aln_format: (str) or None
    :return: returns the biopython alignment object
    """

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


def parse_msa_sequences_from_file(inputfile, excluded=()):
    """
    Reads a Pfam/CATH MSA and returns a pandas table with a
    collection of protein IDs and sequences

    :param inputfile: path to the MSA file
    :param excluded: option to exclude some columns
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
        entry = parse_sequence_info_from_description(desc, entry)
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
    msa_types = {key: str for key in list(table)}
    for col in table:
        try:
            table[col] = table[col].astype(msa_types[col])
        except ValueError:
            # there are some NaNs in there
            pass

    if table.empty:
        raise ValueError('{} resulted in an empty DataFrame...'.format(inputfile))

    return table


def parse_sequence_info_from_description(desc, entry):
    """
    Parses the Biopython alignment sequence description and tries to guess
    the content. (only works for known formats).

    :param desc: Biopython's 'description' field
    :param entry: Mutable dictionary
    :return: (updated) Dictionary
    """

    # trying the UniProt fasta seq description
    pattern = re.compile("([a-zA-Z])+\|([A-Z0-9])+\|([A-Z0-9])+_([A-Z0-9])+")
    uniprot_fasta = re.search(pattern, desc, flags=0)

    # trying the Pfam Stockholm seq description
    pattern = re.compile("([A-Z0-9])+_([A-Z0-9])+\/[\-]?([0-9])+-[\-]?([0-9])+")
    pfam_sth = re.search(pattern, desc, flags=0)

    # trying the CATH fasta seq description
    pattern = re.compile("([a-zA-Z])+\|([0-9])(.|-)([0-9])(.|-)([0-9])\|"
                         "([a-zA-Z0-9])+\/[\-]?([0-9])+-[\-]?([0-9])+")
    # cath_fasta = re.match(pattern, desc, flags=0)
    cath_fasta = re.search(pattern, desc, flags=0)

    # trying a generic sequence description
    pattern = re.compile("([A-Z0-9])+[\_]?([a-zA-Z0-9])+\/[\-]?([0-9])+-[\-]?([0-9])+")
    gen_fasta = re.search(pattern, desc, flags=0)

    if uniprot_fasta:
        return parse_uniprot_fasta_seq_description(uniprot_fasta.group(), desc, entry)
    elif pfam_sth:
        return parse_pfam_sth_seq_description(pfam_sth.group(), desc, entry)
    elif cath_fasta:
        return parse_cath_fasta_seq_description(cath_fasta.group(), desc, entry)
    elif gen_fasta:
        return parse_generic_seq_description(gen_fasta.group(), desc, entry)
    else:
        logger.debug("Nothing parsed from the MSA sequence description...")
        return entry


def parse_uniprot_fasta_seq_description(match, desc, entry):
    """
    Pattern: <source>|<Accession_ID>|<Accession_Name> ++
    Example: sp|P00439|PH4H_HUMAN Phenylalanine-4-hydroxylase (...)
        OS=Homo sapiens GN=PAH PE=1 SV=1

    :param match: matched from the regex pattern
    :param desc: Biopython's 'description' field
    :param entry: Mutable dictionary
    :return: (updated) Dictionary
    """

    # pattern matching
    pattern = re.compile("([a-zA-Z])+\|")
    source = re.search(pattern, match, flags=0)
    if source:
        source = source.group().rstrip('|')
        entry['UniProt_Source'] = source

    pattern = re.compile("\|([a-zA-Z0-9])+\|")
    identifier = re.search(pattern, match, flags=0)
    if identifier:
        identifier = identifier.group().lstrip('|').rstrip('|')
        entry['UniProt_ID'] = identifier

    pattern = re.compile("\|([A-Z0-9])+_([A-Z0-9])+")
    name = re.search(pattern, match, flags=0)
    if name:
        name = name.group().lstrip('|')
        entry['UniProt_Name'] = name

    entry['Source'] = 'UniProt'

    # remaining description
    if desc != match:
        desc = desc.replace(match, "")
        entry['Description'] = desc.strip()

    return entry


def parse_pfam_sth_seq_description(match, desc, entry):
    """
    Pattern: <Accession_Name>/<Start>-<End>
    Example: C7P4T5_HALMD/44-372

    :param match: matched from the regex pattern
    :param desc: Biopython's 'description' field
    :param entry: Mutable dictionary
    :return: (updated) Dictionary
    """

    # pattern matching
    pattern = re.compile("([A-Z0-9])+[\_]?([a-zA-Z0-9])+\/[\-]?([0-9])+-[\-]?([0-9])+")
    # same as the generic seq description
    nmatch = re.search(pattern, match, flags=0)
    if nmatch:
        entry = parse_generic_seq_description(nmatch.group(), nmatch.group(), entry)

    entry['Source'] = 'Pfam'

    # remaining description
    if desc != match:
        desc = desc.replace(match, "")
        entry['Description'] = desc.strip()

    return entry


def parse_cath_fasta_seq_description(match, desc, entry):
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

    :param match: matched from the regex pattern
    :param desc: Biopython's 'description' field
    :param entry: Mutable dictionary
    :return: (updated) Dictionary
    """

    # pattern matching
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
                                              nmatch.group(), entry)

    entry['Source'] = 'CATH'

    # remaining description
    if desc != match:
        desc = desc.replace(match, "")
        entry['Description'] = desc.strip()

    return entry


def parse_generic_seq_description(match, desc, entry):
    """
    Pattern: <Accession_ID>/<Start>-<End> or <Accession_Name>/<Start>-<End>
    Example: P00439/24-145

    :param match: matched from the regex pattern
    :param desc: Biopython's 'description' field
    :param entry: Mutable dictionary
    :return: (updated) Dictionary
    """

    # pattern matching
    pattern = re.compile("([A-Z0-9])+_([a-zA-Z0-9])+\/")
    name = re.search(pattern, match, flags=0)
    if name:
        name = name.group().rstrip('/')
        entry['UniProt_Name'] = name

    pattern = re.compile("([a-zA-Z0-9])+\/")
    identifier = re.search(pattern, match, flags=0)
    if identifier:
        identifier = identifier.group().rstrip('/')
        entry['UniProt_ID'] = identifier

    pattern = re.compile("\/[\-]?([0-9])+")
    start = re.search(pattern, match, flags=0)
    if start:
        start = int(start.group().lstrip('/'))
        entry['Start'] = start

    pattern = re.compile("-[\-]?([0-9])+")
    end = re.search(pattern, match, flags=0)
    if end:
        end = int((end.group()).replace('--', '-'))
        entry['End'] = end

    # UniProt ID missing
    if name:
        r = fetch_uniprot_id_from_name(name, cached=True)
        try:
            identifier = r.json()[0]["id"]
        except:
            # new - to match the update in the UniProt API endpoint
            identifier = str(r.content, encoding='utf-8').strip()
        entry['UniProt_ID'] = identifier

    entry['Source'] = 'GenericParser'

    # remaining description
    if desc != match:
        desc = desc.replace(match, "")
        entry['Description'] = desc.strip()

    return entry


if __name__ == '__main__':
    pass
