#!/usr/bin/env python
# -*- coding: utf-8 -*-

import click
import click_log
import logging

from prointvar.fetchers import fetch_best_structures_pdbe
from prointvar.fetchers import download_structure_from_pdbe
from prointvar.fetchers import download_sifts_from_ebi
from prointvar.fetchers import download_data_from_uniprot
from prointvar.fetchers import download_alignment_from_cath
from prointvar.fetchers import download_alignment_from_pfam

from prointvar.utils import flash


logger = logging.getLogger("prointvar")


def read_alignment(inputfile):
    """
    Reads an input multiple sequence alignment (MSA).
    Formats recognised by Biopython:
        "clustal", "emboss", "nexus", "fasta", "phylip" and "stockholm"

    :param inputfile: Input MSA (read by biopython)
    :return: returns the biopython alignment object
    """

    from Bio import AlignIO

    if inputfile.endswith('.fasta') or inputfile.endswith('.fa'):
        align_format = 'fasta'
    elif inputfile.endswith('.sto') or inputfile.endswith('.sth'):
        align_format = 'stockholm'
    elif inputfile.endswith('.aln') or inputfile.endswith('.clw'):
        align_format = 'clustal'
    else:
        raise ValueError("Alignment format unrecognised...")

    alignment = AlignIO.read(inputfile, align_format)
    return alignment


# main application
@click.group(chain=True,
             context_settings={'help_option_names': ['-h', '--help']})
def cli():
    """Main script that process all the available sub-commands and options.
    """
    pass


@cli.command('alignment')
@click.option('-i', '--input', 'inputfile', type=str,
              multiple=False, help='The input file to open.',
              required=True)
@click.option('-o', '--output', 'outputfile', type=str,
              multiple=False, help='The output file to open.',
              required=False)
def alignment(inputfile, outputfile=None):

    from prointvar.fetchers import fetch_uniprot_id_from_name
    from prointvar.fetchers import fetch_best_structures_pdbe
    from prointvar.fetchers import fetch_uniprot_variants_ebi
    from prointvar.merger import TableMerger

    align = read_alignment(inputfile)
    var_rows = []
    for record in align:
        # seq = str(record.seq)
        desc = str(record.description)
        uniprot_name = desc.split('/')[0]
        r = fetch_uniprot_id_from_name(uniprot_name, cached=True)
        uniprot_id = str(r.content, encoding='utf-8').strip()
        print(uniprot_id)
        # Fetch/download structures
        r = fetch_best_structures_pdbe(uniprot_id, cached=True)
        pdb_ids = []
        if r is not None:
            for entry in r.json()[uniprot_id]:
                pdb_id = entry["pdb_id"]
                if pdb_id not in pdb_ids:
                    pdb_ids.append(pdb_id)

        for pdb_id in pdb_ids:
            file_downloader([pdb_id], mmcif=True, bio=True, sifts=True)

        # Fetch/download variants
        r = fetch_uniprot_variants_ebi(uniprot_id, cached=True)
        if r is not None:
            from prointvar.variants import flatten_uniprot_variants_ebi
            var_table = flatten_uniprot_variants_ebi(r)


@cli.command('mmcif2pdb')
@click.option('-i', '--input', 'inputfile', type=str,
              multiple=False, help='The input file to open.',
              required=True)
@click.option('-o', '--output', 'outputfile', type=str,
              multiple=False, help='The output file to open.',
              required=False)
@click.option('--override', 'override',
              multiple=False, help='Overrides any existing file, if available.',
              default=False, is_flag=True, required=False)
def mmcif2pdb(inputfile, outputfile=None, override=False):
    """
    Converts mmCIF to PDB format.

    :param inputfile: (str) path to the file
    :param outputfile: (str) path to the file
    :param override: boolean
    :return: (side-effects)
    """
    from prointvar.pdbx import PDBXwriter
    w = PDBXwriter(inputfile=inputfile, outputfile=outputfile)
    w.run(format_type='pdb', override=override)


@cli.command('clean_mmcif')
@click.option('-i', '--input', 'inputfile', type=str,
              multiple=False, help='The input file to open.',
              required=True)
@click.option('-o', '--output', 'outputfile', type=str,
              multiple=False, help='The output file to open.',
              required=False)
@click.option('--remove_hydrogens', 'rm_hydrogens',
              multiple=False, help='Removes hydrogens.',
              default=False, is_flag=True, required=False)
@click.option('--remove_altlocs', 'rm_altlocs',
              multiple=False, help='Removes alternative locations.',
              default=False, is_flag=True, required=False)
@click.option('--remove_partial_res', 'rm_partial_res',
              multiple=False, help='Removes incomplete residues.',
              default=False, is_flag=True, required=False)
@click.option('--override', 'override',
              multiple=False, help='Overrides any existing file, if available.',
              default=False, is_flag=True, required=False)
def clean_mmcif(inputfile, outputfile=None, rm_hydrogens=False, rm_altlocs=False,
              rm_partial_res=False, override=False):
    """
    Cleans the mmCIF file.

    :param inputfile: (str) path to the file
    :param outputfile: (str) path to the file
    :param rm_hydrogens: (bool) removes explicit hydrogens
    :param rm_altlocs: (bool) removes alternative locations
    :param rm_partial_res: (bool) removes incomplete residues
    :param override: boolean
    :return: (side-effects)
    """
    from prointvar.pdbx import PDBXwriter, PDBXreader
    r = PDBXreader(inputfile=inputfile)
    data = r.atoms(format_type='mmcif',
                   remove_hydrogens=rm_hydrogens, remove_altloc=rm_altlocs,
                   remove_partial_res=rm_partial_res)
    w = PDBXwriter(inputfile=None, outputfile=outputfile)
    w.run(data, format_type='mmcif', override=override)


@cli.command('clean_pdb')
@click.option('-i', '--input', 'inputfile', type=str,
              multiple=False, help='The input file to open.',
              required=True)
@click.option('-o', '--output', 'outputfile', type=str,
              multiple=False, help='The output file to open.',
              required=False)
@click.option('--remove_hydrogens', 'rm_hydrogens',
              multiple=False, help='Removes hydrogens.',
              default=False, is_flag=True, required=False)
@click.option('--remove_altlocs', 'rm_altlocs',
              multiple=False, help='Removes alternative locations.',
              default=False, is_flag=True, required=False)
@click.option('--remove_partial_res', 'rm_partial_res',
              multiple=False, help='Removes incomplete residues.',
              default=False, is_flag=True, required=False)
@click.option('--override', 'override',
              multiple=False, help='Overrides any existing file, if available.',
              default=False, is_flag=True, required=False)
@click.option('--pro_format', 'pro_format',
              multiple=False, help='',
              default=False, is_flag=True, required=False)
def clean_pdb(inputfile, outputfile=None, rm_hydrogens=False, rm_altlocs=False,
              rm_partial_res=False, override=False, pro_format=False):
    """
    Cleans the PDB file.

    :param inputfile: (str) path to the file
    :param outputfile: (str) path to the file
    :param rm_hydrogens: (bool) removes explicit hydrogens
    :param rm_altlocs: (bool) removes alternative locations
    :param rm_partial_res: (bool) removes incomplete residues
    :param override: boolean
    :return: (side-effects)
    """
    from prointvar.pdbx import PDBXwriter, PDBXreader
    r = PDBXreader(inputfile=inputfile)
    data = r.atoms(format_type='pdb',
                   remove_hydrogens=rm_hydrogens, remove_altloc=rm_altlocs,
                   remove_partial_res=rm_partial_res)
    w = PDBXwriter(inputfile=None, outputfile=outputfile)
    w.run(data, format_type='pdb', override=override, pro_format=pro_format)


def file_downloader(ids, pdb=False, mmcif=False, bio=False, sifts=False,
                    fasta=False, gff=False, txt=False, cath=False,
                    pfam=False, best_structures=False, override=False):
    for pid in ids:
        # best_structures takes precedence
        if best_structures:
            data = fetch_best_structures_pdbe(pid)
            if data is not None:
                data = data.json()
                for i, entry in enumerate(data[pid]):
                    pdb_id = entry['pdb_id']
                    chain_id = entry['chain_id']
                    flash("{}\t{}\t{}".format(i + 1, pdb_id, chain_id))
                    file_downloader([pdb_id], pdb=pdb, mmcif=mmcif,
                                    bio=bio, sifts=sifts)
                # once finish the first loop, skips over for the next pid
                continue
            else:
                flash('Best structures not available from the PDBe API for '
                      '{}'.format(pid))
        if pdb:
            download_structure_from_pdbe(pid, pdb=True,
                                         override=override)
        if mmcif:
            download_structure_from_pdbe(pid, pdb=False, bio=False,
                                         override=override)
        if bio:
            download_structure_from_pdbe(pid, pdb=False, bio=True,
                                         override=override)
        if sifts:
            download_sifts_from_ebi(pid, override=override)

        if fasta:
            download_data_from_uniprot(pid, file_format="fasta",
                                       override=override)
        if gff:
            download_data_from_uniprot(pid, file_format="gff",
                                       override=override)
        if txt:
            download_data_from_uniprot(pid, file_format="txt",
                                       override=override)
        if cath:
            download_alignment_from_cath(pid, max_sequences=20000,
                                         override=override)
        if pfam:
            download_alignment_from_pfam(pid, override=override)


@cli.command('download')
@click.argument('ids', nargs=-1, required=True)
@click_log.init("biodownloader")
@click_log.simple_verbosity_option()
@click.option('--pdb', 'pdb', multiple=False,
              help='Structure in PDB format (expects PDB ID).',
              default=False, is_flag=True, required=False)
@click.option('--mmcif', 'mmcif', multiple=False,
              help='Structure in mmCIF format (expects PDB ID).',
              default=False, is_flag=True, required=False)
@click.option('--bio', 'bio', multiple=False,
              help=('Preferred BioUnit instead of the asymmetric unit. '
                    'This option only works paired with --mmcif'),
              default=False, is_flag=True, required=False)
@click.option('--sifts', 'sifts', multiple=False,
              help='SIFTS xml format (expects PDB ID).',
              default=False, is_flag=True, required=False)
@click.option('--fasta', 'fasta', multiple=False,
              help='UniProt sequence in fasta format (expects UniProt ID).',
              default=False, is_flag=True, required=False)
@click.option('--gff', 'gff', multiple=False,
              help='UniProt record in gff format (expects UniProt ID).',
              default=False, is_flag=True, required=False)
@click.option('--txt', 'txt', multiple=False,
              help='UniProt record in txt format (expects UniProt ID).',
              default=False, is_flag=True, required=False)
@click.option('--cath', 'cath', multiple=False,
              help=('CATH Funfam alignment in fasta format '
                    '(expects a CATH <Superfamily>_<Funfam> ID).'),
              default=False, is_flag=True, required=False)
@click.option('--pfam', 'pfam', multiple=False,
              help=('Pfam alignment in Stockholm format '
                    '(expects a Pfam ID).'),
              default=False, is_flag=True, required=False)
@click.option('--best_structures', 'best_structures', multiple=False,
              help=('Structures that map to the UniProt ID (expects UniProt ID). '
                    'You should still pass which file types you need to download!'),
              default=False, is_flag=True, required=False)
@click.option('--override', 'override',
              multiple=False, help='Overrides any existing file, if available.',
              default=False, is_flag=True, required=False)
# TODO variants in vcf, json, etc.
def downloads(ids, pdb=False, mmcif=False, bio=False, sifts=False,
              fasta=False, gff=False, txt=False, cath=False,
              pfam=False, best_structures=False, override=False):
    """
    Downloads a variety of files from main repositories.

    Pass 1 or more accession IDs (e.g. '2pah' or '2pah 3kic').\n
    Currently accepted: PDB, UniProt, Pfam and CATH.
    """

    file_downloader(ids, pdb=pdb, mmcif=mmcif, bio=bio, sifts=sifts,
                    fasta=fasta, gff=gff, txt=txt, cath=cath, pfam=pfam,
                    best_structures=best_structures, override=override)


if __name__ == '__main__':
    cli()
