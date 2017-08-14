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
from prointvar.msas import parse_msa_sequences_from_file

from prointvar.utils import flash


logger = logging.getLogger("prointvar")


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

    # print(inputfile)
    # print(outputfile)
    align = parse_msa_sequences_from_file(inputfile,
                                          get_uniprot_id=True, cached=True)
    var_rows = []
    for ix in align.index:
        # get descriptions
        uniprot_id = align.loc[ix, 'Accession']
        print(uniprot_id)
        # print(align.loc[ix, :])
        # print(align.loc[ix, 'Description'])
        # print(align.loc[ix, 'Source'])
        # print(align.loc[ix, 'Desc'])

        # # Fetch/download structures
        # r = fetch_best_structures_pdbe(uniprot_id, cached=True)
        # pdb_ids = []
        # if r is not None:
        #     # print(r.json())
        #     for entry in r.json()[uniprot_id]:
        #         pdb_id = entry["pdb_id"]
        #         if pdb_id not in pdb_ids:
        #             pdb_ids.append(pdb_id)
        #
        # for pdb_id in pdb_ids:
        #     file_downloader([pdb_id], mmcif=True, bio=True, sifts=True)

        # # # UniProt Variants
        # # Fetch/download variants
        # r = fetch_uniprot_variants_ebi(uniprot_id, cached=True)
        # if r is not None:
        #     from prointvar.variants import flatten_uniprot_variants_ebi
        #     var_table = flatten_uniprot_variants_ebi(r)
        #     print(var_table.head())
        #     print(var_table.loc[0, :])
        #     print(list(var_table))
        #     break

        # # # Ensembl Variants
        from prointvar.fetchers import fetch_uniprot_species_from_id
        from prointvar.fetchers import fetch_ensembl_uniprot_ensembl_mapping
        from prointvar.fetchers import fetch_ensembl_transcript_variants
        from prointvar.fetchers import fetch_ensembl_somatic_variants
        from prointvar.fetchers import fetch_ensembl_variants_by_id
        info = fetch_uniprot_species_from_id(uniprot_id)
        organism = str(info.content, encoding='utf-8').split('\n')[1]
        species = '_'.join(organism.split()[0:2]).lower()
        info = fetch_ensembl_uniprot_ensembl_mapping(uniprot_id, cached=True,
                                                     species=species).json()
        ensps = []
        for entry in info:
            if 'type' in entry and 'id' in entry:
                if entry['type'] == 'translation':
                    if entry['id'] not in ensps:
                        ensps.append(entry['id'])
        for ensp in ensps:
            print(ensp)
            ens_vars = fetch_ensembl_transcript_variants(ensp, cached=False).json()
            # print(ens_vars)

            som_vars = fetch_ensembl_somatic_variants(ensp, cached=False).json()
            # print(som_vars)

            var_ids = []
            for entry in ens_vars:
                print(entry)
                if 'id' in entry:
                    if entry['id'] not in var_ids:
                        var_ids.append(entry['id'])
            for entry in som_vars:
                print(entry)
                if 'id' in entry:
                    if entry['id'] not in var_ids:
                        var_ids.append(entry['id'])
            print(var_ids, len(var_ids))

            # maybe not needed?
            vars = fetch_ensembl_variants_by_id(var_ids, cached=False, species=species).json()
            print([vars])




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
