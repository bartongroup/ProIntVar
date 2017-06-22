#!/usr/bin/env python
# -*- coding: utf-8 -*-

import click


# main application
@click.group(chain=True,
             context_settings={'help_option_names': ['-h', '--help']})
def cli():
    """Main script that process all the available sub-commands and options.
    """
    pass


@cli.command('mmcif2pdb')
@click.option('-i', '--input', 'inputfile', type=str,
              multiple=False, help='The input file to open.',
              required=True)
@click.option('-o', '--output', 'outputfile', type=str,
              multiple=False, help='The output file to open.',
              required=False)
def mmcif2pdb(inputfile, outputfile=None):
    from prointvar.pdbx import PDBXwriter
    w = PDBXwriter(inputfile=inputfile, outputfile=outputfile)
    w.run(format_type='pdb')


@cli.command('mmcif_clean')
@click.option('-i', '--input', 'inputfile', type=str,
              multiple=False, help='The input file to open.',
              required=True)
@click.option('-o', '--output', 'outputfile', type=str,
              multiple=False, help='The output file to open.',
              required=False)
@click.option('--remove_hydrogens', 'rm_hydrogens',
              multiple=False, help='', default=False, is_flag=True,
              required=False)
@click.option('--remove_altlocs', 'rm_altlocs',
              multiple=False, help='', default=False, is_flag=True,
              required=False)
def mmcif_clean(inputfile, outputfile=None, rm_hydrogens=False, rm_altlocs=False):
    from prointvar.pdbx import PDBXwriter, PDBXreader
    r = PDBXreader(inputfile=inputfile)
    data = r.atoms(format_type='mmcif',
                   remove_hydrogens=rm_hydrogens, remove_altloc=rm_altlocs)
    w = PDBXwriter(inputfile=None, outputfile=outputfile)
    w.run(data, format_type='mmcif')


@cli.command('pdb_clean')
@click.option('-i', '--input', 'inputfile', type=str,
              multiple=False, help='The input file to open.',
              required=True)
@click.option('-o', '--output', 'outputfile', type=str,
              multiple=False, help='The output file to open.',
              required=False)
@click.option('--remove_hydrogens', 'rm_hydrogens',
              multiple=False, help='',
              default=False, is_flag=True, required=False)
@click.option('--remove_altlocs', 'rm_altlocs',
              multiple=False, help='',
              default=False, is_flag=True, required=False)
@click.option('--remove_partial_res', 'rm_partial_res',
              multiple=False, help='',
              default=False, is_flag=True, required=False)
@click.option('--override', 'override',
              multiple=False, help='',
              default=False, is_flag=True, required=False)
@click.option('--pro_format', 'pro_format',
              multiple=False, help='',
              default=False, is_flag=True, required=False)
def pdb_clean(inputfile, outputfile=None, rm_hydrogens=False, rm_altlocs=False,
              rm_partial_res=False, override=False, pro_format=False):
    from prointvar.pdbx import PDBXwriter, PDBXreader
    r = PDBXreader(inputfile=inputfile)
    data = r.atoms(format_type='pdb',
                   remove_hydrogens=rm_hydrogens, remove_altloc=rm_altlocs,
                   remove_partial_res=rm_partial_res)
    w = PDBXwriter(inputfile=None, outputfile=outputfile)
    w.run(data, format_type='pdb', override=override, pro_format=pro_format)


if __name__ == '__main__':
    cli()
