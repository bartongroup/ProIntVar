#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
from setuptools import setup
from setuptools import find_packages

from prointvar import __version__


def gather_dependencies():
    with open('requirements.txt', 'r') as f_in:
        return [l for l in f_in.read().rsplit(os.linesep)
                if l and not l.startswith("#")]
DEPENDENCIES = gather_dependencies()


setup(
    name="ProIntVar",
    version=__version__,
    packages=find_packages(exclude=["tests", 'tests.*']),
    # should always match the entries in requirements.txt
    install_requires=DEPENDENCIES,
    package_data={'prointvar': ['config_template.ini']},
    include_package_data=True,

    entry_points={
        "console_scripts": [
            "ProIntVar-config-setup=prointvar.config:config_setup",
            "ProIntVar-config-load=prointvar.config:config_load",
            "ProIntVar=cli.main:cli",
        ]
    },

    author="Fábio Madeira",
    author_email="fabiomadeira@me.com",
    url="https://github.com/bartongroup/ProIntVar",
    download_url="https://github.com/bartongroup/ProIntVar/archive/master.zip",
    license='MIT',
    keywords='python pdb structures pandas dssp sifts ensembl uniprot alignments',
    description=('Python module that implements methods for working with Protein Structures '
                 'and Genetic Variation'),
    long_description=open('README.md').read(),
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License"
    ]
)
