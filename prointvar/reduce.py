#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This defines the methods that work with REDUCE.

Fábio Madeira, 2017+

"""

import os

from prointvar.mmcif import MMCIFwriter

from prointvar.utils import flash

from prointvar.config import config


class REDUCEgenerator(object):
    def __init__(self, inputfile, outputfile=None, verbose=False):
        """
        :param inputfile: Needs to point to a valid PDB or mmCIF file.
        :param outputfile: if not provided will use the same file name and <*.h> extension
        :param verbose: boolean
        """
        self.inputfile = inputfile
        self.outputfile = outputfile
        self.verbose = verbose
        self.data = None

        # generate outputfile if missing
        self._generate_output()

        if not os.path.isfile(inputfile):
            raise IOError("{} not available or could not be read...".format(inputfile))

        # inputfile needs to be in PDB format
        filename, extension = os.path.splitext(inputfile)
        if extension in ['.pdb', '.ent']:
            pass
        elif extension in ['.cif']:
            self._generate_pdb_from_mmcif()
        else:
            raise ValueError("{} is expected to be in mmCIF or PDB format..."
                             "".format(inputfile))

    def _generate_output(self):
        if not self.outputfile:
            filename, extension = os.path.splitext(self.inputfile)
            self.outputfile = filename + ".h"

    def _generate_pdb_from_mmcif(self):
        filename, extension = os.path.splitext(self.inputfile)
        w = MMCIFwriter(inputfile=self.inputfile, outputfile=filename + ".pdb")
        w.run(format_type="pdb")
        self.inputfile = filename + ".pdb"

    def _run(self, reduce_bin):

        # reduce generates new coordinates adding Hydrogen atoms
        cmd = "{} -noflip -quiet {} > {}".format(reduce_bin, self.inputfile,
                                                 self.outputfile)
        os.system(cmd)
        if not os.path.isfile(self.outputfile):
            raise IOError("Reduce output not generated for {}"
                          "".format(self.outputfile))

    def run(self, override=False):
        if not os.path.exists(self.outputfile) or override:
            if os.path.isfile(config.reduce_bin):
                reduce_bin = config.reduce_bin
            elif os.path.isfile(config.reduce_bin_local):
                reduce_bin = config.reduce_bin_local
            else:
                raise IOError('REDUCE executables are not available...')

            # run probe and generate output
            self._run(reduce_bin)

        else:
            flash('REDUCE for {} already available...'.format(self.outputfile))
        return


if __name__ == '__main__':
    pass
