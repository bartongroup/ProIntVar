#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This defines the methods that work with REDUCE.

Fábio Madeira, 2017+

"""

import os

from prointvar.mmcif import MMCIFwriter

from prointvar.utils import flash
from prointvar.utils import lazy_file_remover

from prointvar.config import config


class REDUCEgenerator(object):
    def __init__(self, inputfile, outputfile=None, verbose=False):
        """
        :param inputfile: Needs to point to a valid PDB or mmCIF file.
        :param outputfile: if not provided will use the same file name and
          <*.h> extension
        :param verbose: boolean
        """
        self.inputfile = inputfile
        self.inputfile_back = inputfile
        self.outputfile = outputfile
        self.verbose = verbose
        self.data = None

        if not os.path.isfile(self.inputfile):
            raise IOError("{} not available or could not be read..."
                          "".format(self.inputfile))

        # inputfile needs to be in PDB or mmCIF format
        filename, extension = os.path.splitext(self.inputfile)
        if extension not in ['.pdb', '.ent', '.cif']:
            raise ValueError("{} is expected to be in mmCIF or PDB format..."
                             "".format(self.inputfile))

    def _generate_output(self):
        filename, extension = os.path.splitext(self.inputfile)
        self.outputfile = filename + ".h"

    def _generate_pdb(self, override=False):
        filename, extension = os.path.splitext(self.inputfile)
        self.inputfile = filename + "_new.pdb"
        w = MMCIFwriter(inputfile=self.inputfile_back, outputfile=self.inputfile)
        w.run(format_type="pdb", override=override)

    def _run(self, reduce_bin, clean_output=False):

        # reduce generates new coordinates adding Hydrogen atoms
        cmd = "{} -noflip -quiet {} > {}".format(reduce_bin, self.inputfile,
                                                 self.outputfile)
        os.system(cmd)
        if not os.path.isfile(self.outputfile):
            raise IOError("Reduce output not generated for {}"
                          "".format(self.outputfile))

    def run(self, override=False, clean_output=True, save_new_input=False):

        # generate outputfile if missing
        if not self.outputfile:
            self._generate_output()

        if not os.path.exists(self.outputfile) or override:
            if os.path.isfile(config.reduce_bin):
                reduce_bin = config.reduce_bin
            else:
                raise IOError('REDUCE executables are not available...')

            # inputfile needs to be in PDB format
            filename, extension = os.path.splitext(self.inputfile)
            if extension == '.cif':
                self._generate_pdb(override=override)

            # run probe and generate output - also clean unnecessary output
            self._run(reduce_bin, clean_output=clean_output)

            # clean the new PDB input file generated
            if not save_new_input:
                if self.inputfile != self.inputfile_back:
                    lazy_file_remover(self.inputfile)

        else:
            flash('REDUCE for {} already available...'.format(self.outputfile))
        return


if __name__ == '__main__':
    pass
