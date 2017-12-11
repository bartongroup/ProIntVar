# -*- coding: utf-8 -*-

"""

This defines the methods that work with REDUCE.

Fábio Madeira, 2017+

"""

import os
import logging

from proteofav.structures import PDB, mmCIF, filter_structures
from proteofav.utils import GenericInputs
from proteofav.utils import InputFileHandler

from prointvar.utils import lazy_file_remover

from prointvar.config import config

logger = logging.getLogger("prointvar")


def run_reduce(filename_input, filename_output=None,
               save_new_input=False, overwrite=False):
    """
    Runs REDUCE to add explicit Hydrogen atoms to a PDB
    structure.

    :param filename_input: path to input file.
      Needs to point to a valid PDB or mmCIF file.
    :param filename_output: path to output file
      if not provided will use the same file name and <*.h> extension
    :param save_new_input: boolean
    :param overwrite: boolean
    :return: Runs REDUCE on the provided PDB structure
    """

    InputFileHandler(filename_input)

    # inputfile needs to be in PDB or mmCIF format
    filename, extension = os.path.splitext(filename_input)
    if extension not in ['.pdb', '.ent', '.cif']:
        raise ValueError("{} is expected to be in mmCIF or PDB format..."
                         "".format(filename_input))

    if not filename_output:
        filename, extension = os.path.splitext(filename_input)
        filename_output = filename + ".h"

    if not os.path.exists(filename_output) or overwrite:
        if os.path.isfile(config.reduce_bin):
            reduce_bin = config.reduce_bin
        else:
            raise IOError('REDUCE executables are not available...')

        # inputfile needs to be in PDB format
        filename_input_back = filename_input
        filename, extension = os.path.splitext(filename_input)
        if extension == '.cif':
            filename, extension = os.path.splitext(filename_input)
            filename_input = filename + "_new.pdb"
            r = mmCIF.read(filename=filename_input_back)
            table = filter_structures(r, add_res_full=False, add_contacts=False,
                                      lines='ATOM', category='auth')
            PDB.write(table=table, filename=filename_input,
                      output_format="pdb", overwrite=overwrite)

        # run probe and generate output - also clean unnecessary output
        # reduce generates new coordinates adding Hydrogen atoms
        cmd = "{} -noflip -quiet {} > {}".format(reduce_bin, filename_input,
                                                 filename_output)
        os.system(cmd)
        if not os.path.isfile(filename_output):
            raise IOError("Reduce output not generated for {}"
                          "".format(filename_output))

        # clean the new PDB input file generated
        if not save_new_input:
            if filename_input != filename_input_back:
                lazy_file_remover(filename_input)
    else:
        logger.info("REDUCE for %s already available...", filename_output)


class _REDUCE(GenericInputs):
    def run(self, filename_input=None, filename_output=None, **kwargs):
        self.table = run_reduce(filename_input=filename_input,
                                filename_output=filename_output, **kwargs)
        return self.table


REDUCE = _REDUCE()
