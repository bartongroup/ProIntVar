# -*- coding: utf-8 -*-

import os
import sys
import logging
import unittest

try:
    from mock import patch
except ImportError:
    from unittest.mock import patch

from prointvar.config import config
from prointvar.dssp import run_dssp, dssp_generate_output_filename, DSSP


@patch("prointvar.config.config", config)
class TestDSSP(unittest.TestCase):
    """Test the DSSP parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.pdbid = '2pah'
        self.inputcif = os.path.join(os.path.dirname(__file__), "testdata",
                                     "mmcif", "{}.cif".format(self.pdbid))
        self.inputdssp = os.path.join(os.path.dirname(__file__), "testdata",
                                      "dssp", "{}.dssp".format(self.pdbid))

        self.inputbiocif = os.path.join(os.path.dirname(__file__), "testdata",
                                        "mmcif", "{}_bio.cif".format(self.pdbid))
        self.inputbiodssp = os.path.join(os.path.dirname(__file__), "testdata",
                                         "dssp", "{}_bio.dssp".format(self.pdbid))

        self.notfound = ""
        self.DSSP = DSSP

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.inputcif = None
        self.inputdssp = None
        self.inputbiocif = None
        self.inputbiodssp = None
        self.notfound = None
        self.DSSP = None

        logging.disable(logging.NOTSET)

    def test_file_not_found_generator(self):
        with self.assertRaises(OSError):
            self.DSSP.generate(self.notfound)

    def test_generator_cif_exec(self):
        if os.path.isfile(self.inputcif):
            self.DSSP.generate(filename_input=self.inputcif,
                               filename_output=self.inputdssp + '.test',
                               overwrite=True)
            msg = ("DSSP execution failed: make sure the settings "
                   "are set properly in config.ini!")
            self.assertTrue(os.path.isfile(self.inputdssp + '.test'), msg)
            os.remove(self.inputdssp + '.test')
        else:
            raise IOError("%s" % self.inputcif)

    def test_generator_cif(self):
        if os.path.isfile(self.inputcif):
            self.DSSP.generate(self.inputcif, self.inputdssp)
            self.assertTrue(os.path.isfile(self.inputdssp))
        else:
            raise IOError("%s" % self.inputcif)

    def test_generator_biocif(self):
        if os.path.isfile(self.inputbiocif):
            self.DSSP.generate(self.inputbiocif, self.inputbiodssp)
            self.assertTrue(os.path.isfile(self.inputbiodssp))
        else:
            raise IOError("%s" % self.inputbiocif)

    def test_generator_run_unbound(self):
        if os.path.isfile(self.inputcif):
            filename, extension = os.path.splitext(self.inputdssp)
            self.DSSP.generate(self.inputcif, filename + '_unbound.dssp',
                               run_unbound=True, overwrite=True, category='auth')
            self.assertTrue(os.path.isfile(filename + '_unbound.dssp'))
        else:
            raise IOError("%s" % self.inputcif)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("prointvar").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDSSP)
    unittest.TextTestRunner(verbosity=2).run(suite)
