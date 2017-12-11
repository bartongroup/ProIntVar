# -*- coding: utf-8 -*-


import os
import sys
import logging
import unittest

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
try:
    from mock import patch
except ImportError:
    from unittest.mock import patch

from prointvar.reduce import REDUCE

from prointvar.config import config


@patch("prointvar.config.config", config)
class TestREDUCE(unittest.TestCase):
    """Test the REDUCE parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.pdbid = '2pah'
        self.inputpdb = os.path.join(os.path.dirname(__file__), "testdata",
                                     config.db_pdb, "{}.pdb".format(self.pdbid))
        self.inputcif = os.path.join(os.path.dirname(__file__), "testdata",
                                     config.db_mmcif, "{}.cif".format(self.pdbid))
        self.outputred = os.path.join(os.path.dirname(__file__), "testdata",
                                      config.db_pdb, "{}.reduce.pdb".format(self.pdbid))
        self.excluded = ()
        self.REDUCE = REDUCE

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.inputpdb = None
        self.inputcif = None
        self.outputred = None
        self.excluded = None
        self.REDUCE = None

        logging.disable(logging.NOTSET)

    def test_generator_pdb(self):
        if os.path.isfile(self.inputpdb):
            self.REDUCE.run(self.inputpdb, self.outputred,
                            overwrite=True, save_new_input=False)
            self.assertTrue(os.path.isfile(self.outputred))
        else:
            raise IOError("%s" % self.inputpdb)

    def test_generator_cif(self):
        if os.path.isfile(self.inputcif):
            self.REDUCE.run(self.inputcif, self.outputred,
                            overwrite=True, save_new_input=False)
            self.assertTrue(os.path.isfile(self.outputred))
        else:
            raise IOError("%s" % self.inputcif)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("prointvar").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestREDUCE)
    unittest.TextTestRunner(verbosity=2).run(suite)
