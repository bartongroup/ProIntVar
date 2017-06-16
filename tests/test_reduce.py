#!/local/bin/python
# -*- coding: utf-8 -*-


import os
import unittest

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
try:
    from mock import patch
except ImportError:
    from unittest.mock import patch

from prointvar.reduce import REDUCEgenerator

from prointvar.config import config as c

root = os.path.abspath(os.path.dirname(__file__))
c.db_root = "{}/testdata/".format(root)


@patch("prointvar.config.config.db_root", c.db_root)
class TestREDUCE(unittest.TestCase):
    """Test the REDUCE parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.pdbid = '2pah'
        self.inputpdb = "{}{}{}.pdb".format(c.db_root, c.db_cif, self.pdbid)
        self.inputcif = "{}{}{}.cif".format(c.db_root, c.db_cif, self.pdbid)
        self.outputred = "{}{}{}.reduce.pdb".format(c.db_root, c.db_cif, self.pdbid)
        self.emptyfile = "{}{}{}.tmp".format(c.db_root, c.tmp_dir_local, self.pdbid)
        self.notfound = ""
        self.excluded = ()

        self.generator = REDUCEgenerator

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.inputpdb = None
        self.inputcif = None
        self.outputred = None

        self.emptyfile = None
        self.notfound = None
        self.excluded = None

        self.generator = None

    def test_file_not_found_generator(self):
        with self.assertRaises(IOError):
            self.generator(self.notfound)

    def test_generator_pdb(self):
        if os.path.isfile(self.inputpdb):
            self.generator(self.inputpdb, self.outputred).run(override=True)
            self.assertTrue(os.path.isfile(self.outputred))
        else:
            raise IOError("%s" % self.inputpdb)

    def test_generator_cif(self):
        if os.path.isfile(self.inputcif):
            self.generator(self.inputcif, self.outputred).run(override=True)
            self.assertTrue(os.path.isfile(self.outputred))
        else:
            raise IOError("%s" % self.inputcif)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestREDUCE)
    unittest.TextTestRunner(verbosity=2).run(suite)
