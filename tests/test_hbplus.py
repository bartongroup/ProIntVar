#!/local/bin/python
# -*- coding: utf-8 -*-


import os
import sys
import json
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

from prointvar.hbplus import (HBPLUSreader, HBPLUSrunner,
                              parse_hb2_from_file,
                              get_hbplus_selected_from_table)

from prointvar.config import config as c

root = os.path.abspath(os.path.dirname(__file__))
c.db_root = "{}/testdata/".format(root)


@patch("prointvar.config.config.db_root", c.db_root)
class TestHBPLUS(unittest.TestCase):
    """Test the HBPLUS parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.pdbid = '2pah'
        self.inputpdb = os.path.join(c.db_root, c.db_pdbx, "{}.pdb".format(self.pdbid))
        self.inputcif = os.path.join(c.db_root, c.db_pdbx, "{}.cif".format(self.pdbid))
        self.outputhbplus = os.path.join(c.db_root, c.db_contacts, "{}.h2b".format(self.pdbid))
        self.outputhbplus_h = os.path.join(c.db_root, c.db_pdbx, "{}.hbplus.pdb".format(self.pdbid))
        self.emptyfile = os.path.join(c.db_root, c.db_tmp, "{}.tmp".format(self.pdbid))
        self.notfound = ""
        self.excluded = ()

        self.parser = parse_hb2_from_file
        self.reader = HBPLUSreader
        self.generator = HBPLUSrunner
        self.filter = get_hbplus_selected_from_table

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.inputpdb = None
        self.inputcif = None
        self.outputhbplus = None

        self.emptyfile = None
        self.notfound = None
        self.excluded = None
        self.parser = None
        self.reader = None
        self.generator = None
        self.filter = None

        logging.disable(logging.NOTSET)

    def test_file_not_found_reader(self):
        with self.assertRaises(IOError):
            self.reader(self.notfound)

    def test_file_not_found_generator(self):
        with self.assertRaises(IOError):
            self.generator(self.notfound)

    def test_file_not_found_parser(self):
        with self.assertRaises(IOError):
            self.parser(self.notfound)

    def test_empty_file_reader(self):
        with self.assertRaises(ValueError):
            open(self.emptyfile, 'w').close()
            self.reader(self.emptyfile).read()
            os.remove(self.emptyfile)

    def test_generator_pdb_exec(self):
        if os.path.isfile(self.inputpdb):
            self.generator(self.inputpdb,
                           self.outputhbplus + '.test').run(clean_output=True,
                                                            override=True)
            msg = ("HBPLUS execution failed: make sure the settings "
                   "are set properly in config.ini!")
            self.assertTrue(os.path.isfile(self.outputhbplus + '.test'), msg)
            os.remove(self.outputhbplus + '.test')
        else:
            raise IOError("%s" % self.inputpdb)

    def test_generator_pdb(self):
        if os.path.isfile(self.inputpdb):
            self.generator(self.inputpdb, self.outputhbplus).run()
            self.assertTrue(os.path.isfile(self.outputhbplus))
        else:
            raise IOError("%s" % self.inputpdb)

    def test_generator_cif(self):
        if os.path.isfile(self.inputcif):
            self.generator(self.inputcif, self.outputhbplus).run()
            self.assertTrue(os.path.isfile(self.outputhbplus))
        else:
            raise IOError("%s" % self.inputcif)

    def test_parser_keys(self):
        self.assertListEqual([k for k in self.parser(self.outputhbplus).CHAIN_A.unique()],
                             ['A', 'B'])
        self.assertListEqual([k for k in self.parser(self.outputhbplus).CHAIN_D.unique()],
                             ['A', 'B'])

    def test_reader_data(self):
        reader = self.reader(self.outputhbplus)
        data = reader.read()
        self.assertEqual(data.loc[0, 'CHAIN_D'], 'A')
        self.assertEqual(data.loc[0, 'CHAIN_A'], 'A')
        self.assertEqual(data.loc[0, 'RES_D'], '123')
        self.assertEqual(data.loc[0, 'RES_A'], '127')
        self.assertEqual(data.loc[0, 'COMP_D'], 'ARG')
        self.assertEqual(data.loc[0, 'COMP_A'], 'GLU')
        self.assertEqual(data.loc[0, 'INSCODE_D'], '?')
        self.assertEqual(data.loc[0, 'INSCODE_A'], '?')
        self.assertEqual(data.loc[0, 'ATOM_D'], 'N')
        self.assertEqual(data.loc[0, 'ATOM_A'], 'OE1')
        self.assertEqual(data.loc[0, 'DIST_DA'], 2.96)

    def test_reader_to_json_pretty(self):
        reader = self.reader(self.outputhbplus)
        reader.read()
        data = reader.to_json()
        self.assertEqual(json.loads(data)[0]['CHAIN_D'], 'A')
        self.assertEqual(json.loads(data)[0]['DIST_DA'], 2.96)

    def test_reader_to_json(self):
        reader = self.reader(self.outputhbplus)
        reader.read()
        data = reader.to_json(pretty=False)
        self.assertEqual(json.loads(data)[0]['CHAIN_D'], 'A')
        self.assertEqual(json.loads(data)[0]['DIST_DA'], 2.96)

    def test_reader_default_excluded(self):
        reader = self.reader(self.outputhbplus)
        keys = reader.read()
        self.assertNotIn("NUM_AAS", keys)
        self.assertNotIn("DIST_CA-CA", keys)
        self.assertNotIn("DIST_H-A", keys)
        self.assertNotIn("ANGLE_D-H-A", keys)
        self.assertNotIn("ANGLE_H-A-AA", keys)
        self.assertNotIn("ANGLE_D-A-AA", keys)

    def test_reader_new_excluded(self):
        reader = self.reader(self.outputhbplus)
        keys = reader.read(excluded=self.excluded)
        self.assertIn("NUM_AAS", keys)
        self.assertIn("DIST_CA-CA", keys)
        self.assertIn("DIST_H-A", keys)
        self.assertIn("ANGLE_D-H-A", keys)
        self.assertIn("ANGLE_H-A-AA", keys)
        self.assertIn("ANGLE_D-A-AA", keys)

    def test_filter_chain_donor(self):
        reader = self.reader(self.outputhbplus)
        reader.read(excluded=self.excluded)
        data = self.filter(reader.data, chain_D=('A',))
        self.assertNotIn("B", data.CHAIN_D.unique())

    def test_filter_chain_acceptor(self):
        reader = self.reader(self.outputhbplus)
        reader.read(excluded=self.excluded)
        data = self.filter(reader.data, chain_A=('B',))
        self.assertNotIn("A", data.CHAIN_A.unique())

    def test_filter_res_donor(self):
        reader = self.reader(self.outputhbplus)
        reader.read(excluded=self.excluded)
        data = self.filter(reader.data, res_D=('123',))
        self.assertNotIn('119', data.RES_D.unique())

    def test_filter_res_acceptor(self):
        reader = self.reader(self.outputhbplus)
        reader.read(excluded=self.excluded)
        data = self.filter(reader.data, res_A=('127',))
        self.assertNotIn('119', data.RES_A.unique())

    def test_generator_pdb_hydrogen(self):
        if os.path.isfile(self.inputpdb):
            self.generator(self.inputpdb, self.outputhbplus_h).run(hydro_pdb_out=True,
                                                                   run_clean=False,
                                                                   override=True)
            self.assertTrue(os.path.isfile(self.outputhbplus_h))
        else:
            raise IOError("%s" % self.inputpdb)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("prointvar").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestHBPLUS)
    unittest.TextTestRunner(verbosity=2).run(suite)
