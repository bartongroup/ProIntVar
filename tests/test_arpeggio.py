#!/local/bin/python
# -*- coding: utf-8 -*-


import os
import sys
import json
import unittest
from contextlib import contextmanager

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
try:
    from mock import patch
except ImportError:
    from unittest.mock import patch

from prointvar.arpeggio import (ARPEGGIOreader, ARPEGGIOgenerator,
                                parse_arpeggio_from_file,
                                get_arpeggio_selected_from_table)

from prointvar.config import config as c
root = os.path.abspath(os.path.dirname(__file__))
c.db_root = "{}/testdata/".format(root)


@contextmanager
def captured_output():
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err


@patch("prointvar.config.config.db_root", c.db_root)
class TestARPEGGIO(unittest.TestCase):
    """Test the ARPEGGIO parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.pdbid = '2pah'
        self.pdbid_small = '2rea'
        self.inputpdb = "{}{}{}.pdb".format(c.db_root, c.db_cif, self.pdbid)
        self.inputpdb_fast = "{}{}{}.pdb".format(c.db_root, c.db_cif, self.pdbid_small)
        self.inputcif = "{}{}{}.cif".format(c.db_root, c.db_cif, self.pdbid)
        self.inputarpeggio = "{}{}{}.contacts".format(c.db_root,
                                                      c.db_contacts_generated, self.pdbid)
        self.inputarpeggio_fast = "{}{}{}.contacts".format(c.db_root,
                                                           c.db_contacts_generated,
                                                           self.pdbid_small)
        self.emptyfile = "{}{}{}.tmp".format(c.db_root, c.tmp_dir_local, self.pdbid)
        self.notfound = ""
        self.excluded = ()

        self.parser = parse_arpeggio_from_file
        self.reader = ARPEGGIOreader
        self.generator = ARPEGGIOgenerator
        self.filter = get_arpeggio_selected_from_table

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.pdbid_small = None
        self.inputpdb = None
        self.inputpdb_fast = None
        self.inputcif = None
        self.inputarpeggio = None
        self.inputarpeggio_fast = None

        self.emptyfile = None
        self.notfound = None
        self.excluded = None
        self.parser = None
        self.reader = None
        self.generator = None
        self.filter = None

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
        if os.path.isfile(self.inputpdb_fast):
            self.generator(self.inputpdb_fast,
                           self.inputarpeggio_fast).run(clean_output=True,
                                                        override=True)
            msg = ("Arpeggio execution failed: make sure the settings "
                   "are set properly in config.ini!")
            self.assertTrue(os.path.isfile(self.inputarpeggio_fast), msg)
            os.remove(self.inputarpeggio_fast)
        else:
            raise IOError("%s" % self.inputpdb_fast)

    def test_generator_pdb(self):
        if os.path.isfile(self.inputpdb):
            self.generator(self.inputpdb, self.inputarpeggio).run()
            self.assertTrue(os.path.isfile(self.inputarpeggio))
        else:
            raise IOError("%s" % self.inputpdb)

    def test_generator_cif(self):
        if os.path.isfile(self.inputcif):
            self.generator(self.inputcif, self.inputarpeggio).run()
            self.assertTrue(os.path.isfile(self.inputarpeggio))
        else:
            raise IOError("%s" % self.inputcif)

    def test_reader_arpeggio_verbose(self):
        with captured_output() as (out, err):
            self.reader(self.inputarpeggio, verbose=True).read()
        # This can go inside or outside the `with` block
        output = out.getvalue().strip()
        self.assertEqual(output, 'Parsing ARPEGGIO from lines...')

    def test_parser_keys(self):
        self.assertListEqual(sorted([k for k in
                                     self.parser(self.inputarpeggio).CHAIN_A.unique()]),
                             ['A', 'B'])
        self.assertListEqual(sorted([k for k in
                                     self.parser(self.inputarpeggio).CHAIN_B.unique()]),
                             ['A', 'B'])

    def test_reader_data(self):
        reader = self.reader(self.inputarpeggio)
        data = reader.read()
        self.assertEqual(data.loc[0, 'CHAIN_A'], 'B')
        self.assertEqual(data.loc[0, 'CHAIN_B'], 'B')
        self.assertEqual(data.loc[0, 'RES_A'], '376')
        self.assertEqual(data.loc[0, 'RES_B'], '374')
        self.assertEqual(data.loc[0, 'INSCODE_A'], '?')
        self.assertEqual(data.loc[0, 'INSCODE_B'], '?')
        self.assertEqual(data.loc[0, 'ATOM_A'], 'ND2')
        self.assertEqual(data.loc[0, 'ATOM_B'], 'O')
        self.assertEqual(data.loc[0, 'DIST'], 4.971)
        self.assertEqual(data.loc[0, 'VDW_DIST'], 1.901)

    def test_reader_to_json_pretty(self):
        reader = self.reader(self.inputarpeggio)
        reader.read()
        data = reader.to_json()
        self.assertEqual(json.loads(data)[0]['CHAIN_A'], 'B')
        self.assertEqual(json.loads(data)[0]['RES_A'], '376')

    def test_reader_to_json(self):
        reader = self.reader(self.inputarpeggio)
        reader.read()
        data = reader.to_json(pretty=False)
        self.assertEqual(json.loads(data)[0]['DIST'], 4.971)
        self.assertEqual(json.loads(data)[0]['VDW_DIST'], 1.901)

    def test_reader_default_excluded(self):
        reader = self.reader(self.inputarpeggio)
        keys = reader.read()
        self.assertNotIn("ENTRY_A", keys)
        self.assertNotIn("ENTRY_B", keys)

    def test_reader_new_excluded(self):
        reader = self.reader(self.inputarpeggio)
        keys = reader.read(excluded=self.excluded)
        self.assertIn("ENTRY_A", keys)
        self.assertIn("ENTRY_B", keys)

    def test_filter_chain(self):
        reader = self.reader(self.inputarpeggio)
        reader.read(excluded=self.excluded)
        data = self.filter(reader.data, chain_B=('A',))
        self.assertNotIn("B", data.CHAIN_B.unique())

    def test_filter_res(self):
        reader = self.reader(self.inputarpeggio)
        reader.read(excluded=self.excluded)
        data = self.filter(reader.data, res_A=('374',))
        self.assertNotIn('119', data.RES_A.unique())


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestARPEGGIO)
    unittest.TextTestRunner(verbosity=2).run(suite)
