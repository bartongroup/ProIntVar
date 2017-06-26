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

from prointvar.dssp import (DSSPreader, DSSPrunner,
                            parse_dssp_from_file, get_dssp_selected_from_table,
                            add_dssp_full_chain, add_dssp_rsa, add_dssp_rsa_class,
                            add_dssp_ss_reduced)

from prointvar.config import config as c
root = os.path.abspath(os.path.dirname(__file__))
c.db_root = "{}/testdata/".format(root)


@patch("prointvar.config.config.db_root", c.db_root)
class TestDSSP(unittest.TestCase):
    """Test the DSSP parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.pdbid = '2pah'
        self.inputcif = os.path.join(c.db_root, c.db_pdbx, "{}.cif".format(self.pdbid))
        self.inputdssp = os.path.join(c.db_root, c.db_dssp, "{}.dssp".format(self.pdbid))

        self.inputbiocif = os.path.join(c.db_root, c.db_pdbx, "{}.cif".format(self.pdbid))
        self.inputbiodssp = os.path.join(c.db_root, c.db_dssp, "{}_bio.dssp".format(self.pdbid))

        self.emptyfile = os.path.join(c.db_root, c.db_tmp, "{}.tmp".format(self.pdbid))
        self.notfound = ""
        self.excluded = ()

        self.parser = parse_dssp_from_file
        self.reader = DSSPreader
        self.generator = DSSPrunner
        self.filter = get_dssp_selected_from_table
        self.add_full_chain = add_dssp_full_chain
        self.add_rsa = add_dssp_rsa
        self.add_rsa_class = add_dssp_rsa_class
        self.add_ss_reduced = add_dssp_ss_reduced

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.inputcif = None
        self.inputdssp = None
        self.inputbiocif = None
        self.inputbiodssp = None
        self.emptyfile = None
        self.notfound = None
        self.excluded = None
        self.parser = None
        self.reader = None
        self.generator = None
        self.filter = None
        self.add_full_chain = None
        self.add_rsa = None
        self.add_rsa_class = None
        self.add_ss_reduced = None

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

    def test_generator_cif_exec(self):
        if os.path.isfile(self.inputcif):
            self.generator(self.inputcif, self.inputdssp + '.test').run(override=True)
            msg = ("DSSP execution failed: make sure the settings "
                   "are set properly in config.ini!")
            self.assertTrue(os.path.isfile(self.inputdssp + '.test'), msg)
            os.remove(self.inputdssp + '.test')
        else:
            raise IOError("%s" % self.inputcif)

    def test_generator_cif(self):
        if os.path.isfile(self.inputcif):
            self.generator(self.inputcif, self.inputdssp).run()
            self.assertTrue(os.path.isfile(self.inputdssp))
        else:
            raise IOError("%s" % self.inputcif)

    def test_generator_biocif(self):
        if os.path.isfile(self.inputbiocif):
            self.generator(self.inputbiocif, self.inputbiodssp).run()
            self.assertTrue(os.path.isfile(self.inputbiodssp))
        else:
            raise IOError("%s" % self.inputbiocif)

    def test_parser_keys(self):
        self.assertListEqual([k for k in self.parser(self.inputdssp).CHAIN.unique()],
                             ['A', 'B'])

    def test_reader_cif_keys(self):
        self.assertListEqual([k for k in self.reader(self.inputdssp).read().CHAIN.unique()],
                             ['A', 'B'])

    def test_reader_biocif_keys(self):
        self.assertListEqual([k for k in self.reader(self.inputbiodssp).read().CHAIN.unique()],
                             ['A', 'B'])

    def test_reader_cif_data(self):
        reader = self.reader(self.inputdssp)
        data = reader.read()
        self.assertEqual(data.loc[0, 'CHAIN'], 'A')
        self.assertEqual(data.loc[0, 'RES'], '118')
        self.assertEqual(data.loc[331, 'CHAIN'], 'B')
        self.assertEqual(data.loc[331, 'RES'], '118')

    def test_reader_biocif_data(self):
        reader = self.reader(self.inputbiodssp)
        data = reader.read()
        self.assertEqual(data.loc[0, 'CHAIN'], 'A')
        self.assertEqual(data.loc[0, 'RES'], '118')
        self.assertEqual(data.loc[662, 'CHAIN'], 'B')
        self.assertEqual(data.loc[662, 'RES'], '118')

    def test_reader_residues_cif_data(self):
        reader = self.reader(self.inputdssp)
        data = reader.residues()
        self.assertEqual(data.loc[0, 'CHAIN'], 'A')
        self.assertEqual(data.loc[0, 'RES'], '118')
        self.assertEqual(data.loc[331, 'CHAIN'], 'B')
        self.assertEqual(data.loc[331, 'RES'], '118')

    def test_reader_cif_to_json_pretty(self):
        reader = self.reader(self.inputdssp)
        reader.read()
        data = reader.to_json()
        self.assertEqual(json.loads(data)[0]['CHAIN'], 'A')
        self.assertEqual(json.loads(data)[0]['RES'], '118')

    def test_reader_cif_to_json(self):
        reader = self.reader(self.inputdssp)
        reader.read()
        data = reader.to_json(pretty=False)
        self.assertEqual(json.loads(data)[0]['CHAIN'], 'A')
        self.assertEqual(json.loads(data)[0]['RES'], '118')

    def test_reader_default_excluded(self):
        reader = self.reader(self.inputdssp)
        keys = reader.read()
        self.assertNotIn("LINE", keys)
        self.assertNotIn("STRUCTURE", keys)
        self.assertNotIn("BP1", keys)
        self.assertNotIn("BP2", keys)
        self.assertNotIn("BP2_CHAIN", keys)
        self.assertNotIn("X-CA", keys)
        self.assertNotIn("Y-CA", keys)
        self.assertNotIn("Z-CA", keys)

    def test_reader_new_excluded(self):
        reader = self.reader(self.inputdssp)
        keys = reader.read(excluded=self.excluded)
        self.assertIn("LINE", keys)
        self.assertIn("STRUCTURE", keys)
        self.assertIn("BP1", keys)
        self.assertIn("BP2", keys)
        self.assertIn("BP2_CHAIN", keys)
        self.assertIn("X-CA", keys)
        self.assertIn("Y-CA", keys)
        self.assertIn("Z-CA", keys)

    def test_filter_chain(self):
        reader = self.reader(self.inputdssp)
        reader.read(excluded=self.excluded)
        data = self.filter(reader.data, chain=('A',))
        self.assertNotIn("B", data.CHAIN.unique())

    def test_filter_chain_full(self):
        reader = self.reader(self.inputbiodssp)
        reader.read(excluded=self.excluded)
        data = self.filter(reader.data, chain_full=('BA',))
        self.assertIn("B", data.CHAIN.unique())
        self.assertNotIn("B", data.CHAIN_FULL.unique())

    def test_filter_res(self):
        reader = self.reader(self.inputdssp)
        reader.read(excluded=self.excluded)
        data = self.filter(reader.data, res=('118',))
        self.assertNotIn('119', data.RES.unique())

    def test_reader_add_chain_full(self):
        reader = self.reader(self.inputbiodssp)
        data = reader.read(add_full_chain=True)
        self.assertIn("BA", data.CHAIN_FULL.unique())

    # @unittest.expectedFailure
    def test_add_chain_full(self):
        reader = self.reader(self.inputbiodssp)
        data = reader.read(add_full_chain=True)
        data = self.add_full_chain(data)
        # because the chain breaks "!*" have been removed already
        # FIXME add an boolean option to remove chain breaks?
        self.assertNotIn("AA", data.CHAIN_FULL.unique())

    def test_reader_add_rsa(self):
        reader = self.reader(self.inputbiodssp)
        data = reader.read(add_rsa=True)
        self.assertEqual(52.863, data.loc[2, 'RSA'])

    def test_add_rsa(self):
        reader = self.reader(self.inputbiodssp)
        data = reader.read(add_rsa=True)
        data = self.add_rsa(data)
        self.assertEqual(52.863, data.loc[2, 'RSA'])

    def test_reader_add_ss_reduced(self):
        reader = self.reader(self.inputbiodssp)
        data = reader.read(add_ss_reduced=True)
        self.assertEqual('H', data.loc[29, 'SS_CLASS'])

    def test_ss_reduced(self):
        reader = self.reader(self.inputbiodssp)
        data = reader.read(add_ss_reduced=True)
        data = self.add_ss_reduced(data)
        self.assertEqual('H', data.loc[29, 'SS_CLASS'])

    def test_reader_add_rsa_class(self):
        reader = self.reader(self.inputbiodssp)
        data = reader.read(add_rsa=True, add_rsa_class=True)
        self.assertEqual('Surface', data.loc[2, 'RSA_CLASS'])

    def test_add_rsa_class(self):
        reader = self.reader(self.inputbiodssp)
        data = reader.read(add_rsa=True)
        data = self.add_rsa_class(data)
        self.assertEqual('Surface', data.loc[2, 'RSA_CLASS'])

    def test_generator_run_unbound(self):
        if os.path.isfile(self.inputcif):
            filename, extension = os.path.splitext(self.inputdssp)
            g = self.generator(self.inputcif, filename + '_unbound.dssp')
            g.run(run_unbound=True, override=True)
            self.assertTrue(os.path.isfile(filename + '_unbound.dssp'))
        else:
            raise IOError("%s" % self.inputcif)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("prointvar").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDSSP)
    unittest.TextTestRunner(verbosity=2).run(suite)
