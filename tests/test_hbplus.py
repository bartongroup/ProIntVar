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

from prointvar.hbplus import (parse_hb2_from_file,
                              filter_hbplus, hbplus_generate_output_filename,
                              run_hbplus, HBPLUS)

from prointvar.config import config


@patch("prointvar.config.config", config)
class TestHBPLUS(unittest.TestCase):
    """Test the HBPLUS parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.pdbid = '2pah'
        self.inputpdb = os.path.join(os.path.dirname(__file__), "testdata",
                                     config.db_pdb, "{}.pdb".format(self.pdbid))
        self.inputcif = os.path.join(os.path.dirname(__file__), "testdata",
                                     config.db_mmcif, "{}.cif".format(self.pdbid))
        self.outputhbplus = os.path.join(os.path.dirname(__file__), "testdata",
                                         config.db_contacts, "{}.h2b".format(self.pdbid))
        self.outputhbplus_h = os.path.join(os.path.dirname(__file__), "testdata",
                                           config.db_pdb, "{}.hbplus.pdb".format(self.pdbid))
        self.excluded = ()
        self.parser = parse_hb2_from_file
        self.filter_hbplus = filter_hbplus
        self.run_hbplus = run_hbplus
        self.generate_output_filename = hbplus_generate_output_filename
        self.HBPLUS = HBPLUS

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.inputpdb = None
        self.inputcif = None
        self.outputhbplus = None
        self.excluded = None
        self.parser = None
        self.parser = None
        self.filter_hbplus = None
        self.run_hbplus = None
        self.generate_output_filename = None
        self.HBPLUS = None

        logging.disable(logging.NOTSET)

    def test_generator_pdb_exec(self):
        if os.path.isfile(self.inputpdb):
            self.HBPLUS.run(self.inputpdb,
                            self.outputhbplus + '.test',
                            clean_output=True, overwrite=True)
            msg = ("HBPLUS execution failed: make sure the settings "
                   "are set properly in config.ini!")
            self.assertTrue(os.path.isfile(self.outputhbplus + '.test'), msg)
            os.remove(self.outputhbplus + '.test')
        else:
            raise IOError("%s" % self.inputpdb)

    def test_generator_pdb(self):
        if os.path.isfile(self.inputpdb):
            self.HBPLUS.run(self.inputpdb, self.outputhbplus)
            self.assertTrue(os.path.isfile(self.outputhbplus))
        else:
            raise IOError("%s" % self.inputpdb)

    def test_generator_cif(self):
        if os.path.isfile(self.inputcif):
            self.HBPLUS.run(self.inputcif, self.outputhbplus)
            self.assertTrue(os.path.isfile(self.outputhbplus))
        else:
            raise IOError("%s" % self.inputcif)

    def test_run_hbplus_pdb(self):
        if os.path.isfile(self.inputpdb):
            self.run_hbplus(self.inputpdb, self.outputhbplus)
            self.assertTrue(os.path.isfile(self.outputhbplus))
        else:
            raise IOError("%s" % self.inputpdb)

    def test_run_hbplus_cif(self):
        if os.path.isfile(self.inputcif):
            self.run_hbplus(self.inputcif, self.outputhbplus)
            self.assertTrue(os.path.isfile(self.outputhbplus))
        else:
            raise IOError("%s" % self.inputcif)

    def test_parser_keys(self):
        self.assertListEqual([k for k in self.parser(self.outputhbplus).CHAIN_A.unique()],
                             ['A', 'B'])
        self.assertListEqual([k for k in self.parser(self.outputhbplus).CHAIN_D.unique()],
                             ['A', 'B'])

    def test_reader_data(self):
        data = self.HBPLUS.read(self.outputhbplus)
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

    def test_reader_default_excluded(self):
        data = self.HBPLUS.read(self.outputhbplus)
        self.assertIn("NUM_AAS", data)
        self.assertIn("DIST_CA-CA", data)
        self.assertIn("DIST_H-A", data)
        self.assertIn("ANGLE_D-H-A", data)
        self.assertIn("ANGLE_H-A-AA", data)
        self.assertIn("ANGLE_D-A-AA", data)

    def test_reader_new_excluded(self):
        data = self.HBPLUS.read(self.outputhbplus,
                                excluded_cols=self.excluded)
        self.assertIn("NUM_AAS", data)
        self.assertIn("DIST_CA-CA", data)
        self.assertIn("DIST_H-A", data)
        self.assertIn("ANGLE_D-H-A", data)
        self.assertIn("ANGLE_H-A-AA", data)
        self.assertIn("ANGLE_D-A-AA", data)
        excluded = ("NUM_AAS", "DIST_CA-CA", "ANGLE_D-H-A",
                    "DIST_H-A", "ANGLE_H-A-AA", "ANGLE_D-A-AA")
        data = self.HBPLUS.read(self.outputhbplus, excluded_cols=excluded)
        self.assertNotIn("NUM_AAS", data)
        self.assertNotIn("DIST_CA-CA", data)
        self.assertNotIn("DIST_H-A", data)
        self.assertNotIn("ANGLE_D-H-A", data)
        self.assertNotIn("ANGLE_H-A-AA", data)
        self.assertNotIn("ANGLE_D-A-AA", data)

    def test_filter_chain_donor(self):
        data = self.HBPLUS.read(self.outputhbplus,
                                excluded_cols=self.excluded)
        data = self.filter_hbplus(data, chain_D=('A',))
        self.assertNotIn("B", data.CHAIN_D.unique())

    def test_filter_chain_acceptor(self):
        data = self.HBPLUS.read(self.outputhbplus,
                                excluded_cols=self.excluded)
        data = self.filter_hbplus(data, chain_A=('B',))
        self.assertNotIn("A", data.CHAIN_A.unique())

    def test_filter_res_donor(self):
        data = self.HBPLUS.read(self.outputhbplus,
                                excluded_cols=self.excluded)
        data = self.filter_hbplus(data, res_D=('123',))
        self.assertNotIn('119', data.RES_D.unique())

    def test_filter_res_acceptor(self):
        data = self.HBPLUS.read(self.outputhbplus,
                                excluded_cols=self.excluded)
        data = self.filter_hbplus(data, res_A=('127',))
        self.assertNotIn('119', data.RES_A.unique())

    def test_generate_output_filename(self):
        output = self.generate_output_filename(self.inputcif)
        self.assertEqual(output, self.inputcif.replace('.cif', '.hb2'))

    def test_generate_output_filename_hydro_pdb_out(self):
        output = self.generate_output_filename(self.inputcif, hydro_pdb_out=True)
        self.assertEqual(output, self.inputcif.replace('.cif', '.h.pdb'))

    def test_generator_pdb_hydrogen(self):
        if os.path.isfile(self.inputpdb):
            self.HBPLUS.run(self.inputcif, self.outputhbplus_h,
                            hydro_pdb_out=True, run_clean=False, overwrite=True)
            self.assertTrue(os.path.isfile(self.outputhbplus_h))
        else:
            raise IOError("%s" % self.inputpdb)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("prointvar").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestHBPLUS)
    unittest.TextTestRunner(verbosity=2).run(suite)
