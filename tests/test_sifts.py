#!/local/bin/python
# -*- coding: utf-8 -*-


import json
import os
import sys
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

from prointvar.sifts import (SIFTSreader, parse_sifts_residues_from_file,
                             get_sifts_selected_from_table,
                             parse_sifts_regions_from_file,
                             parse_sifts_dbs_from_file)

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
class TestSIFTS(unittest.TestCase):
    """Test the DSSP parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.pdbid = '2pah'
        self.inputsifts = "{}{}{}.xml".format(c.db_root, c.db_sifts, self.pdbid)
        self.emptyfile = "{}{}{}.tmp".format(c.db_root, c.db_tmp, self.pdbid)
        self.notfound = ""
        self.excluded = ()

        self.parser = parse_sifts_residues_from_file
        self.parser_regions = parse_sifts_regions_from_file
        self.parser_dbs = parse_sifts_dbs_from_file
        self.reader = SIFTSreader
        self.filter = get_sifts_selected_from_table

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.inputsifts = None
        self.emptyfile = None
        self.notfound = None
        self.excluded = None
        self.parser = None
        self.parser_regions = None
        self.parser_dbs = None
        self.reader = None
        self.filter = None

    def test_file_not_found_reader(self):
        with self.assertRaises(IOError):
            self.reader(self.notfound)

    def test_file_not_found_parser(self):
        with self.assertRaises(IOError):
            self.parser(self.notfound)

    def test_empty_file_reader(self):
        with self.assertRaises(IOError):
            open(self.emptyfile, 'w').close()
            self.reader(self.emptyfile).read()
            os.remove(self.emptyfile)

    def test_reader_cif_verbose(self):
        with captured_output() as (out, err):
            self.reader(self.inputsifts, verbose=True).read()
        # This can go inside or outside the `with` block
        output = out.getvalue().strip()
        self.assertEqual(output, 'Parsing SIFTS residues from lines...')

    def test_parser_keys(self):
        self.assertListEqual([k for k in self.parser(self.inputsifts).PDB_entityId.unique()],
                             ['A', 'B'])

    def test_reader_keys(self):
        self.assertListEqual([k for k in self.reader(self.inputsifts).read().PDB_entityId.unique()],
                             ['A', 'B'])

    def test_reader_sifts_data(self):
        reader = self.reader(self.inputsifts)
        data = reader.read()
        # RESIDUES
        self.assertEqual(data.loc[0, 'PDB_Annotation'], 'Observed')

    def test_reader_residues_sifts_data(self):
        reader = self.reader(self.inputsifts)
        data = reader.residues()
        # RESIDUES
        self.assertEqual(data.loc[0, 'PDB_Annotation'], 'Observed')

    def test_reader_regions_sifts_data(self):
        reader = self.reader(self.inputsifts)
        data = reader.regions()
        # REGIONS
        self.assertEqual(data['B']['PDB']['1']['dbAccessionId'], self.pdbid)
        self.assertEqual(data['B']['PDB']['1']['start'], 1)
        self.assertEqual(data['B']['PDB']['1']['end'], 335)
        self.assertEqual(data['B']['PDB']['1']['dbCoordSys'], 'PDBresnum')

    def test_reader_dbs_sifts_data(self):
        reader = self.reader(self.inputsifts)
        data = reader.dbs()
        # DB versions
        self.assertEqual(data['UniProt']['dbSource'], 'UniProt')
        self.assertEqual(data['UniProt']['dbCoordSys'], 'UniProt')
        self.assertEqual(data['UniProt']['dbVersion'], '2017.03')

    def test_reader_sifts_to_json_pretty(self):
        reader = self.reader(self.inputsifts)
        reader.read()
        data = reader.to_json()
        self.assertEqual(json.loads(data)[0]['PDB_dbChainId'], 'A')

    def test_reader_sifts_to_json(self):
        reader = self.reader(self.inputsifts)
        reader.regions()
        data = reader.to_json(pretty=False)
        self.assertEqual(json.loads(data)['B']['Pfam']['1']['dbCoordSys'], 'UniProt')

    def test_reader_default_excluded(self):
        reader = self.reader(self.inputsifts)
        data = reader.regions()
        keys = [k for k in data['A']]
        self.assertNotIn("InterPro", keys)
        self.assertNotIn("GO", keys)
        self.assertNotIn("EC", keys)

    def test_reader_new_excluded(self):
        reader = self.reader(self.inputsifts)
        data = reader.regions(excluded=self.excluded)
        keys = [k for k in data['A']]
        self.assertIn("InterPro", keys)
        self.assertIn("GO", keys)
        self.assertIn("EC", keys)

    def test_reader_sifts_add_regions(self):
        reader = self.reader(self.inputsifts)
        data = reader.read(add_regions=True)
        self.assertEqual(data.loc[0, 'PDB_regionId'], '1')
        self.assertEqual(data.loc[0, 'PDB_regionStart'], 1)
        self.assertEqual(data.loc[0, 'PDB_regionEnd'], 335)

    def test_reader_sifts_add_dbs(self):
        reader = self.reader(self.inputsifts)
        data = reader.read(add_dbs=True)
        self.assertIn('PDB_dbVersion', data)
        self.assertEqual(data.loc[0, 'PDB_dbVersion'], '10.17')

    def test_filter_uniprot_id(self):
        reader = self.reader(self.inputsifts)
        reader.read(add_regions=False, add_dbs=False)
        data = self.filter(reader.data, uniprot=('P00439', ))
        self.assertEqual("P00439", data.loc[0, 'UniProt_dbAccessionId'])

    def test_filter_chain(self):
        reader = self.reader(self.inputsifts)
        reader.read(add_regions=False, add_dbs=False)
        data = self.filter(reader.data, chain=('A', ))
        self.assertIn("A", data.PDB_entityId.unique())
        self.assertNotIn("B", data.PDB_entityId.unique())

    def test_parser_regions(self):
        data = self.parser_regions(self.inputsifts)
        # REGIONS
        self.assertEqual(data['B']['PDB']['1']['dbAccessionId'], self.pdbid)
        self.assertEqual(data['B']['PDB']['1']['start'], 1)
        self.assertEqual(data['B']['PDB']['1']['end'], 335)
        self.assertEqual(data['B']['PDB']['1']['dbCoordSys'], 'PDBresnum')

    def test_parser_dbs(self):
        data = self.parser_dbs(self.inputsifts)
        # DB versions
        self.assertEqual(data['UniProt']['dbSource'], 'UniProt')
        self.assertEqual(data['UniProt']['dbCoordSys'], 'UniProt')
        self.assertEqual(data['UniProt']['dbVersion'], '2017.03')


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSIFTS)
    unittest.TextTestRunner(verbosity=2).run(suite)
