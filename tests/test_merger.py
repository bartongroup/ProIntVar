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

from prointvar.dssp import DSSPreader, DSSPgenerator
from prointvar.sifts import SIFTSreader
from prointvar.mmcif import MMCIFreader, MMCIFwriter, get_mmcif_selected_from_table

from prointvar.merger import (TableMerger, table_merger,
                              mmcif_dssp_table_merger, mmcif_sifts_table_merger,
                              dssp_sifts_table_merger, table_generator, dssp_dssp_table_merger)

from prointvar.config import config as c
root = os.path.abspath(os.path.dirname(__file__))
c.db_root = "{}/testdata/".format(root)


@patch("prointvar.config.config.db_root", c.db_root)
class TestMerger(unittest.TestCase):
    """Test the Merger methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.pdbid = TestMerger.pdbid
        self.inputcif = TestMerger.inputcif
        self.inputdssp = TestMerger.inputdssp

        self.inputbiocif = TestMerger.inputbiocif
        self.inputbiodssp = TestMerger.inputbiodssp

        self.inputsifts = TestMerger.inputsifts

        self.mmcif = TestMerger.mmcif
        self.mmcif_bio = TestMerger.mmcif_bio
        self.dssp = TestMerger.dssp
        self.dssp_bio = TestMerger.dssp_bio
        self.dssp_unbound = TestMerger.dssp_unbound
        self.sifts = TestMerger.sifts

        self.mmcif_sifts = mmcif_sifts_table_merger
        self.mmcif_dssp = mmcif_dssp_table_merger
        self.dssp_sifts = dssp_sifts_table_merger
        self.dssp_dssp = dssp_dssp_table_merger
        self.merger = TableMerger
        self.table_merger = table_merger
        self.generator = table_generator

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.inputcif = None
        self.inputdssp = None
        self.inputbiocif = None
        self.inputbiodssp = None
        self.inputsifts = None

        self.mmcif = None
        self.mmcif_bio = None
        self.dssp = None
        self.dssp_bio = None
        self.dssp_unbound = None
        self.sifts = None

        self.mmcif_sifts = None
        self.mmcif_dssp = None
        self.dssp_sifts = None
        self.dssp_dssp = None
        self.merger = None
        self.table_merger = None
        self.generator = None

    @classmethod
    def setUpClass(cls):
        # to be run only once
        super(TestMerger, cls).setUpClass()

        cls.pdbid = '2pah'
        cls.inputcif = "{}{}{}.cif".format(c.db_root, c.db_cif, cls.pdbid)
        cls.inputdssp = "{}{}{}.dssp".format(c.db_root, c.db_dssp_generated, cls.pdbid)

        cls.inputbiocif = "{}{}{}.cif".format(c.db_root, c.db_cif_biounit, cls.pdbid)
        cls.inputbiodssp = "{}{}{}_bio.dssp".format(c.db_root, c.db_dssp_generated, cls.pdbid)

        cls.inputsifts = "{}{}{}.xml".format(c.db_root, c.db_sifts_xml, cls.pdbid)

        d = MMCIFreader(cls.inputcif)
        cls.mmcif = d.atoms(add_res_full=True)
        cls.mmcif = get_mmcif_selected_from_table(cls.mmcif, atom=('CA',))
        cls.mmcif = TestMerger.mmcif

        d = MMCIFreader(cls.inputbiocif)
        cls.mmcif_bio = d.atoms(add_res_full=True)
        cls.mmcif_bio = get_mmcif_selected_from_table(cls.mmcif_bio, atom=('CA',))

        d = DSSPreader(cls.inputdssp)
        cls.dssp = d.residues(add_rsa_class=True, add_ss_reduced=True)
        d = DSSPreader(cls.inputbiodssp)
        cls.dssp_bio = d.residues(add_rsa_class=True, add_ss_reduced=True)

        # chain A only: unbound
        cls.outputcif_A = "{}{}{}_A.cif".format(c.db_root, c.db_cif, cls.pdbid)
        w = MMCIFwriter(inputfile=cls.inputcif,
                        outputfile=cls.outputcif_A)
        w.run(chain=('A',))
        cls.outputdssp_A = "{}{}{}_A.dssp".format(c.db_root, c.db_dssp_generated, cls.pdbid)
        d = DSSPgenerator(inputfile=cls.outputcif_A,
                          outputfile=cls.outputdssp_A)
        d.run()
        os.remove(cls.outputcif_A)
        d = DSSPreader(cls.outputdssp_A)
        cls.dssp_unbound = d.residues(add_full_chain=True, add_ss_reduced=True,
                                      add_rsa=True, add_rsa_class=True)
        os.remove(cls.outputdssp_A)

        d = SIFTSreader(cls.inputsifts)
        cls.sifts = d.read(add_regions=True, add_dbs=False)

    @classmethod
    def tearDownClass(cls):

        cls.pdbid = None
        cls.inputcif = None
        cls.inputdssp = None
        cls.inputbiocif = None
        cls.inputbiodssp = None
        cls.inputsifts = None

        cls.mmcif = None
        cls.mmcif_bio = None
        cls.dssp = None
        cls.dssp_bio = None
        cls.sifts = None

    def test_mmcif_dssp_merger(self):
        table = self.mmcif_dssp(self.mmcif, self.dssp)
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertNotIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertNotIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[0, 'label_atom_id'])
        self.assertEqual('A', table.loc[0, 'label_asym_id'])
        self.assertEqual('118', table.loc[0, 'RES'])
        self.assertEqual('V', table.loc[0, 'AA'])

    def test_mmcif_dssp_bio_merger(self):
        table = self.mmcif_dssp(self.mmcif_bio, self.dssp_bio)
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertNotIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertNotIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[329, 'label_atom_id'])
        self.assertEqual('AA', table.loc[329, 'label_asym_id'])
        self.assertEqual('118', table.loc[329, 'RES'])
        self.assertEqual('V', table.loc[329, 'AA'])

    def test_mmcif_sifts_merger(self):
        table = self.mmcif_sifts(self.mmcif, self.sifts)
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertNotIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertNotIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[0, 'label_atom_id'])
        self.assertEqual('A', table.loc[0, 'label_asym_id'])
        self.assertEqual('118', table.loc[0, 'PDB_dbResNum'])
        self.assertEqual('VAL', table.loc[0, 'PDB_dbResName'])

    def test_mmcif_sifts_bio_merger(self):
        table = self.mmcif_sifts(self.mmcif_bio, self.sifts)
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertNotIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertNotIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[329, 'label_atom_id'])
        self.assertEqual('AA', table.loc[329, 'label_asym_id'])
        self.assertEqual('118', table.loc[329, 'PDB_dbResNum'])
        self.assertEqual('VAL', table.loc[329, 'PDB_dbResName'])

    def test_dssp_sifts_merger(self):
        table = self.dssp_sifts(self.dssp, self.sifts)
        # Chain level
        self.assertNotIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertNotIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('A', table.loc[0, 'PDB_entityId'])
        self.assertEqual('118', table.loc[0, 'RES'])
        self.assertEqual('VAL', table.loc[0, 'PDB_dbResName'])

    def test_dssp_sifts_bio_merger(self):
        table = self.dssp_sifts(self.dssp, self.sifts)
        # Chain level
        self.assertNotIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertNotIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('B', table.loc[329, 'PDB_entityId'])
        self.assertEqual('118', table.loc[329, 'RES'])
        self.assertEqual('VAL', table.loc[329, 'PDB_dbResName'])

    def test_table_merger(self):
        table = self.merger(self.mmcif, self.dssp, self.sifts).merge()
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[0, 'label_atom_id'])
        self.assertEqual('A', table.loc[0, 'label_asym_id'])
        self.assertEqual('118', table.loc[0, 'RES'])
        self.assertEqual('VAL', table.loc[0, 'PDB_dbResName'])

    def test_table_merger_method(self):
        table = self.table_merger(self.mmcif, self.dssp, self.sifts)
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[0, 'label_atom_id'])
        self.assertEqual('A', table.loc[0, 'label_asym_id'])
        self.assertEqual('118', table.loc[0, 'RES'])
        self.assertEqual('VAL', table.loc[0, 'PDB_dbResName'])
        self.assertEqual('V', table.loc[0, 'UniProt_dbResName'])

    def test_table_merger_method_bio(self):
        table = self.table_merger(self.mmcif_bio, self.dssp_bio, self.sifts)
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[329, 'label_atom_id'])
        self.assertEqual('AA', table.loc[329, 'label_asym_id'])
        self.assertEqual('118', table.loc[329, 'RES'])
        self.assertEqual('VAL', table.loc[329, 'PDB_dbResName'])
        self.assertEqual('V', table.loc[329, 'UniProt_dbResName'])

    def test_table_generator(self):
        mmcif_table, dssp_table, sifts_table = \
            self.generator(uniprot_id=None, pdb_id=self.pdbid, chain=None,
                           res=None, site=None, atom=('CA',), lines=None,
                           bio=False, add_dssp=True)

        table = self.table_merger(mmcif_table, dssp_table, sifts_table)
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[0, 'label_atom_id'])
        self.assertEqual('A', table.loc[0, 'label_asym_id'])
        self.assertEqual('118', table.loc[0, 'RES'])
        self.assertEqual('VAL', table.loc[0, 'PDB_dbResName'])
        self.assertEqual('V', table.loc[0, 'UniProt_dbResName'])

    def test_table_generator_bio(self):
        mmcif_table, dssp_table, sifts_table = \
            self.generator(uniprot_id=None, pdb_id=self.pdbid, chain=None,
                           res=None, site=None, atom=('CA',), lines=None,
                           bio=True, add_dssp=True)

        table = self.table_merger(mmcif_table, dssp_table, sifts_table)
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[329, 'label_atom_id'])
        self.assertEqual('AA', table.loc[329, 'label_asym_id'])
        self.assertEqual('118', table.loc[329, 'RES'])
        self.assertEqual('VAL', table.loc[329, 'PDB_dbResName'])
        self.assertEqual('V', table.loc[329, 'UniProt_dbResName'])

    def test_table_generator_full_dssp(self):
        mmcif_table, dssp_table, sifts_table = \
            self.generator(uniprot_id=None, pdb_id=self.pdbid, chain=None,
                           res=None, site=None, atom=('CA',), lines=None,
                           bio=False, add_dssp=True, dssp_unbound=True)

        table = self.table_merger(mmcif_table, dssp_table, sifts_table)
        self.assertIn('CHAIN_FULL', table)
        self.assertIn('CHAIN_FULL_UNB', table)
        self.assertIn('RSA', table)
        self.assertIn('RSA_UNB', table)
        # values
        self.assertEqual(80.488, table.loc[648, 'RSA'])
        self.assertEqual(80.488, table.loc[648, 'RSA_UNB'])
        self.assertEqual(63.314, table.loc[649, 'RSA'])
        self.assertEqual(74.556, table.loc[649, 'RSA_UNB'])

    def test_table_merger_get_filename(self):
        filename = self.merger()._get_filename(pdb_id=self.pdbid, bio=True)
        self.assertTrue(os.path.exists(filename))

    def test_table_merger_private_dump(self):
        filename = self.merger()._get_filename(pdb_id=self.pdbid, bio=True)
        if os.path.exists(filename):
            os.remove(filename)
        t = self.merger(self.mmcif_bio, self.dssp_bio, self.sifts)
        t.merge()
        # alternatively
        # t = self.merger(self.mmcif_bio, self.dssp_bio, self.sifts).merge(outputfile=filename)
        t._dump_merged_table(outputfile=filename)

    def test_table_merger_run(self):
        table = self.merger().run(pdb_id=self.pdbid, atom=('CA',), bio=True)
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[329, 'label_atom_id'])
        self.assertEqual('AA', table.loc[329, 'label_asym_id'])
        self.assertEqual('118', table.loc[329, 'RES'])
        self.assertEqual('VAL', table.loc[329, 'PDB_dbResName'])
        self.assertEqual('V', table.loc[329, 'UniProt_dbResName'])

    def test_table_merger_private_load(self):
        filename = self.merger()._get_filename(pdb_id=self.pdbid, bio=True)
        if not os.path.exists(filename):
            t = self.merger(self.mmcif_bio, self.dssp_bio, self.sifts)
            t.merge()
            t._dump_merged_table(outputfile=filename)
        table = self.merger()._load_merged_table(filename)
        self.assertIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)

    def test_table_merger_load(self):
        filename = self.merger()._get_filename(pdb_id=self.pdbid, bio=True)
        if not os.path.exists(filename):
            t = self.merger(self.mmcif_bio, self.dssp_bio, self.sifts)
            t.merge()
            t._dump_merged_table(outputfile=filename)
        t = self.merger()
        table = t.load(pdb_id=self.pdbid, bio=True)
        # Chain level
        self.assertIn('label_asym_id', table)
        self.assertIn('CHAIN_FULL', table)
        self.assertIn('PDB_entityId', table)
        # Res level
        self.assertIn('label_seq_id_full', table)
        self.assertIn('RES', table)
        self.assertIn('PDB_dbResNum', table)
        # values
        self.assertEqual('CA', table.loc[329, 'label_atom_id'])
        self.assertEqual('AA', table.loc[329, 'label_asym_id'])
        self.assertEqual('118', table.loc[329, 'RES'])
        self.assertEqual('VAL', table.loc[329, 'PDB_dbResName'])
        self.assertEqual('V', table.loc[329, 'UniProt_dbResName'])

    def test_table_merger_dssp_dssp(self):
        table = self.dssp_dssp(self.dssp, self.dssp_unbound)
        table = table.loc[table['RES'] == '438']
        self.assertIn('RES_UNB', table)
        self.assertIn('SS_UNB', table)
        self.assertIn('ACC_UNB', table)
        self.assertEqual(36.943, table.loc[314, 'RSA'])
        self.assertEqual(79.618, table.loc[314, 'RSA_UNB'])
        self.assertEqual(58.0, table.loc[314, 'ACC'])
        self.assertEqual(125.0, table.loc[314, 'ACC_UNB'])

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMerger)
    unittest.TextTestRunner(verbosity=2).run(suite)
