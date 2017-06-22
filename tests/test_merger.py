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
from prointvar.pdbx import PDBXreader, PDBXwriter, get_mmcif_selected_from_table
from prointvar.arpeggio import ARPEGGIOreader

from prointvar.merger import (TableMerger, table_merger,
                              mmcif_dssp_table_merger, mmcif_sifts_table_merger,
                              dssp_sifts_table_merger, table_generator,
                              dssp_dssp_table_merger, contacts_mmcif_table_merger,
                              load_merged_table, dump_merged_table)

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
        self.contacts = TestMerger.contacts

        self.mmcif_sifts = mmcif_sifts_table_merger
        self.mmcif_dssp = mmcif_dssp_table_merger
        self.dssp_sifts = dssp_sifts_table_merger
        self.dssp_dssp = dssp_dssp_table_merger
        self.contacts_mmcif = contacts_mmcif_table_merger
        self.merger = TableMerger
        self.table_merger = table_merger
        self.generator = table_generator
        self.dump_merged_table = dump_merged_table
        self.load_merged_table = load_merged_table

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
        self.contacts= None

        self.mmcif_sifts = None
        self.mmcif_dssp = None
        self.dssp_sifts = None
        self.dssp_dssp = None
        self.contacts_mmcif = None
        self.merger = None
        self.table_merger = None
        self.generator = None
        self.dump_merged_table = None
        self.load_merged_table = None

    @classmethod
    def setUpClass(cls):
        # to be run only once
        super(TestMerger, cls).setUpClass()

        cls.pdbid = '2pah'
        cls.inputcif = "{}{}{}.cif".format(c.db_root, c.db_pdbx, cls.pdbid)
        cls.inputdssp = "{}{}{}.dssp".format(c.db_root, c.db_dssp, cls.pdbid)

        cls.inputbiocif = "{}{}{}_bio.cif".format(c.db_root, c.db_pdbx, cls.pdbid)
        cls.inputbiodssp = "{}{}{}_bio.dssp".format(c.db_root, c.db_dssp, cls.pdbid)

        cls.inputsifts = "{}{}{}.xml".format(c.db_root, c.db_sifts, cls.pdbid)

        cls.inputcontacts = "{}{}{}.contacts" \
                            "".format(c.db_root, c.db_contacts, cls.pdbid)

        d = PDBXreader(cls.inputcif)
        cls.mmcif = d.atoms(add_res_full=True)
        cls.mmcif = get_mmcif_selected_from_table(cls.mmcif, atom=('CA',))
        cls.mmcif = TestMerger.mmcif

        d = PDBXreader(cls.inputbiocif)
        cls.mmcif_bio = d.atoms(add_res_full=True)
        cls.mmcif_bio = get_mmcif_selected_from_table(cls.mmcif_bio, atom=('CA',))

        d = DSSPreader(cls.inputdssp)
        cls.dssp = d.residues(add_rsa_class=True, add_ss_reduced=True)
        d = DSSPreader(cls.inputbiodssp)
        cls.dssp_bio = d.residues(add_rsa_class=True, add_ss_reduced=True)

        # chain A only: unbound
        cls.outputcif_A = "{}{}{}_A.cif".format(c.db_root, c.db_pdbx, cls.pdbid)
        w = PDBXwriter(inputfile=cls.inputcif,
                       outputfile=cls.outputcif_A)
        w.run(chain=('A',))
        cls.outputdssp_A = "{}{}{}_A.dssp".format(c.db_root, c.db_dssp, cls.pdbid)
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

        r = ARPEGGIOreader(cls.inputcontacts)
        cls.contacts = r.contacts(residue_agg=True, agg_method="minimum",
                                  collapsed_cont=True, col_method="full",
                                  int_filter=True, int_mode='inter-chain',
                                  parse_special=True)

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
        mmcif_table, dssp_table, sifts_table, contacts_table = \
            self.generator(uniprot_id=None, pdb_id=self.pdbid, chain=None,
                           res=None, site=None, atom=('CA',), lines=None,
                           bio=False, dssp=True)

        table = self.table_merger(mmcif_table, dssp_table, sifts_table, contacts_table)
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
        mmcif_table, dssp_table, sifts_table, contacts_table = \
            self.generator(uniprot_id=None, pdb_id=self.pdbid, chain=None,
                           res=None, site=None, atom=('CA',), lines=None,
                           bio=True, dssp=True)

        table = self.table_merger(mmcif_table, dssp_table, sifts_table, contacts_table)
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
        mmcif_table, dssp_table, sifts_table, contacts_table = \
            self.generator(uniprot_id=None, pdb_id=self.pdbid, chain=None,
                           res=None, site=None, atom=('CA',), lines=None,
                           bio=False, dssp=True, dssp_unbound=True)

        table = self.table_merger(mmcif_table, dssp_table, sifts_table, contacts_table)
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

    def test_table_merger_dump(self):
        filename = self.merger()._get_filename(pdb_id=self.pdbid, bio=True,
                                               lines=('ATOM', ))
        t = self.merger(self.mmcif_bio, self.dssp_bio, self.sifts)
        t.merge()
        self.dump_merged_table(t.merged_table, outputfile=filename)
        self.assertTrue(os.path.exists(filename))
        os.remove(filename)

    def test_table_merger_run(self):
        table = self.merger().run(pdb_id=self.pdbid, atom=('CA',), bio=True, dssp=True)
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
        filename = self.merger()._get_filename(pdb_id=self.pdbid, atom=('CA',),
                                               bio=True, dssp=True)
        os.remove(filename)

    def test_table_merger_load(self):
        filename = self.merger()._get_filename(pdb_id=self.pdbid, bio=True)
        if not os.path.exists(filename):
            t = self.merger(self.mmcif_bio, self.dssp_bio, self.sifts)
            t.merge()
            self.dump_merged_table(t.merged_table, outputfile=filename)
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

    def test_contacts_mmcif_merger(self):
        table = self.contacts_mmcif(self.contacts, self.mmcif, suffix='A')
        self.assertIn('CHAIN_A', list(table))
        self.assertIn('label_asym_id_A', list(table))
        self.assertIn('label_seq_id_full_A', list(table))
        self.assertTrue('366', table.loc[0, 'RES_FULL_A'])
        self.assertTrue(1.694, table.loc[0, 'VDW_DIST'])

    def test_contacts_mmcif_merger_both_sides(self):
        table = self.contacts_mmcif(self.contacts, self.mmcif, suffix='A')
        table = self.contacts_mmcif(table, self.mmcif, suffix='B')
        self.assertIn('CHAIN_A', list(table))
        self.assertIn('CHAIN_B', list(table))
        self.assertIn('label_asym_id_A', list(table))
        self.assertIn('label_asym_id_B', list(table))
        self.assertIn('label_seq_id_full_A', list(table))
        self.assertIn('label_seq_id_full_B', list(table))
        self.assertTrue('366', table.loc[0, 'RES_FULL_A'])
        self.assertTrue(1.694, table.loc[0, 'VDW_DIST'])

    def test_table_generator_contacts_mmcif(self):
        mmcif_table, dssp_table, sifts_table, contacts_table = \
            self.generator(uniprot_id=None, pdb_id=self.pdbid, chain=None,
                           res=None, site=None, atom=('CA',), lines=None,
                           bio=False, sifts=False, dssp=False, dssp_unbound=False,
                           contacts=True)
        table = self.contacts_mmcif(contacts_table, mmcif_table, suffix='A')
        self.assertIn('CHAIN_A', list(table))
        self.assertIn('label_asym_id_A', list(table))
        self.assertIn('label_seq_id_full_A', list(table))
        self.assertTrue('366', table.loc[0, 'RES_FULL_A'])
        self.assertTrue(1.694, table.loc[0, 'VDW_DIST'])

    def test_table_merger_contacts_mmcif(self):
        mmcif_table, dssp_table, sifts_table, contacts_table = \
            self.generator(uniprot_id=None, pdb_id=self.pdbid, chain=None,
                           res=None, site=None, atom=('CA',), lines=None,
                           bio=False, sifts=False, dssp=False, dssp_unbound=False,
                           contacts=True)

        table = self.table_merger(mmcif_table, dssp_table, sifts_table, contacts_table)
        self.assertIn('CHAIN_A', list(table))
        self.assertIn('CHAIN_B', list(table))
        self.assertIn('label_asym_id_A', list(table))
        self.assertIn('label_asym_id_B', list(table))
        self.assertIn('label_seq_id_full_A', list(table))
        self.assertIn('label_seq_id_full_B', list(table))
        self.assertTrue('366', table.loc[0, 'RES_FULL_A'])
        self.assertTrue(1.694, table.loc[0, 'VDW_DIST'])

    def test_table_merger_contacts_mmcif_bio_sifts_dssp_no_residue_agg(self):
        mmcif_table, dssp_table, sifts_table, contacts_table = \
            self.generator(uniprot_id=None, pdb_id=self.pdbid, chain=None,
                           res=None, site=None, atom=None, lines=None,
                           bio=True, sifts=True, dssp=True, dssp_unbound=True,
                           contacts=True)

        table = self.table_merger(mmcif_table, dssp_table, sifts_table, contacts_table)
        self.assertIn('CHAIN_A', list(table))
        self.assertIn('CHAIN_B', list(table))
        self.assertIn('new_asym_id_A', list(table))
        self.assertIn('new_asym_id_B', list(table))
        self.assertIn('label_asym_id_A', list(table))
        self.assertIn('label_asym_id_B', list(table))
        self.assertIn('label_seq_id_full_A', list(table))
        self.assertIn('label_seq_id_full_B', list(table))
        self.assertTrue('844', table.loc[0, 'RES_FULL_A'])
        self.assertTrue('844', table.loc[0, 'new_seq_id_A'])
        self.assertTrue('203', table.loc[0, 'label_seq_id_full_A'])
        self.assertTrue(3.109, table.loc[0, 'DIST'])

    def test_table_merger_contacts_mmcif_bio_sifts_dssp_residue_agg(self):
        mmcif_table, dssp_table, sifts_table, contacts_table = \
            self.generator(uniprot_id=None, pdb_id=self.pdbid, chain=None,
                           res=None, site=None, atom=None, lines=None,
                           bio=True, sifts=True, dssp=True, dssp_unbound=True,
                           contacts=True, residue_agg=True)

        table = self.table_merger(mmcif_table, dssp_table, sifts_table, contacts_table)
        self.assertIn('CHAIN_A', list(table))
        self.assertIn('CHAIN_B', list(table))
        self.assertIn('label_asym_id_A', list(table))
        self.assertIn('label_asym_id_B', list(table))
        self.assertIn('label_seq_id_full_A', list(table))
        self.assertIn('label_seq_id_full_B', list(table))
        self.assertIn('Cartn_x_A', list(table))
        self.assertIn('Cartn_x_A', list(table))
        self.assertIn('RSA_A', list(table))
        self.assertIn('RSA_B', list(table))
        self.assertIn('UniProt_dbAccessionId_A', list(table))
        self.assertIn('UniProt_dbAccessionId_B', list(table))
        self.assertTrue('366', table.loc[0, 'RES_FULL_A'])
        self.assertTrue(1.694, table.loc[0, 'VDW_DIST'])
        self.assertIn('Amide-Amide', table.loc[0, 'Int_Types'])

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMerger)
    unittest.TextTestRunner(verbosity=2).run(suite)
