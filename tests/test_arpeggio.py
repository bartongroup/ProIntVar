# -*- coding: utf-8 -*-


import os
import sys
import logging
import unittest
from unittest.mock import patch

from proteofav.structures import PDB, filter_structures

from prointvar.arpeggio import (parse_arpeggio_from_file,
                                add_arpeggio_res_split,
                                interaction_modes, residues_aggregation,
                                collapsed_contacts, ignore_consecutive_residues,
                                parse_arpeggio_spec_from_file, add_contact_info,
                                add_special_cont_types,
                                filter_arpeggio, run_arpeggio, ARPEGGIO)

from prointvar.config import config


@patch("prointvar.config.config", config)
class TestARPEGGIO(unittest.TestCase):
    """Test the ARPEGGIO parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.pdbid = '2pah'
        self.pdbid_small = '2rea'
        self.pdbid2 = "1ejg"
        self.inputpdb = os.path.join(os.path.dirname(__file__), "testdata",
                                     config.db_pdb, "{}.pdb".format(self.pdbid))
        self.inputpdb_fast = os.path.join(os.path.dirname(__file__), "testdata",
                                          config.db_pdb, "{}.pdb".format(self.pdbid_small))
        self.inputcif = os.path.join(os.path.dirname(__file__), "testdata",
                                     config.db_mmcif, "{}.cif".format(self.pdbid))
        self.inputarpeggio = os.path.join(os.path.dirname(__file__), "testdata",
                                          config.db_contacts, "{}.contacts".format(self.pdbid))
        self.input_amam = os.path.join(os.path.dirname(__file__), "testdata",
                                       config.db_contacts, "{}.amam".format(self.pdbid))
        self.input_amri = os.path.join(os.path.dirname(__file__), "testdata",
                                       config.db_contacts, "{}.amri".format(self.pdbid))
        self.input_ari = os.path.join(os.path.dirname(__file__), "testdata",
                                      config.db_contacts, "{}.ari".format(self.pdbid))
        self.input_ri = os.path.join(os.path.dirname(__file__), "testdata",
                                     config.db_contacts, "{}.ri".format(self.pdbid))
        self.inputarpeggio_fast = os.path.join(os.path.dirname(__file__), "testdata",
                                               config.db_contacts,
                                               "{}.contacts".format(self.pdbid_small))
        self.excluded = ("ENTRY_A", "ENTRY_B", "ENTITIES")
        self.parser = parse_arpeggio_from_file
        self.filter_arpeggio = filter_arpeggio
        self.add_arpeggio_res_split = add_arpeggio_res_split
        self.interaction_modes = interaction_modes
        self.residues_aggregation = residues_aggregation
        self.collapsed_contacts = collapsed_contacts
        self.ignore_consecutive = ignore_consecutive_residues
        self.parser_spec = parse_arpeggio_spec_from_file
        self.add_contact_info = add_contact_info
        self.add_special_cont_types = add_special_cont_types
        self.run_arpeggio = run_arpeggio
        self.ARPEGGIO = ARPEGGIO

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.pdbid_small = None
        self.pdbid2 = None
        self.inputpdb = None
        self.inputpdb_fast = None
        self.inputcif = None
        self.inputarpeggio = None
        self.input_amam = None
        self.input_amri = None
        self.input_ari = None
        self.input_ri = None
        self.inputarpeggio_fast = None

        self.emptyfile = None
        self.notfound = None
        self.excluded = None
        self.parser = None
        self.filter_arpeggio = None
        self.add_arpeggio_res_split = None
        self.interaction_modes = None
        self.residues_aggregation = None
        self.collapsed_contacts = None
        self.ignore_consecutive = None
        self.parser_spec = None
        self.add_contact_info = None
        self.add_special_cont_types = None
        self.run_arpeggio = None
        self.ARPEGGIO = None

        logging.disable(logging.NOTSET)

    def test_generator_pdb_exec(self):
        if os.path.isfile(self.inputpdb_fast):
            self.ARPEGGIO.run(self.inputpdb_fast, self.inputarpeggio_fast,
                              clean_output=True, overwrite=True)
            msg = ("Arpeggio execution failed: make sure the settings "
                   "are set properly in config.ini!")
            self.assertTrue(os.path.isfile(self.inputarpeggio_fast), msg)
            os.remove(self.inputarpeggio_fast)
        else:
            raise IOError("%s" % self.inputpdb_fast)

    # @unittest.expectedFailure
    def test_generator_pdb_exec_fail(self):

        inputpdb = os.path.join(os.path.dirname(__file__), "testdata",
                                config.db_pdb, "{}.pdb".format(self.pdbid2))
        inputarpeggio = os.path.join(os.path.dirname(__file__), "testdata",
                                     config.db_pdb, "{}.contacts".format(self.pdbid2))
        try:
            self.ARPEGGIO.run(inputpdb, inputarpeggio,
                              clean_output=True, overwrite=True)
        except (FileNotFoundError, OSError):
            # expected failure
            msg = "PDB with residues have missing atoms..."
            self.assertFalse(os.path.isfile(inputarpeggio), msg)

        inputpdb_new = os.path.join(os.path.dirname(__file__), "testdata",
                                    config.db_pdb, "{}_new.pdb".format(self.pdbid2))
        data = PDB.read(filename=inputpdb)
        data = filter_structures(data, remove_altloc=True, remove_hydrogens=True,
                                 reset_atom_id=True, remove_partial_res=True)

        PDB.write(table=data, filename=inputpdb_new, category="auth")

        self.ARPEGGIO.run(inputpdb_new, inputarpeggio,
                          clean_output=True, overwrite=True)
        self.assertTrue(os.path.isfile(inputarpeggio))
        os.remove(inputpdb_new)
        os.remove(inputarpeggio)
        os.remove(os.path.join(os.path.dirname(__file__), "testdata",
                               config.db_pdb, "{}.amam".format(self.pdbid2)))
        os.remove(os.path.join(os.path.dirname(__file__), "testdata",
                               config.db_pdb, "{}.amri".format(self.pdbid2)))
        os.remove(os.path.join(os.path.dirname(__file__), "testdata",
                               config.db_pdb, "{}.ari".format(self.pdbid2)))
        os.remove(os.path.join(os.path.dirname(__file__), "testdata",
                               config.db_pdb, "{}.ri".format(self.pdbid2)))

    def test_generator_pdb(self):
        if os.path.isfile(self.inputpdb):
            self.ARPEGGIO.run(self.inputpdb, self.inputarpeggio)
            self.assertTrue(os.path.isfile(self.inputarpeggio))
        else:
            raise IOError("%s" % self.inputpdb)

    def test_generator_cif(self):
        if os.path.isfile(self.inputcif):
            self.ARPEGGIO.run(self.inputcif, self.inputarpeggio)
            self.assertTrue(os.path.isfile(self.inputarpeggio))
        else:
            raise IOError("%s" % self.inputcif)

    def test_run_arpeggio_pdb(self):
        if os.path.isfile(self.inputpdb):
            self.run_arpeggio(self.inputpdb, self.inputarpeggio)
            self.assertTrue(os.path.isfile(self.inputarpeggio))
        else:
            raise IOError("%s" % self.inputpdb)

    def test_run_arpeggio_cif(self):
        if os.path.isfile(self.inputcif):
            self.run_arpeggio(self.inputcif, self.inputarpeggio)
            self.assertTrue(os.path.isfile(self.inputarpeggio))
        else:
            raise IOError("%s" % self.inputcif)

    def test_parser_keys(self):
        self.assertListEqual(sorted([k for k in
                                     self.parser(self.inputarpeggio).CHAIN_A.unique()]),
                             ['A', 'B'])
        self.assertListEqual(sorted([k for k in
                                     self.parser(self.inputarpeggio).CHAIN_B.unique()]),
                             ['A', 'B'])

    def test_reader_data(self):
        data = self.ARPEGGIO.read(self.inputarpeggio)
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

    def test_reader_default_excluded(self):
        data = self.ARPEGGIO.read(self.inputarpeggio)
        self.assertIn("ENTRY_A", data)
        self.assertIn("ENTRY_B", data)

    def test_reader_new_excluded(self):
        data = self.ARPEGGIO.read(self.inputarpeggio, excluded_cols=self.excluded)
        self.assertNotIn("ENTRY_A", data)
        self.assertNotIn("ENTRY_B", data)

    def test_filter_chain(self):
        data = self.ARPEGGIO.read(self.inputarpeggio, excluded_cols=self.excluded)
        data = self.filter_arpeggio(data, chain_B=('A',))
        self.assertNotIn("B", data.CHAIN_B.unique())

    def test_filter_res(self):
        data = self.ARPEGGIO.read(self.inputarpeggio, excluded_cols=self.excluded)
        data = self.filter_arpeggio(data, res_A=('374',))
        self.assertNotIn('119', data.RES_A.unique())

    def test_add_split_res(self):
        data = self.ARPEGGIO.read(self.inputarpeggio, add_res_split=False)
        self.assertNotIn('CHAIN_A', data)
        self.assertNotIn('CHAIN_B', data)
        data = self.add_arpeggio_res_split(data)
        self.assertIn('CHAIN_A', data)
        self.assertIn('CHAIN_B', data)
        self.assertEqual(data.loc[0, 'CHAIN_A'], 'B')

    def test_interaction_modes(self):
        data = self.ARPEGGIO.read(self.inputarpeggio, excluded_cols=self.excluded)
        data = self.filter_arpeggio(data, int_filter=True, int_mode='inter-chain')
        self.assertEqual(285, len(data))
        self.assertNotEqual(data.loc[0, 'CHAIN_A'], data.loc[0, 'CHAIN_B'])
        data = self.ARPEGGIO.read(self.inputarpeggio, excluded_cols=self.excluded)
        data = self.filter_arpeggio(data, int_filter=True, int_mode='intra-chain')
        self.assertEqual(27135, len(data))
        self.assertEqual(data.loc[0, 'CHAIN_A'], data.loc[0, 'CHAIN_B'])
        data = self.ARPEGGIO.read(self.inputarpeggio, excluded_cols=self.excluded)
        data = self.interaction_modes(data, int_mode='inter-chain')
        self.assertEqual(285, len(data))
        self.assertEqual(data.loc[0, 'RES_FULL_A'], '368')
        self.assertEqual(data.loc[0, 'RES_FULL_B'], '368')

    def test_collapsed_contacts(self):
        data = self.ARPEGGIO.read(self.inputarpeggio, excluded_cols=self.excluded)
        data = self.filter_arpeggio(data, collapsed_cont=True, col_method='full')
        self.assertNotIn('IONIC', list(data))
        self.assertIn('Int_Types', list(data))
        self.assertEqual('Polar-Bond, VDW-Proximal',
                         ', '.join(sorted(list(data.loc[3, 'Int_Types']))))
        data = self.ARPEGGIO.read(self.inputarpeggio, excluded_cols=self.excluded)
        data = self.collapsed_contacts(data, col_method='minimal')
        self.assertNotIn('IONIC', list(data))
        self.assertIn('Int_Types', list(data))
        self.assertEqual('Polar-Bond', ', '.join(list(data.loc[3, 'Int_Types'])))

    def test_residues_agg(self):
        data = self.ARPEGGIO.read(self.inputarpeggio, excluded_cols=self.excluded)
        data = self.filter_arpeggio(data, residue_agg=True, agg_method='minimum')
        self.assertEqual(data.loc[0, 'RES_FULL_A'], '118')
        self.assertEqual(data.loc[0, 'RES_FULL_B'], '312')
        self.assertEqual(data.loc[0, 'ATOM_A'], 'CG1')
        self.assertEqual(data.loc[0, 'ATOM_B'], 'O')
        self.assertEqual(data.loc[0, 'DIST'], 3.758)
        self.assertEqual(data.loc[0, 'VDW_DIST'], 0.538)
        self.assertEqual(data.loc[1, 'RES_FULL_A'], '118')
        self.assertEqual(data.loc[1, 'RES_FULL_B'], '409')
        self.assertEqual(data.loc[1, 'ATOM_A'], 'CG2')
        self.assertEqual(data.loc[1, 'ATOM_B'], 'CG')
        self.assertEqual(data.loc[1, 'DIST'], 4.370)
        self.assertEqual(data.loc[1, 'VDW_DIST'], 0.970)

    def test_residues_agg_method(self):
        data = self.ARPEGGIO.read(self.inputarpeggio, excluded_cols=self.excluded)
        data = self.residues_aggregation(data, agg_method='first')
        self.assertEqual(data.loc[0, 'RES_FULL_A'], '118')
        self.assertEqual(data.loc[0, 'RES_FULL_B'], '312')
        self.assertEqual(data.loc[0, 'ATOM_A'], 'CG1')
        self.assertEqual(data.loc[0, 'ATOM_B'], 'C')
        self.assertEqual(data.loc[0, 'DIST'], 4.458)
        self.assertEqual(data.loc[0, 'VDW_DIST'], 1.058)
        self.assertEqual(data.loc[1, 'RES_FULL_A'], '118')
        self.assertEqual(data.loc[1, 'RES_FULL_B'], '409')
        self.assertEqual(data.loc[1, 'ATOM_A'], 'CG2')
        self.assertEqual(data.loc[1, 'ATOM_B'], 'CD')
        self.assertEqual(data.loc[1, 'DIST'], 4.533)
        self.assertEqual(data.loc[1, 'VDW_DIST'], 1.133)

    def test_residues_agg_collapsed(self):
        data = self.ARPEGGIO.read(self.inputarpeggio, excluded_cols=self.excluded)
        data = self.filter_arpeggio(data, collapsed_cont=True, col_method='minimal',
                                    residue_agg=True, agg_method='unique')
        self.assertEqual('Hydrophobic-Bond', ', '.join(list(data.loc[12, 'Int_Types'])))
        self.assertEqual(len(data.loc[12, 'DIST']), 2)
        data = self.ARPEGGIO.read(self.inputarpeggio, excluded_cols=self.excluded)
        data = self.filter_arpeggio(data, collapsed_cont=True, col_method='minimal',
                                    residue_agg=True, agg_method='minimum')
        self.assertEqual('Hydrophobic-Bond', ', '.join(list(data.loc[12, 'Int_Types'])))
        self.assertAlmostEqual(data.loc[12, 'DIST'], 4.259, places=2)

    def test_ignore_consecutive(self):
        data = self.ARPEGGIO.read(self.inputarpeggio, excluded_cols=self.excluded)
        self.assertEqual(len(data.index), 27420)
        data = self.ignore_consecutive(data, numb_res=5)
        self.assertEqual(len(data.index), 12105)
        data = self.ARPEGGIO.read(self.inputarpeggio, excluded_cols=self.excluded)
        data = self.filter_arpeggio(data, residue_agg=True, agg_method='minimum')
        self.assertEqual(len(data.index), 4626)
        data = self.ignore_consecutive(data, numb_res=5)
        self.assertEqual(len(data.index), 2269)

    def test_parser_spec_amam(self):
        data = self.parser_spec(self.input_amam, int_type="res-res")
        self.assertEqual(data.loc[0, 'CHAIN_A'], 'A')
        self.assertEqual(data.loc[0, 'CHAIN_B'], 'A')
        self.assertEqual(data.loc[0, 'RES_FULL_A'], '159')
        self.assertEqual(data.loc[0, 'RES_FULL_B'], '161')

    def test_parser_spec_amri(self):
        data = self.parser_spec(self.input_amri, int_type="res-res")
        self.assertEqual(data.loc[0, 'CHAIN_A'], 'A')
        self.assertEqual(data.loc[0, 'CHAIN_B'], 'A')
        self.assertEqual(data.loc[0, 'RES_FULL_A'], '121')
        self.assertEqual(data.loc[0, 'RES_FULL_B'], '120')

    def test_parser_spec_ari(self):
        data = self.parser_spec(self.input_ari, int_type="atom-res")
        self.assertEqual(data.loc[0, 'CHAIN_A'], 'B')
        self.assertEqual(data.loc[0, 'CHAIN_B'], 'B')
        self.assertEqual(data.loc[0, 'RES_FULL_A'], '222')
        self.assertEqual(data.loc[0, 'RES_FULL_B'], '220')
        self.assertEqual(data.loc[0, 'ATOM_A'], 'N')

    def test_parser_spec_ri(self):
        data = self.parser_spec(self.input_ri, int_type="res-res")
        self.assertEqual(data.loc[0, 'CHAIN_A'], 'B')
        self.assertEqual(data.loc[0, 'CHAIN_B'], 'B')
        self.assertEqual(data.loc[0, 'RES_FULL_A'], '208')
        self.assertEqual(data.loc[0, 'RES_FULL_B'], '204')

    def test_add_contact_info_res_res(self):
        data = self.parser(self.inputarpeggio,
                           excluded_cols=("ENTRY_A", "ENTRY_B", "ENTITIES"))
        info = self.parser_spec(self.input_ri, int_type="res-res",
                                excluded_cols=("ID_A", "ENTRY_A", "COORDS_A", "ID_B",
                                               "ENTRY_B", "COORDS_B",
                                               "CONT_TYPE", "INT_TYPE", "SELECT"))
        table = self.add_contact_info(data, info, int_type="res-res",
                                      col_name="Aromatic-Aromatic")
        self.assertIn('Aromatic-Aromatic', list(table))

    def test_add_contact_info_atom_res(self):
        data = self.parser(self.inputarpeggio,
                           excluded_cols=("ENTRY_A", "ENTRY_B", "ENTITIES"))
        info = self.parser_spec(self.input_ari, int_type="atom-res",
                                excluded_cols=("ID_A", "ENTRY_A", "COORDS_A", "ID_B", "ENTRY_B",
                                               "COORDS_B", "CONT_TYPE", "INT_TYPE", "SELECT"))
        table = self.add_contact_info(data, info, int_type="atom-res",
                                      col_name="Atom-Ring")
        self.assertIn('Atom-Ring', list(table))

    def test_parse_special_reader(self):
        data = self.ARPEGGIO.read(self.inputarpeggio,
                                  excluded_cols=self.excluded,
                                  parse_special=True)
        self.assertIn('Amide-Amide', list(data))
        self.assertIn('Aromatic-Aromatic', list(data))
        self.assertIn('Amide-Aromatic', list(data))
        self.assertIn('Atom-Ring', list(data))

    def test_add_special_cont_types(self):
        data = self.ARPEGGIO.read(self.inputarpeggio,
                                  excluded_cols=self.excluded,
                                  parse_special=False)
        data = self.add_special_cont_types(self.inputarpeggio, data)
        self.assertIn('Amide-Amide', list(data))
        self.assertIn('Aromatic-Aromatic', list(data))
        self.assertIn('Amide-Aromatic', list(data))
        self.assertIn('Atom-Ring', list(data))


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("prointvar").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestARPEGGIO)
    unittest.TextTestRunner(verbosity=2).run(suite)
