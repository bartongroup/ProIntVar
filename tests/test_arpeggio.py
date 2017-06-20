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

from prointvar.mmcif import MMCIFreader, MMCIFwriter

from prointvar.arpeggio import (ARPEGGIOreader, ARPEGGIOgenerator,
                                parse_arpeggio_from_file,
                                get_arpeggio_selected_from_table,
                                add_arpeggio_res_split,
                                interaction_modes, residues_aggregation,
                                collapsed_contacts, ignore_consecutive_residues,
                                parse_arpeggio_spec_from_file, add_contact_info,
                                add_special_cont_types)

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
        self.input_amam = "{}{}{}.amam".format(c.db_root, c.db_contacts_generated, self.pdbid)
        self.input_amri = "{}{}{}.amri".format(c.db_root, c.db_contacts_generated, self.pdbid)
        self.input_ari = "{}{}{}.ari".format(c.db_root, c.db_contacts_generated, self.pdbid)
        self.input_ri = "{}{}{}.ri".format(c.db_root, c.db_contacts_generated, self.pdbid)
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
        self.add_arpeggio_res_split = add_arpeggio_res_split
        self.interaction_modes = interaction_modes
        self.residues_aggregation = residues_aggregation
        self.collapsed_contacts = collapsed_contacts
        self.ignore_consecutive = ignore_consecutive_residues
        self.parser_spec = parse_arpeggio_spec_from_file
        self.add_contact_info = add_contact_info
        self.add_special_cont_types = add_special_cont_types

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.pdbid_small = None
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
        self.reader = None
        self.generator = None
        self.filter = None
        self.add_arpeggio_res_split = None
        self.interaction_modes = None
        self.residues_aggregation = None
        self.collapsed_contacts = None
        self.ignore_consecutive = None
        self.parser_spec = None
        self.add_contact_info = None
        self.add_special_cont_types = None

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

    # @unittest.expectedFailure
    def test_generator_pdb_exec_fail(self):
        pdbid = "1ejg"
        inputpdb = "{}{}{}.pdb".format(c.db_root, c.db_cif, pdbid)
        inputarpeggio = "{}{}{}.contacts".format(c.db_root, c.db_cif, pdbid)
        try:
            self.generator(inputpdb,
                           inputarpeggio).run(clean_output=True,
                                              override=True)
        except (FileNotFoundError, OSError):
            # expected failure
            msg = "PDB with residues have missing atoms..."
            self.assertFalse(os.path.isfile(inputarpeggio), msg)

        inputpdb_new = "{}{}{}_new.pdb".format(c.db_root, c.db_cif, pdbid)
        r = MMCIFreader(inputpdb)
        data = r.atoms(format_type="pdb", remove_altloc=True,
                       remove_hydrogens=True, reset_atom_id=True,
                       remove_partial_res=True)
        w = MMCIFwriter(inputfile=None, outputfile=inputpdb_new)
        w.run(data, format_type='pdb')

        self.generator(inputpdb_new,
                       inputarpeggio).run(clean_output=True,
                                          override=True)
        self.assertTrue(os.path.isfile(inputarpeggio))
        os.remove(inputpdb_new)
        os.remove(inputarpeggio)
        os.remove("{}{}{}.amam".format(c.db_root, c.db_cif, pdbid))
        os.remove("{}{}{}.amri".format(c.db_root, c.db_cif, pdbid))
        os.remove("{}{}{}.ari".format(c.db_root, c.db_cif, pdbid))
        os.remove("{}{}{}.ri".format(c.db_root, c.db_cif, pdbid))

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

    def test_add_split_res(self):
        reader = self.reader(self.inputarpeggio)
        data = reader.contacts(excluded=self.excluded, add_res_split=False)
        self.assertNotIn('CHAIN_A', data)
        self.assertNotIn('CHAIN_B', data)
        data = self.add_arpeggio_res_split(data)
        self.assertIn('CHAIN_A', data)
        self.assertIn('CHAIN_B', data)
        self.assertEqual(data.loc[0, 'CHAIN_A'], 'B')

    def test_interaction_modes(self):
        reader = self.reader(self.inputarpeggio)
        data = reader.contacts(int_filter=True, int_mode='inter-chain')
        self.assertEqual(285, len(data))
        self.assertNotEqual(data.loc[0, 'CHAIN_A'], data.loc[0, 'CHAIN_B'])
        data = reader.contacts(int_filter=True, int_mode='intra-chain')
        self.assertEqual(27135, len(data))
        self.assertEqual(data.loc[0, 'CHAIN_A'], data.loc[0, 'CHAIN_B'])
        data = reader.contacts()
        data = self.interaction_modes(data, int_mode='inter-chain')
        self.assertEqual(285, len(data))
        self.assertEqual(data.loc[0, 'RES_FULL_A'], '368')
        self.assertEqual(data.loc[0, 'RES_FULL_B'], '368')

    def test_collapsed_contacts(self):
        reader = self.reader(self.inputarpeggio)
        data = reader.contacts(collapsed_cont=True, col_method='full')
        self.assertNotIn('IONIC', list(data))
        self.assertIn('Int_Types', list(data))
        self.assertEqual('VDW-Proximal, Polar-Bond',
                         ', '.join(list(data.loc[3, 'Int_Types'])))
        reader = self.reader(self.inputarpeggio)
        data = reader.contacts()
        data = self.collapsed_contacts(data, col_method='minimal')
        self.assertNotIn('IONIC', list(data))
        self.assertIn('Int_Types', list(data))
        self.assertEqual('Polar-Bond', ', '.join(list(data.loc[3, 'Int_Types'])))

    def test_residues_agg(self):
        reader = self.reader(self.inputarpeggio)
        data = reader.contacts(residue_agg=True, agg_method='minimum')
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
        reader = self.reader(self.inputarpeggio)
        data = reader.contacts()
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
        reader = self.reader(self.inputarpeggio)
        data = reader.contacts(collapsed_cont=True, col_method='minimal',
                               residue_agg=True, agg_method='unique')
        self.assertEqual('Hydrophobic-Bond', ', '.join(list(data.loc[1, 'Int_Types'])))
        self.assertEqual(len(data.loc[1, 'DIST']), 2)
        data = reader.contacts(collapsed_cont=True, col_method='minimal',
                               residue_agg=True, agg_method='minimum')
        self.assertEqual('Hydrophobic-Bond', ', '.join(list(data.loc[1, 'Int_Types'])))
        self.assertEqual(data.loc[1, 'DIST'], 4.370)

    def test_ignore_consecutive(self):
        reader = self.reader(self.inputarpeggio)
        data = reader.contacts()
        self.assertEqual(len(data.index), 27420)
        data = self.ignore_consecutive(data, numb_res=5)
        self.assertEqual(len(data.index), 12105)
        data = reader.contacts(residue_agg=True, agg_method='minimum')
        self.assertEqual(len(data.index), 4609)
        data = self.ignore_consecutive(data, numb_res=5)
        self.assertEqual(len(data.index), 2252)

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
        data = self.parser(self.inputarpeggio, excluded=("ENTRY_A", "ENTRY_B", "ENTITIES"))
        info = self.parser_spec(self.input_ri, int_type="res-res", excluded=("ID_A", "ENTRY_A",
                                                                             "COORDS_A", "ID_B",
                                                                             "ENTRY_B", "COORDS_B",
                                                                             "CONT_TYPE", "INT_TYPE",
                                                                             "SELECT"))
        table = self.add_contact_info(data, info, int_type="res-res",
                                      col_name="Aromatic-Aromatic")
        self.assertIn('Aromatic-Aromatic', list(table))

    def test_add_contact_info_atom_res(self):
        data = self.parser(self.inputarpeggio, excluded=("ENTRY_A", "ENTRY_B", "ENTITIES"))
        info = self.parser_spec(self.input_ari, int_type="atom-res", excluded=("ID_A", "ENTRY_A",
                                                                               "COORDS_A", "ID_B",
                                                                               "ENTRY_B", "COORDS_B",
                                                                               "CONT_TYPE", "INT_TYPE",
                                                                               "SELECT"))
        table = self.add_contact_info(data, info, int_type="atom-res",
                                      col_name="Atom-Ring")
        self.assertIn('Atom-Ring', list(table))

    def test_parse_special_reader(self):
        reader = self.reader(self.inputarpeggio)
        data = reader.contacts(parse_special=True)
        self.assertIn('Amide-Amide', list(data))
        self.assertIn('Aromatic-Aromatic', list(data))
        self.assertIn('Amide-Aromatic', list(data))
        self.assertIn('Atom-Ring', list(data))

    def test_add_special_cont_types(self):
        reader = self.reader(self.inputarpeggio)
        data = reader.contacts(parse_special=False)
        data = self.add_special_cont_types(self.inputarpeggio, data)
        self.assertIn('Amide-Amide', list(data))
        self.assertIn('Aromatic-Aromatic', list(data))
        self.assertIn('Amide-Aromatic', list(data))
        self.assertIn('Atom-Ring', list(data))


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestARPEGGIO)
    unittest.TextTestRunner(verbosity=2).run(suite)
