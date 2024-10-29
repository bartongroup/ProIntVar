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

from prointvar.pdbx import PDBXreader, PDBXwriter

from prointvar.arpeggio import (ARPEGGIOreader, ARPEGGIOrunner,
                                parse_arpeggio_from_file,
                                get_arpeggio_selected_from_table,
                                interaction_modes, residues_aggregation,
                                ignore_consecutive_residues)

from prointvar.config import config as c

root = os.path.abspath(os.path.dirname(__file__))
c.db_root = "{}/testdata/".format(root)


@patch("prointvar.config.config.db_root", c.db_root)
class TestARPEGGIO(unittest.TestCase):
    """Test the ARPEGGIO parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.pdbid = '2pah'
        self.pdbid_small = '2rea'
        self.inputpdb = os.path.join(c.db_root, c.db_pdbx, "{}.pdb".format(self.pdbid))
        self.inputpdb_fast = os.path.join(c.db_root, c.db_pdbx,
                                          "{}.pdb".format(self.pdbid_small))
        self.inputcif = os.path.join(c.db_root, c.db_pdbx, "{}.cif".format(self.pdbid))
        self.inputcif_fast = os.path.join(c.db_root, c.db_pdbx, "{}.cif".format(self.pdbid_small))
        self.inputarpeggio = os.path.join(c.db_root, c.db_contacts,
                                          "{}.json".format(self.pdbid))
        self.inputarpeggio_fast = os.path.join(c.db_root, c.db_contacts,
                                               "{}.json".format(self.pdbid_small))
        self.inputarpeggio_demo = os.path.join(c.db_root, c.db_contacts, "pdbe-arpeggio-demo.json")
        self.emptyfile = os.path.join(c.db_root, c.db_tmp, "{}.tmp".format(self.pdbid))
        self.notfound = ""

        self.parser = parse_arpeggio_from_file
        self.reader = ARPEGGIOreader
        self.generator = ARPEGGIOrunner
        self.filter = get_arpeggio_selected_from_table
        self.interaction_modes = interaction_modes
        self.residues_aggregation = residues_aggregation
        self.ignore_consecutive = ignore_consecutive_residues

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.pdbid_small = None
        self.inputpdb = None
        self.inputpdb_fast = None
        self.inputcif = None
        self.inputcif_fast = None
        self.inputarpeggio = None
        self.inputarpeggio_fast = None

        self.emptyfile = None
        self.notfound = None
        self.parser = None
        self.reader = None
        self.generator = None
        self.filter = None
        self.interaction_modes = None
        self.residues_aggregation = None
        self.ignore_consecutive = None
        self.parser_spec = None

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
        if os.path.isfile(self.inputpdb_fast):
            try:
                self.generator(self.inputpdb_fast,
                               self.inputarpeggio_fast).run(clean_output=True,
                                                            override=True)
            except OSError as e:
                if str(e) == "ARPEGGIO executable is not available...":
                    self.skipTest("ARPEGGIO executable is not available")
                else:
                    raise e
            msg = ("Arpeggio execution failed: make sure the settings "
                   "are set properly in config.ini!")
            self.assertTrue(os.path.isfile(self.inputarpeggio_fast), msg)
            os.remove(self.inputarpeggio_fast)
        else:
            raise IOError("%s" % self.inputpdb_fast)

    def test_generator_pdb_exec_fail(self):
        pdbid = "1ejg"
        inputpdb = os.path.join(c.db_root, c.db_pdbx, "{}.pdb".format(pdbid))
        inputarpeggio = os.path.join(c.db_root, c.db_pdbx, "{}.json".format(pdbid))
        try:
            self.generator(inputpdb,
                           inputarpeggio).run(clean_output=True,
                                              override=True)
        except (FileNotFoundError, OSError):
            # expected failure
            msg = "PDB with residues have missing atoms..."
            self.assertFalse(os.path.isfile(inputarpeggio), msg)

        inputpdb_new = os.path.join(c.db_root, c.db_pdbx, "{}_new.pdb".format(pdbid))
        r = PDBXreader(inputpdb)
        data = r.atoms(format_type="pdb", remove_altloc=True,
                       remove_hydrogens=True, reset_atom_id=True,
                       remove_partial_res=True)
        w = PDBXwriter(inputfile=None, outputfile=inputpdb_new)
        w.run(data, format_type='pdb')

        try:
            self.generator(inputpdb_new,
                           inputarpeggio).run(clean_output=True,
                                              override=True)
        except OSError as e:
            if str(e) == "ARPEGGIO executable is not available...":
                self.skipTest("ARPEGGIO executable is not available")
            else:
                raise e
        self.assertTrue(os.path.isfile(inputarpeggio))
        os.remove(inputpdb_new)
        os.remove(inputarpeggio)

    # TODO: This test might not test the generator as 2pah.json is already generated. Review.
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

    def test_parser_keys(self):
        self.assertListEqual(sorted([k for k in
                                     self.parser(self.inputarpeggio)['bgn.auth_asym_id'].unique()]),
                             ['A', 'B'])
        self.assertListEqual(sorted([k for k in
                                     self.parser(self.inputarpeggio)['end.auth_asym_id'].unique()]),
                             ['A', 'B'])

    def test_reader_data(self):
        reader = self.reader(self.inputarpeggio)
        data = reader.read()
        self.assertEqual(data.loc[27419, 'bgn.auth_asym_id'], 'B')
        self.assertEqual(data.loc[27419, 'end.auth_asym_id'], 'B')
        self.assertEqual(data.loc[27419, 'end.auth_seq_id'], 376)
        self.assertEqual(data.loc[27419, 'bgn.auth_seq_id'], 374)
        self.assertEqual(data.loc[27419, 'bgn.pdbx_PDB_ins_code'], '?')
        self.assertEqual(data.loc[27419, 'end.pdbx_PDB_ins_code'], '?')
        self.assertEqual(data.loc[27419, 'end.auth_atom_id'], 'ND2')
        self.assertEqual(data.loc[27419, 'bgn.auth_atom_id'], 'O')
        self.assertEqual(data.loc[27419, 'distance'], 4.97)

    def test_reader_to_json_pretty(self):
        reader = self.reader(self.inputarpeggio)
        reader.read()
        data = reader.to_json()
        self.assertEqual(json.loads(data)[13851]['end.auth_asym_id'], 'B')
        self.assertEqual(json.loads(data)[13851]['end.auth_seq_id'], 376)

    def test_reader_to_json(self):
        reader = self.reader(self.inputarpeggio)
        reader.read()
        data = reader.to_json(pretty=False)
        self.assertEqual(json.loads(data)[13851]['distance'], 4.97)

    def test_filter_chain(self):
        reader = self.reader(self.inputarpeggio)
        reader.read()
        data = self.filter(reader.data, chain_B=('A',))
        self.assertNotIn("B", data['end.auth_asym_id'].unique())

    def test_filter_res(self):
        reader = self.reader(self.inputarpeggio)
        reader.read()
        data = self.filter(reader.data, res_A=(374,))
        self.assertNotIn('119', data['bgn.auth_seq_id'].unique())

    def test_interaction_modes(self):
        reader = self.reader(self.inputarpeggio)
        data = reader.contacts(int_filter=True, int_mode='inter-chain')
        self.assertEqual(285, len(data))
        self.assertNotEqual(data.loc[0, 'bgn.auth_asym_id'], data.loc[0, 'end.auth_asym_id'])
        data = reader.contacts(int_filter=True, int_mode='intra-chain')
        self.assertEqual(27135, len(data))
        self.assertEqual(data.loc[0, 'bgn.auth_asym_id'], data.loc[0, 'end.auth_asym_id'])
        data = reader.contacts()
        data = self.interaction_modes(data, int_mode='inter-chain')
        self.assertEqual(285, len(data))
        self.assertEqual(data.loc[0, 'bgn.auth_seq_id'], 431)
        self.assertEqual(data.loc[0, 'end.auth_seq_id'], 452)
        self.assertEqual(data.loc[124, 'bgn.auth_seq_id'], 368)
        self.assertEqual(data.loc[124, 'end.auth_seq_id'], 368)

    def test_residues_agg(self):
        reader = self.reader(self.inputarpeggio)
        data = reader.contacts(residue_agg=True, agg_method='minimum')
        self.assertEqual(data.loc[1, 'bgn.auth_seq_id'], 118)
        self.assertEqual(data.loc[1, 'end.auth_seq_id'], 312)
        self.assertEqual(data.loc[1, 'bgn.auth_atom_id'], 'CG1')
        self.assertEqual(data.loc[1, 'end.auth_atom_id'], 'O')
        self.assertEqual(data.loc[1, 'distance'], 3.76)
        self.assertEqual(data.loc[2, 'bgn.auth_seq_id'], 118)
        self.assertEqual(data.loc[2, 'end.auth_seq_id'], 409)
        self.assertEqual(data.loc[2, 'bgn.auth_atom_id'], 'CG2')
        self.assertEqual(data.loc[2, 'end.auth_atom_id'], 'CG')
        self.assertEqual(data.loc[2, 'distance'], 4.37)

    def test_residues_agg_method(self):
        reader = self.reader(self.inputarpeggio)
        data = reader.contacts()
        data = self.residues_aggregation(data, agg_method='first')
        self.assertEqual(data.loc[1, 'bgn.auth_seq_id'], 118)
        self.assertEqual(data.loc[1, 'end.auth_seq_id'], 312)
        self.assertEqual(data.loc[1, 'bgn.auth_atom_id'], 'CB')
        self.assertEqual(data.loc[1, 'end.auth_atom_id'], 'O')
        self.assertEqual(data.loc[1, 'distance'], 4.98)
        self.assertEqual(data.loc[2, 'bgn.auth_seq_id'], 118)
        self.assertEqual(data.loc[2, 'end.auth_seq_id'], 409)
        self.assertEqual(data.loc[2, 'bgn.auth_atom_id'], 'CG2')
        self.assertEqual(data.loc[2, 'end.auth_atom_id'], 'CG')
        self.assertEqual(data.loc[2, 'distance'], 4.37)
        self.assertEqual(data.loc[1444, 'bgn.auth_seq_id'], 118)
        self.assertEqual(data.loc[1444, 'end.auth_seq_id'], 409)
        self.assertEqual(data.loc[1444, 'bgn.auth_atom_id'], 'CG2')
        self.assertEqual(data.loc[1444, 'end.auth_atom_id'], 'CD')
        self.assertEqual(data.loc[1444, 'distance'], 4.58)

    def test_ignore_consecutive(self):
        reader = self.reader(self.inputarpeggio)
        data = reader.contacts()
        self.assertEqual(len(data.index), 27420)
        data = self.ignore_consecutive(data, numb_res=5)
        self.assertEqual(len(data.index), 12105)
        data = reader.contacts(residue_agg=True, agg_method='minimum')
        self.assertEqual(len(data.index), 2856)  # reduced since aggregation after reordering bgn/end
        data = self.ignore_consecutive(data, numb_res=5)
        self.assertEqual(len(data.index), 1485)

    def test_parse_special_reader(self):
        reader = self.reader(self.inputarpeggio_demo)
        data = reader.contacts(parse_special=True)['type']
        self.assertIn('atom-atom', list(data))
        self.assertIn('group-group', list(data))
        self.assertIn('atom-plane', list(data))
        self.assertIn('plane-plane', list(data))


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("prointvar").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestARPEGGIO)
    unittest.TextTestRunner(verbosity=2).run(suite)
