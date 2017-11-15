# -*- coding: utf-8 -*-

import os
import sys
import logging
import unittest

try:
    from mock import patch
except ImportError:
    # python 3.5
    from unittest.mock import patch

from proteofav.structures import mmCIF

from prointvar.config import config
from prointvar.structures import (parse_mmcif_categories_from_file,
                                  get_contact_indexes_from_table, add_mmcif_contacts,
                                  add_mmcif_new_pro_ids, get_coordinates)


@patch("prointvar.config.config", config)
class TestStructures(unittest.TestCase):
    """Test the PDBx parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.pdbid = '2pah'
        self.inputcif = os.path.join(os.path.dirname(__file__), "testdata",
                                     "mmcif", "{}.cif".format(self.pdbid))
        self.inputbiocif = os.path.join(os.path.dirname(__file__), "testdata",
                                        "mmcif", "{}_bio.cif".format(self.pdbid))

        self.parser_categories = parse_mmcif_categories_from_file
        self.get_contacts = get_contact_indexes_from_table
        self.add_contacts = add_mmcif_contacts
        self.add_mmcif_new_pro_ids = add_mmcif_new_pro_ids
        self.get_coordinates = get_coordinates

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.inputcif = None
        self.inputbiocif = None

        self.parser_categories = None
        self.get_contacts = None
        self.add_contacts = None
        self.add_mmcif_new_pro_ids = None
        self.get_coordinates = None

        logging.disable(logging.NOTSET)

    def test_parser_categories_cif_data(self):
        data = self.parser_categories(self.inputcif)
        # CATEGORIES
        self.assertEqual(data['_entry']['id'], [self.pdbid.upper()])
        self.assertEqual(data['_cell']['entry_id'], [self.pdbid.upper()])
        self.assertEqual(data['_symmetry']['entry_id'], [self.pdbid.upper()])

    def test_parser_categories_selected_cif_data(self):
        # CATEGORIES
        data = self.parser_categories(self.inputcif, category='_entry')
        self.assertEqual(data['_entry']['id'], [self.pdbid.upper()])
        data = self.parser_categories(self.inputcif, category='_cell')
        self.assertEqual(data['_cell']['entry_id'], [self.pdbid.upper()])
        data = self.parser_categories(self.inputcif, category='_symmetry')
        self.assertEqual(data['_symmetry']['entry_id'], [self.pdbid.upper()])

    def test_parser_categories_excluded_cif_data(self):
        # CATEGORIES
        data = self.parser_categories(self.inputcif, excluded=('_atom_site',))
        keys = [k for k in data.keys()]
        self.assertNotIn('_atom_site', keys)

    def test_get_contact_indexes(self):
        data = mmCIF.read(self.inputcif)
        indexes = self.get_contacts(data)
        self.assertEqual(indexes[0], [0, 1, 2, 3, 4, 5, 6, 7, 8, 13])

    def test_reader_add_contacts(self):
        data = mmCIF.read(self.inputcif)
        data = self.add_contacts(data)
        keys = [k for k in data]
        self.assertIn('contact_indexes', keys)

    def test_add_contacts_value(self):
        data = mmCIF.read(self.inputcif)
        data = self.add_contacts(data)
        self.assertEqual(data.loc[0, 'contact_indexes'], '0,1,2,3,4,5,6,7,8,13')

    def test_add_new_pro_id(self):
        data = mmCIF.read(self.inputbiocif)
        keys = [k for k in data]
        self.assertNotIn("new_asym_id", keys)
        self.assertNotIn("new_seq_id", keys)
        self.assertListEqual([k for k in data.label_asym_id.unique()],
                             ['A', 'AA', 'B', 'BA', 'C', 'CA', 'D', 'DA'])
        data = self.add_mmcif_new_pro_ids(data)
        keys = [k for k in data]
        self.assertIn("new_asym_id", keys)
        self.assertIn("new_seq_id", keys)
        self.assertListEqual([k for k in data.new_asym_id.unique()], ['A'])
        self.assertEqual(1306, len([k for k in data.new_seq_id.unique()]))

    def test_get_coordinates(self):
        data = mmCIF.read(self.inputcif)
        coords = get_coordinates(data)
        self.assertEqual(coords[0].tolist(), [-7.069, 21.943, 18.77])


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("prointvar").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestStructures)
    unittest.TextTestRunner(verbosity=2).run(suite)
