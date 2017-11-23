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
from prointvar.structures import add_mmcif_new_pro_ids


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

        self.add_mmcif_new_pro_ids = add_mmcif_new_pro_ids

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.inputcif = None
        self.inputbiocif = None

        self.add_mmcif_new_pro_ids = None

        logging.disable(logging.NOTSET)

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


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("prointvar").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestStructures)
    unittest.TextTestRunner(verbosity=2).run(suite)
