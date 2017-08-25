#!/local/bin/python
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

from prointvar.stamp import parse_stamp_scan_scores_from_file

from prointvar.config import config as c

root = os.path.abspath(os.path.dirname(__file__))
c.db_root = "{}/testdata/".format(root)


@patch("prointvar.config.config.db_root", c.db_root)
class TestSTAMP(unittest.TestCase):
    """Test the STAMP parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.input_scan = os.path.join(c.db_root, c.db_stamp, "scan.scores")
        self.excluded = ()

        self.parse_stamp_scan_scores_from_file = parse_stamp_scan_scores_from_file

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.input_scan = None
        self.excluded = None

        self.parse_stamp_scan_scores_from_file = None

        logging.disable(logging.NOTSET)

    def test_parse_stamp_scan_scores_from_file(self):
        table = self.parse_stamp_scan_scores_from_file(self.input_scan,
                                                       self.excluded)
        # print(table.head())
        self.assertIn('Domain1', list(table))
        self.assertIn('Domain2', list(table))
        self.assertIn('RMS', list(table))
        self.assertIn('Sc', list(table))
        self.assertIn('Pm', list(table))
        self.assertEqual(9.795, table.loc[0, 'Sc'])
        self.assertEqual(0.084, table.loc[0, 'RMS'])
        self.assertEqual(280, table.loc[0, 'A_Len'])
        self.assertEqual(98.57, table.loc[0, 'PID'])


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("prointvar").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSTAMP)
    unittest.TextTestRunner(verbosity=2).run(suite)
