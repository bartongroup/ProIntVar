#!/local/bin/python
# -*- coding: utf-8 -*-

"""
This suit doesn't check the data fetched because it will change overtime.
That is tested in the main fetching method implemented in test_utils.py!

Here, we only test whether the endpoints still exist or not!
"""


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

from prointvar.fetchers import fetch_best_structures_pdbe

from prointvar.config import config as c
root = os.path.abspath(os.path.dirname(__file__))
c.db_root = "{}/testdata/".format(root)


@patch("prointvar.config.config.db_root", c.db_root)
class TestFetchers(unittest.TestCase):
    """Test the DSSP parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.uniprotid = "P00439"
        self.fetch_best_structures_pdbe = fetch_best_structures_pdbe

    def tearDown(self):
        """Remove testing framework."""

        self.uniprotid = None
        self.fetch_best_structures_pdbe = None

    def test_best_structure_pdbe(self):
        r = self.fetch_best_structures_pdbe(self.uniprotid)
        self.assertTrue(r.ok)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFetchers)
    unittest.TextTestRunner(verbosity=2).run(suite)
