#!/local/bin/python
# -*- coding: utf-8 -*-


import os
import unittest
from unittest.mock import patch

from prointvar.config import Defaults
from prointvar.config import config


mock_config = """\
[Global]
stamp_bin = /stamp_dir/

[Other]
test = /test/value
"""


@patch("prointvar.config.config.db_root", "/root/")
@patch("prointvar.config.config.tmp_dir", "/tmp/")
class TestUTILS(unittest.TestCase):
    """Test the Config methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.config = config
        self.defaults = Defaults

    def tearDown(self):
        """Remove testing framework."""

        self.config = None
        self.defaults = None

    def test_loading_config_defaults(self):
        config = self.config
        self.assertTrue(hasattr(config, 'db_root'))
        self.assertEqual(config.db_root, "/root/")
        self.assertTrue(hasattr(config, 'tmp_dir'))
        self.assertEqual(config.tmp_dir, "/tmp/")
        self.assertFalse(hasattr(config, 'test'))

    def test_updating_config_defaults(self):
        config = self.config
        config.db_root = '/new_value/'
        self.assertNotEqual(config.db_root, "/stamp_dir/")
        self.assertEqual(config.db_root, "/new_value/")

    def test_deleting_entry_config(self):
        config = self.config
        del config.db_root
        self.assertFalse(hasattr(config, 'db_root'))

    def test_loading_config_from_file(self):
        # using mocking config in this case
        didnt_exist = False
        if not os.path.exists(self.config.tmp_dir):
            didnt_exist = True
            os.makedirs(self.config.tmp_dir)
        new_config_file = "{}mock_config.ini".format(self.config.tmp_dir)
        with open(new_config_file, 'w') as out:
            out.write(mock_config)
        config = self.defaults(config_file=new_config_file)
        self.assertFalse(hasattr(config, 'some_name'))
        self.assertTrue(hasattr(config, 'test'))
        self.assertEqual(config.test, "/test/value")
        os.remove("{}mock_config.ini".format(self.config.tmp_dir))
        if didnt_exist:
            os.removedirs(self.config.tmp_dir)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestUTILS)
    unittest.TextTestRunner(verbosity=2).run(suite)
