#!/local/bin/python
# -*- coding: utf-8 -*-


import os
import sys
import json
import math
import unittest
import numpy as np
import pandas as pd

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
try:
    from mock import patch
except ImportError:
    from unittest.mock import patch

from prointvar.contacts import get_distance_between_atoms

from prointvar.contacts import (get_unit_vector, get_angle_between_vectors, get_vectors_from_points,
                                get_vector_subtraction, get_vector_scaled)

from prointvar.config import config as c
root = os.path.abspath(os.path.dirname(__file__))
c.db_root = "{}/testdata/".format(root)


@patch("prointvar.config.config.db_root", c.db_root)
class TestContacts(unittest.TestCase):
    """Test the DSSP parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.get_unit_vector = get_unit_vector
        self.get_angle_between = get_angle_between_vectors
        self.get_distance_between_atoms = get_distance_between_atoms
        self.get_vectors_from_atom_coords = get_vectors_from_points
        self.get_vector_subtraction = get_vector_subtraction
        self.get_vector_scaled = get_vector_scaled

    def tearDown(self):
        """Remove testing framework."""

        self.get_unit_vector = None
        self.get_angle_between = None
        self.get_distance_between_atoms = None
        self.get_vectors_from_atom_coords = None
        self.get_vector_subtraction = None
        self.get_vector_scaled = None

    def assertListAlmostEqual(self, list1, list2, places):
        self.assertEqual(len(list1), len(list2))
        for a, b in zip(list1, list2):
            self.assertAlmostEqual(a, b, places=places)

    def test_unit_vector(self):

        uv = self.get_unit_vector(np.array((1, 0, 0)))
        self.assertEqual([1, 0, 0], uv.tolist())
        # 2pah: A ASN 167 N
        uv = self.get_unit_vector(np.array((15.679, 9.840, 31.623)))
        self.assertListAlmostEqual([0.42789, 0.26854, 0.86301], uv.tolist(), places=4)

    def test_vectors_from_atom_coords(self):

        # 2pah: Donor A ASN 167 N --- Acceptor A ASP 163 O
        a1, a2, a3 = [15.696, 8.382, 31.566], [15.679, 9.840, 31.623], [18.116, 11.104, 32.820]
        v1, v2 = self.get_vectors_from_atom_coords(a1, a2, a3)
        self.assertListAlmostEqual([0.01699, -1.45800, -0.05700], v1, places=4)
        self.assertListAlmostEqual([-2.42000, -2.72199, -1.25400], v2, places=4)

        # 2pah: Donor A GLN 449 NE2 --- Acceptor A CYS 445 O
        a1, a2, a3 = [-36.909, 7.937, 54.987], [-38.831, 9.910, 52.686], [-35.435, 9.462, 51.388]
        v1, v2 = self.get_vectors_from_atom_coords(a1, a2, a3)
        self.assertListAlmostEqual([1.92200, -1.97299, 2.30100], v1, places=4)
        self.assertListAlmostEqual([-1.47399, -1.52499, 3.59900], v2, places=4)

    def test_vector_subtraction(self):

        v1, v2 = [-38.831, 9.910, 52.686], [-35.435, 9.462, 51.388]
        v = self.get_vector_subtraction(v1, v2)
        self.assertListAlmostEqual([-3.39600, 0.44800, 1.29800], v, places=4)

    def test_vector_scaled(self):

        v1, v2 = np.array((1, 0, 0)), np.array((1, 1, 0))
        nv1 = self.get_vector_scaled(v1, 0.5)
        self.assertListAlmostEqual([0.5, 0.0, 0.0], nv1, places=4)
        nv2 = self.get_vector_scaled(v2, 10.0)
        self.assertListAlmostEqual([10.0, 10.0, 0.0], nv2, places=4)

    def test_angle_between_vectors(self):

        v1, v2 = np.array((1, 0, 0)), np.array((1, 0, 0))
        angle = self.get_angle_between(v1, v2)
        self.assertEqual(0.0, angle)
        self.assertAlmostEqual(0.0, math.degrees(angle), places=4)

        v1, v2 = np.array((1, 0, 0)), np.array((0, 1, 0))
        angle = self.get_angle_between(v1, v2)
        self.assertEqual(math.pi / 2, angle)
        self.assertAlmostEqual(1.57079, angle, places=4)
        self.assertAlmostEqual(90.0, math.degrees(angle), places=4)

        v1, v2 = np.array((1, 0, 0)), np.array((-1, 0, 0))
        angle = self.get_angle_between(v1, v2)
        self.assertEqual(math.pi, angle)
        self.assertAlmostEqual(3.14159, angle, places=4)
        self.assertAlmostEqual(180.0, math.degrees(angle), places=4)

        # 2pah: Donor A ASN 167 N --- Acceptor A ASP 163 O
        v1, v2 = [15.679, 9.840, 31.623], [18.116, 11.104, 32.820]
        v1, v2 = np.array(v1), np.array(v2)
        angle = self.get_angle_between(v1, v2)
        self.assertAlmostEqual(0.04530, angle, places=4)
        self.assertAlmostEqual(2.59605, math.degrees(angle), places=4)

        # 2pah: Donor A GLN 449 NE2 --- Acceptor A CYS 445 O
        v1, v2 = [-38.831, 9.910, 52.686], [-35.435, 9.462, 51.388]
        v1, v2 = np.array(v1), np.array(v2)
        angle = self.get_angle_between(v1, v2)
        self.assertAlmostEqual(0.03110, angle, places=4)
        self.assertAlmostEqual(1.78238, math.degrees(angle), places=4)

    def test_distance_between_vectors(self):

        # 2pah: Donor A ASN 167 N --- Acceptor A ASP 163 O
        v1, v2 = [15.679, 9.840, 31.623], [18.116, 11.104, 32.820]
        data = {'Cartn_x': v1[0], 'Cartn_y': v1[1], 'Cartn_z': v1[2],
                'Cartn_x_2': v2[0], 'Cartn_y_2': v2[1], 'Cartn_z_2': v2[2]}
        t = pd.DataFrame(data, index=[0])
        dist = self.get_distance_between_atoms(t)
        self.assertEqual(2.99, dist)

        # 2pah: Donor A GLN 449 NE2 --- Acceptor A CYS 445 O
        v1, v2 = [-38.831, 9.910, 52.686], [-35.435, 9.462, 51.388]
        data = {'Cartn_x': v1[0], 'Cartn_y': v1[1], 'Cartn_z': v1[2],
                'Cartn_x_2': v2[0], 'Cartn_y_2': v2[1], 'Cartn_z_2': v2[2]}
        t = pd.DataFrame(data, index=[0])
        dist = self.get_distance_between_atoms(t)
        self.assertEqual(3.66, dist)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestContacts)
    unittest.TextTestRunner(verbosity=2).run(suite)
