#!/local/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import logging
import unittest
import numpy as np

try:
    from mock import patch
except ImportError:
    from unittest.mock import patch

from prointvar.pdbx import (PDBXreader, get_mmcif_selected_from_table,
                            get_coordinates, get_sequence)

from prointvar.rmsd import (rmsd, centroid,
                            kabsch_rmsd, kabsch_rotate, kabsch,
                            quaternion_rmsd, quaternion_rotate, quaternion_transform,
                            makeQ, makeW)

from prointvar.utils import get_pairwise_alignment, get_pairwise_indexes

from prointvar.config import config as c

root = os.path.abspath(os.path.dirname(__file__))
c.db_root = "{}/testdata/".format(root)


@patch("prointvar.config.config.db_root", c.db_root)
class TestSupRMSD(unittest.TestCase):
    """Test the RMSD methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.inputcif = TestSupRMSD.inputcif
        self.table1 = TestSupRMSD.table1
        self.table2 = TestSupRMSD.table2
        self.coords1 = TestSupRMSD.coords1
        self.coords2 = TestSupRMSD.coords2
        self.centroid = centroid
        self.rmsd = rmsd
        self.get_coordinates = get_coordinates

        self.kabsch_rmsd = kabsch_rmsd
        self.kabsch_rotate = kabsch_rotate
        self.kabsch_algo = kabsch

        self.quaternion_rmsd = quaternion_rmsd
        self.quaternion_rotate = quaternion_rotate
        self.quaternion_transform = quaternion_transform
        self.makeQ = makeQ
        self.makeW = makeW

    def tearDown(self):
        """Remove testing framework."""

        self.inputcif = None
        self.table1 = None
        self.table2 = None
        self.coords1 = None
        self.coords2 = None
        self.centroid = None
        self.rmsd = None
        self.get_coordinates = None

        self.kabsch_rmsd = None
        self.kabsch_rotate = None
        self.kabsch_algo = None

        self.quaternion_rmsd = None
        self.quaternion_rotate = None
        self.quaternion_transform = None
        self.makeQ = None
        self.makeW = None

    @classmethod
    def setUpClass(cls):
        # to be run only once
        super(TestSupRMSD, cls).setUpClass()

        cls.inputcif = os.path.join(c.db_root, c.db_pdbx, "2pah.cif")

        r = PDBXreader(inputfile=cls.inputcif)
        table = r.atoms(format_type="mmcif")
        table = get_mmcif_selected_from_table(table, chain=('A',),
                                              atom=('CA',), lines=('ATOM',))
        table1 = table.reset_index()
        seq1 = get_sequence(table1)

        table = r.atoms(format_type="mmcif")
        table = get_mmcif_selected_from_table(table, chain=('B',),
                                              atom=('CA',), lines=('ATOM',))
        table2 = table.reset_index()
        seq2 = get_sequence(table2)

        seq1, seq2 = get_pairwise_alignment(seq1, seq2)
        match, drop = get_pairwise_indexes(seq1, seq2)

        cls.table1 = table1.drop(table1.index[drop.ix1])
        cls.table2 = table2.drop(table2.index[drop.ix2])

        cls.coords1 = get_coordinates(cls.table1)
        cls.coords2 = get_coordinates(cls.table2)

    @classmethod
    def tearDownClass(cls):
        cls.inputcif = None
        cls.table1 = None
        cls.table2 = None
        cls.coords1 = None
        cls.coords2 = None

    def assertListAlmostEqual(self, list1, list2, places):
        self.assertEqual(len(list1), len(list2))
        for a, b in zip(list1, list2):
            self.assertAlmostEqual(a, b, places=places)

    def test_centroid(self):
        a1 = np.array([-19.658, 17.18, 25.163], dtype=float)
        a2 = np.array([-20.573, 18.059, 25.88], dtype=float)
        a3 = np.array([-22.018, 17.551, 26.0], dtype=float)
        atms = np.asarray([a1, a2, a3])
        centroid = self.centroid(atms)
        self.assertEqual(3, len(centroid))
        self.assertListAlmostEqual([-20.7496, 17.5966, 25.6810],
                                   centroid, places=3)

    def test_rmsd(self):
        rmsd = self.rmsd(self.coords1, self.coords2)
        # self.assertAlmostEqual(32.5860, rmsd, places=3)
        self.assertAlmostEqual(57.4916, rmsd, places=3)

    def test_kabash_algorith(self):
        U = self.kabsch_algo(self.coords1, self.coords2)
        self.assertListAlmostEqual([-0.2333, 0.5816, -0.7792],
                                   U[0].tolist(), places=3)

    def test_kabash_rotate(self):
        nP = self.kabsch_rotate(self.coords1, self.coords2)
        # self.assertListAlmostEqual([1.6261, -9.0014, 23.2175],
        self.assertListAlmostEqual([-21.3124, 6.9790, 20.6722],
                                   nP[0].tolist(), places=3)

    def test_kabash_rmsd(self):
        Pc = self.centroid(self.coords1)
        Qc = self.centroid(self.coords2)
        self.coords1 -= Pc
        self.coords2 -= Qc
        rmsd = self.kabsch_rmsd(self.coords1, self.coords2)
        self.assertAlmostEqual(2.1780, rmsd, places=3)

    def test_quaternion_rmsd(self):
        Pc = self.centroid(self.coords1)
        Qc = self.centroid(self.coords2)
        self.coords1 -= Pc
        self.coords2 -= Qc
        rmsd = self.quaternion_rmsd(self.coords1, self.coords2)
        self.assertAlmostEqual(2.1780, rmsd, places=3)

    def test_quaternion_rotate(self):
        nP = self.quaternion_rotate(self.coords1, self.coords2)
        # self.assertListAlmostEqual([-0.5346, 0.8450, 0.0044],
        self.assertListAlmostEqual([-0.2333, 0.5816, -0.7792],
                                   nP[0].tolist(), places=3)

    def test_quaternion_transform(self):
        r = [-0.31019, -0.59291, 0.63612, -0.38415]
        U = self.quaternion_transform(r)
        self.assertListAlmostEqual([-0.5124, 0.8565, 0.0608],
                                   U[0].tolist(), places=3)

    def test_makeQ(self):
        r = [-0.31019, -0.59291, 0.63612, -0.38415]
        Q_r = self.makeQ(*r)
        self.assertListAlmostEqual([-0.3841, -0.6361, -0.5929, -0.3101],
                                   Q_r[0].tolist(), places=3)

    def test_makeW(self):
        r = [-0.31019, -0.59291, 0.63612, -0.38415]
        Wt_r = self.makeW(*r)
        self.assertListAlmostEqual([-0.3841, 0.6361, 0.5929, -0.3101],
                                   Wt_r[0].tolist(), places=3)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("prointvar").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSupRMSD)
    unittest.TextTestRunner(verbosity=2).run(suite)
