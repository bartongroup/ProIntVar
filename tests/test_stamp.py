#!/local/bin/python
# -*- coding: utf-8 -*-


import os
import sys
import json
import logging
import unittest
import pandas as pd

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
try:
    from mock import patch
except ImportError:
    from unittest.mock import patch

from prointvar.sifts import SIFTSreader, get_sifts_selected_from_table
from prointvar.pdbx import PDBXreader, PDBXwriter, get_mmcif_selected_from_table
from prointvar.merger import TableMerger

from prointvar.stamp import (parse_stamp_scan_scores_from_file,
                             get_stamp_domain_line,
                             parse_stamp_domain_definitions_from_line,
                             parse_stamp_domain_definitions_from_from_file,
                             write_stamp_domain_definitions_from_table)

from prointvar.config import config as c

root = os.path.abspath(os.path.dirname(__file__))
c.db_root = "{}/testdata/".format(root)


@patch("prointvar.config.config.db_root", c.db_root)
class TestSTAMP(unittest.TestCase):
    """Test the STAMP parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.pdb_id = TestSTAMP.pdb_id
        self.chain_id = TestSTAMP.chain_id
        self.domains_info = TestSTAMP.domains_info

        self.input_scan = os.path.join(c.db_root, c.db_stamp, "scan.scores")
        self.domainfile = os.path.join(c.db_root, c.db_stamp, "domain.file")

        self.excluded = ()

        self.parse_stamp_scan_scores_from_file = parse_stamp_scan_scores_from_file
        self.get_stamp_domain_line = get_stamp_domain_line
        self.parse_stamp_domain_definition = parse_stamp_domain_definitions_from_line
        self.parse_stamp_domain_definitions = parse_stamp_domain_definitions_from_from_file
        self.write_stamp_domain_definitions = write_stamp_domain_definitions_from_table

        self.domain_def = "testdata/pdbx/2pah_A_1.pdb 2pah_A_1 { A 118 _ TO A 331 _ }"

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.pdb_id = None
        self.chain_id = None
        self.domains_info = None

        self.input_scan = None
        self.domainfile = None

        self.excluded = None

        self.parse_stamp_scan_scores_from_file = None
        self.get_stamp_domain_line = None
        self.parse_stamp_domain_definition = None
        self.parse_stamp_domain_definitions = None
        self.write_stamp_domain_definitions = None

        self.domain_def = None

        logging.disable(logging.NOTSET)

    @classmethod
    def setUpClass(cls):
        # to be run only once
        super(TestSTAMP, cls).setUpClass()

        cls.pdb_id = "2pah"
        cls.chain_id = "A"
        domains_info = []
        inputsifts = os.path.join(c.db_root, c.db_sifts, "%s.xml" % cls.pdb_id)
        s = SIFTSreader(inputfile=inputsifts)
        r = s.regions()
        if 'Pfam' in r[cls.chain_id]:
            for reg in r[cls.chain_id]['Pfam']:
                # UniProt Coordinates
                start = r[cls.chain_id]['Pfam'][reg]['start']
                end = r[cls.chain_id]['Pfam'][reg]['end']

                inputcif = os.path.join(c.db_root, c.db_pdbx, "%s.cif" % cls.pdb_id)
                p = PDBXreader(inputfile=inputcif)
                cif_table = p.atoms(category='auth', format_type="mmcif")
                cif_table = get_mmcif_selected_from_table(cif_table, chain=(cls.chain_id,),
                                                          category='auth',
                                                          lines=("ATOM",))
                inputsifts = os.path.join(c.db_root, c.db_sifts, "%s.xml" % cls.pdb_id)
                s = SIFTSreader(inputfile=inputsifts)
                r = s.regions()
                sifts_table = s.residues()
                sifts_table = get_sifts_selected_from_table(sifts_table,
                                                            chain_auth=(cls.chain_id,))

                t = TableMerger(pdbx_table=cif_table, sifts_table=sifts_table, store=False)
                table = t.merge()
                uni_sites = tuple([str(i) for i in range(start, end)])
                table = get_sifts_selected_from_table(table, chain_auth=(cls.chain_id,),
                                                      site=uni_sites)

                outputpdb = os.path.join(c.db_root, c.db_pdbx, "%s_%s_%s.pdb" % (cls.pdb_id,
                                                                                 cls.chain_id,
                                                                                 reg))
                p = PDBXwriter(outputfile=outputpdb)
                p.run(data=cif_table, format_type="pdb")

                # TODO loop over residues for getting the actual res ranges in the table
                # this is needed for STAMP
                start = table.loc[list(table.index)[0], "auth_seq_id"]
                start_inscode = table.loc[list(table.index)[0], "pdbx_PDB_ins_code"]
                end = table.loc[list(table.index)[-1], "auth_seq_id"]
                end_inscode = table.loc[list(table.index)[-1], "pdbx_PDB_ins_code"]
                info = {'pdb_id': cls.pdb_id,
                        'start_chain': (cls.chain_id,),
                        'end_chain': (cls.chain_id,),
                        'start': (start,), 'end': (end,),
                        'start_inscode': (start_inscode,),
                        'end_inscode': (end_inscode,),
                        'path': outputpdb,
                        'domain_id': '%s_%s_%s' % (cls.pdb_id, cls.chain_id, reg)}
                domains_info.append(info)
        cls.domains_info = pd.DataFrame(domains_info)

    @classmethod
    def tearDownClass(cls):

        cls.pdb_id = None
        cls.chain_id = None
        cls.domains_info = None

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

    def test_parse_stamp_domain_definition_from_line(self):
        r = self.parse_stamp_domain_definition(self.domain_def)
        self.assertEqual(r["domain_id"], '2pah_A_1')
        self.assertEqual(r["start_chain"], ('A',))
        self.assertEqual(r["start"], ('118',))
        self.assertEqual(r["start_inscode"], ('_',))
        self.assertEqual(r["end_chain"], ('A',))
        self.assertEqual(r["end"], ('331',))
        self.assertEqual(r["end_inscode"], ('_',))

    def test_get_stamp_domain_line(self):
        data = self.get_stamp_domain_line(self.domains_info)
        r = self.parse_stamp_domain_definition(data)
        self.assertEqual(r["domain_id"], '2pah_A_1')
        self.assertEqual(r["start_chain"], ('A',))
        self.assertEqual(r["start"], ('118',))
        self.assertEqual(r["start_inscode"], ('_',))
        self.assertEqual(r["end_chain"], ('A',))
        self.assertEqual(r["end"], ('331',))
        self.assertEqual(r["end_inscode"], ('_',))

    def test_write_stamp_domain_definitions(self):
        self.write_stamp_domain_definitions(self.domainfile,
                                            self.domains_info,
                                            override=True)
        self.assertTrue(os.path.isfile(self.domainfile))

    def test_parse_stamp_domain_definitions(self):
        r = self.parse_stamp_domain_definitions(self.domainfile)
        self.assertEqual(r.loc[0, "domain_id"], '2pah_A_1')
        self.assertEqual(r.loc[0, "start_chain"], ('A',))
        self.assertEqual(r.loc[0, "start"], ('118',))
        self.assertEqual(r.loc[0, "start_inscode"], ('_',))
        self.assertEqual(r.loc[0, "end_chain"], ('A',))
        self.assertEqual(r.loc[0, "end"], ('331',))
        self.assertEqual(r.loc[0, "end_inscode"], ('_',))


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("prointvar").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSTAMP)
    unittest.TextTestRunner(verbosity=2).run(suite)
