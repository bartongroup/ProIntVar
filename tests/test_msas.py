#!/local/bin/python
# -*- coding: utf-8 -*-


import os
import sys
import json
import logging
import unittest

from Bio.Align import MultipleSeqAlignment

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
try:
    from mock import patch
except ImportError:
    from unittest.mock import patch

from prointvar.msas import (read_alignment, parse_msa_sequences_from_file,
                            parse_sequence_info_from_description,
                            parse_uniprot_fasta_seq_description,
                            parse_pfam_sth_seq_description,
                            parse_cath_fasta_seq_description,
                            parse_generic_seq_description,
                            MSAreader)

from prointvar.config import config as c

root = os.path.abspath(os.path.dirname(__file__))
c.db_root = "{}/testdata/".format(root)


@patch("prointvar.config.config.db_root", c.db_root)
class TestMSAS(unittest.TestCase):
    """Test the MSAs parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.cathid = '1.50.10.100_1318'
        self.pfamid = 'PF00118'

        self.uniprot_fasta = ("tr|A0A067NRR9|A0A067NRR9_PLEOS Polysaccharide lyase "
                              "family 8 protein OS=Pleurotus ostreatus PC15 "
                              "GN=PLEOSDRAFT_1036297 PE=4 SV=1")
        self.pfam_sto = "A0A067NRR9_PLEOS/22-313"
        self.cath_fasta = ("cath|4.1.0|1rwhA01/4-372 CATH_S35=1.50.10.100.1;"
                           "GENE=P84141_ARTAU;GO=GO:0005576,GO:0005975,GO:0016837,"
                           "GO:0030246;MDA=1.50.10.100;ORG=Arthrobacter_aurescens;"
                           "TAXON=43663;UNIPROT=P84141")
        self.cath_fasta2 = ("biomap|4.1.0|08901c114dd3dfa7eb9d39cf55143681/343-713 "
                            "EC=4.2.2.1;GENE=Q8VLQ7_STRSU;GO=GO:0005576,GO:0005975,GO:"
                            "0016020,GO:0030246,GO:0030340;MDA=2.60.120.260_2.60.40.1380_"
                            "1.50.10.100_2.70.98.10_2.60.220.10;ORG=Streptococcus_suis;"
                            "TAXON=1307;UNIPROT=Q8VLQ7")
        self.generic_id = "A0A067NRR9/22-313"

        self.inputcath = os.path.join(c.db_root, c.db_cath,
                                      "{}.fasta".format(self.cathid))
        self.inputpfam = os.path.join(c.db_root, c.db_pfam,
                                      "{}.sth".format(self.pfamid))
        self.emptyfile = os.path.join(c.db_root, c.db_tmp, "2pah.tmp")
        self.notfound = ""
        self.excluded = ()
        self.entry = {}

        self.reader = MSAreader
        self.read_alignment = read_alignment
        self.parser = parse_msa_sequences_from_file

        self.uniprot_fasta_seq_description = parse_uniprot_fasta_seq_description
        self.pfam_sth_seq_description = parse_pfam_sth_seq_description
        self.cath_fasta_seq_description = parse_cath_fasta_seq_description
        self.generic_seq_description = parse_generic_seq_description
        self.sequence_info_from_description = parse_sequence_info_from_description

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.cathid = None
        self.pfamid = None

        self.uniprot_fasta = None
        self.pfam_sto = None
        self.cath_fasta = None
        self.generic_id = None

        self.inputcath = None
        self.inputpfam = None

        self.emptyfile = None
        self.notfound = None
        self.excluded = None
        self.entry = None

        self.reader = None
        self.read_alignment = None
        self.parser = None

        self.uniprot_fasta_seq_description = None
        self.pfam_sth_seq_description = None
        self.cath_fasta_seq_description = None
        self.generic_seq_description = None
        self.sequence_info_from_description = None

        logging.disable(logging.NOTSET)

    def test_file_not_found_reader(self):
        with self.assertRaises(IOError):
            self.read_alignment(self.notfound)

    def test_file_not_found_parser(self):
        with self.assertRaises(IOError):
            self.parser(self.notfound)

    def test_empty_file_reader(self):
        with self.assertRaises(ValueError):
            open(self.emptyfile, 'w').close()
            self.read_alignment(self.emptyfile).read()
            os.remove(self.emptyfile)

    def test_empty_file_parser(self):
        with self.assertRaises(ValueError):
            open(self.emptyfile, 'w').close()
            self.parser(self.emptyfile)
            os.remove(self.emptyfile)

    def test_read_alignment(self):
        align = self.read_alignment(self.inputcath)
        self.assertTrue(type(align), MultipleSeqAlignment)
        self.assertEqual(len(align), 42)
        align = self.read_alignment(self.inputpfam)
        self.assertTrue(type(align), MultipleSeqAlignment)
        self.assertEqual(len(align), 57)

    def test_parse_uniprot_fasta_seq_description(self):
        self.uniprot_fasta_seq_description(self.uniprot_fasta,
                                           self.entry)
        self.assertEqual(self.entry['Source'], "UniProt")
        self.assertEqual(self.entry['Collection'], "tr")
        self.assertEqual(self.entry['Accession'], "A0A067NRR9")
        self.assertEqual(self.entry['Name'], "A0A067NRR9_PLEOS")

    def test_parse_pfam_sth_seq_description(self):
        self.pfam_sth_seq_description(self.pfam_sto, self.entry,
                                      get_uniprot_id=True, cached=False)
        self.assertEqual(self.entry['Source'], "Pfam")
        self.assertEqual(self.entry['Accession'], "A0A067NRR9")
        self.assertEqual(self.entry['Name'], "A0A067NRR9_PLEOS")
        self.assertEqual(self.entry['Start'], 22)
        self.assertEqual(self.entry['End'], 313)

    def test_parse_cath_fasta_seq_description(self):
        self.cath_fasta_seq_description(self.cath_fasta,
                                        self.entry)
        self.assertEqual(self.entry['Source'], "CATH")
        self.assertEqual(self.entry['Collection'], "cath")
        self.assertEqual(self.entry['Version'], "4.1.0")
        self.assertEqual(self.entry['Accession'], "1rwhA01")
        self.assertEqual(self.entry['Start'], 4)
        self.assertEqual(self.entry['End'], 372)

    def test_parse_cath_fasta_seq_description_2(self):
        self.cath_fasta_seq_description(self.cath_fasta2,
                                        self.entry)
        self.assertEqual(self.entry['Source'], "CATH")
        self.assertEqual(self.entry['Collection'], "biomap")
        self.assertEqual(self.entry['Version'], "4.1.0")
        self.assertEqual(self.entry['Accession'],
                         "08901c114dd3dfa7eb9d39cf55143681")
        self.assertEqual(self.entry['Start'], 343)
        self.assertEqual(self.entry['End'], 713)

    def test_parse_generic_seq_description(self):
        self.generic_seq_description(self.generic_id,
                                     self.entry)
        self.assertEqual(self.entry['Source'], "GenericParser")
        self.assertEqual(self.entry['Accession'], "A0A067NRR9")
        self.assertEqual(self.entry['Start'], 22)
        self.assertEqual(self.entry['End'], 313)

    def test_parse_sequence_info_from_description(self):
        self.sequence_info_from_description(self.uniprot_fasta,
                                            self.entry)
        self.assertEqual(self.entry['Source'], "UniProt")
        self.assertEqual(self.entry['Collection'], "tr")
        self.assertEqual(self.entry['Accession'], "A0A067NRR9")
        self.assertEqual(self.entry['Name'], "A0A067NRR9_PLEOS")

        self.sequence_info_from_description(self.pfam_sto, self.entry,
                                            get_uniprot_id=True, cached=False)
        self.assertEqual(self.entry['Source'], "Pfam")
        self.assertEqual(self.entry['Accession'], "A0A067NRR9")
        self.assertEqual(self.entry['Name'], "A0A067NRR9_PLEOS")
        self.assertEqual(self.entry['Start'], 22)
        self.assertEqual(self.entry['End'], 313)

        self.sequence_info_from_description(self.cath_fasta,
                                            self.entry)
        self.assertEqual(self.entry['Source'], "CATH")
        self.assertEqual(self.entry['Collection'], "cath")
        self.assertEqual(self.entry['Version'], "4.1.0")
        self.assertEqual(self.entry['Accession'], "1rwhA01")
        self.assertEqual(self.entry['Start'], 4)
        self.assertEqual(self.entry['End'], 372)

        self.sequence_info_from_description(self.generic_id,
                                            self.entry)
        self.assertEqual(self.entry['Source'], "GenericParser")
        self.assertEqual(self.entry['Accession'], "A0A067NRR9")
        self.assertEqual(self.entry['Start'], 22)
        self.assertEqual(self.entry['End'], 313)

    def test_read_alignment_to_table_cath(self):
        data = self.parser(inputfile=self.inputcath,
                           excluded=self.excluded)
        self.assertIn('Sequence', list(data))
        self.assertIn('Source', list(data))
        self.assertIn('Description', list(data))
        self.assertIn('Collection', list(data))
        self.assertEqual(data.loc[0, 'Accession'], '1hm3A01')
        self.assertEqual(data.loc[0, 'Collection'], 'cath')
        self.assertEqual(data.loc[0, 'Start'], 27)
        self.assertEqual(data.loc[0, 'End'], 338)

    def test_read_alignment_to_table_pfam(self):
        data = self.parser(inputfile=self.inputpfam,
                           excluded=self.excluded,
                           get_uniprot_id=True, cached=False)
        self.assertIn('Sequence', list(data))
        self.assertIn('Source', list(data))
        self.assertIn('Name', list(data))
        self.assertIn('Accession', list(data))
        self.assertEqual(data.loc[0, 'Accession'], 'B9LRY6')
        self.assertEqual(data.loc[0, 'Source'], 'Pfam')
        self.assertEqual(data.loc[0, 'Start'], 27)
        self.assertEqual(data.loc[0, 'End'], 514)

    def test_reader_class_pfam(self):
        r = self.reader(inputfile=self.inputcath)
        data = r.msas(excluded=self.excluded, cached=False,
                      get_uniprot_id=True)
        self.assertIn('Sequence', list(data))
        self.assertIn('Source', list(data))
        self.assertIn('Description', list(data))
        self.assertIn('Collection', list(data))
        self.assertEqual(data.loc[0, 'Accession'], '1hm3A01')
        self.assertEqual(data.loc[0, 'Collection'], 'cath')
        self.assertEqual(data.loc[0, 'Start'], 27)
        self.assertEqual(data.loc[0, 'End'], 338)

    def test_reader_class_pfam_to_json(self):
        r = self.reader(inputfile=self.inputcath)
        r.read(excluded=self.excluded, cached=False,
               get_uniprot_id=True)
        data = r.to_json()
        self.assertEqual(json.loads(data)[0]['Accession'], '1hm3A01')
        self.assertEqual(json.loads(data)[0]['Collection'], 'cath')
        self.assertEqual(json.loads(data)[0]['Start'], 27)
        self.assertEqual(json.loads(data)[0]['End'], 338)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("prointvar").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMSAS)
    unittest.TextTestRunner(verbosity=2).run(suite)
