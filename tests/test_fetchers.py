#!/local/bin/python
# -*- coding: utf-8 -*-

"""
This suit doesn't check the data fetched because it will change overtime.
That is tested in the main fetching method implemented in test_utils.py!

Here, we only test whether the endpoints still exist or not!
"""

import os
import sys
import logging
import unittest
import requests_cache

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
try:
    from mock import patch
except ImportError:
    from unittest.mock import patch

from prointvar.fetchers import (fetch_best_structures_pdbe,
                                fetch_uniprot_variants_ebi,
                                fetch_summary_properties_pdbe,
                                get_preferred_assembly_id,
                                download_structure_from_pdbe,
                                download_sifts_from_ebi,
                                download_data_from_uniprot,
                                download_alignment_from_cath,
                                download_alignment_from_pfam,
                                fetch_uniprot_fasta,
                                fetch_uniprot_id_from_name,
                                BioFetcher, BioDownloader,
                                fetch_uniprot_species_from_id,
                                fetch_ensembl_uniprot_ensembl_mapping,
                                fetch_ensembl_ensembl_uniprot_mapping,
                                fetch_ensembl_transcript_variants,
                                fetch_ensembl_somatic_variants,
                                fetch_ensembl_variants_by_id,
                                fetch_ensembl_sequence_from_id)

from prointvar.variants import (get_ensembl_protein_id_from_mapping,
                                get_uniprot_id_from_mapping)

from prointvar.config import config as c

root = os.path.abspath(os.path.dirname(__file__))
c.db_root = "{}/testdata/".format(root)


@patch("prointvar.config.config.db_root", c.db_root)
@patch("prointvar.config.config.db_pdbx", 'tmp/')
@patch("prointvar.config.config.db_sifts", 'tmp/')
@patch("prointvar.config.config.db_uniprot", 'tmp/')
@patch("prointvar.config.config.db_cath", 'tmp/')
@patch("prointvar.config.config.db_pfam", 'tmp/')
class TestFetchers(unittest.TestCase):
    """Test the DSSP parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        requests_cache.uninstall_cache()

        self.uniprotid = "P00439"
        self.pdbid = "2pah"
        self.cathid = "1.50.10.100_1318"
        self.pfamid = "PF08124"
        self.ensemblid = "ENSP00000448059"
        self.varid = "rs750420403"
        self.fetch_best_structures_pdbe = fetch_best_structures_pdbe
        self.fetch_summary_properties_pdbe = fetch_summary_properties_pdbe
        self.fetch_preferred_assembly_id = get_preferred_assembly_id
        self.fetch_uniprot_variants_ebi = fetch_uniprot_variants_ebi
        self.download_structure_from_pdbe = download_structure_from_pdbe
        self.download_sifts_from_ebi = download_sifts_from_ebi
        self.download_data_from_uniprot = download_data_from_uniprot
        self.download_alignment_from_cath = download_alignment_from_cath
        self.download_alignment_from_pfam = download_alignment_from_pfam
        self.fetch_uniprot_fasta = fetch_uniprot_fasta
        self.fetch_uniprot_id_from_name = fetch_uniprot_id_from_name
        self.fetch_uniprot_species_from_id = fetch_uniprot_species_from_id
        self.fetch_ensembl_uniprot_ensembl_mapping = fetch_ensembl_uniprot_ensembl_mapping
        self.fetch_ensembl_ensembl_uniprot_mapping = fetch_ensembl_ensembl_uniprot_mapping
        self.fetch_ensembl_transcript_variants = fetch_ensembl_transcript_variants
        self.fetch_ensembl_somatic_variants = fetch_ensembl_somatic_variants
        self.fetch_ensembl_variants_by_id = fetch_ensembl_variants_by_id
        self.fetch_ensembl_sequence_from_id = fetch_ensembl_sequence_from_id
        self.BioFetcher = BioFetcher
        self.BioDownloader = BioDownloader
        self.get_ensembl_protein_id_from_mapping = get_ensembl_protein_id_from_mapping
        self.get_uniprot_id_from_mapping = get_uniprot_id_from_mapping

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.uniprotid = None
        self.pdbid = None
        self.cathid = None
        self.pfamid = None
        self.ensemblid = None
        self.varid = None
        self.fetch_best_structures_pdbe = None
        self.fetch_summary_properties_pdbe = None
        self.fetch_preferred_assembly_id = None
        self.fetch_uniprot_variants_ebi = None
        self.download_structure_from_pdbe = None
        self.download_sifts_from_ebi = None
        self.download_data_from_uniprot = None
        self.download_alignment_from_cath = None
        self.download_alignment_from_pfam = None
        self.fetch_uniprot_fasta = None
        self.fetch_uniprot_id_from_name = None
        self.fetch_uniprot_species_from_id = None
        self.fetch_ensembl_uniprot_ensembl_mapping = None
        self.fetch_ensembl_ensembl_uniprot_mapping = None
        self.fetch_ensembl_transcript_variants = None
        self.fetch_ensembl_somatic_variants = None
        self.fetch_ensembl_variants_by_id = None
        self.fetch_ensembl_sequence_from_id = None
        self.BioFetcher = None
        self.BioDownloader = None
        self.get_ensembl_protein_id_from_mapping = None
        self.get_uniprot_id_from_mapping = None
        logging.disable(logging.NOTSET)

    def test_best_structure_pdbe(self):
        r = self.fetch_best_structures_pdbe(self.uniprotid)
        self.assertTrue(r.ok)

    def test_best_structure_pdbe_cached(self):
        pickled = os.path.join(c.db_root, c.db_pickled, "{}_bs.pkl".format(self.uniprotid))
        self.assertFalse(os.path.isfile(pickled))
        r = self.fetch_best_structures_pdbe(self.uniprotid, cached=True)
        self.assertTrue(r.ok)
        self.assertTrue(os.path.isfile(pickled))
        os.remove(pickled)

    def test_summary_properties_pdbe(self):
        r = self.fetch_summary_properties_pdbe(self.pdbid)
        self.assertTrue(r.ok)

    def test_summary_properties_cached(self):
        pickled = os.path.join(c.db_root, c.db_pickled, "{}_sp.pkl".format(self.pdbid))
        self.assertFalse(os.path.isfile(pickled))
        r = self.fetch_summary_properties_pdbe(self.pdbid, cached=True)
        self.assertTrue(r.ok)
        self.assertTrue(os.path.isfile(pickled))
        os.remove(pickled)

    def test_preferred_assembly_pdbe(self):
        r = self.fetch_preferred_assembly_id(self.pdbid)
        self.assertEqual("1", r)

    def test_uniprot_variants_ebi(self):
        r = self.fetch_uniprot_variants_ebi(self.uniprotid)
        self.assertTrue(r.ok)

    def test_uniprot_variants_ebi_cached(self):
        pickled = os.path.join(c.db_root, c.db_pickled,
                               "{}_vars_uni.pkl".format(self.uniprotid))
        self.assertFalse(os.path.isfile(pickled))
        r = self.fetch_uniprot_variants_ebi(self.uniprotid, cached=True)
        self.assertTrue(r.ok)
        self.assertTrue(os.path.isfile(pickled))
        os.remove(pickled)

    def test_download_structure_from_pdbe_pdb(self):
        self.download_structure_from_pdbe(self.pdbid, pdb=True)
        os.remove(os.path.join(c.db_root, c.db_pdbx, "{}.pdb".format(self.pdbid)))

    def test_download_structure_from_pdbe_mmcif(self):
        self.download_structure_from_pdbe(self.pdbid, pdb=False)
        os.remove(os.path.join(c.db_root, c.db_pdbx, "{}.cif".format(self.pdbid)))

    def test_download_structure_from_pdbe_mmcif_bio(self):
        self.download_structure_from_pdbe(self.pdbid, pdb=False, bio=True)
        os.remove(os.path.join(c.db_root, c.db_pdbx, "{}_bio.cif".format(self.pdbid)))

    def test_download_sifts_from_ebi(self):
        self.download_sifts_from_ebi(self.pdbid)
        os.remove(os.path.join(c.db_root, c.db_sifts, "{}.xml".format(self.pdbid)))

    def test_download_data_from_uniprot_fasta(self):
        self.download_data_from_uniprot(self.uniprotid, file_format="fasta")
        os.remove(os.path.join(c.db_root, c.db_uniprot, "{}.fasta".format(self.uniprotid)))

    def test_download_data_from_uniprot_gff(self):
        self.download_data_from_uniprot(self.uniprotid, file_format="gff")
        os.remove(os.path.join(c.db_root, c.db_uniprot, "{}.gff".format(self.uniprotid)))

    def test_download_data_from_uniprot_txt(self):
        self.download_data_from_uniprot(self.uniprotid, file_format="txt")
        os.remove(os.path.join(c.db_root, c.db_uniprot, "{}.txt".format(self.uniprotid)))

    def test_download_alignment_from_cath(self):
        self.download_alignment_from_cath(self.cathid)
        os.remove(os.path.join(c.db_root, c.db_cath, "{}.fasta".format(self.cathid)))

    def test_download_alignment_from_pfam(self):
        self.download_alignment_from_pfam(self.pfamid)
        os.remove(os.path.join(c.db_root, c.db_pfam, "{}.sth".format(self.pfamid)))

    def test_fetch_uniprot_fasta(self):
        r = self.fetch_uniprot_fasta(self.uniprotid)
        self.assertTrue(r.ok)

    def test_fetch_uniprot_fasta_cached(self):
        pickled = os.path.join(c.db_root, c.db_pickled, "{}_fasta.pkl".format(self.uniprotid))
        self.assertFalse(os.path.isfile(pickled))
        r = self.fetch_uniprot_fasta(self.uniprotid, cached=True)
        self.assertTrue(r.ok)
        self.assertTrue(os.path.isfile(pickled))
        os.remove(pickled)

    def test_fetch_uniprot_id_from_name(self):
        r = self.fetch_uniprot_id_from_name("PH4H_HUMAN")
        self.assertTrue(r.ok)
        self.assertEqual("P00439", str(r.content, encoding='utf-8').strip())

    def test_fetch_uniprot_id_from_name_cached(self):
        pickled = os.path.join(c.db_root, c.db_pickled, "{}_id.pkl".format("PH4H_HUMAN"))
        self.assertFalse(os.path.isfile(pickled))
        r = self.fetch_uniprot_id_from_name("PH4H_HUMAN", cached=True)
        self.assertTrue(r.ok)
        self.assertTrue(os.path.isfile(pickled))
        os.remove(pickled)

    def test_fetch_uniprot_species_from_id(self):
        r = self.fetch_uniprot_species_from_id(self.uniprotid)
        self.assertTrue(r.ok)
        organism = str(r.content, encoding='utf-8').split('\n')[1]
        species = '_'.join(organism.split()[0:2]).lower()
        self.assertEqual(species, "homo_sapiens")

    def test_fetch_uniprot_species_from_id_cached(self):
        pickled = os.path.join(c.db_root, c.db_pickled,
                               "{}_org.pkl".format(self.uniprotid))
        self.assertFalse(os.path.isfile(pickled))
        r = self.fetch_uniprot_species_from_id(self.uniprotid, cached=True)
        self.assertTrue(r.ok)
        self.assertTrue(os.path.isfile(pickled))
        os.remove(pickled)

    def test_fetch_ensembl_uniprot_ensembl_mapping(self):
        r = self.fetch_ensembl_uniprot_ensembl_mapping(self.uniprotid)
        self.assertTrue(r.ok)
        ensps = self.get_ensembl_protein_id_from_mapping(r.json())
        self.assertEqual(ensps, [self.ensemblid])

    def test_fetch_ensembl_ensembl_uniprot_mapping(self):
        r = self.fetch_ensembl_ensembl_uniprot_mapping(self.ensemblid)
        self.assertTrue(r.ok)
        uniprots = self.get_uniprot_id_from_mapping(r.json())
        self.assertEqual(uniprots, ['A0A024RBG4', self.uniprotid])

    def test_fetch_ensembl_transcript_variants(self):
        r = self.fetch_ensembl_transcript_variants(self.ensemblid)
        self.assertTrue(r.ok)

    def fetch_fetch_ensembl_somatic_variants(self):
        r = self.fetch_ensembl_somatic_variants(self.ensemblid)
        self.assertTrue(r.ok)

    def fetch_fetch_ensembl_variants_by_id(self):
        r = self.fetch_ensembl_variants_by_id(self.varid)
        self.assertTrue(r.ok)

    def fetch_fetch_ensembl_sequence_from_id(self):
        r = self.fetch_ensembl_sequence_from_id(self.ensemblid)
        self.assertTrue(r.ok)

    def test_biofetcher_uniprot_fasta(self):
        url_root = c.http_uniprot
        url_endpoint = self.uniprotid + ".fasta"
        url = url_root + url_endpoint
        b = self.BioFetcher(url, cached=False)
        r = b.response
        self.assertTrue(r.ok)

    def test_biodownloader_uniprot_fasta(self):
        filename = self.uniprotid + ".fasta"
        outputfile = os.path.join(c.db_root, c.db_uniprot, filename)
        url = c.http_uniprot + filename
        self.BioDownloader(url=url, outputfile=outputfile,
                           decompress=True, override=False)
        os.remove(os.path.join(c.db_root, c.db_uniprot, filename))


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("prointvar").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFetchers)
    unittest.TextTestRunner(verbosity=2).run(suite)
