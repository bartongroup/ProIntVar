#!/local/bin/python
# -*- coding: utf-8 -*-


import os
import sys
import logging
import unittest
import numpy as np

from prointvar.variants import (flatten_uniprot_variants_ebi,
                                collapse_unique_values, flatten_json)

from prointvar.config import config as c

root = os.path.abspath(os.path.dirname(__file__))
c.db_root = "{}/testdata/".format(root)

example_uniprot_variants = {
    'accession': 'P40227',
    'entryName': 'TCPZ_HUMAN',
    'sequence': 'MAAVKTLNPKAEVARAQAAMKATNILLVDEIMRAGMSSLKG',
    'sequenceChecksum': '4876218983604824961',
    'taxid': 9606,
    'features': [
        {'type': 'VARIANT',
         'alternativeSequence': 'F',
         'begin': '231',
         'end': '231',
         'xrefs': [{'name': '1000Genomes',
                    'id': 'rs148616984',
                    'url': 'http://www.ensembl.org/Homo_sapiens/Variation/Explore?v=rs148616984;vdb=variation'},
                   {'name': 'ESP',
                    'id': 'rs148616984',
                    'url': 'http://evs.gs.washington.edu/EVS/PopStatsServlet?searchBy=rsID&target=rs148616984&x=0&y=0'},
                   {'name': 'ExAC',
                    'id': 'rs148616984',
                    'url': 'http://exac.broadinstitute.org/awesome?query=rs148616984'}
                   ],
         'wildType': 'L',
         'frequency': 0.000399361,
         'polyphenPrediction': 'benign',
         'polyphenScore': 0.41,
         'siftPrediction': 'deleterious',
         'siftScore': 0.01,
         'somaticStatus': 0,
         'cytogeneticBand': '7p11.2',
         'consequenceType': 'missense',
         'genomicLocation': 'NC_000007.14:g.56058069C>T',
         'sourceType': 'large_scale_study'},
        {'type': 'VARIANT',
         'alternativeSequence': 'S',
         'begin': '83',
         'end': '83',
         'xrefs': [{'name': 'ExAC',
                    'id': 'rs779741978',
                    'url': 'http://exac.broadinstitute.org/awesome?query=rs779741978'}
                   ],
         'wildType': 'A',
         'frequency': 0.0,
         'polyphenPrediction': 'benign',
         'polyphenScore': 0.305,
         'siftPrediction': 'tolerated',
         'siftScore': 0.06,
         'somaticStatus': 0,
         'cytogeneticBand': '7p11.2',
         'consequenceType': 'missense',
         'genomicLocation': 'NC_000007.14:g.56054414G>T',
         'sourceType': 'large_scale_study'},
        {'type': 'VARIANT',
         'description': '[LSS_COSMIC]: primary tissue(s): lung',
         'alternativeSequence': 'V',
         'begin': '292',
         'end': '292',
         'xrefs': [{'name': 'cosmic curated',
                    'id': 'COSM1549991',
                    'url': 'http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=1549991'}
                   ],
         'evidences': [{'code': 'ECO:0000313',
                        'source': {'name': 'cosmic_study',
                                   'id': 'COSU:417',
                                   'url': 'http://cancer.sanger.ac.uk/cosmic/study/overview?study_id=417'}
                        }
                       ],
         'wildType': 'I',
         'polyphenPrediction': 'benign',
         'polyphenScore': 0.021,
         'siftPrediction': 'tolerated',
         'siftScore': 0.1, 'somaticStatus': 0,
         'cytogeneticBand': '7p11.2',
         'consequenceType': 'missense',
         'genomicLocation': '7:g.56058510A>G',
         'sourceType': 'large_scale_study'}
    ]
}


class TestVariants(unittest.TestCase):
    """Test the VAR parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.data = example_uniprot_variants
        self.flatten_uniprot_variants_ebi = flatten_uniprot_variants_ebi
        self.flatten_json = flatten_json
        self.collapse_unique_values = collapse_unique_values

        logging.disable(logging.DEBUG)

    def tearDown(self):
        """Remove testing framework."""

        self.data = None
        self.flatten_uniprot_variants_ebi = None
        self.flatten_json = None
        self.collapse_unique_values = None

        logging.disable(logging.NOTSET)

    def test_flatten_uniprot_variants_ebi(self):
        r = self.flatten_uniprot_variants_ebi(self.data)

        # flattening the accession
        self.assertEqual('P40227', self.data['accession'])
        self.assertEqual('P40227', r.loc[0, 'accession'])

        # flattening the 'xrefs'
        self.assertEqual('1000Genomes', self.data['features'][0]['xrefs'][0]['name'])
        self.assertEqual('1000Genomes', r.loc[0, 'xrefs_name'][0])
        self.assertEqual(['1000Genomes', 'ESP', 'ExAC'], r.loc[0, 'xrefs_name'])
        self.assertTrue(np.isnan(r.loc[1, 'evidences_source_name']))

        # flattening the 'evidences'
        self.assertEqual('cosmic_study', self.data['features'][2]['evidences'][0]['source']['name'])
        self.assertEqual('cosmic_study', r.loc[2, 'evidences_source_name'])

    def test_flatten_json(self):
        modified = {}
        self.flatten_json(self.data['features'][0], modified)
        self.assertIn('1_xrefs_name', modified)
        self.assertIn('2_xrefs_name', modified)
        self.assertIn('3_xrefs_name', modified)
        self.assertIn('1_xrefs_id', modified)
        self.assertIn('2_xrefs_id', modified)
        self.assertIn('3_xrefs_id', modified)
        self.assertEqual('1000Genomes', modified['1_xrefs_name'])
        self.assertEqual('rs148616984', modified['1_xrefs_id'])

    def test_collapse_unique_values(self):
        modified = {}
        self.flatten_json(self.data['features'][0], modified)
        r = self.collapse_unique_values(modified)
        self.assertIn('type', r)
        self.assertIn('xrefs_name', r)
        self.assertIn('xrefs_id', r)
        self.assertIn('xrefs_url', r)
        self.assertIn('1000Genomes', r['xrefs_name'])
        self.assertIn('rs148616984', r['xrefs_id'])
        self.assertEqual('VARIANT', r['type'])


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("prointvar").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestVariants)
    unittest.TextTestRunner(verbosity=2).run(suite)
