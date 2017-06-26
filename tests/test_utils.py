#!/local/bin/python
# -*- coding: utf-8 -*-

import re
import sys
import json
import logging
import requests
import responses
import unittest
import datetime
import pandas as pd
from io import StringIO
from datetime import datetime
from contextlib import contextmanager

try:
    from mock import patch, MagicMock
except ImportError:
    from unittest.mock import patch, MagicMock

from prointvar.utils import check_sequence
from prointvar.utils import compare_sequences
from prointvar.utils import convert_str_to_bool
from prointvar.utils import count_mismatches
from prointvar.utils import current_date
from prointvar.utils import current_time
from prointvar.utils import fetch_from_url_or_retry
from prointvar.utils import flash
from prointvar.utils import get_pairwise_alignment
from prointvar.utils import string_split
from prointvar.utils import get_new_pro_ids
from prointvar.utils import row_selector

from prointvar.config import config as c


@contextmanager
def captured_output():
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err


class FakeDatetime(datetime):
    """A manipulable date replacement"""
    def __new__(cls, *args, **kwargs):
        return datetime.__new__(datetime, *args, **kwargs)


def response_mocker(kwargs, base_url, endpoint_url, status=200,
                    content_type='application/json', post=False, data=None):
    """
    Generates a mocked requests response for a given set of
    kwargs, base url and endpoint url
    """

    url = re.sub('\{\{(?P<m>[a-zA-Z_]+)\}\}', lambda m: "%s" % kwargs.get(m.group(1)),
                      base_url + endpoint_url)
    with responses.RequestsMock() as rsps:
        if post:
            rsps.add(responses.POST, url,
                     body=b'{"data": "some json formatted output"}',
                     status=status, content_type='application/json')
            response = requests.post(url, data=data)

        elif content_type == 'application/json':
            rsps.add(responses.GET, url,
                     body=b'{"data": "some json formatted output"}',
                     status=status, content_type='application/json')
            response = requests.get(url)
        elif content_type == 'text/plain':
            rsps.add(responses.GET, url,
                     body="Some text-based content\n spanning multiple lines",
                     status=status, content_type='text/plain')
            response = requests.get(url)
        else:
            rsps.add(responses.GET, url,
                     body=b"Some other binary stuff...",
                     status=status, content_type='application/octet-stream')
            response = requests.get(url)
    return response


class TestUTILS(unittest.TestCase):
    """Test the utility methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.flash = flash
        self.string_split = string_split
        self.fetch_from_url_or_retry = fetch_from_url_or_retry
        self.seq1 = 'ADEK--XX*PTSV'
        self.seq2 = 'GEDL--XX*PSSV'
        self.seq3 = 'GEDL----*PSSV'
        self.check_sequence = check_sequence
        self.count_mismatches = count_mismatches
        self.compare_sequences = compare_sequences
        self.pairwise_alignment = get_pairwise_alignment
        self.current_time = current_time
        self.current_date = current_date
        self.convert_str_to_bool = convert_str_to_bool
        self.get_new_pro_ids = get_new_pro_ids
        self.row_selector = row_selector

    def tearDown(self):
        """Remove testing framework."""

        self.flash = None
        self.string_split = None
        self.fetch_from_url_or_retry = None
        self.seq1 = None
        self.seq2 = None
        self.seq3 = None
        self.check_sequence = None
        self.count_mismatches = None
        self.compare_sequences = None
        self.pairwise_alignment = None
        self.current_time = None
        self.current_date = None
        self.convert_str_to_bool = None
        self.get_new_pro_ids = None
        self.row_selector = None

    def test_flash(self):
        with captured_output() as (out, err):
            self.flash('Testing the flash method...')
        output = out.getvalue().strip()
        self.assertEqual(output, 'Testing the flash method...')

    def test_string_split(self):
        r1 = self.string_split('118A')
        r2 = self.string_split('118A12')
        r3 = self.string_split('A121')
        r4 = self.string_split('p.ALA53THR')
        r5 = self.string_split('-1')
        self.assertEqual(r1, ['118', 'A'])
        self.assertEqual(r2, ['118', 'A', '12'])
        self.assertEqual(r3, ['A', '121'])
        self.assertEqual(r4, ['p.ALA', '53', 'THR'])
        self.assertEqual(r5, ['-', '1'])

    def test_fetch_from_url_or_retry_get_text(self):
        # mocked requests
        identifier = "P00439"
        base_url = c.http_uniprot
        endpoint_url = "{}.fasta".format(identifier)
        response = response_mocker(kwargs={"P00439.fasta"}, base_url=base_url,
                                   endpoint_url="",
                                   content_type='text/plain')
        self.fetch_from_url_or_retry = MagicMock(return_value=response)
        url = base_url + endpoint_url
        r = self.fetch_from_url_or_retry(url, json=True,
                                         header={'text/plain'}).content
        self.assertEqual(str(r, 'utf-8'),
                         "Some text-based content\n spanning multiple lines")

    def test_fetch_from_url_or_retry_get_json(self):
        # mocked requests
        identifier = "2pah"
        base_url = c.api_pdbe
        endpoint_url = "pdb/entry/summary/"
        response = response_mocker(kwargs={identifier}, base_url=base_url,
                                   endpoint_url=endpoint_url,
                                   content_type='application/json')
        self.fetch_from_url_or_retry = MagicMock(return_value=response)
        url = base_url + endpoint_url + identifier
        r = self.fetch_from_url_or_retry(url, json=True,
                                         header={'application/json'}).json()
        self.assertEqual(r, json.loads('{"data": "some json formatted output"}'))

    def test_fetch_from_url_or_retry_get_binary(self):
                # mocked requests
        identifier = "P00439"
        base_url = c.http_uniprot
        endpoint_url = "{}.fasta".format(identifier)
        response = response_mocker(kwargs={"P00439.fasta"}, base_url=base_url,
                                   endpoint_url="",
                                   content_type='application/octet-stream')
        self.fetch_from_url_or_retry = MagicMock(return_value=response)
        url = base_url + endpoint_url
        r = self.fetch_from_url_or_retry(url, json=True,
                                         header={'application/octet-stream'},
                                         retry_in=None, wait=0,
                                         n_retries=10, stream=False).content
        self.assertEqual(r, b"Some other binary stuff...")

    def test_fetch_from_url_or_retry_post_json(self):
        # mocked requests
        identifier = "1csb, 2pah"
        base_url = c.api_pdbe
        endpoint_url = "pdb/entry/summary/"
        response = response_mocker(kwargs={}, base_url=base_url,
                                   endpoint_url=endpoint_url,
                                   content_type='application/octet-stream',
                                   post=True, data=identifier)
        self.fetch_from_url_or_retry = MagicMock(return_value=response)
        url = base_url + endpoint_url + identifier
        r = self.fetch_from_url_or_retry(url, json=True, post=True, data=identifier,
                                         header={'application/octet-stream'},
                                         retry_in=None, wait=0,
                                         n_retries=10, stream=False).json()
        self.assertEqual(r, json.loads('{"data": "some json formatted output"}'))

    def test_fetch_from_url_or_retry_get_404(self):
        # mocked requests
        identifier = "P00439"
        base_url = c.http_uniprot
        endpoint_url = "{}.fasta".format(identifier)
        response = response_mocker(kwargs={"P00439.fasta"}, base_url=base_url,
                                   endpoint_url="",
                                   content_type='text/plain', status=404)
        self.fetch_from_url_or_retry = MagicMock(return_value=response)
        url = base_url + endpoint_url
        r = self.fetch_from_url_or_retry(url, json=True, header={'text/plain'},
                                         retry_in=None, wait=0,
                                         n_retries=10, stream=False)
        self.assertEqual(r.status_code, 404)
        self.assertFalse(r.ok)

    def test_fetch_from_url_or_retry_get_500(self):
        # mocked requests
        identifier = "P00439"
        base_url = c.http_uniprot
        endpoint_url = "{}.fasta".format(identifier)
        response = response_mocker(kwargs={"P00439.fasta"}, base_url=base_url,
                                   endpoint_url="",
                                   content_type='text/plain', status=500)
        self.fetch_from_url_or_retry = MagicMock(return_value=response)
        url = base_url + endpoint_url
        r = self.fetch_from_url_or_retry(url, json=True, header={'text/plain'},
                                         retry_in=(500, ), wait=1,
                                         n_retries=10, stream=False)
        self.assertEqual(r.status_code, 500)
        self.assertFalse(r.ok)

    def test_check_sequence(self):
        seq = self.check_sequence(self.seq1, gap_symbol='-', new_gap_symbol='X', ambiguous='X')
        self.assertNotEqual(self.seq1, seq)
        self.assertEqual(seq, 'ADEKXXXXXPTSV')
        seq = self.check_sequence(self.seq2, gap_symbol='X', new_gap_symbol='-', ambiguous='-')
        self.assertNotEqual(self.seq2, seq)
        self.assertEqual(seq, 'GEDL-----PSSV')
        seq = self.check_sequence(self.seq3, gap_symbol='-', new_gap_symbol='-', ambiguous='X')
        self.assertNotEqual(self.seq3, seq)
        self.assertEqual(seq, 'GEDL----XPSSV')

    def test_count_mismatches(self):
        r = self.count_mismatches(self.seq1, self.seq2)
        self.assertEqual(r, 5)
        r = self.count_mismatches(self.seq1, self.seq3)
        self.assertEqual(r, 7)
        r = self.count_mismatches(self.seq2, self.seq3)
        self.assertEqual(r, 2)
        seq2 = self.check_sequence(self.seq2, gap_symbol='X', new_gap_symbol='-', ambiguous='-')
        seq3 = self.check_sequence(self.seq3, gap_symbol='-', new_gap_symbol='-', ambiguous='-')
        r = self.count_mismatches(seq2, seq3)
        self.assertEqual(r, 0)

    def test_compare_sequences(self):
        r = self.compare_sequences(self.seq1, self.seq2, permissive=False, n_mismatches=0)
        self.assertEqual(r, False)
        r = self.compare_sequences(self.seq1, self.seq3, permissive=False, n_mismatches=0)
        self.assertEqual(r, False)
        r = self.compare_sequences(self.seq1, self.seq2, permissive=True, n_mismatches=0)
        self.assertEqual(r, False)
        r = self.compare_sequences(self.seq2, self.seq3, permissive=True, n_mismatches=2)
        self.assertEqual(r, True)

    def test_pairwise_alignment(self):
        r1, r2 = self.pairwise_alignment(self.seq1, self.seq1)
        self.assertEqual(r1, self.seq1)
        self.assertEqual(r2, self.seq1)
        self.assertEqual(r1, r2)
        r1, r2 = self.pairwise_alignment(self.seq1, self.seq2, method='global')
        self.assertEqual(r1, 'ADEKXXXXXPTSV')
        self.assertEqual(r2, 'GEDLXXXXXPSSV')
        r1, r2 = self.pairwise_alignment(self.seq1, self.seq2, method='local')
        self.assertEqual(r1, 'ADEKXXXXXPTSV')
        self.assertEqual(r2, 'GEDLXXXXXPSSV')
        nr1 = self.check_sequence(r1, gap_symbol='X', new_gap_symbol='-', ambiguous='-')
        nr2 = self.check_sequence(r2, gap_symbol='X', new_gap_symbol='-', ambiguous='-')
        self.assertEqual(nr1, 'ADEK-----PTSV')
        self.assertEqual(nr2, 'GEDL-----PSSV')

    @patch('datetime.datetime', FakeDatetime)
    def test_current_time(self):
        FakeDatetime.now = classmethod(lambda cls: datetime(2017, 2, 1, 5, 0, 0))
        mocked_time = FakeDatetime.now()
        self.assertEqual(self.current_time(input_time=mocked_time), '01/02/2017 05:00:00')

    @patch('datetime.datetime', FakeDatetime)
    def test_current_date(self):
        FakeDatetime.now = classmethod(lambda cls: datetime(2017, 2, 1, 0, 0, 0))
        mocked_time = FakeDatetime.now()
        self.assertEqual(self.current_date(input_date=mocked_time), '01/02/2017')

    def test_convert_str_to_bool(self):
        self.assertTrue(self.convert_str_to_bool('y'))
        self.assertTrue(self.convert_str_to_bool('Yes'))
        self.assertTrue(self.convert_str_to_bool('YES'))
        self.assertTrue(self.convert_str_to_bool('t'))
        self.assertTrue(self.convert_str_to_bool('True'))
        self.assertTrue(self.convert_str_to_bool('true'))
        self.assertTrue(self.convert_str_to_bool('TRUE'))
        self.assertTrue(self.convert_str_to_bool('1'))
        self.assertFalse(self.convert_str_to_bool('False'))
        self.assertFalse(self.convert_str_to_bool('F'))
        self.assertFalse(self.convert_str_to_bool('NO'))
        self.assertFalse(self.convert_str_to_bool('No'))
        self.assertFalse(self.convert_str_to_bool('0'))

    def test_get_new_pro_ids(self):
        pros = self.get_new_pro_ids()
        asym_id_1, seq_id_1 = next(pros)
        self.assertEqual(asym_id_1, 'A')
        self.assertEqual(seq_id_1, '1')
        asym_id_2, seq_id_2 = next(pros)
        self.assertEqual(asym_id_2, 'A')
        self.assertEqual(seq_id_2, '2')

    def test_row_selector(self):
        data = pd.DataFrame([{'label': '1', 'value': 1},
                             {'label': '2', 'value': 1},
                             {'label': '3', 'value': 2},
                             {'label': '4', 'value': 3},
                             {'label': '5', 'value': 5}])
        d = self.row_selector(data, key='value', value=3, method='equals')
        self.assertEqual(len(d.index), 1)
        d = self.row_selector(data, key='value', value=3, method='diffs')
        self.assertEqual(len(d.index), 4)
        d = self.row_selector(data, key='value', value=None, method='first')
        self.assertEqual(len(d.index), 2)
        d = self.row_selector(data, key='value', value=(2, 3), method='isin')
        self.assertEqual(len(d.index), 2)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("prointvar").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestUTILS)
    unittest.TextTestRunner(verbosity=2).run(suite)
