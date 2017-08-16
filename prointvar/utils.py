#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""

Defines a number of general purpose routines. Some need to be removed or
moved to more appropriated locations.

Fábio Madeira, 2015+

"""

import os
import re
import sys
import time
import json
import logging
import requests
import pandas as pd
from string import digits
from string import ascii_lowercase
from string import ascii_uppercase
from datetime import datetime
from contextlib import suppress
from collections import OrderedDict

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

from prointvar.library import (ASA_Miller, ASA_Wilke, ASA_Sander)
from prointvar.library import aa_codes_1to3_common
from prointvar.library import aa_codes_1to3_extended

logger = logging.getLogger("prointvar")


def string_split(s):
    """
    Splits string on numbers. Example:
    'p.ALA53THR' gets split into ['p.ALA', '53', 'THR']
    """
    return list(filter(None, re.split(r'(\d+)', s)))


def current_date(input_date=None):
    """
    Gets the current date.

    :return: outputs formatted date as 'Day/Month/Year'
    """

    if input_date is None:
        input_date = datetime.now()
    # month abbreviation would be "%b"
    return input_date.strftime("%d/%m/%Y")


def current_time(input_time=None):
    """
    Gets current datetime and time.

    :return: outputs formatted time as 'Day/Month/Year H:M:S'
    """

    if input_time is None:
        input_time = datetime.now()
    return input_time.strftime("%d/%m/%Y %H:%M:%S")


# TODO python 3.5 no longer need this with option 'exist_ok=True'
def create_directory(directory):
    """
    Creates a directory structure if it does not exist.

    :param directory: directory name (expects full path)
    :return: creates a directory if it does not exist yet
    """

    if not os.path.exists(directory):
        os.makedirs(directory)
    return


def lazy_file_remover(filename):
    with suppress(FileNotFoundError):
        os.remove(filename)


def flash(message):
    """
    Flashes a message out.
    :param message: input message str()
    """
    print(str(message))
    sys.stdout.flush()


def logging_out(message, log='', logger=None, verbose=False):
    """
    Generic method to log messages out.

    :param message: string
    :param log: string
    :param logger: base logger or None
    :param verbose: boolean
    """

    if verbose:
        flash(message)
    if logger is not None:
        if log == 'debug':
            logger.debug(message)
        elif log == 'info':
            logger.info(message)
        elif log == 'warning':
            logger.warning(message)
        elif log == 'error':
            logger.error(message)
        elif log == 'critical':
            logger.critical(message)
    return


def convert_str_to_bool(v):
    """
    Converts a string input on boolean.
    """
    return v.lower() in ("yes", "y", "true", "t", "1")


def fetch_from_url_or_retry(url, json=True, header=None, post=False, data=None,
                            retry_in=None, wait=1, n_retries=10, stream=False, **params):
    """
    Fetch an url using Requests or retry fetching it if the server is
    complaining with retry_in error. There is a limit to the number of retries.

    Retry code examples: 429, 500 and 503
    
    :param url: url to be fetched as a string
    :param json: json output
    :param header: dictionary
    :param post: boolean
    :param data: dictionary: only if post is True
    :param retry_in: http codes for retrying
    :param wait: sleeping between tries in seconds
    :param n_retries: number of retry attempts
    :param stream: boolean
    :param params: request.get kwargs.
    :return: url content
    """

    if retry_in is None:
        retry_in = ()
    else:
        assert type(retry_in) is tuple or type(retry_in) is list

    if header is None:
        header = {}
    else:
        assert type(header) is dict

    if json:
        header.update({"Content-Type": "application/json"})
    else:
        if "Content-Type" not in header:
            header.update({"Content-Type": "text/plain"})

    if post:
        if data is not None:
            assert type(data) is dict or type(data) is str
            response = requests.post(url, headers=header, data=data)
        else:
            return None
    else:
        response = requests.get(url, headers=header, params=params, stream=stream)

    if response.ok:
        return response
    elif response.status_code in retry_in and n_retries >= 0:
        time.sleep(wait)
        return fetch_from_url_or_retry(url, json, header, post, data, retry_in, wait,
                                       (n_retries - 1), stream, **params)
    else:
        try:
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            logger.debug('%s: Unable to retrieve %s for %s',
                         response.status_code, url, e)


def store_data(data, path_root, directory, filename, extension="json"):
    """
    Generic routine to store information to the filesystem.

    :param data: data to be printed down
    :param path_root: root of the file path
    :param directory: name of the leading directory name (method)
    :param filename: filename or same identifier
    :param extension: file extension
    """

    local_dir = "{}{}".format(path_root, directory)
    create_directory(local_dir)

    output_file = "{}{}{}.{}".format(path_root, directory, filename, extension)
    if extension == 'json':
        try:
            # indent=4
            string = json.dumps(data, sort_keys=False)
        except ValueError:
            string = data
    elif type(data) is list or type(data) is tuple:
        # assumes that each item is supposed to be printed in a new line
        string = "\n".join(data)
    else:
        string = str(data)

    # write down the content
    outfile = open(output_file, "w")
    outfile.write(string)
    outfile.close()
    return


def load_data(path_root, directory, filename, extension="json"):
    """
    Generic routine to load information from local filesystem.

    :param path_root: root of the file path
    :param directory: name of the leading directory name (method)
    :param filename: filename or same identifier
    :param extension: file extension
    """

    data = None
    input_file = "{}{}{}.{}".format(path_root, directory, filename, extension)
    if os.path.isfile(input_file):
        with open(input_file, "rb") as lines:
            inlines = lines.read()
        if extension == "json":
            try:
                data = json.loads(inlines,
                                  object_pairs_hook=OrderedDict)
            except ValueError:
                data = inlines
        else:
            data = inlines
    return data


def check_sequence(sequence, gap_symbol='-', new_gap_symbol='-', ambiguous='X'):
    """
    Checks an input sequence for uncommon residue symbols.

    :param sequence: (str) protein sequence
    :param gap_symbol: (str) 1-letter symbol for gaps
    :param new_gap_symbol: (str) 1-letter symbol for gaps
    :param ambiguous: (str) 1-letter symbol for ambiguous residues
    :return: returns modified sequence
    """

    new_sequence = "".join([ambiguous if aa not in list(aa_codes_1to3_extended.keys()) else aa
                            for aa in sequence])
    if gap_symbol != new_gap_symbol:
        new_sequence = new_sequence.replace(gap_symbol, new_gap_symbol)

    return new_sequence


def count_mismatches(sequence1, sequence2):
    """
    Counts the number of mismatches between two sequences
    of the same length.

    :param sequence1: sequence 1
    :param sequence2: sequence 2
    :return: The number of mismatches between sequences 1 and 2.
    """
    return sum(i != j for i, j in zip(sequence1, sequence2))


def compare_sequences(sequence1, sequence2, permissive=True, n_mismatches=0):
    """Compares two given sequences in terms of length and sequence content.

    :param sequence1: First sequence
    :param sequence2: Second sequence
    :param permissive: if True it allow sequences to have mismatches
    :param n_mismatches: number of allowed mismatches
    :return: simply a true or false
    :rtype: boolean
    """
    if (sequence1 == sequence2 or
            (permissive and count_mismatches(sequence1, sequence2) <= n_mismatches)):
        return True
    return False


def get_pairwise_alignment(sequence1, sequence2, method="global",
                           gap_open=-10.0, gap_extend=-0.5,
                           gap_beg=-10.0, gap_end=-10.0,
                           matrix=matlist.blosum62):
    """
    Gets a pairwise alignment using a simple method that uses default
    distance matrix (BLOSUM62) and gap extension and opening penalties.

    :param sequence1: input sequence 1
    :param sequence2: input sequence 2
    :param method: pairwise alignment method
    :param gap_open: (float)
    :param gap_extend: (float)
    :param gap_beg: (float)
    :param gap_end: (float)
    :param matrix: biopython matlist object
    :return: returns the two aligned sequences
    """

    # if sequences are exactly the same
    if sequence1 == sequence2:
        return sequence1, sequence2

    sequence1 = check_sequence(sequence1, gap_symbol='-', new_gap_symbol='X')
    sequence2 = check_sequence(sequence2, gap_symbol='-', new_gap_symbol='X')

    if method == "global":
        aligns = pairwise2.align.globalds(sequence1, sequence2, matrix,
                                          gap_open, gap_extend)
    elif method == "local":
        # aligns = pairwise2.align.localds(sequence1, sequence2, blosum,
        #                                  gap_open, gap_extend)
        aligns = pairwise2.align.localdd(sequence1, sequence2, matrix,
                                         gap_open, gap_extend, gap_beg, gap_end)
    else:
        raise Exception("Error: method {method} not currently implemented...")

    top_align = aligns[0]
    sequence1, sequence2, score, begin, end = top_align

    try:
        assert len(sequence1) == len(sequence2)
    except:
        raise AssertionError("Error: Something went wrong when trying to align sequences...")

    return sequence1, sequence2


def get_rsa_class(rsa):
    """
    Gets a class based on the RSA value

    :param rsa: (float) RSA score or string
    :return: RSA
    """
    rsa_class = '-'
    try:
        rsa = float(rsa)
        # surface is rsa >= 25% (value = 2)
        # exposed is rsa >= 5% and rsa < 25% (value = 1)
        # core is rsa < 5% (value = 0)
        if rsa >= 25.0:
            rsa_class = 'Surface'
        elif 5.0 <= rsa < 25.0:
            rsa_class = 'Part. Exposed'
        else:
            rsa_class = 'Core'
    except ValueError:
        # returns a string
        pass
    return rsa_class


def compute_rsa(acc, resname, method="Sander"):
    """
    Computes Relative Solvent Accessibility (RSA) from an input
    DSSP ACC value, and according to ASA standard values.

    :param acc: (int) DSSP ACC
    :param resname: single letter residue name
    :param method: name of the method
    :return: (float) RSA value
    """

    rsa = ""
    try:
        acc = float(acc)
    except ValueError:
        return rsa

    try:
        assert len(resname) == 1
    except AssertionError:
        return rsa

    if method == "Miller":
        sasa = ASA_Miller
    elif method == "Wilke":
        sasa = ASA_Wilke
    elif method == "Sander":
        sasa = ASA_Sander
    else:
        raise ValueError("Method {} is not implemented...".format(method))

    try:
        rsa = round((acc / sasa[aa_codes_1to3_extended[resname]] * 100), 3)
    except KeyError:
        return rsa

    return rsa


def row_selector(data, key=None, value=None, method="isin"):
    """
    Generic method to filter columns
    :param data: pandas DataFrame
    :param key: pandas DataFrame column name
    :param value: value(s) to be looked for
    :param method: operator method
    :return:
    """

    table = data
    assert type(table) is pd.core.frame.DataFrame
    if ((key is not None and value is not None) or
            (key is not None and method == 'first')):
        assert type(key) is str
        assert type(method) is str
        if key in table:
            if method == "isin":
                assert type(value) is tuple
                table = table.loc[table[key].isin(value)]
            elif method == "equals":
                # assert type(values) is str
                table = table.loc[table[key] == value]
            elif method == "diffs":
                table = table.loc[table[key] != value]
            elif method == "first":
                value = table[key].iloc[0]
                table = table.loc[table[key] == value]
        else:
            logger.debug("%s not in the DataFrame...", key)

    if table.empty:
        message = 'Your filters resulted in an empty DataFrame...'
        logger.debug(message)
        raise ValueError(message)

    return table


def get_new_pro_ids():
    """
    Both DSSP and arpeggio work with single-letter characters such
        as [A-Z], [a-z], [0-9]; which gives a total of 62 chars
    Now this generator returns a new asym_id + seq_id:
        asym_id: 62 chars
        seq_id: from '1' to '9999' [up to 4 characters]
    :return: (str) chainid and (str) seqid
    """
    chain_ids = ascii_uppercase + ascii_lowercase + digits
    for char in chain_ids:
        for seqid in range(1, 10000, 1):
            yield char, str(seqid)


if __name__ == '__main__':
    pass
