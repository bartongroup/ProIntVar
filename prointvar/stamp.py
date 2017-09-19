#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This defines the methods that work with
STAMP - Structural Alignment of Proteins.

Fábio Madeira, 2017+

"""

# STAMP file types:
# blocks format
# Domain description file

import os
import re
import logging
import pandas as pd
from io import StringIO

from prointvar.pdbx import PDBXwriter

from prointvar.utils import constrain_column_types
from prointvar.library import stamp_types

from prointvar.config import config

logger = logging.getLogger("prointvar")



def parse_stamp_domain_definitions_from_line(string):
    """
    Parses STAMP domain definitions from a line String.

    STAMP domain definition format:
      <file_path> <domain_id> { <residue_ranges> }
        where <residue_ranges> =>
            <CHAIN> <RES> <INSCODE> TO <CHAIN> <RES> <INSCODE>

        <INSCODE> == '_' if None

        e.g.
            1hww.pdb 1hwwA { A 648 _ TO A 927 _ }

            1hww.pdb 1hwwB { A 648 b TO A 650 _  A 655 _ TO A 927 _ }

    :param string: input string
    :return: returns a dictionary with parsed items
    """

    string = string.replace('}', '')
    dom_info = string.split('{')[0]
    ranges = string.split('{')[1]
    start = []
    start_inscode = []
    start_chain = []
    end = []
    end_inscode = []
    end_chain = []
    pattern = re.compile("([A-Z]) ([0-9])+ ([A-Za-z0-9_?]) TO "
                         "([A-Z]) ([0-9])+ ([A-Za-z0-9_?])")
    match = re.finditer(pattern, ranges, flags=0)
    for m in match:
        m = m.group()
        start_chain.append(m.split()[0])
        start.append(m.split()[1])
        start_inscode.append(m.split()[2])
        end_chain.append(m.split()[4])
        end.append(m.split()[5])
        end_inscode.append(m.split()[6])

    info = {'path': dom_info.split()[0],
            'domain_id': dom_info.split()[1],
            'start_chain': tuple(start_chain),
            'end_chain': tuple(end_chain),
            'start': tuple(start),
            'end': tuple(end),
            'start_inscode': tuple(start_inscode),
            'end_inscode': tuple(end_inscode)}
    return info


def parse_stamp_domain_definitions_from_file(inputfile):
    """
    Parses a STAMP domain definitions file.

    :param inputfile: path to input file
    :return: pandas DataFrame
    """

    lines = []
    with open(inputfile, 'r') as inlines:
        for line in inlines:
            if not line.startswith('% STAMP'):
                info = parse_stamp_domain_definitions_from_line(line)
                lines.append(info)
    return pd.DataFrame(lines)


def get_stamp_domain_line(data, index=0):
    """
    Returns a STAMP domain definition line.

    STAMP domain definition format:
      <file_path> <domain_id> { <residue_ranges> }
        where <residue_ranges> =>
            <CHAIN> <RES> <INSCODE> TO <CHAIN> <RES> <INSCODE>

        <INSCODE> == '_' if None

        e.g.
            1hww.pdb 1hwwA { A 648 _ TO A 927 _ }

            1hww.pdb 1hwwB { A 648 b TO A 650 _  A 655 _ TO A 927 _ }

    :param data: pandas DataFrame object
    :param index: (int) index
    :return: returns a STAMP domain-formatted line
    """

    table = data
    ix = index
    keys = ["path", "domain_id", "start", "end", "start_inscode", "end_inscode",
            "start_chain", "end_chain"]

    if set(keys).issubset(table.columns):

        path = table.loc[ix, "path"]
        domain = table.loc[ix, "domain_id"]
        start_chain = table.loc[ix, "start_chain"]
        start = table.loc[ix, "start"]
        start_inscode = table.loc[ix, "start_inscode"]
        end_chain = table.loc[ix, "end_chain"]
        end = table.loc[ix, "end"]
        end_inscode = table.loc[ix, "end_inscode"]

        assert type(start) is tuple
        assert type(start_chain) is tuple
        assert type(start_inscode) is tuple
        assert type(end) is tuple
        assert type(end_chain) is tuple
        assert type(end_inscode) is tuple

        string = []
        for i, j, k, l, a, b in zip(start, start_inscode,
                                    end, end_inscode,
                                    start_chain, end_chain):
            if j == " " or j == "?":
                j = "_"
            if l == " " or l == "?":
                l = "_"
            string.append("%s %s %s TO %s %s %s" % (a, i, j, b, k, l))
        string = "{ %s }" % (" ".join(string))
        domain_definition = " ".join([path, domain, string])

    else:
        message = "Expected columns not found in the Pandas Table..."
        logger.debug(message)
        raise ValueError(message)

    return domain_definition


def write_stamp_domain_definitions_from_table(outputfile, data, override=False):
    """
    Generic method that writes STAMP domain definitions from a pre-formatted
    Pandas DataFrame.

    :param outputfile: path to the PDB file
    :param data: pandas DataFrame object
    :param override: boolean
    :return: (side effects) writes to file
    """

    domain_lines = ['% STAMP domains file generated by ProIntVar']
    for ix in data.index:
        domain_lines.append(get_stamp_domain_line(data=data, index=ix))

    # write the final output
    if not os.path.exists(outputfile) or override:
        with open(outputfile, 'w') as outlines:
            outlines.write("\n".join(domain_lines) + "\n")
    else:
        logger.info("Domain definitions for %s already available...", outputfile)
    return


def parse_stamp_scan_scores_from_file(inputfile, excluded=()):
    """
    Parse STAMP SCAN mode scores.

    :param inputfile: path to the mmCIF file
    :param excluded: option to exclude mmCIF columns
    :return: returns a pandas DataFrame
    """

    logger.info("Parsing STAMP SCAN scores from lines...")

    # example lines with some problems
    """
    (...)
         Domain1         Domain2          Fits  Sc      RMS   Len1 Len2 Align Fit   Eq. Secs    %I    %S     P(m)
    Scan 2uuar_          1ekcc_             1   4.782   1.101   73   50   65   45   36    1  41.67  77.78 6.09e-07
    Scan 2uuar_          1ekch_             1   4.684   1.091   73   50   64   44   37    1  40.54  72.97 5.41e-06
    Scan 2uuar_          1fjfr_             1   9.648   0.526   73   73   73   73   73    7  98.63  93.15 6.57e-71
    Scan 2uuar_          1fjgr_             1   9.596   0.511   73   73   73   73   73    6  98.63  90.41 6.57e-71
    Scan 2uuar_          1fkar_             1   4.028   1.082   73   50   62   39   33    1  45.45  72.73 1.56e-07
    """

    if not os.path.isfile(inputfile):
        raise IOError("{} not available or could not be read...".format(inputfile))

    # parsing atom lines
    header = ["Domain1", "Domain2", "Fits", "Sc", "RMS", "A_Len", "B_Len",
              "Align_Len", "N_Fit", "N_Equiv", "N_SS_Equiv", "PID", "SS_PID", "Pm"]
    lines = []
    with open(inputfile) as inlines:
        for line in inlines:
            if line.startswith("Scan"):
                lines.append(line.lstrip("Scan "))
    lines = "".join(lines)

    all_str = {key: str for key in header}
    table = pd.read_table(StringIO(lines), delim_whitespace=True, low_memory=False,
                          names=header, compression=None, converters=all_str,
                          keep_default_na=False)

    # excluding columns
    if excluded is not None:
        assert type(excluded) is tuple
        try:
            table = table.drop(list(excluded), axis=1)
        except ValueError:
            # most likely theses are not in there
            pass

    # enforce some specific column types
    table = constrain_column_types(table, stamp_types)

    if table.empty:
        raise ValueError('{} resulted in an empty DataFrame...'.format(inputfile))

    return table


if __name__ == '__main__':
    pass
