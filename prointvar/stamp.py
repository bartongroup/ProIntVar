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
import logging
import pandas as pd
from io import StringIO

from prointvar.library import stamp_types

logger = logging.getLogger("prointvar")


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
    for col in table:
        if col in stamp_types:
            try:
                table[col] = table[col].astype(stamp_types[col])
            except ValueError:
                # there are some NaNs in there
                pass

    if table.empty:
        raise ValueError('{} resulted in an empty DataFrame...'.format(inputfile))

    return table


if __name__ == '__main__':
    pass
