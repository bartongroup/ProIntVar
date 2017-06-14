#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""

This is where methods that fetch data from 'bio' resources are
implemented.

Fábio Madeira, 2015+

"""

from prointvar.utils import fetch_from_url_or_retry

from prointvar.config import config


def fetch_best_structures_pdbe(identifier, retry_in=(429,)):
    """
    Queries the PDBe API SIFTS mappings best_structures endpoint.

    :param identifier: UniProt ID
    :param retry_in: http code for retrying connections
    :return: dictionary
    """

    pdbe_endpoint = "mappings/best_structures/"
    url = config.api_pdbe + pdbe_endpoint + identifier

    return fetch_from_url_or_retry(url, retry_in=retry_in, json=True)


if __name__ == '__main__':
    pass
