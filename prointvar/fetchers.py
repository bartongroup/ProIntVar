#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""

This is where methods that fetch data from 'bio' resources are
implemented.

Fábio Madeira, 2015+

"""

import os
import gzip
import shutil
import pickle
import logging

from prointvar.utils import fetch_from_url_or_retry

from prointvar.config import config

logger = logging.getLogger("prointvar")


class BioFetcher(object):
    def __init__(self, url, cached=False, cache_output=None, **kwargs):
        """
        :param url: (str) Full web-address
        :param cached: (boolean) if True, stores a pickle file locally
        :param cache_output: (str) file name if 'cached=True'
        """

        self.url = url
        self.cached = cached
        self.cache_output = cache_output
        self.pickled = os.path.join(config.db_root, config.db_pickled,
                                    self.cache_output)
        self.kwargs = kwargs
        self.response = None
        self._fetch()

    def _fetch(self):

        if self.cached and os.path.isfile(self.pickled):
            self.response = pickle.load(open(self.pickled, 'rb'))
        else:
            self.response = fetch_from_url_or_retry(self.url, **self.kwargs)
            if self.cached and self.response is not None:
                with open(self.pickled, 'wb') as output:
                    pickle.dump(self.response, output, -1)
        return self.response


def fetch_uniprot_variants_ebi(identifier, cached=False, retry_in=(429,)):
    """
    Queries the EBI Proteins API for variants.
    based on UniProt identifiers (e.g. O15294).

    :param identifier: UniProt ID
    :param cached: (boolean) if True, stores a pickle file locally
    :param retry_in: http code for retrying connections
    :return: pandas table dataframe
    """

    url_root = config.api_proteins
    url_enpoint = "variation/"
    url = url_root + url_enpoint + identifier
    b = BioFetcher(url=url, cached=cached,
                   cache_output="{}_vars.pkl".format(identifier),
                   json=True, retry_in=retry_in)
    return b.response


def fetch_best_structures_pdbe(identifier, cached=False, retry_in=(429,)):
    """
    Queries the PDBe API SIFTS mappings best_structures endpoint.

    :param identifier: UniProt ID
    :param cached: (boolean) if True, stores a pickle file locally
    :param retry_in: http code for retrying connections
    :return: response object
    """

    url_root = config.api_pdbe
    url_enpoint = "mappings/best_structures/"
    url = url_root + url_enpoint + identifier
    b = BioFetcher(url=url, cached=cached,
                   cache_output="{}_bs.pkl".format(identifier),
                   json=True, retry_in=retry_in)
    return b.response


def fetch_summary_properties_pdbe(identifier, cached=False, retry_in=(429,)):
    """
    Queries the PDBe API to get summary properties.

    :param identifier: PDB ID
    :param cached: (boolean) if True, stores a pickle file locally
    :param retry_in: http code for retrying connections
    :return: response object
    """

    url_root = config.api_pdbe
    url_enpoint = "pdb/entry/summary/"
    url = url_root + url_enpoint + identifier
    b = BioFetcher(url=url, cached=cached,
                   cache_output="{}_sp.pkl".format(identifier),
                   json=True, retry_in=retry_in)
    return b.response


def get_preferred_assembly_id(identifier):
    """
    Gets the preferred assembly id for the given PDB ID, from the PDBe API.

    :param identifier: PDB ID
    :return: (str)
    """

    # getting the preferred biological assembly from the PDBe API
    pref_assembly = "1"
    try:
        data = fetch_summary_properties_pdbe(identifier)
    except Exception as e:
        logger.debug("Something went wrong for %s... %s", identifier, e)
    try:
        if data is not None:
            data = data.json()
            nassemblies = data[identifier][0]["assemblies"]
            if len(nassemblies) > 1:
                for entry in nassemblies:
                    if entry["preferred"]:
                        pref_assembly = entry["assembly_id"]
                        break
            else:
                pref_assembly = data[identifier][0]["assemblies"][0]["assembly_id"]
    except Exception as e:
        pref_assembly = "1"
        logger.debug("Something went wrong for %s... %s", identifier, e)

    bio_best = str(pref_assembly)
    return bio_best


class BioDownloader(object):
    def __init__(self, url, outputfile, decompress=True, override=False):
        """
        :param url: (str) Full web-address
        :param outputfile: (str) Output filename
        :param decompress: (boolean) Decompresses the file
        :param override: (boolean) Overrides any existing file, if available
        """

        self.url = url
        self.outputfile = outputfile
        self.outputfile_origin = outputfile
        self.decompress = decompress
        self.override = override

        if self.decompress:
            if self.outputfile_origin.endswith('.gz'):
                self.outputfile = self.outputfile_origin.rstrip('.gz')

        if not os.path.exists(self.outputfile) or self.override:
            self._download()
            if self.outputfile_origin.endswith('.gz') and self.decompress:
                self._decompress()
        else:
            logger.info("%s already available...", self.outputfile)

    def _download(self):
        try:
            try:
                import urllib.request
                from urllib.error import URLError, HTTPError
                with urllib.request.urlopen(self.url) as response, \
                        open(self.outputfile_origin, 'wb') as outfile:
                    shutil.copyfileobj(response, outfile)
            except (AttributeError, ImportError):
                import urllib
                urllib.urlretrieve(self.url, self.outputfile_origin)
        except (URLError, HTTPError, IOError, Exception) as e:
            logger.debug("Unable to retrieve %s for %s", self.url, e)

    def _decompress(self):
        with gzip.open(self.outputfile_origin, 'rb') as infile, \
                open(self.outputfile, 'wb') as outfile:
            shutil.copyfileobj(infile, outfile)
            os.remove(self.outputfile_origin)
            logger.info("Decompressed %s to %s",
                        self.outputfile_origin, self.outputfile)


def download_structure_from_pdbe(identifier, pdb=False, bio=False, override=False):
    """
    Downloads a structure from the PDBe to the filesystem.

    :param identifier: (str) PDB ID
    :param pdb: (boolean) PDB formatted if True, otherwise mmCIF format
    :param bio: (boolean) if true downloads the preferred Biological Assembly
    :param override: (boolean)
    :return: (side effects)
    """

    if pdb:
        filename = "{}.pdb".format(identifier)
    else:
        if bio:
            filename = "{}_bio.cif.gz".format(identifier)
        else:
            filename = "{}.cif".format(identifier)

    outputfile = os.path.join(config.db_root, config.db_pdbx, filename)
    os.makedirs(os.path.join(config.db_root, config.db_pdbx), exist_ok=True)

    if pdb:
        url_endpoint = "entry-files/download/pdb{}.ent".format(identifier)
    else:
        if bio:
            # atom lines only?
            # url_endpoint = ("static/entry/download/"
            #                "{}-assembly-{}_atom_site.cif.gz".format(identifier, pref))
            pref = get_preferred_assembly_id(identifier=identifier)
            url_endpoint = ("static/entry/download/"
                            "{}-assembly-{}.cif.gz".format(identifier, pref))
        else:
            # original mmCIF?
            # url_endpoint = "entry-files/download/{}.cif".format(pdbid)
            url_endpoint = "entry-files/download/{}_updated.cif".format(identifier)

    url_root = config.http_pdbe
    url = url_root + url_endpoint
    BioDownloader(url=url, outputfile=outputfile,
                  decompress=True, override=override)
    return


def download_sifts_from_ebi(identifier, override=False):
    """
    Downloads a SIFTS xml from the EBI FTP to the filesystem.

    :param identifier: (str) PDB ID
    :param override: (boolean)
    :return: (side effects)
    """

    filename = "{}.xml.gz".format(identifier)
    outputfile = os.path.join(config.db_root, config.db_sifts, filename)
    os.makedirs(os.path.join(config.db_root, config.db_sifts), exist_ok=True)

    url_root = config.ftp_sifts
    url_endpoint = "{}.xml.gz".format(identifier)
    url = url_root + url_endpoint
    BioDownloader(url=url, outputfile=outputfile,
                  decompress=True, override=override)
    return


def download_data_from_uniprot(identifier, file_format="fasta", override=False):
    """
    Downloads a UniProt fasta, gff or txt to the filesystem.

    :param identifier: (str) UniProt ID
    :param file_format: (str) endpoint
    :param override: (boolean)
    :return: (side effects)
    """

    file_format = file_format.lstrip('.')
    if file_format in ['txt', 'fasta', 'gff']:
        filename = "{}.{}".format(identifier, file_format)
        outputfile = os.path.join(config.db_root, config.db_uniprot, filename)
        os.makedirs(os.path.join(config.db_root, config.db_uniprot), exist_ok=True)

        url_root = config.http_uniprot
        url_endpoint = "{}.{}".format(identifier, file_format)
        url = url_root + url_endpoint
        BioDownloader(url=url, outputfile=outputfile,
                      decompress=True, override=override)
    else:
        raise ValueError("File format {} is not currently implemented..."
                         "".format(file_format))
    return


def download_alignment_from_cath(identifier, max_sequences=200, override=False):
    """
    Downloads a MSA in fasta format from CATH to the filesystem.

    :param identifier: (str) CATH ID (<Superfamily>_<Funfam>)
    :param max_sequences: (str) Maximum number of sequences (default = 200)
    :param override: (boolean)
    :return: (side effects)
    """

    if '_' in identifier:
        filename = "{}.fasta".format(identifier)
        superfamily, funfam = identifier.split('_')[0], identifier.split('_')[1]
        outputfile = os.path.join(config.db_root, config.db_cath, filename)
        os.makedirs(os.path.join(config.db_root, config.db_cath), exist_ok=True)

        url_root = config.http_cath
        url_endpoint = ("superfamily/{}/funfam/{}/files/seed_alignment.fasta"
                        "?max_sequences={}".format(superfamily, funfam,
                                                   max_sequences))
        url = url_root + url_endpoint
        BioDownloader(url=url, outputfile=outputfile,
                      decompress=True, override=override)
    else:
        raise ValueError("Expected CATH  ID but got {}..."
                         "".format(identifier))
    return


def download_alignment_from_pfam(identifier, alignment_size="seed",
                                 override=False):
    """
    Downloads a MSA in Stockholm format from Pfam to the filesystem.

    :param identifier: (str) PFam ID
    :param alignment_size: (str) either "seed" or "full"
    :param override: (boolean)
    :return: (side effects)
    """

    filename = "{}.sth".format(identifier)
    outputfile = os.path.join(config.db_root, config.db_pfam, filename)
    os.makedirs(os.path.join(config.db_root, config.db_pfam), exist_ok=True)

    url_root = config.http_pfam
    url_endpoint = ("family/{}/alignment/{}"
                    "".format(identifier, alignment_size))
    url = url_root + url_endpoint
    BioDownloader(url=url, outputfile=outputfile,
                  decompress=True, override=override)
    return


if __name__ == '__main__':
    pass
