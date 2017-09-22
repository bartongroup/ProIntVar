#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""

This defines the methods that work with SIFTS files.

Fábio Madeira, 2017+

"""

import os
import json
import logging
import pandas as pd
from lxml import etree
from collections import OrderedDict

from prointvar.utils import row_selector
from prointvar.utils import constrain_column_types
from prointvar.library import sifts_types

logger = logging.getLogger("prointvar")


def parse_sifts_residues_from_file(inputfile, excluded=(),
                                   add_regions=True, add_dbs=False):
    """
    Parses the residue fields of a SIFTS XML file.

    :param inputfile: path to the SIFTS file
    :param excluded: option to exclude SIFTS dbSources
    :param add_regions: boolean
    :param add_dbs: boolean
    :return: returns a pandas DataFrame
    """

    logger.info("Parsing SIFTS residues from lines...")

    # example lines with some problems
    """
    <?xml version='1.0' encoding='UTF-8' standalone='yes'?>
    <entry xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:align="http://www.ebi.ac.uk/pdbe/docs/sifts/alignment.xsd" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:data="http://www.ebi.ac.uk/pdbe/docs/sifts/dataTypes.xsd" xmlns="http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd" dbSource="PDBe" dbVersion="1.3" dbCoordSys="PDBe" dbAccessionId="2pah" dbEntryVersion="2011-07-13" date="2014-07-21" xsi:schemaLocation="http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd">
      <rdf:RDF>
        <rdf:Description rdf:about="self">
          <dc:rights rdf:resource="http://pdbe.org/sifts">
            Copyright notice: (c) 2004-2013, EMBL-EBI, PDBe-UniProt
            Jose M. Dana, Paul Gane, Jie Luo, Glen van Ginkel, Claire O'Donovan, Maria J. Martin, Sameer Velankar.
            The information included is supplied as-is, under the terms and conditions of the licence agreement.
        </dc:rights>
        </rdf:Description>
      </rdf:RDF>
      <listDB>
        <db dbSource="Pfam" dbCoordSys="UniProt" dbVersion="27.0"/>
        <db dbSource="InterPro" dbCoordSys="UniProt" dbVersion="48.0"/>
        <db dbSource="CATH" dbCoordSys="PDBresnum" dbVersion="3.5.0"/>
        <db dbSource="EC" dbCoordSys="UniProt" dbVersion="30.14"/>
        <db dbSource="UniProt" dbCoordSys="UniProt" dbVersion="2014.08"/>
        <db dbSource="SCOP" dbCoordSys="PDBresnum" dbVersion="1.75"/>
        <db dbSource="GO" dbCoordSys="UniProt" dbVersion="20140708"/>
        <db dbSource="PDB" dbCoordSys="PDBresnum" dbVersion="30.14"/>
      </listDB>
      <entity type="protein" entityId="A">
        <segment segId="2pah_A_1_335" start="1" end="335">
          <listResidue>
            <residue dbSource="PDBe" dbCoordSys="PDBe" dbResNum="1" dbResName="VAL">
              <crossRefDb dbSource="PDB" dbCoordSys="PDBresnum" dbAccessionId="2pah" dbResNum="118" dbResName="VAL" dbChainId="A"/>
              <crossRefDb dbSource="UniProt" dbCoordSys="UniProt" dbAccessionId="P00439" dbResNum="118" dbResName="V"/>
              <crossRefDb dbSource="CATH" dbCoordSys="PDBresnum" dbAccessionId="1.10.800.10" dbResNum="118" dbResName="VAL" dbChainId="A"/>
              <crossRefDb dbSource="SCOP" dbCoordSys="PDBresnum" dbAccessionId="42581" dbResNum="118" dbResName="VAL" dbChainId="A"/>
              <crossRefDb dbSource="NCBI" dbCoordSys="UniProt" dbAccessionId="9606" dbResNum="118" dbResName="V"/>
              <crossRefDb dbSource="InterPro" dbCoordSys="UniProt" dbAccessionId="IPR001273" dbResNum="118" dbResName="V" dbEvidence="G3DSA:1.10.800.10"/>
              <crossRefDb dbSource="InterPro" dbCoordSys="UniProt" dbAccessionId="IPR019774" dbResNum="118" dbResName="V" dbEvidence="PS51410"/>
              <crossRefDb dbSource="InterPro" dbCoordSys="UniProt" dbAccessionId="IPR019774" dbResNum="118" dbResName="V" dbEvidence="SSF56534"/>
              <crossRefDb dbSource="InterPro" dbCoordSys="UniProt" dbAccessionId="IPR019774" dbResNum="118" dbResName="V" dbEvidence="SSF56534"/>
              <crossRefDb dbSource="InterPro" dbCoordSys="UniProt" dbAccessionId="IPR005961" dbResNum="118" dbResName="V" dbEvidence="TIGR01268"/>
              <crossRefDb dbSource="InterPro" dbCoordSys="UniProt" dbAccessionId="IPR001273" dbResNum="118" dbResName="V" dbEvidence="PTHR11473"/>
              <residueDetail dbSource="PDBe" property="codeSecondaryStructure">T</residueDetail>
              <residueDetail dbSource="PDBe" property="nameSecondaryStructure">loop</residueDetail>
            </residue>
            (...)
    """

    if not os.path.isfile(inputfile):
        raise IOError("{} not available or could not be read...".format(inputfile))

    # parse regions first
    if add_regions:
        regions = parse_sifts_regions_from_file(inputfile=inputfile, excluded=excluded)
    if add_dbs:
        dbs = parse_sifts_dbs_from_file(inputfile=inputfile, excluded=excluded)

    # parsing sifts segments
    try:
        tree = etree.parse(inputfile)
    except etree.XMLSyntaxError:
        raise IOError("{} not available or could not be read...".format(inputfile))
    root = tree.getroot()
    namespace = root.nsmap[None]
    namespace_map = {'ns': namespace}
    cross_reference = "{{{}}}crossRefDb".format(namespace)
    residue_detail = "{{{}}}residueDetail".format(namespace)

    rows = []
    # Entities
    for entity_list in root.iterfind('.//ns:entity[@type="protein"]',
                                     namespaces=namespace_map):

        entity_id = entity_list.attrib['entityId']
        # Entities : Segments
        for segment in entity_list:

            # Entities : Segments : Residues
            for list_residue in segment.iterfind('.//ns:listResidue',
                                                 namespaces=namespace_map):
                for residue in list_residue:
                    # get residue annotations
                    residue_annotation = OrderedDict()
                    # key, value pairs
                    for k, v in residue.attrib.items():
                        # skipping dbSource
                        if k == 'dbResNum':
                            # adding to the dictionary
                            # residue_annotation[k] = v
                            resnum = int(v)
                        else:
                            continue

                    # parse extra annotations for each residue
                    for annotation in residue:
                        for k, v in annotation.attrib.items():
                            # crossRefDb entries
                            if annotation.tag == cross_reference:
                                if annotation.attrib["dbSource"] not in excluded:
                                    # skipping various fields
                                    if k == 'dbSource' or k == 'dbCoordSys':
                                        continue
                                    elif (annotation.attrib["dbSource"] != "PDB" and
                                                annotation.attrib["dbSource"] != "UniProt"):
                                        if k == 'dbResName' or k == 'dbResNum' or k == 'dbChainId':
                                            continue
                                    # elif annotation.attrib["dbSource"] == "PDB" and k == "dbAccessionId":
                                    #     continue

                                    # adding a new column with the regionId from the 'regions'
                                    if add_regions:
                                        if k == "dbAccessionId":
                                            source = annotation.attrib["dbSource"]
                                            if source in regions[entity_id]:
                                                keys = regions[entity_id][source]
                                                for key in [k for k in keys]:
                                                    entry = regions[entity_id][source][key]
                                                    if v == entry["dbAccessionId"]:
                                                        start = int(entry["start"])
                                                        end = int(entry["end"])
                                                        if resnum in range(start, end + 1, 1):
                                                            nk = "{}_regionId".format(source)
                                                            residue_annotation[nk] = key
                                                            nk = "{}_regionStart".format(source)
                                                            residue_annotation[nk] = start
                                                            nk = "{}_regionEnd".format(source)
                                                            residue_annotation[nk] = end
                                                            nk = "{}_regionResNum".format(source)
                                                            residue_annotation[nk] = resnum

                                    if add_dbs:
                                        source = annotation.attrib["dbSource"]
                                        if source in dbs:
                                            nk = "{}_dbVersion".format(source)
                                            residue_annotation[nk] = dbs[source]['dbVersion']

                                    # renaming all keys with dbSource prefix
                                    k = "{}_{}".format(
                                        annotation.attrib["dbSource"], k)
                                    if k == "{}_1".format(annotation.attrib["dbSource"]):
                                        continue

                            if annotation.tag == residue_detail:
                                if 'PDB' not in excluded:
                                    k = "PDB_{}".format(annotation.attrib["property"])
                                    # value is the text field in the XML
                                    v = annotation.text

                            # adding to the dictionary
                            if "_" in k:
                                try:
                                    if v in residue_annotation[k]:
                                        continue
                                    if k != "PDB_Annotation":
                                        residue_annotation[k].append(v)
                                    else:
                                        residue_annotation[k] = v
                                except KeyError:
                                    residue_annotation[k] = v
                                except AttributeError:
                                    residue_annotation[k] = [residue_annotation[k]]
                                    residue_annotation[k].append(v)
                                except TypeError:
                                    # bool column for annotation
                                    residue_annotation[k] = v
                        if 'PDB' not in excluded and "PDB_Annotation" not in residue_annotation:
                            residue_annotation["PDB_Annotation"] = "Observed"

                        # adding the entity_id
                        if 'PDB_entityId' not in residue_annotation:
                            residue_annotation['PDB_entityId'] = entity_id

                    rows.append(residue_annotation)

    table = pd.DataFrame(rows)

    # enforce some specific column types
    table = constrain_column_types(table, sifts_types)

    for c in list(table):
        if '_regionId' in c:
            table[c] = table[c].fillna('-').astype(str)
        elif '_regionStart' in c or '_regionEnd' in c:
            table[c] = table[c].fillna(0).astype(int)

    if table.empty:
        raise ValueError('{} resulted in an empty DataFrame...'.format(inputfile))

    return table


def parse_sifts_dbs_from_file(inputfile, excluded=()):
    """
    Parses the residue fields of a SIFTS XML file.

    :param inputfile: path to the SIFTS file
    :param excluded: option to exclude SIFTS dbSources
    :return: returns a nested dictionary
    """

    logger.info("Parsing SIFTS dbs from lines...")

    if not os.path.isfile(inputfile):
        raise IOError("{} not available or could not be read...".format(inputfile))

    # parsing sifts segments
    try:
        tree = etree.parse(inputfile)
    except etree.XMLSyntaxError:
        raise IOError("{} not available or could not be read...".format(inputfile))
    root = tree.getroot()
    namespace = root.nsmap[None]
    namespace_map = {'ns': namespace}

    # DB entries and versions
    sifts_object = OrderedDict()
    for db_list in root.iterfind('.//ns:listDB', namespaces=namespace_map):
        for db in db_list:
            source = ''
            db_entries = OrderedDict()
            for k, v in db.attrib.items():
                db_entries[k] = v
                if k == 'dbSource' and v not in excluded:
                    source = v
            if source not in excluded and source != '':
                sifts_object[source] = db_entries

    return sifts_object


def parse_sifts_regions_from_file(inputfile, excluded=()):
    """
    Parses the residue fields of a SIFTS XML file.

    :param inputfile: path to the SIFTS file
    :param excluded: option to exclude SIFTS dbSources
    :return: returns a nested dictionary
    """

    logger.info("Parsing SIFTS regions from lines...")

    if not os.path.isfile(inputfile):
        raise IOError("{} not available or could not be read...".format(inputfile))

    # parsing sifts segments
    try:
        tree = etree.parse(inputfile)
    except etree.XMLSyntaxError:
        raise IOError("{} not available or could not be read...".format(inputfile))
    root = tree.getroot()
    namespace = root.nsmap[None]
    namespace_map = {'ns': namespace}
    db_reference = "{{{}}}db".format(namespace)

    sifts_object = OrderedDict()

    # Entities
    for entity_list in root.iterfind('.//ns:entity[@type="protein"]',
                                     namespaces=namespace_map):
        entity_id = entity_list.attrib['entityId']
        regions_full = OrderedDict()
        # Entities : Segments
        for segment in entity_list:

            # parse the regions found for this segment
            # Entities : Segments : Regions
            for region_list in segment.iterfind('.//ns:listMapRegion',
                                                namespaces=namespace_map):
                for region in region_list:
                    # get region annotations
                    # parse extra annotations for each region
                    for annotation in region:
                        for k, v in annotation.attrib.items():
                            # db entries
                            if annotation.tag == db_reference:
                                if k == 'dbSource' and v not in excluded:
                                    source = v
                                else:
                                    continue
                                if source not in regions_full:
                                    counter = 1
                                    region_object = OrderedDict()
                                else:
                                    region_object = regions_full[source]
                                    keys = regions_full[source]
                                    last_counter = [int(k) for k in keys][-1]
                                    counter = last_counter + 1
                                try:
                                    coord = annotation.attrib['dbCoordSys']
                                except KeyError:
                                    coord = '-'
                                region_annotation = OrderedDict()
                                region_annotation['dbAccessionId'] = annotation.attrib['dbAccessionId']
                                region_annotation['start'] = int(region.attrib['start'])
                                region_annotation['end'] = int(region.attrib['end'])
                                region_annotation['dbCoordSys'] = coord
                                region_object["{}".format(counter)] = region_annotation
                                regions_full[source] = region_object

        sifts_object[entity_id] = regions_full

    return sifts_object


def get_sifts_selected_from_table(data, chain=None, chain_auth=None, res=None,
                                  uniprot=None, site=None):
    """
    Utility that filters a pandas DataFrame by the input tuples.

    :param data: pandas DataFrame object
    :param chain: (tuple) chain IDs or None
    :param chain_auth: (tuple) chain IDs or None
    :param res: (tuple) res IDs or None
    :param uniprot: (tuple) UniProt IDs or None
    :param site: (tuple) UniProt (positional) sites or None
    :return: returns a modified pandas DataFrame
    """

    # excluding rows
    table = data
    if chain is not None:
        table = row_selector(table, 'PDB_entityId', chain)
        logger.info("SIFTS table filtered by PDB_entityId...")

    if chain_auth is not None:
        table = row_selector(table, 'PDB_dbChainId', chain_auth)
        logger.info("SIFTS table filtered by PDB_dbChainId...")

    if res is not None:
        table = row_selector(table, 'PDB_dbResNum', res)
        logger.info("SIFTS table filtered by PDB_dbResNum...")

    if uniprot is not None:
        table = row_selector(table, 'UniProt_dbAccessionId', uniprot)
        logger.info("SIFTS table filtered by UniProt_dbAccessionId...")

    if site is not None:
        table = row_selector(table, 'UniProt_dbResNum', site)
        logger.info("SIFTS table filtered by UniProt_dbResNum...")

    return table


class SIFTSreader(object):
    def __init__(self, inputfile):
        """
        :param inputfile: Needs to point to a valid SIFTS file.
        """
        self.inputfile = inputfile
        self.data = None
        self.excluded = ("InterPro", "GO", "EC", "NCBI")

        if not os.path.isfile(inputfile):
            raise IOError("{} not available or could not be read...".format(inputfile))

    def read(self, **kwargs):
        return self.residues(**kwargs)

    def residues(self, excluded=None, add_regions=True, add_dbs=False):
        if excluded is None:
            excluded = self.excluded
        self.data = parse_sifts_residues_from_file(self.inputfile, excluded=excluded,
                                                   add_regions=add_regions, add_dbs=add_dbs)
        return self.data

    def regions(self, excluded=None):
        if excluded is None:
            excluded = self.excluded
        self.data = parse_sifts_regions_from_file(self.inputfile, excluded=excluded)
        return self.data

    def dbs(self, excluded=None):
        if excluded is None:
            excluded = self.excluded
        self.data = parse_sifts_dbs_from_file(self.inputfile, excluded=excluded)
        return self.data

    def to_json(self, pretty=True):
        if self.data is not None:
            if type(self.data) is pd.core.frame.DataFrame:
                data = self.data.to_dict(orient='records')
            else:
                data = self.data
            if pretty:
                return json.dumps(data, sort_keys=False, indent=4)
            else:
                return json.dumps(data)
        else:
            logger.info("No SIFTS data parsed...")


if __name__ == '__main__':
    pass
