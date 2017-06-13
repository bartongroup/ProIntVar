#!/local/bin/python
# -*- coding: utf-8 -*-


import os
import sys
import json
import unittest
import pandas as pd
from contextlib import contextmanager

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
try:
    from mock import patch
except ImportError:
    from unittest.mock import patch

from prointvar.mmcif import (MMCIFreader, MMCIFwriter, parse_mmcif_atoms_from_file,
                             parse_mmcif_categories_from_file,
                             write_mmcif_from_table, get_mmcif_selected_from_table,
                             get_contact_indexes_from_table, add_mmcif_contacts,
                             add_mmcif_res_full, parse_pdb_atoms_from_file,
                             add_mmcif_res_split, get_atom_line, write_pdb_from_table,
                             add_mmcif_atom_altloc, residues_aggregation,
                             remove_multiple_altlocs)

from prointvar.config import config as c
root = os.path.abspath(os.path.dirname(__file__))
c.db_root = "{}/testdata/".format(root)


@contextmanager
def captured_output():
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err


@patch("prointvar.config.config.db_root", c.db_root)
class TestMMCIF(unittest.TestCase):
    """Test the DSSP parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.pdbid = '2pah'
        self.pdbid2 = '1ejg'
        self.inputcif = "{}{}{}.cif".format(c.db_root, c.db_cif, self.pdbid)
        self.inputbiocif = "{}{}{}.cif".format(c.db_root, c.db_cif_biounit, self.pdbid)
        self.outputcif = "{}{}{}.cif".format(c.db_root, c.tmp_dir_local, self.pdbid)
        self.inputpdb = "{}{}{}.pdb".format(c.db_root, c.db_cif, self.pdbid)
        self.inputpdb2 = "{}{}{}.pdb".format(c.db_root, c.db_cif, self.pdbid2)
        self.outputpdb = "{}{}{}.pdb".format(c.db_root, c.tmp_dir_local, self.pdbid)
        self.emptyfile = "{}{}{}.tmp".format(c.db_root, c.tmp_dir_local, self.pdbid)
        self.notfound = ""
        self.excluded = ()

        self.parser = parse_mmcif_atoms_from_file
        self.parser_pdb = parse_pdb_atoms_from_file
        self.reader = MMCIFreader
        self.writer = MMCIFwriter
        self.writer_method = write_mmcif_from_table
        self.filter = get_mmcif_selected_from_table
        self.contacts = get_contact_indexes_from_table
        self.parser_categories = parse_mmcif_categories_from_file
        self.add_contacts = add_mmcif_contacts
        self.add_res_full = add_mmcif_res_full
        self.add_mmcif_res_split = add_mmcif_res_split
        self.write_pdb_from_table = write_pdb_from_table
        self.get_atom_line = get_atom_line
        self.add_mmcif_atom_altloc = add_mmcif_atom_altloc
        self.residues_aggregation = residues_aggregation
        self.remove_altloc = remove_multiple_altlocs

    def tearDown(self):
        """Remove testing framework."""

        self.pdbid = None
        self.pdbid2 = None
        self.inputcif = None
        self.inputbiocif = None
        self.outputcif = None
        self.inputpdb = None
        self.inputpdb2 = None
        self.outputpdb = None
        self.emptyfile = None
        self.notfound = None
        self.excluded = None
        self.parser = None
        self.parser_pdb = None
        self.reader = None
        self.writer = None
        self.writer_method = None
        self.filter = None
        self.contacts = None
        self.parser_categories = None
        self.add_contacts = None
        self.add_res_full = None
        self.add_mmcif_res_split = None
        self.write_pdb_from_table = None
        self.get_atom_line = None
        self.add_mmcif_atom_altloc = None
        self.residues_aggregation = None
        self.remove_altloc = None

    def test_file_not_found_reader(self):
        with self.assertRaises(IOError):
            self.reader(self.notfound)

    def test_file_not_found_parser(self):
        with self.assertRaises(IOError):
            self.parser(self.notfound)

    def test_empty_file_reader(self):
        with self.assertRaises(ValueError):
            open(self.emptyfile, 'w').close()
            self.reader(self.emptyfile).read()
            os.remove(self.emptyfile)

    def test_reader_cif_verbose(self):
        with captured_output() as (out, err):
            self.reader(self.inputcif, verbose=True).read()
        # This can go inside or outside the `with` block
        output = out.getvalue().strip()
        self.assertEqual(output, 'Parsing mmCIF atoms from lines...')

    def test_parser_keys(self):
        self.assertListEqual([k for k in self.parser(self.inputcif).label_asym_id.unique()],
                             ['A', 'B', 'C', 'D'])

    def test_reader_keys(self):
        self.assertListEqual([k for k in self.reader(self.inputcif).read().label_asym_id.unique()],
                             ['A', 'B', 'C', 'D'])

    def test_reader_biocif_keys(self):
        self.assertListEqual([k for k in self.reader(self.inputbiocif).read().label_asym_id.unique()],
                             ['A', 'AA', 'B', 'BA', 'C', 'CA', 'D', 'DA'])

    def test_reader_cif_data(self):
        reader = self.reader(self.inputcif)
        data = reader.read()
        # ATOMS
        self.assertEqual(data.loc[0, 'label_asym_id'], 'A')
        self.assertEqual(data.loc[2686, 'label_asym_id'], 'B')
        self.assertEqual(data.loc[5315, 'label_asym_id'], 'C')
        self.assertEqual(data.loc[5316, 'label_asym_id'], 'D')
        self.assertEqual(data.loc[1, 'label_atom_id'], 'CA')
        self.assertEqual(data.loc[2687, 'label_atom_id'], 'CA')
        self.assertEqual(data.loc[5315, 'label_atom_id'], 'FE')
        self.assertEqual(data.loc[5316, 'label_atom_id'], 'FE')

    def test_reader_biocif_data(self):
        reader = self.reader(self.inputbiocif)
        data = reader.read()
        self.assertEqual(data.loc[0, 'label_asym_id'], 'A')
        self.assertEqual(data.loc[5372, 'label_asym_id'], 'B')
        self.assertEqual(data.loc[10630, 'label_asym_id'], 'C')
        self.assertEqual(data.loc[10632, 'label_asym_id'], 'D')
        self.assertEqual(data.loc[2686, 'label_asym_id'], 'AA')
        self.assertEqual(data.loc[8001, 'label_asym_id'], 'BA')
        self.assertEqual(data.loc[10631, 'label_asym_id'], 'CA')
        self.assertEqual(data.loc[10633, 'label_asym_id'], 'DA')

    def test_reader_atoms_cif_data(self):
        reader = self.reader(self.inputcif)
        data = reader.atoms()
        # ATOMS
        self.assertEqual(data.loc[0, 'label_asym_id'], 'A')
        self.assertEqual(data.loc[2686, 'label_asym_id'], 'B')
        self.assertEqual(data.loc[5315, 'label_asym_id'], 'C')
        self.assertEqual(data.loc[5316, 'label_asym_id'], 'D')
        self.assertEqual(data.loc[1, 'label_atom_id'], 'CA')
        self.assertEqual(data.loc[2687, 'label_atom_id'], 'CA')
        self.assertEqual(data.loc[5315, 'label_atom_id'], 'FE')
        self.assertEqual(data.loc[5316, 'label_atom_id'], 'FE')

    def test_parser_categories_cif_data(self):
        data = self.parser_categories(self.inputcif)
        # CATEGORIES
        self.assertEqual(data['_entry']['id'], [self.pdbid.upper()])
        self.assertEqual(data['_cell']['entry_id'], [self.pdbid.upper()])
        self.assertEqual(data['_symmetry']['entry_id'], [self.pdbid.upper()])

    def test_reader_categories_cif_data(self):
        reader = self.reader(self.inputcif)
        data = reader.categories()
        # CATEGORIES
        self.assertEqual(data['_entry']['id'], [self.pdbid.upper()])
        self.assertEqual(data['_cell']['entry_id'], [self.pdbid.upper()])
        self.assertEqual(data['_symmetry']['entry_id'], [self.pdbid.upper()])

    def test_reader_categories_selected_cif_data(self):
        reader = self.reader(self.inputcif)
        # CATEGORIES
        data = reader.categories(category='_entry')
        self.assertEqual(data['_entry']['id'], [self.pdbid.upper()])
        data = reader.categories(category='_cell')
        self.assertEqual(data['_cell']['entry_id'], [self.pdbid.upper()])
        data = reader.categories(category='_symmetry')
        self.assertEqual(data['_symmetry']['entry_id'], [self.pdbid.upper()])

    def test_reader_categories_excluded_cif_data(self):
        reader = self.reader(self.inputcif)
        # CATEGORIES
        data = reader.categories(excluded=('_atom_site',))
        keys = [k for k in data.keys()]
        self.assertNotIn('_atom_site', keys)

    def test_reader_cif_to_json_pretty(self):
        reader = self.reader(self.inputcif)
        reader.read()
        data = reader.to_json()
        self.assertEqual(json.loads(data)[0]['label_asym_id'], 'A')

    def test_reader_cif_to_json(self):
        reader = self.reader(self.inputbiocif)
        reader.read()
        data = reader.to_json(pretty=False)
        self.assertEqual(json.loads(data)[8001]['label_asym_id'], 'BA')

    def test_reader_default_excluded(self):
        reader = self.reader(self.inputcif)
        data = reader.atoms()
        keys = [k for k in data]
        self.assertNotIn("Cartn_x_esd", keys)
        self.assertNotIn("Cartn_y_esd", keys)
        self.assertNotIn("Cartn_z_esd", keys)
        self.assertNotIn("occupancy_esd", keys)
        self.assertNotIn("B_iso_or_equiv_esd", keys)

    def test_reader_new_excluded(self):
        reader = self.reader(self.inputcif)
        data = reader.atoms(excluded=self.excluded)
        keys = [k for k in data]
        self.assertIn("Cartn_x_esd", keys)
        self.assertIn("Cartn_y_esd", keys)
        self.assertIn("Cartn_z_esd", keys)
        self.assertIn("occupancy_esd", keys)
        self.assertIn("B_iso_or_equiv_esd", keys)

    def test_reader_default_category(self):
        reader = self.reader(self.inputcif)
        data = reader.atoms()
        keys = [k for k in data.label_asym_id.unique()]
        self.assertEqual(keys, ['A', 'B', 'C', 'D'])

    def test_reader_new_category(self):
        reader = self.reader(self.inputcif)
        data = reader.atoms()
        keys = [k for k in data.auth_asym_id.unique()]
        self.assertEqual(keys, ['A', 'B'])

    def test_filter_chain_id(self):
        reader = self.reader(self.inputcif)
        reader.atoms()
        data = self.filter(reader.data, chain=('A', ))
        self.assertIn("A", data.label_asym_id.unique())
        self.assertNotIn("B", data.label_asym_id.unique())

    def test_filter_atom_lines(self):
        reader = self.reader(self.inputcif)
        reader.atoms()
        data = self.filter(reader.data, lines=('ATOM', ))
        self.assertIn("ATOM", data.group_PDB.unique())
        self.assertNotIn("HETATM", data.group_PDB.unique())

    def test_writer_cif(self):
        if os.path.isfile(self.inputcif):
            self.writer(self.inputcif, self.outputcif).run()
            self.assertTrue(os.path.isfile(self.outputcif))
            os.remove(self.outputcif)
        else:
            raise IOError("%s" % self.inputcif)

    def test_writer_biocif(self):
        if os.path.isfile(self.inputbiocif):
            self.writer(self.inputbiocif, self.outputcif).run()
            self.assertTrue(os.path.isfile(self.outputcif))
            os.remove(self.outputcif)
        else:
            raise IOError("%s" % self.inputbiocif)

    def test_writer_method_cif(self):
        if os.path.isfile(self.inputcif):
            data = self.reader(self.inputcif).read()
            data = self.filter(data, chain=None, res=None, atom=None,
                               lines=None)
            self.writer_method(self.outputcif, data, override=True)
            self.assertTrue(os.path.isfile(self.outputcif))
            os.remove(self.outputcif)
        else:
            raise IOError("%s" % self.inputcif)

    def test_writer_cif_chain(self):
        if os.path.isfile(self.inputcif):
            self.writer(self.inputcif, self.outputcif).run(chain=('A',))
            self.assertTrue(os.path.isfile(self.outputcif))
            # reading the output
            reader = self.reader(self.outputcif)
            data = reader.read()
            self.assertNotIn('B', [k for k in data.keys()])
            self.assertNotIn('C', [k for k in data.keys()])
            os.remove(self.outputcif)
        else:
            raise IOError("%s" % self.inputcif)

    def test_writer_cif_res(self):
        if os.path.isfile(self.inputcif):
            self.writer(self.inputcif, self.outputcif).run(res=('1', '2', '3'))
            self.assertTrue(os.path.isfile(self.outputcif))
            # reading the output
            reader = self.reader(self.outputcif)
            data = reader.read()
            self.assertNotIn('4', [k for k in data.loc[:, 'label_seq_id'].unique()])
            self.assertNotIn('5', [k for k in data.loc[:, 'label_seq_id'].unique()])
            os.remove(self.outputcif)
        else:
            raise IOError("%s" % self.inputcif)

    def test_writer_cif_atom(self):
        if os.path.isfile(self.inputcif):
            self.writer(self.inputcif, self.outputcif).run(atom=('CA',))
            self.assertTrue(os.path.isfile(self.outputcif))
            # reading the output
            reader = self.reader(self.outputcif)
            data = reader.read()
            self.assertNotIn('N', [k for k in data.loc[:, 'label_atom_id'].unique()])
            self.assertNotIn('CB', [k for k in data.loc[:, 'label_atom_id'].unique()])
            os.remove(self.outputcif)
        else:
            raise IOError("%s" % self.inputcif)

    def test_writer_pdb_atom(self):
        if os.path.isfile(self.inputcif):
            self.writer(self.inputcif, self.outputpdb).run(atom=('CA',),
                                                           format_type="pdb")
            self.assertTrue(os.path.isfile(self.outputpdb))
            # reading the output
            reader = self.reader(self.outputpdb)
            data = reader.read(format_type="pdb")
            self.assertEqual(data.loc[0, 'label_asym_id'], 'A')
            self.assertEqual(data.loc[0, 'label_seq_id'], '118')
            self.assertEqual(data.loc[0, 'label_atom_id'], 'CA')
            self.assertEqual(data.loc[0, 'group_PDB'], 'ATOM')
            os.remove(self.outputpdb)
        else:
            raise IOError("%s" % self.inputcif)

    def test_writer_method_pdb(self):
        if os.path.isfile(self.inputcif):
            reader = self.reader(self.inputcif)
            data = self.filter(reader.read(), atom=('CA',))
            self.write_pdb_from_table(self.outputpdb, data, override=True)
            self.assertTrue(os.path.isfile(self.outputpdb))
            # reading the output
            reader = self.reader(self.outputpdb)
            data = reader.read(format_type="pdb")
            self.assertEqual(data.loc[0, 'label_asym_id'], 'A')
            self.assertEqual(data.loc[0, 'label_seq_id'], '118')
            self.assertEqual(data.loc[0, 'label_atom_id'], 'CA')
            self.assertEqual(data.loc[0, 'group_PDB'], 'ATOM')
            os.remove(self.outputpdb)
        else:
            raise IOError("%s" % self.inputcif)

    def test_get_atom_line(self):
        if os.path.isfile(self.inputcif):
            reader = self.reader(self.inputcif)
            data = reader.read()
            line = self.get_atom_line(data, index=0, atom_number=1,
                                      coords=None, new_chain=None)
            r = 'ATOM      1  N   VAL A 118      -7.069  21.943  18.770  1.00 56.51           N  \n'
            self.assertEqual(str(line), r)
            line = self.get_atom_line(data, index=1000, atom_number=1000,
                                      coords=None, new_chain=None)
            r = 'ATOM   1000  CD  ARG A 241       1.614   7.798  40.365  1.00 23.29           C  \n'
            self.assertEqual(str(line), r)
        else:
            raise IOError("%s" % self.inputcif)

    def test_writer_cif_lines(self):
        if os.path.isfile(self.inputcif):
            try:
                os.remove(self.outputcif)
            except FileNotFoundError:
                pass
            self.writer(self.inputcif, self.outputcif).run(lines=('HETATM', ))
            self.assertTrue(os.path.isfile(self.outputcif))
            # reading the output
            reader = self.reader(self.outputcif)
            data = reader.read()
            self.assertNotIn('A', [k for k in data.loc[:, 'label_asym_id'].unique()])
            self.assertNotIn('B', [k for k in data.loc[:, 'label_asym_id'].unique()])
            self.assertEqual(data.loc[0, 'group_PDB'], 'HETATM')
            self.assertEqual(data.loc[1, 'group_PDB'], 'HETATM')
            os.remove(self.outputcif)
        else:
            raise IOError("%s" % self.inputcif)

    def test_reader_add_contacts(self):
        reader = self.reader(self.inputcif)
        data = reader.read(add_contacts=True)
        keys = [k for k in data]
        self.assertIn('contact_indexes', keys)

    def test_reader_add_contacts_value(self):
        reader = self.reader(self.inputcif)
        data = reader.read(add_contacts=True)
        self.assertEqual(data.loc[0, 'contact_indexes'], '0,1,2,3,4,5,6,7,8,13')

    def test_add_contacts_value(self):
        reader = self.reader(self.inputcif)
        data = reader.read()
        data = self.add_contacts(data)
        self.assertEqual(data.loc[0, 'contact_indexes'], '0,1,2,3,4,5,6,7,8,13')

    def test_get_contact_indexes(self):
        reader = self.reader(self.inputcif)
        data = reader.read()
        indexes = self.contacts(data)
        self.assertEqual(indexes[0], [0, 1, 2, 3, 4, 5, 6, 7, 8, 13])

    def test_reader_add_res_full(self):
        reader = self.reader(self.inputcif)
        data = reader.read(add_res_full=True)
        self.assertIn('label_seq_id_full', data)
        self.assertIn('auth_seq_id_full', data)
        self.assertEqual(data.loc[0, 'auth_seq_id_full'], '118')

    def test_add_res_full(self):
        reader = self.reader(self.inputcif)
        data = reader.read()
        data = self.add_res_full(data)
        self.assertIn('label_seq_id_full', data)
        self.assertIn('auth_seq_id_full', data)
        self.assertEqual(data.loc[0, 'auth_seq_id_full'], '118')

    def test_parse_pdb(self):
        if os.path.isfile(self.inputpdb):
            data = self.parser_pdb(self.inputpdb)
            self.assertListEqual([k for k in data.label_asym_id.unique()], ['A', 'B'])
            self.assertIn('label_seq_id', [k for k in data])
            self.assertIn('label_asym_id', [k for k in data])
            self.assertIn('label_atom_id', [k for k in data])
            self.assertIn('Cartn_x', [k for k in data])
            self.assertIn('Cartn_y', [k for k in data])
            self.assertIn('Cartn_z', [k for k in data])
        else:
            raise IOError("%s" % self.inputpdb)

    def test_read_pdb(self):
        if os.path.isfile(self.inputpdb):
            reader = self.reader(self.inputpdb)
            data = reader.read(format_type="pdb")
            self.assertListEqual([k for k in data.label_asym_id.unique()], ['A', 'B'])
            self.assertIn('label_seq_id', [k for k in data])
            self.assertIn('label_asym_id', [k for k in data])
            self.assertIn('label_atom_id', [k for k in data])
            self.assertIn('Cartn_x', [k for k in data])
            self.assertIn('Cartn_y', [k for k in data])
            self.assertIn('Cartn_z', [k for k in data])
        else:
            raise IOError("%s" % self.inputpdb)

    def test_get_res_split(self):
        if os.path.isfile(self.inputpdb):
            data = self.parser_pdb(self.inputpdb)
            data = self.add_mmcif_res_split(data)
            self.assertIn('label_seq_id_full', [k for k in data])
            self.assertIn('label_seq_id', [k for k in data])
            self.assertIn('auth_seq_id', [k for k in data])
            self.assertIn('pdbx_PDB_ins_code', [k for k in data])
            self.assertEqual('118', data.loc[0, 'label_seq_id_full'])
            self.assertEqual('118', data.loc[0, 'label_seq_id'])
            self.assertEqual('118', data.loc[0, 'auth_seq_id'])
            self.assertEqual('?', data.loc[0, 'pdbx_PDB_ins_code'])
            self.assertEqual('.', data.loc[0, 'label_alt_id'])
        else:
            raise IOError("%s" % self.inputpdb)

    def test_get_res_split_various(self):
        data = pd.DataFrame([{'label_seq_id_full': '1'}])
        data = self.add_mmcif_res_split(data)
        self.assertEqual('1', data.loc[0, 'label_seq_id'])
        self.assertEqual('?', data.loc[0, 'pdbx_PDB_ins_code'])
        data = pd.DataFrame([{'label_seq_id_full': '-10'}])
        data = self.add_mmcif_res_split(data)
        self.assertEqual('-10', data.loc[0, 'label_seq_id'])
        self.assertEqual('?', data.loc[0, 'pdbx_PDB_ins_code'])
        data = pd.DataFrame([{'label_seq_id_full': '-100A'}])
        data = self.add_mmcif_res_split(data)
        self.assertEqual('-100', data.loc[0, 'label_seq_id'])
        self.assertEqual('A', data.loc[0, 'pdbx_PDB_ins_code'])

    def test_parse_pdb_altloc(self):
        if os.path.isfile(self.inputpdb2):
            data = self.parser_pdb(self.inputpdb2, remove_altloc=False)
            data = self.add_mmcif_atom_altloc(data)
            self.assertIn('label_atom_altloc_id', [k for k in data])
            self.assertIn('auth_atom_altloc_id', [k for k in data])
            self.assertEqual('N.A', data.loc[0, 'label_atom_altloc_id'])
            self.assertEqual('N.B', data.loc[1, 'label_atom_altloc_id'])
            self.assertEqual('CA.A', data.loc[2, 'label_atom_altloc_id'])
            self.assertEqual('CA.B', data.loc[3, 'label_atom_altloc_id'])
            self.assertEqual('C', data.loc[4, 'label_atom_altloc_id'])
            self.assertEqual('O', data.loc[5, 'label_atom_altloc_id'])
        else:
            raise IOError("%s" % self.inputpdb2)

    def test_residues_aggregation_first(self):
        reader = self.reader(self.inputcif)
        data = reader.atoms()
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('VAL', data.loc[0, 'label_comp_id'])
        self.assertEqual('118', data.loc[0, 'auth_seq_id'])
        self.assertEqual(-7.069, data.loc[0, 'Cartn_x'])
        data = self.residues_aggregation(data, agg_method='first', category='label')
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('VAL', data.loc[0, 'label_comp_id'])
        self.assertEqual('118', data.loc[0, 'auth_seq_id'])
        self.assertEqual(-7.069, data.loc[0, 'Cartn_x'])

    def test_residues_aggregation_centroid(self):
        reader = self.reader(self.inputcif)
        data = reader.atoms()
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('VAL', data.loc[0, 'label_comp_id'])
        self.assertEqual('118', data.loc[0, 'auth_seq_id'])
        self.assertEqual(-7.069, data.loc[0, 'Cartn_x'])
        data = self.residues_aggregation(data, agg_method='centroid', category='label')
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('VAL', data.loc[0, 'label_comp_id'])
        self.assertEqual('118', data.loc[0, 'auth_seq_id'])
        self.assertAlmostEqual(-7.310, data.loc[0, 'Cartn_x'], places=2)

    def test_reader_remove_altloc(self):
        reader = self.reader(self.inputpdb2)
        data = reader.atoms(format_type='pdb', remove_altloc=True)
        self.assertEqual(1, data.loc[0, 'id'])
        self.assertEqual('N', data.loc[0, 'type_symbol'])
        self.assertEqual('N', data.loc[0, 'label_atom_id'])
        self.assertEqual('.', data.loc[0, 'label_alt_id'])

    def test_remove_altloc(self):
        reader = self.reader(self.inputpdb2)
        data = reader.atoms(format_type='pdb', remove_altloc=False, reset_atom_id=False)
        ndata = self.remove_altloc(data)
        self.assertEqual(1, ndata.loc[0, 'id'])
        self.assertEqual('N', ndata.loc[0, 'type_symbol'])
        self.assertEqual('N', ndata.loc[0, 'label_atom_id'])
        self.assertEqual('.', ndata.loc[0, 'label_alt_id'])


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMMCIF)
    unittest.TextTestRunner(verbosity=2).run(suite)
