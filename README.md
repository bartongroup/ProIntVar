### ProIntVar Core
ProIntVar Core is a Python module that implements methods for working with Protein Structures (handles mmCIF, DSSP, SIFTS, protein interactions, etc.) and genetic variation (via UniProt and Ensembl APIs). 

ProIntVar core is now separated from ProtIntVar Analysis, which contains analysis scripts that use ProIntVar Core components. 

### Dependencies
Using Python 3.5+.

Check [requirements.txt](./requirements.txt) for all dependencies.

### Installing 

Setting up a virtual environment 
```sh
$ virtualenv --python `which python` env
$ source env/bin/activate
```

Installing ProIntVar

```sh
$ wget https://github.com/bartongroup/ProIntVar-Core/archive/master.zip -O ProIntVar-Core.zip
$ unzip ProIntVar-Core.zip
  
# alternatively
$ git clone https://github.com/bartongroup/ProIntVar-Core.git
  
# installing requirements
$ cd ProIntVar-Core
$ pip install -r requirements.txt
  
# then...  
$ python setup.py test
$ python setup.py install
```

### Configuration

Editing the provided template configuration settings
```sh
$ cd /path/to/desired/working/dir/
  
# Get a copy of the template config.ini file shipped with ProIntVar
$ ProIntVar-Core-config-setup new_config.ini
  
# Update the settings according to user preferences and push them
$ ProIntVar-Core-config-load new_config.ini
```

Testing that the new values are correctly loaded by ProIntVar
```sh
$ python
>>> from prointvar.config import config
>>> config.tmp_dir
'/new_config_path_to_tmp_dir/'
```


### How to use

TODO.

### Additional Information

#### Project Structure

TODO

#### Guidelines on file names and extensions
**PDB/PDBx/mmCIF Macromolecular structures**
* PDB and mmCIF formatted files are read and written from `db_pdbx` folder, as defined in the configuration file `config.ini`
    - PDB/mmCIF files are written as `<pdb_id>.pdb` or `<pdb_id>.cif`
    - BioUnits from PDBe are written as `<pdb_id>_bio.cif` 
    - New structure files written for running DSSP, Reduce, HBPLUS or Arpeggio are generally written as `<4char>_new.pdb` format
    - By-chain/entity structures are written as `<pdb_id>_<chain_id>.pdb`

**DSSP Secondary Structure**
* DSSP files are read and written from `db_dssp` folder  
    - DSSP files are generally written as `<pdb_id>.dssp`
    - By-chain/entity DSSP outputs are written as `<pdb_id>_<chain_id>.dssp`
    - Unbound-state DSSP are written as `<pdb_id>_unbound.dssp`

**SIFTS Structure-Sequence (PDB-UniProt) cross-reference**
* SIFTS files are read and written from `db_sifts` folder
    - SIFTS files are written as `<pdb_id>.xml`

**Arpeggio Interface Contacts**
* Arpeggio files are read and written from `db_contacts` folder
    - Arpeggio files are written as `<pdb_id>.contacts`, `<pdb_id>.amam`, `<pdb_id>.amri`, `<pdb_id>.ari` and `<pdb_id>.ri`

**HBPLUS Hydrogen-Bond Contacts**
* HBPLUS files are read and written from `db_contacts` folder
    - HBPLUS files are written as `<pdb_id>.h2b`
    - HBPLUS Hydrogen-filled PDBs are written as `<pdb_id>.h.pdb` in `db_pdbx`

**Reduce PDBs filled with Hydrogen**
* Reduce files are read and written from `db_pdbx` folder
    - Reduce Hydrogen-filled PDBs are written as `<pdb_id>.h.pdb` in `db_pdbx`


#### Key features

* PDBx/mmCIF support in both reading and writing


### Licensing
[GNU GPL3](LICENSE.md)
