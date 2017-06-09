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

### Project Structure

TODO.

### Licensing
[GNU GPL3](LICENSE.md)
