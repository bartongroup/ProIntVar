#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""

This defines the methods that load and validate user defined
parameters.

Fábio Madeira, 2015+

"""

import os
import sys
import click
import logging
import pkg_resources
from configparser import ConfigParser

log = logging.getLogger(__name__)
logging.captureWarnings(True)
logging.basicConfig(stream=sys.stderr, level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s ')

CONFIG_FILE = "config.ini"
CONFIG_FILE_TEMPLATE = "config_template.ini"


class Defaults(object):
    def __init__(self, config_file=None):
        if config_file is not None:
            if not os.path.isfile(config_file):
                raise IOError("{config_file} not available!")
        elif os.path.isfile(os.path.join(os.path.dirname(__file__), CONFIG_FILE)):
            config_file = os.path.join(os.path.dirname(__file__), CONFIG_FILE)
        elif os.path.isfile(os.path.join(os.path.dirname(__file__), CONFIG_FILE_TEMPLATE)):
            config_file = os.path.join(os.path.dirname(__file__), CONFIG_FILE_TEMPLATE)
        _config = ConfigParser()
        _config.read(config_file)
        self.__config = _config
        self.__populate_attributes()
        # log.info("Loaded values from {}".format(config_file))

    def __populate_attributes(self):
        for section in self.__config.sections():
            for var_name, var_par in self.__config.items(section):
                setattr(self, var_name, var_par)

    def __iter__(self):
        for section in self.__config.sections():
            for var_name, var_par in self.__config.items(section):
                yield var_name, var_par


config = Defaults()


@click.command()
@click.argument("filename")
def config_setup(filename):

    if os.path.isfile(filename):
        click.confirm("Config file already exist. Do you want to override it?",
                      abort=True)

    with open(filename, 'wb') as f:
        f.write(pkg_resources.resource_string(
            'prointvar', 'config_template.ini'))

    log.info("Wrote a template config file at {}".format(filename))


@click.command()
@click.confirmation_option()
@click.argument("filename")
def config_load(filename):

    if not os.path.isfile(filename):
        raise IOError("{filename} not available!")

    with open(filename, 'rb') as r:
        r = r.read()

    setup_dir_file = os.path.join(os.path.dirname(__file__), 'config.ini')
    with open(setup_dir_file, 'wb') as f:
        f.write(r)
    log.info("Config file loaded!")


if __name__ == "__main__":
    pass
