#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging

logging.getLogger("prointvar").addHandler(logging.NullHandler())
logging.captureWarnings(True)
logging.basicConfig(stream=sys.stderr, level=logging.INFO,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
