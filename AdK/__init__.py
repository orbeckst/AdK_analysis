# $Id: __init__.py 4096 2010-04-26 15:51:19Z root $
"""
AdK python scripts
==================

The AdK module collects scripts that are very specific for the AdK
project, such as the Modeller script or the analysis fo the angles.


Installation
------------

In order to use this module you will have to modify your
PYTHONPATH. The simplest thing to do is to setup your shell by running

  source bin/setup.bash

NOTE: You MUST run the script from the base directory, i.e. the one
containing bin, lib, atom_files etc.
"""

# config.py is auto generated
__all__ = []

try:
    import config
except ImportError:
    raise ImportError("Run 'bin/install.sh' first to setup AdK.config.")

import warnings
try:
    from Model import ModelAdK
except ImportError:
    warnings.warn("WARNING (ImportError): Modeller is not available; ModelAdK() will not work.")



