import os
import sys
# add root dir in the path so current is reachable
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
# import importlib
import imp
import utils

# importing base and current init as init
from init_base import Settings

init = None
settings = None
# init.__dict__.update(settings)
# del _

def meta_init(datadir: str):
    # make sure datadir has trailing slash
    datadir = utils.add_trailing_slash(datadir)
    global settings
    settings = Settings(datadir)
    global init
    # init = importlib.import_module(module_name, package=None)
    # print("datadir is", datadir+'init.py')
    init = imp.load_source('init', datadir+'init.py')
    # print("after loading init", dir(init))
