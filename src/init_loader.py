import os
import sys
# add root dir in the path so current is reachable
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# importing base and current init as init
import init_base as init
import current.init as _
init.__dict__.update(_.__dict__)
del _