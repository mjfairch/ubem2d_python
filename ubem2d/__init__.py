# Package imports
from ubem2d.aerodynamics import *
from ubem2d.fluids import *
from ubem2d.geometry import *
from ubem2d.math import *
from ubem2d.motion import *
from ubem2d.panel import *
from ubem2d.solvers import *
from ubem2d.unsteady import *
from ubem2d.util import *

# Module imports
from ubem2d.Errors import *

# Attributes
import os
__version__ = '1.0'
__base_dir = os.path.dirname(os.path.dirname(__file__))
__plot_dir = os.path.join(__base_dir, 'plots')
__data_dir = os.path.join(__base_dir, 'data')
__movie_dir = os.path.join(__base_dir, 'movies')
