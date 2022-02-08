# -*- coding: utf-8 -*-
#!/usr/bin/python

__version__="0.2"

__doc__='''
####################################################################################################
#  author: Pedro H. A. Hasselmann, Rio de Janeiro, Brasil (LESIA-OBSPM/ON-MCTI)
#
#  Shapeimager
#  ===========
#
#  Shapeimager is a FOV rendering tool for scientific spectral and BRDF analysis for Solar System Bodies.
#
#  Python 3.8
#
#  requirements:
#      essential: cython 0.29, spiceypy 4.0, pandas 1.2.4, numpy 1.20, scipy 1.6, astropy 4.2, matplotlib 3.3, PIL
#      support:   mayavi, gdal, osgeo, plyfile
#
#
#  Path Files :
#
# "folder" : Directory of the calibrated and aligned images.
# "core"   : Directory of the renderings and HDF5 databases.
# "aux"    : Directory of auxiliary .dat or .txt files containing important values and 3D Shape Model.
# "kern"   : Directory of NAIF/SPICE Kernels
#
#
######################################################################################################
'''

# Garbage Collector disabled.
#import gc
#gc.disable()
#gc.set_debug(gc.DEBUG_LEAK)

# System
from os import path, walk, listdir, mkdir, remove, rename, system, name, stat, environ, sep, getcwd
import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter(action="ignore")

# Load File path
home = path.expanduser('~')
global folder, core, aux, kern
folder = getcwd()
core   = getcwd()
aux    = getcwd() 
kern   = getcwd()

# cythonization
from . import cythonize


####################
## Physical Units ##
####################
#import astropy.units as u
au_km = 149597870.659999996424 #u.au.to(u.km)



#END
