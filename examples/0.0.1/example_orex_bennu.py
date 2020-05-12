#!/usr/bin/python
# -*- coding: utf-8 -*-
# orex image simulator
# Author: Pedro H. A. Hasselmann
from __future__ import print_function, absolute_import, division

import numpy as np
from shapeimager import *

def ls_disk(x):
  from numpy import cos
  return 2e0*cos(x.inc_)/(cos(x.inc_)+cos(x.emi_))

###
### OSIRIS-REX GEO - BENNU
###
label='occ2_shaded2'
ccd = (1024, 1024)

# NAIF SPICE KERNELS
spc = pos.from_spice(body=2101955,obs=-64,ins=-64361) # OCAMS/MAPCAM/POLYCAM
spc.furnish('orex_2019-10-02_EQ.mk')

# SHAPE MODEL
sha = "g_01590mm_spc_obj_0000n00000_v032.obj"
S = ShapeModel(sha, comments=70)
S.normal_vector(True)

utctime = '2019-06-06T20:07:41.837'
spc.load_time(utctime)
sun = spc.solar_coord(spc.body_frame)[0]
sc = spc.sc_coord(spc.body_frame)[0]

print('sun',sun)
print('sc',sc)

# FRAMING ONTO FIEL-OF-VIEW OF THE INSTRUMENT
FOV, CamMatrix, boresight = spc.instrument_frame()

##
## IMAGER
##
Im = Imager(S, CamMatrix, boresight,
             sun, sc, 
             visible=True, 
             illuminated=True,
             raytrace=False,
             shaded=2, 
             occ=4)

#Im.plot_v(FOV, ccd, 'test', 1, False)
Im.imaging(FOV, ccd)
#Im.irradiance(Sw=1e0)
XYZ = Im.onto_target_frame()

# Image, 2d-array
property_image = Im.broadcast1(ls_disk(Im))
to_fits('test.fit', property_image)
to_fits('test_XYZ.fit', XYZ)

to_npz('test', utctime, Im)



# END
