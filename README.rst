Shapeimager
===========

Shapeimager is a Python 3 package for scientific FOV simulation from any digital terrain models, NAIF SPICE kernels, or detector settings.

+-------------------+--------+-----------------+------------+--------------+
| Code Coverage     | Docs   | Chat            |  Citation  |  Code Style  |
+===================+========+=================+============+==============+
|   Coverage Status |   Docs |   Join the chat | Citation   |              |
+-------------------+--------+-----------------+------------+--------------+
IN PROGRESS.

Introduction
------------

Shapeimager is a Python Package dedicated to the imaging and FOV simulations of digital terrain models of bodies of the Solar System
for scientific purposes. The tool simulates surfaces at varied observational settings, mimicking the conditions in which the data
were originally obtained. 

The main goal is to provide a detailed and precise list of facets that fall inside a given pixel or spectral acquisition, 
and then, to be able to recover the facet's incidence, emergence, azimuth and phase angles from inside that spot. Shapeimager is bridging an extremely crucial procedure for Bi-directional reflectance distribution analysis of different terrains ans soils. When a disk-resolved data are available for the surface of an asteroid or comet, it can be important to study all pixels or spectral acquisitions available, together or separately, to understand how it scatters light. Such information leads to constrains on grain size distribution and mineralogical composition.

Shapeimager is mainly an image simulator, but also works as a "mesh texturing" where every measurement is match with their corresponding facets,
thus can be used for reproducing RGB images of the body in varied perpectives.


Benchmarks
----------

Benchmarks with an image of 1024 x 1024 pixels using NASA/OSIRIS-REx SPICE Kernels 
and a shape model of asteroid Bennu of 786432 facets perform in less than 20 seconds in general.
Time variations can arise depending on the asteroid distance, shadowing and shapemodel resolution.



Requirements
------------

Shapeimager 1.0.0 is built on Python 3.8 using the following packages:

  - Numpy >1.15
  - Pandas >1.0
  - Scipy >1.2
  - Cython >0.29
  - Astropy >4.0
  - SpiceyPy >2.0 written by `AndrewAnnex <https://github.com/AndrewAnnex/SpiceyPy>`__  (NASA NAIF/SPICE)
  - Shapeimager comes with adapted and incorporated `PyPDS <https://github.com/RyanBalfanz/PyPDS>`__ for Python 3.8.

The module "Viewer" is optional and will require:

  - Mayavi 4.6
  - PyQt 4 or 5
  - wxPython 3.0


Publications with Shapeimager
----------------------------

 - Hasselmann et al. (2019). Pronounced morphological changes in a southern active zone on comet 67P/Churyumov-Gerasimenko. Astronomy & Astrophysics, Volume 630, id.A8, 19 pp.
 - Hasselmann et al. (2020). Modeling optical roughness and first-order scattering processes from OSIRIS-REx color images of the rough surface of asteroid (101955) Bennu. Icarus 357, 114106. 


Installation
------------
SOON.


Documentation
-------------

In preparation.


First Usage
-----------

First, ensure yourself to have available all required SPICE Kernels, mk files and Digital Terrain Models for your analysis.



In your script call:

``import shapeimager``

To change the global pathnames, in case your data are not located at same folder as your script:
::

  shapeimager.folder = [Directory of calibrated and aligned images.]

  shapeimager.core   = [Directory of renderings and geo files.]

  shapeimager.aux    = [Directory of auxiliary .dat or .txt files and also .obj Shape Model.]

  shapeimager.kern   = [Directory of NAIF/SPICE Kernels.]

Whether the pathnames are changed or not, call again:

``from shapeimager.main import *``



Then, load the SPICE kernels:

::

  spc = position.from_spice(body=[BODY CODE],obs=[S/C CODE],ins=[INSTRUMENT CODE])

  spc.furnish('[FILENAME].mk')

Also, load your DTM or Shape Model:

::

  S = ShapeModel('[DTM NAME].obj', comments=[COMMENTED LINES])

  S.normal_vector(True)

It pre-loads the DTM and pre-calculates the normal vectors for future speed up.

Chose a date:

::

  spc.load_time('YYYY-MM-DDThh:mm:ss.sss')

  sun = spc.solar_coord(spc.body_frame)[0]

  sc = spc.sc_coord(spc.body_frame)[0]

Compute the Camera Matrix and boresight vector:

::

  FOV, CamMatrix, boresight = spc.instrument_frame()

Load the Imager Class to compute the FOV:

::

  Im = Imager(S, CamMatrix, boresight, sun, sc, visible=True, illuminated=True, raytrace=False, shaded=4, occ=4)

============== ========================================================
  flags                       description                            
============== ========================================================
 visible          only visible facets                                   
 illuminated      only illuminated facets                               
 raytrace         higher precision but slower calculation of occlusions 
 shaded           >2, shadowing precision                               
 occ              >2, occlusion precision with raytrace=False           
============== ========================================================

Visualize mesh and check if the FOV is correct:

::

  Im.plot_v(FOV, ccd, 'test', 1, save=False)

ccd :: 2-tuple with the CCD dimensions.

Run the Imaging function:

::

  Im.imaging(FOV, ccd)

What is calculated by Im.imaging?

====================== ========================================================
  properties                       description                            
====================== ========================================================
 d                       S/C Distance to target                        
 inc                     Incidence angle                   
 emi                     Emergence angle  
 pha                     Phase angle      
 facetid                 Active facet index
 facet_area_pix          Portion of facet under a pixel/acquisition
 solid_angle_inc         Incoming solid angle
 solid_angle_emi         Oucoming solid angle
 facet_pix               Link among facets and image pixel
 facet_image             Image with the central facet index
====================== ========================================================

Get the Cartesian coordinates as image cube, for geo-referencing:

::

  XYZ = Im.onto_target_frame()

Make a FOV image applying a scattering law to compute surface brightness:

::

  def ls_disk(x):

    from numpy import cos
  
    return 2e0*cos(x.inc_)/(cos(x.inc_)+cos(x.emi_))



::
  
  property_image = Im.broadcast1(ls_disk(Im))  # Less accurate but faster
  
or

::

  property_image = Im.broadcast2(ls_disk(Im), plot=False) # Accurate but slightly less faster
  
Images can be saved into FITS format using:
  
::
  
  to_fits('test.fit', property_image)
    
And Imager properties can be saved into npz format:

::

  to_npz('[LABEL]', 'YYYY-MM-DDThh:mm:ss.sss', Im)

If the flux-calibrated image or acquisition is available, the pixel-facet matching can be performed and stored:

::

  import numpy as np
  
  image = read_image([IMAGE PATH], channel=0)
  
  HDF = Store(columns=['dist', 'pha', 'emi', 's', 'inc'], label=[LABEL]) 
  
  geo_data = np.load([NPZ FILEPATH], mmap_mode='r')
  
  HDF.image_dataframe([VALUE NAME], image.data, geo_data, offset=(0,0), threshold=1e-5) 
  
  HDF.storing_image()
  
  HDF.close()
 
 
 
