Shapeimager
===========

Shapeimager is a Python package for scientific FOV simulation from digital terrain models, NAIF SPICE kernels and any detector settings.

+------------------------+-------------------+--------+-----------------+------------+--------------+
| Continuous Integration | Code Coverage     | Docs   | Chat            |  Citation  |  Code Style  |
+========================+===================+========+=================+============+==============+
| |Travis Build Status|  | |Coverage Status| | |Docs| | |Join the chat| | |Citation| |  |Black|     |
| |Windows Build Status| |                   |        |                 | |JOSS|     |              |
+------------------------+-------------------+--------+-----------------+------------+--------------+

.. |Travis Build Status| image:: 
   :alt: Travis - Build Status
   :target: 
.. |Windows Build Status| image:: 
   :alt: Appveyor - Build Status
   :target: 
.. |Coverage Status| 
   :alt: Codecov - Test Coverage
   :target: 
.. |Docs| image:: 
   :alt: Readthedocs - Documentation
   :target: 
.. |Join the chat| image:: 
   :alt: Gitter - Chat room
   :target: 
.. |Citation| image:: 
   :alt: Citation Information: Zenodo
   :target: 
.. |JOSS| image:: 
   :alt: Citation Information: Journal of Open Source Software
   :target: 
.. |Black| image:: 
   :alt: Code Style - Black
   :target: 

Introduction
------------

Shapeimager is a Python Package dedicated to the image and FOV simulations of digital terrain models of bodies of the Solar System
for scientific purposes. The tool simulates surfaces at varied observational settings, mimicking the conditions in which the data
were originally obtained. 

The main goal is to provide the list of facets that fall inside a given pixel or spectral acquisition
and get the incidence, emergence, azimuth and phase angles of the facets in that spot.

Shapeimager is simulator but also works as a "mesh texturing" where every measurement is match with their corresponding facets.


Requirements
------------

Shapeimager 0.0.1 is built on Python 2.7 using the following packages:
- Numpy >1.15
- Pandas 0.24
- Scipy >1.2
- Cython 3.0
- Astropy 4.0
- SpiceyPy 2.0 written by `AndrewAnnex <https://github.com/AndrewAnnex/SpiceyPy>`__  (NASA NAIF/SPICE)
- Shapeimager comes with adapted and incorporated `PyPDS <https://github.com/RyanBalfanz/PyPDS>`__.

The module "Viewer" will require:
- Mayavi 4.6
- PyQt 4 or 5
- wxPython 3.0

Soon support for Python 3.


Publications with Shapeimager
----------------------------

Hasselmann et al. (2019). Pronounced morphological changes in a southern active zone on comet 67P/Churyumov-Gerasimenko. Astronomy & Astrophysics, Volume 630, id.A8, 19 pp.
Hasselmann et al. (2020). Modeling first-order scattering processes from OSIRIS-REx color images of the rough surface of asteroid (101955) Bennu. Submitted to Icarus Journal. Soon.


Installation
------------

No installation procedure for now. Just download the source codes and add shapeimager into your work folder.


Documentation
-------------

In preparation.


First Usage
-----------

First, ensure yourself to have all required SPICE Kernels, mk files and DTMs for your analysis.

In spec.py you must edit the path files:
``"obj"    : Object name or label.``

``"folder" : Directory of calibrated and aligned images.```

``"core"   : Directory of renderings and geo files.``

``"aux"    : Directory of auxiliary .dat or .txt files and also .obj Shape Model.``

``"kern"   : Directorty of NAIF/SPICE Kernels.``

``"prod"   : Directory where products and secondary data structures are stored.``

``"filter" : filter name.``


In your script call:
``from shapeimager import *``

Load SPICE kernels:

``spc = pos.from_spice(body=[BODY CODE],obs=[S/C CODE],ins=[INSTRUMENT CODE])``

``spc.furnish('[FILENAME].mk')``

Load the DTM or Shape Model:
``sha = "[DTM NAME].obj"``

``S = ShapeModel(sha, comments=[IF THERE ARE COMMENTS IN THE FILE])``

``S.normal_vector(True)`` 

It pre-loads the DTM and pre-calculate the normal vectors.

Chose a date:

``spc.load_time('YYYY-MM-DDThh:mm:ss.sss')``

``sun = spc.solar_coord(spc.body_frame)[0]``

``sc = spc.sc_coord(spc.body_frame)[0]``

Compute the Camera Matrix and boresight vector:

``FOV, CamMatrix, boresight = spc.orex_instrument_frame()``

Load the Imager Class to compute the FOV:

``Im = Imager(S, CamMatrix, boresight, sun, sc, visible=True, illuminated=True, raytrace=False, shaded=4, occ=4)``

flags:
visible     :: only visible facets
illuminated :: only illuminated facets
raytrace    :: higher precision but slower calculation of occlusion
shaded      :: >2, shadowing precision
occ         :: >2, occlusion precision with raytrace=False

Visualize mesh and check if the FOV is correct:

``Im.plot_v(FOV, ccd, 'test', 1, save=False)``

ccd :: 2-tuple with the CCD dimensions.

Run the Imaging function:

``Im.imaging(FOV, ccd)``

``XYZ = Im.onto_obj_frame()``

Make a 2d-array in the FOV using a scattering law:

  ``def ls_disk(x):``

    ``from numpy import cos``
  
    ``return 2e0*cos(x.inc_)/(cos(x.inc_)+cos(x.emi_))``

``property_image_ = Im.broadcast2(ls_disk(Im))``



