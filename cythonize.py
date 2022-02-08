# -*- coding: utf-8 -*-
#!/usr/bin/python

'''
  Setup.py to externally build a cython code.
  
  python -O cythonize.py build_ext --inplace
'''
from setuptools import setup, find_packages, Extension
#from distutils.core import setup
#from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from numpy import get_include
from os import path, name, remove, getcwd
from pathlib import Path
print(name)

pyxfiles = [path.join(Path(__file__).resolve().parent,'geo', 'imager_utils.pyx'),
            path.join(Path(__file__).resolve().parent,'geo', 'raytracing.pyx'),
            path.join(Path(__file__).resolve().parent,'geo', 'lineclipping.pyx'),
            ]

module_names = ["imager_utils", "raytracing", "lineclipping"]

for n, pyx in enumerate(zip(pyxfiles, module_names)):
  pyxfile, module_name = pyx[0], pyx[1]
  print(pyxfile)
  #print(tree)
  #raw_input([path.join(getcwd(),'source','model',md) for md in module_names[:n]])

  # Linux
  if name != 'nt':
    tree = pyxfile.split('/')

    ext_modules = [Extension('.'.join((tree[-3],tree[-2],module_name)), 
                           [pyxfile],
                           extra_link_args=["-fopenmp"],
                           include_dirs=[get_include()],
                           extra_compile_args=["-O3", "-fno-math-errno", "-funsafe-math-optimizations", "-fcx-limited-range", "-fexcess-precision=fast", "-fopenmp","-w"])] #"-ffast-math"

  # Windows
  if name == 'nt':
    tree = pyxfile.split('\\')

    ext_modules = [Extension('.'.join((tree[-3],tree[-2],module_name)), 
                           [pyxfile],
                           extra_link_args=["/openmp"],
                           include_dirs=[get_include()],
                           extra_compile_args=["/O2", "/openmp","/W3"])]

  ext_modules[0].cython_directives = {"embedsignature": True}

  setup(
  name = module_name,
  #cmdclass = {'build_ext': build_ext},
  include_dirs = [get_include()],
  ext_modules = cythonize(ext_modules, build_dir="build"),
  script_args=['build_ext'],
  options={'build_ext':{'inplace':True}},
  include_package_data=True,
  #package_data={'source.model': [path.join(getcwd(),'source','model','*.pxd'),]}
  )

# END
