#cython: language_level=3, boundscheck=False, wraparound=False, nonecheck=False, cdivision=True
#distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
#distutils: language = c++
# -*- coding: utf-8 -*-

__version__="0.2"

__doc__='''
  ############## CYTHON RASTERIZATION & UTILS #############################
  ### Author: Pedro Henrique A. Hasselmann, Dr. (OBSPM-Lesia)
  ### Last Modified: April, 2021
  ###
  ### requirements: 
  ### Cython, cpp compiler
  ###
  ### Main Functions:
  ###  rebase_c
  ###   
  ###
  ###
  ###
  ###
  ###
  ################################################################################
'''
# Cythonize
import cython
cimport cython
from cython.parallel import prange, parallel
from cython cimport view as cview
from cpython cimport array
from libc.math cimport ceil
from libc.string cimport memcpy, memset
#from libcpp.vector cimport vector

import numpy as np
cimport numpy as np
from numpy import ndarray, uint8, uint16, uint32, uint64, int32, float32, float64, uint8
from numpy cimport ndarray, uint8_t, uint16_t, uint32_t, int32_t, int64_t, float32_t, float64_t
from numpy.math cimport INFINITY
np.import_array()

# NUMBER OF THREADS
cdef short thrn=4

################
### FUNCTIONS ##
################

cpdef rebase_c(double[:,::1] v, # Vertices N columns
               unsigned int[:,::1] fc,   # Facet indexes
               ):

    cdef:
       ssize_t i, j
       ssize_t F = fc.shape[0]
       double[:,:,::1] R = np.empty((F,3,v.shape[1]), order='C', dtype=float64)
       double[:,::1] t = np.empty((3,v.shape[1]), order='C', dtype=float64) 

    with nogil:
      for i in range(F):
            j =  fc[i][0]          
            R[i][0,...] = v[j,...]

            j = fc[i][1]
            R[i][1,...] = v[j,...]

            j = fc[i][2]
            R[i][2,...] = v[j,...]

    return np.asarray(R)


cpdef  weighted_average_c(float[:,::1] fpa,  # Y X Facet Area
                          float[::1] X, # Values to be broadcast
                          ):

    cdef:
       ssize_t i
       unsigned int f=0
       double value=0., a=0.
       ssize_t N = fpa.shape[0]
       float[:,::1] X_pix = np.empty((N,3), dtype=float32, order='C')

    X_pix[:,:] = -999.

    #ith nogil:
    a = fpa[0,3]
    #f = int(fpa[0,2])
    value = a*X[0]
    for i in range(1,N):
        if (fpa[i,0] == fpa[i-1,0]) and (fpa[i,1] == fpa[i-1,1]):
              #f = int(fpa[i,2])
              value += fpa[i,3]*X[i]
              a += fpa[i,3]

        else:
            #if a>0.99:
            #       value = value/a
            
            X_pix[i,2] = value/a
            X_pix[i,1] = fpa[i-1,1]
            X_pix[i,0] = fpa[i-1,0]
            #print(i, N, np.asarray(X_pix[i,:]))

            a = fpa[i,3]
            #f = int(fpa[i,2])
            value = a*X[i]

         

    return np.unique(np.asarray(X_pix), axis=0)
          

####################################################################
#################### RASTERIZATION #################################
####################################################################

cdef float[:,::1] vc = np.empty((2,3), dtype=float32, order='C')

cpdef raster_c((ssize_t, ssize_t) it,
                float[:,:,::1] V,
                unsigned int[::1] f,
                (short, short) frame, 
                unicode mode=u'pixelation'):
    '''
      Rasterization.

      Fill-in the area inside vertices (skimage.draw.polygon) -- DONE
      https://stackoverflow.com/questions/26807177/overlapping-polygons-in-python-pil
      https://stackoverflow.com/questions/42349587/cython-slow-numpy-arrays
      https://www.scratchapixel.com/lessons/3d-basic-rendering/rasterization-practical-implementation

      Parameters
      ==========
      it       : start & stop iteration (mutiprocessing.pool.map input)
      V        : 2x3xN numpy.array, vertices triplet per facet
      facetid  : N numpy.array, original facet id
      frame    : frame size in pixels
      mode     : 'pixelation' or 'occultation'
      
      Output
      ======
      image : 2D int32 numpy.array
      facet_pix : 3xN numpy.array, X and Y pixel-coordinates alongside with facetid.      
    '''
    #from skimage.draw import polygon, polygon_perimeter

    cdef:
       ssize_t j
    
    pa = TriangleAccumulator(frame[1], frame[0], mem_fac=4)
    
    #with nogil:
    if mode == 'pixelation':
          for j in range(it[0],it[1]):
            # Tringular facet area projection (polygon filling)
            vc[...] = V[...,j]
            pa.superimpose_polygon(vc[0,:],vc[1,:], f[j])

            #print('v', j)
            #print('f', f[j])
            #print('n', pa.n)
            #r = np.asarray(v[0,:], dtype=int32)
            #c = np.asarray(v[1,:], dtype=int32)
            #print(rr_area)
            #print(cc_area)
            #print(r)
            #print(c)
            #image[rr_area, cc_area]=n+1
            #image[rr_peri, cc_peri]=n+1
            #image[r,c]=n+1
            #raw_input(image[np.min(r)-3:np.max(r)+3,np.min(c)-3:np.max(c)+3])

    if mode == 'occultation':
          for j in range(it[0],it[1]):
            # Tringular facet area projection (polygon filling)
            vc[...] = V[...,j]
            pa.add_polygon(vc[0,:],vc[1,:], 1)
    
    return np.asarray(pa.image), np.stack((pa.rr, pa.cc, pa.ff), axis=-1).astype(uint32)



##############################
## DRAW.POLYGON ACCUMULATOR ##
##############################
# https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
cdef bint point_in_triangle(float[::1] xp, 
                            float[::1] yp,
                            float x, 
                            float y) nogil:
     
    cdef:
      float dX = x - xp[0]
      float dY = y - yp[0]
      float dX20 = xp[2] - xp[0]
      float dY20 = yp[2] - yp[0]
      float dX10 = xp[1] - xp[0]
      float dY10 = yp[1] - yp[0]

      float s = (dY20*dX) - (dX20*dY)
      float t = (dX10*dY) - (dY10*dX)
      float D = (dX10*dY20) - (dY10*dX20)
 
    if D<0:
      return (s<=0) & (t<=0) & ((s + t)>=D)
    else:
      return (s>=0) & (t>=0) & ((s + t)<=D)


cdef (float, float) minmax(float[::1] arr) nogil:
    cdef float min = INFINITY
    cdef float max = -INFINITY
    cdef ssize_t i, L = arr.shape[0]
    for i in range(L):
        if arr[i] < min:
            min = arr[i]
        if arr[i] > max:
            max = arr[i]
    return min, max

#https://stackoverflow.com/questions/26807177/overlapping-polygons-in-python-pil
cdef class TriangleAccumulator:

    cdef:
       short width, height
       int n
       unsigned int[:, ::1] image
       unsigned short[::1] rr, cc
       unsigned int[::1] ff

    def __cinit__(self, short width, 
                        short height, 
                        short mem_fac):

        self.n = 0
        self.width = width
        self.height = height
        shape = (height, width)
        
        self.image = np.zeros(shape, dtype=uint32)
        self.rr = np.empty(mem_fac*(self.width-1)*(self.height-1), dtype=uint16, order='C')
        self.cc = np.empty(mem_fac*(self.width-1)*(self.height-1), dtype=uint16, order='C')
        self.ff = np.empty(mem_fac*(self.width-1)*(self.height-1), dtype=uint32, order='C')

    def reset(self):
       self.image[:, :] = 0

    # Rasterization of triangles in self.image are **summed** over the previous values
    cdef void add_polygon(self,
                          float[::1] ya, 
                          float[::1] xa, 
                          unsigned int value) nogil:

        cdef float minya, maxya, minxa, maxxa
        minya, maxya = minmax(ya)
        minxa, maxxa = minmax(xa)

        cdef:
           ssize_t minr = int(max(0, minya))
           ssize_t maxr = int(ceil(maxya))
           ssize_t minc = int(max(0, minxa))
           ssize_t maxc = int(ceil(maxxa))
           unsigned short r, c

        maxr = min(self.height -1, maxr)
        maxc = min(self.width  -1, maxc)

        #with nogil:
        for r in range(minr, maxr+1):
            for c in range(minc, maxc+1):
                if (point_in_triangle(xa, ya, c, r+0.5) or 
                   point_in_triangle(xa, ya, c+0.5, r) or 
                   point_in_triangle(xa, ya, c-0.5, r) or 
                   point_in_triangle(xa, ya, c, r-0.5) or 
                   point_in_triangle(xa, ya, c+0.5, r+0.5) or 
                   point_in_triangle(xa, ya, c-0.5, r-0.5) or 
                   point_in_triangle(xa, ya, c-0.5, r+0.5) or
                   point_in_triangle(xa, ya, c+0.5, r-0.5) or  
                   point_in_triangle(xa, ya, c, r)
                   ):

                    self.image[r, c] += value
                    self.rr[self.n] = r
                    self.cc[self.n] = c
                    self.ff[self.n] = value
                    self.n +=1

        #return rr, cc



    # Rasterization of triangles in self.image are **superimposed** over the previous values
    #@cython.boundscheck(False)  # Deactivate bounds checking
    #@cython.wraparound(False)   # Deactivate negative indexing.
    #@cython.nonecheck(False)
    #@cython.cdivision(True)
    cdef void superimpose_polygon(self, 
                                  float[::1] ya, 
                                  float[::1] xa, 
                                  unsigned int value) nogil:

        cdef float minya, maxya, minxa, maxxa
        minya, maxya = minmax(ya)
        minxa, maxxa = minmax(xa)

        cdef:
           ssize_t minr = int(max(0, minya))
           ssize_t maxr = int(ceil(maxya))
           ssize_t minc = int(max(0, minxa))
           ssize_t maxc = int(ceil(maxxa))
           unsigned short r, c

        maxr = min(self.height -1, maxr)
        maxc = min(self.width  -1, maxc)

        #with nogil:
        for r in range(minr, maxr+1):
            for c in range(minc, maxc+1):
                if (point_in_triangle(xa, ya, c, r+0.5) or 
                   point_in_triangle(xa, ya, c+0.5, r) or 
                   point_in_triangle(xa, ya, c-0.5, r) or 
                   point_in_triangle(xa, ya, c, r-0.5) or 
                   point_in_triangle(xa, ya, c+0.5, r+0.5) or 
                   point_in_triangle(xa, ya, c-0.5, r-0.5) or 
                   point_in_triangle(xa, ya, c-0.5, r+0.5) or
                   point_in_triangle(xa, ya, c+0.5, r-0.5) or  
                   point_in_triangle(xa, ya, c, r)
                   ):

                    self.image[r, c] = value
                    self.rr[self.n] = r
                    self.cc[self.n] = c
                    self.ff[self.n] = value
                    self.n +=1

          #print('sizes')
          #print(self.n)
          #print(maxc-minc+1)
          #print(maxr-minr+1)
          #print('pixel')
          #print(value)
          #print(r, self.rr[self.n-1])
          #print(self.rr[self.n])
          #input('raster pixel enter')
        #return rr, cc


# END
