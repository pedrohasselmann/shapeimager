#cython: language_level=3, boundscheck=False
#distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
#distutils: language = c++
# -*- coding: utf-8 -*-

__version__="0.2"

__doc__='''
  ############## CYTHON RAY-TRACING #############################
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
  ### , boundscheck=False, wraparound=False, nonecheck=False, cdivision=True
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

########################################################################
#################### RAY TRACING NOGIL #################################
########################################################################
def plot_hidden(fxy, # v1_v2 (2,3,N)
                fc,  # facet_center (N)
                hidden # uint8, hidden facets
                ):
    import matplotlib.pyplot as plt
    
    #ind = fc[2,:].argsort()
    #fxy = fxy[...,ind]
    #fc = fc[...,ind]
    #hidden = hidden[ind]
    for v in range(fxy.shape[2]):
      if fc[0,v]!=-999.:
            if hidden[v] == 1: 
               c='black'
               print(v, 1, fc[2,v])
            else:
               c='green'
               print(v, 0, fc[2,v])
            t = plt.Polygon(fxy[...,v].T, color=c, alpha=0.3)
            plt.gca().add_patch(t)
            plt.scatter(fc[0,v], fc[1,v], marker='.', s=0.3, color='red')
            plt.text(fc[0,v], fc[1,v], str(v), fontsize=8, color='red')
    print()
    plt.show()
    plt.close()


cpdef raytracing_loop1_c(buff,  # Buffer from pd.Groupby    # facets intercepted in the ray-pixel
                       const double[:,::1] fc,    # facet centroid + distance
                       const double[:,:,::1] fxy, # facet projected coordinates
                        ):
    '''
       Pixel-wise ray tracing (adapted for pd.Groupby.agg buff).
       Loop through pixels containing overlayed facets and identify them.
       
       Input
       =====
       buff : pandas.Dataframe.GroupBy (X, Y, facet)
       fc   : ndarray, facet centroids & distance to source
       fxy  : ndarray, projected facet vertices
       
       Output
       ======
       ndarray : ray-tracing binary flags 

       Binary Flag:
       1 : Visible facet
       0 : Blocked facet
    '''
    #print(buff)
    #from time import time
    cdef:
       unsigned int[::1]   f = np.asarray(buff.index, order='C', dtype=uint32)
       ssize_t             i, j 
       ssize_t             N = f.shape[0]
       unsigned int[:,::1] o = np.empty((2,N), dtype=uint32, order='C')
       uint8_t[::1]        h = np.ones(N, dtype=uint8, order='C')
       double[:,::1]       c = np.empty((3,N), dtype=float64, order='C')
       double[:,:,::1]     v = np.empty((2,3,N), dtype=float64, order='C')

    
    
    with nogil:
       for i in range(N):
          j = f[i]
          c[...,i] = fc[...,j-1]
          v[...,i] = fxy[...,j-1]

    h = raytracing_c(c, v, h, N)

    with nogil:
       for i in range(N):
          o[0,i] = f[i]
          o[1,i] = h[i]

    #print(np.asarray(o))
    #plot_hidden(np.asarray(v), 
    #            np.asarray(c), 
    #            np.asarray(h))
    return o


cpdef raytracing_loop2_c(const unsigned int[:,::1] fp,  # Y X facets intercepted in the ray-pixel
                         const double[:,::1] fc,    # facet centroid + distance
                         const double[:,:,::1] fxy, # facet projected coordinates
                         int buff_size
                         ):
    '''
       Pixel-wise ray tracing.
       Loop through pixels containing overlayed facets and identify them.
       
       Input
       =====
       fp   : Y X facet
       fc   : ndarray, facet centroids & distance to source
       fxy  : ndarray, projected facet vertices
       
       Output
       ======
       ndarray : ray-tracing binary flags 

       Binary Flag:
       1 : Visible facet
       0 : Blocked facet
    '''
    cdef:
       ssize_t i, j, n = 0
       ssize_t N = fp.shape[1]
       unsigned int f
       unsigned int[:,::1] o = np.empty((2,N), dtype=uint32, order='C')
       uint8_t[::1]        h = np.ones(buff_size, dtype=uint8, order='C')
       double[:,::1]       c = np.empty((3,buff_size), dtype=float64, order='C')
       double[:,:,::1]     v = np.empty((2,3,buff_size), dtype=float64, order='C')

    o[:,:] = 0
    c[:,:] = -999.
    v[:,:] = -999.

    with nogil:

      #f = fp[2,0]
      c[...,0] = fc[...,0]
      v[...,0] = fxy[...,0]
      o[0,:] = fp[2,:]

      for i in range(1,N):

            if (fp[0,i] == fp[0,i-1]) and (fp[1,i] == fp[1,i-1]):
               #print(np.asarray(fp[:,i]))

               #f = fp[2,i]
               c[...,n] = fc[...,i]
               v[...,n] = fxy[...,i]
               #o[0,i] = f
               n+=1

            else:

               #with nogil:
               h[:] = 1
               h = raytracing_c(c, v, h, n)

               #plot_hidden(np.asarray(v), 
               #            np.asarray(c), 
               #            np.asarray(h))

               for j in range(n):
                           o[1,i-n+j+1] = h[j]
                           #print(np.asarray(o[:,i-n+j+1]))
               
               #print(i, n)
               #print(j)
               #input('\n')


               c[:,:] = -999.
               v[:,:] = -999.

               n = 0
               #f = fp[2,i]
               c[...,n] = fc[...,i]
               v[...,n] = fxy[...,i]


    return np.asarray(o)





# Static variables in raytracing_c
cdef:
  double[::1]   pt1 = np.empty(3, dtype=float64, order='C'),
  double[::1]   pt2 = np.empty(3, dtype=float64, order='C'),
  double[:,::1] vtc = np.empty((2,3), dtype=float64, order='C'),
  double[:,::1] tri = np.empty((2,3), dtype=float64, order='C')

cdef uint8_t[::1] raytracing_c(const double[:,::1] ray,  # Rays: facet centers + distance
                                const double[:,:,::1] fxy, # Facet projected coordinates
                                uint8_t[::1] interc,
                                unsigned short L,
                                ) nogil:
    '''
       Pixel-wise ray tracing (shadow casting tool).


       Parameters
       ==========
       ray  : Rays: facet centers + distance to facet
       fxy : Facet projected coordinates

       Output
       ======
       bool : a ray-tracing binary flag

       Binary Flag:
       1 : Visible facet
       0 : Blocked facet

    '''
    cdef:
       ssize_t n, m
       #ssize_t L = ray.shape[1]
       bint answer1, answer2

    #with nogil:
    for n in range(L-1):#, num_threads=thrn, schedule='static'):

         #if ray[0,n] != -999.:

            #print(np.asarray(ray[:,n]))
            pt1[...] = ray[...,n]
            vtc[...] = fxy[...,n]

            for m in range(n+1,L):
               #if ray[0,m] != -999.:

                    tri[...] = fxy[...,m]
                    pt2[...] = ray[...,m]
                    answer1 = inside_triangle(pt1, tri[:,2], tri[:,1], tri[:,0])
                    answer2 = inside_triangle(pt2, vtc[:,2], vtc[:,1], vtc[:,0])

                    #print(n, m)
                    #print(answer1, int(answer1))
                    #print(answer2, int(answer2))

                    if (answer1==1)&(answer2==1)&(pt1[2]<pt2[2]):
                            interc[m] = 0

                    elif (answer1==1)&(answer2==1)&(pt1[2]>pt2[2]):
                            interc[n] = 0

                    elif (answer1==1)&(answer2==0)&(pt1[2]>pt2[2]):
                            interc[n] = 0

                    elif (answer1==0)&(answer2==1)&(pt2[2]>pt1[2]):
                            interc[m] = 0

                    elif (answer1==0)&(answer2==0):
                            pass

    return interc


cdef bint inside_triangle(double[::1] pt, 
                           double[:] v0, 
                           double[:] v1, 
                           double[:] v2) nogil:
    # https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
    cdef:
       double dX = pt[0] - v0[0]
       double dY = pt[1] - v0[1]
       double dX20 = v2[0] - v0[0]
       double dY20 = v2[1] - v0[1]
       double dX10 = v1[0] - v0[0]
       double dY10 = v1[1] - v0[1]

       double s = (dY20*dX) - (dX20*dY)
       double t = (dX10*dY) - (dY10*dX)
       double D = (dX10*dY20) - (dY10*dX20)
 
    if D<0:
      return (s<=0) & (t<=0) & ((s + t)>=D)
    else:
      return (s>=0) & (t>=0) & ((s + t)<=D)



# END