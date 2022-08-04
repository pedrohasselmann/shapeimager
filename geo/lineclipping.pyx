#cython: language_level=3, boundscheck=False, nonecheck=False, wraparound=False, cdivision=True
#distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
#distutils: language = c++
# -*- coding: utf-8 -*-

__version__="0.2"

__doc__='''
  ############## CYTHON LINE CLIP COHEN-SUTHERLAND ###################
  ### Author: Pedro Henrique A. Hasselmann, Dr. (OBSPM-Lesia)
  ### Last Modified: April, 2021
  ###
  ### requirements: 
  ### 
  ###
  ### Main Functions:
  ###
  ###
  ###
  ####################################################################
'''
# Cythonize
import cython
cimport cython
from cython.parallel import prange, parallel
from cython cimport view as cview
from cpython cimport array
from libc.math cimport ceil, fabs, sqrt, acos, pi
from libc.string cimport memcpy, memset
#from libcpp.vector cimport vector

import numpy as np
cimport numpy as np
from numpy import ndarray, uint8, uint16, uint32, uint64, int32, float32, float64, uint8
from numpy cimport ndarray, uint8_t, uint16_t, uint32_t, int32_t, int64_t, float32_t, float64_t
from numpy.math cimport INFINITY
np.import_array()

# Local CIMPORT
from .raytracing cimport inside_triangle


# NUMBER OF THREADS
cdef short thrn=4

################
### FUNCTIONS ##
################

cpdef plot_clipping(fr, v1_v2, yx0, w, h, clip):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    # create the figure and the axis in one shot
    #fig, fr = plt.subplots(1,figsize=(5,5))
    r = mpatches.Rectangle(yx0, w, h, edgecolor='black',facecolor='black', fill=False, alpha=0.8)
    t = mpatches.Polygon(v1_v2, fill=True, alpha=0.2)
    fr.add_patch(r)
    fr.add_patch(t)
    plt.scatter(yx0[0], yx0[1], marker='.', s=1, color='black')
    plt.scatter(clip[::2], clip[1::2], marker='o', s=20, color='red')

    #print(fr.patches)
    plt.draw()


##### COHEN-SUTHERLAND ######
cdef short INSIDE, LEFT, RIGHT, LOWER, UPPER
INSIDE, LEFT, RIGHT, LOWER, UPPER = 0, 1, 2, 4, 8

cdef short _getclip(double xa, 
                    double ya,
                    double xmin, 
                    double ymax, 
                    double xmax,
                    double ymin, 
                    ) nogil:
        #if dbglvl>1: print('point: '),; print(xa,ya)
        cdef short p = INSIDE  # default is inside

        # consider x
        if xa < xmin:
            p |= LEFT
        elif xa > xmax:
            p |= RIGHT

        # consider y
        if ya < ymin:
            p |= LOWER  # bitwise OR
        elif ya > ymax:
            p |= UPPER  # bitwise OR
        return p

cdef (double, double, double, double) cohensutherland(double xmin, 
                                                  double ymax, 
                                                  double xmax,
                                                  double ymin, 
                                                  double x1, 
                                                  double y1, 
                                                  double x2, 
                                                  double y2) nogil:
    """Clips a line to a rectangular area.
    This implements the Cohen-Sutherland line clipping algorithm.  xmin,
    ymax, xmax and ymin denote the clipping area, into which the line
    defined by x1, y1 (start point) and x2, y2 (end point) will be
    clipped.
    If the line does not intersect with the rectangular clipping area,
    four None values will be returned as tuple. Otherwise a tuple of the
    clipped line points will be returned in the form (cx1, cy1, cx2, cy2).
    """

# check for trivially outside lines
    cdef short k1 = _getclip(x1, y1, xmin, ymax, xmax, ymin)
    cdef short k2 = _getclip(x2, y2, xmin, ymax, xmax, ymin)
    cdef short opt
    cdef double x, y  

# %% examine non-trivially outside points
    # bitwise OR |
    while (k1 | k2) != 0:  # if both points are inside box (0000) , ACCEPT trivial whole line in box

        # if line trivially outside window, REJECT
        if (k1 & k2) != 0:  # bitwise AND &
            #if dbglvl>1: print('  REJECT trivially outside box')
            # return nan, nan, nan, nan
            return -999., -999., -999., -999.

        # non-trivial case, at least one point outside window
        # this is not a bitwise or, it's the word "or"
        opt = k1 or k2  # take first non-zero point, short circuit logic
        if opt & UPPER:  # these are bitwise ANDS
            x = x1 + (x2 - x1) * (ymax - y1) / (y2 - y1)
            y = ymax
        elif opt & LOWER:
            x = x1 + (x2 - x1) * (ymin - y1) / (y2 - y1)
            y = ymin
        elif opt & RIGHT:
            y = y1 + (y2 - y1) * (xmax - x1) / (x2 - x1)
            x = xmax
        elif opt & LEFT:
            y = y1 + (y2 - y1) * (xmin - x1) / (x2 - x1)
            x = xmin
        else:
            raise RuntimeError('Undefined clipping state')

        if opt == k1:
            x1, y1 = x, y
            k1 = _getclip(x1, y1, xmin, ymax, xmax, ymin)
            #if dbglvl>1: print('checking k1: ' + str(x) + ',' + str(y) + '    ' + str(k1))
        elif opt == k2:
            #if dbglvl>1: print('checking k2: ' + str(x) + ',' + str(y) + '    ' + str(k2))
            x2, y2 = x, y
            k2 = _getclip(x2, y2, xmin, ymax, xmax, ymin)

    return x1, y1, x2, y2



# Calculate Polygon Area using Shoelace formula (order of the corners matter! clockwise or counterclockwise)
cdef float PolygonArea(short n, double[:,::1] corners) nogil:
    cdef:
        ssize_t i=0, j=0
        float area = 0.0

    while (corners[i,0] != -999.):
          j = (i + 1) % n
          area += corners[i][1] * corners[j][2]
          area -= corners[j][1] * corners[i][2]
          i+=1
    area = fabs(area) / 2.0
    return area



cdef double vec_angle(double[::1] v1, double[::1] v2, double[::1] ct) nogil:
    # Angle between first vector & other vectors. 
    # Sort vectors in counterclockwise.
    cdef:
       double v11 =  v1[1]-ct[0]
       double v21 = -v2[1]-ct[0]
       double v12 =  v1[2]-ct[1]
       double v22 = -v2[2]-ct[1]
       double cos0 = v11*v21 + v12*v22
       double    a = sqrt(v11*v11 + v12*v12)
       double    b = sqrt(v21*v21 + v22*v22)

    if(v11*v22 - v12*v21 < 0.):
       return -acos(cos0/fabs(a*b))
    else:
       return acos(cos0/fabs(a*b))




cdef unsigned short[:,::1]        c = np.array(((0,1),(1,2),(2,0)), dtype=uint16, order='C')
cdef double[:,::1]               cc = np.empty((10,3), dtype=float64, order='C')
cdef double[::1]              acos0 = np.empty(10, dtype=float64, order='C')
cdef double[:,::1]              rcc = np.empty((8,3), dtype=float64, order='C')
cdef double[:,::1]               vc = np.empty((3,2), dtype=float64, order='C')
cdef double[::1]                 xy = np.empty(2, dtype=float64, order='C')
cdef double[:,::1]           corner = np.empty((4,2), dtype=float64, order='C')
cdef double[::1]                 ct = np.zeros(2, dtype=float64, order='C')

cpdef clipping_loop_c(unsigned int[:,::1] fp, 
                      double[:,:,::1]     v1v2,
                      (double, double)    px,
                      float[::1]          gx,
                      float[::1]          gy
                      ):
    """
       Looping Choen-Sutherland clipping algorithm through all pixeled facets (fp) in the FOV. 

       Parameters
       ==========
       fp   : facet-pixel link
       v1v2 : vertices-facet coordinates
       gx   : grid on X
       gy   : grid on Y

       Output
       ======

    """
    import matplotlib.pyplot as plt

    cdef:
      ssize_t i, j, k, u, t, z, h
      ssize_t N = fp.shape[0]
      bint answer
      double aa, a, b, acos0u
      (double, double, double, double)  ipoints = (0., 0., 0., 0.)
      #float[:,::1] intersects = np.empty((6*v1v2.shape[0], 3), dtype=float32, order='C')
      float[:,::1] px_area = np.zeros((N,2), dtype=float32, order='C')

    #intersects[:] = -999.

    with nogil:#, parallel(num_threads=thrn):
        #j = 0
        for i in range(N):#, schedule='static'):

            cc[:,:]   = -999.
            rcc[:,:]  = -999.
            acos0[:]  = -999.

            vc[...] = v1v2[i,...]

            px_area[i,0] = fp[i,2]

            xy[0] = gx[fp[i,1]]
            xy[1] = gy[fp[i,0]]

            #if (i<N):
            #  if (fp[i,2] != fp[i+1,2]):
            #        j=j+1

            #print('loop')
            #print(i, j, np.array(fp[i,:]))
            #print(xy[0], xy[1])
            #print(np.array(vc))

            #fig, fr = plt.subplots(1,figsize=(5,5))
            for k in range(3):
                ipoints = cohensutherland(xmin=xy[0]-px[0],    ymax=xy[1], 
                                          xmax=xy[0],          ymin=xy[1]-px[0], 
                                          x1=vc[c[k,0],0],  y1=vc[c[k,0],1], 
                                          x2=vc[c[k,1],0],  y2=vc[c[k,1],1])

                if (ipoints[0] != -999.) and (vc[k,0] != ipoints[0]) and (vc[k,1] != ipoints[1]):
                    cc[2*k,0] = fp[i,2]
                    cc[2*k,1] = ipoints[0]
                    cc[2*k,2] = ipoints[1]
                    z = 2*k

                if (ipoints[2] != -999.) and (vc[k,0] != ipoints[2]) and (vc[k,1] != ipoints[3]):
                    cc[2*k+1,0] = fp[i,2]
                    cc[2*k+1,1] = ipoints[2]
                    cc[2*k+1,2] = ipoints[3]
                    z = 2*k+1

                #if (ipoints[0] != -999.) and (ipoints[2] != -999.):
                #    plot_clipping(fr, vc, (xy[0], xy[1]), -px[0], -px[1], clip)

            # Pixel corners are inside or outside the triangle?
            corner[0,0] = xy[0] -px[0]
            corner[0,1] = xy[1] -px[1]
            corner[1,0] = xy[0] -px[0]
            corner[1,1] = xy[1]
            corner[2,0] = xy[0]
            corner[2,1] = xy[1]
            corner[3,0] = xy[0]
            corner[3,1] = xy[1] -px[1]
            # Loop through corners
            t=0
            for u in range(4):
                answer = inside_triangle(corner[u,:], vc[0,:], vc[1,:], vc[2,:])
                if answer == 1:
                    cc[2*k+2+u,0] = fp[i,2]
                    cc[2*k+2+u,1:] = corner[u,:]
                    t+=1

            #print('first ', z)
            #print(np.asarray(cc))

            # return if the pixel is embedded inside the facet
            if t == 4:
                px_area[i,1] = 1.
                #print(1, t, N, 'area ')
                #plt.close()
                continue 

            # Barycenter (average)
            h=0
            ct[:] = 0.
            for u in range(10):
                if cc[u,0] != -999.:
                    ct[0] += cc[u,1]
                    ct[1] += cc[u,2]
                    h=h+1
            # None inside the pixel
            if h == 0:
                #print(2, h, N, 'area ')
                #plt.close()
                continue

            ct[0] = ct[0]/h
            ct[1] = ct[1]/h
            #print('barycenter ', ct[0], ct[1])

            # Re-arrange corner vector - Clockwise/counterwise & all -999. vectors at the end.
            for u in range(10):
                if cc[u,0] != -999.:
                    # Angle between first vector & other vectors. 
                    # Sort vectors in counterclockwise.
                    acos0[u] = vec_angle(cc[u,:], cc[z,:], ct) + pi
                    #print(u, acos0[u], np.degrees(acos0[u]).astype(int32))

            t=0
            for u in range(10): # roll through every side of the pixel/square.
                if (cc[u,0] != -999.):
                    t=1
                    for k in range(10):
                        if (cc[k,0] != -999.) and (acos0[u]>acos0[k]):
                            t+=1

                    for k in range(h):
                         if t == h-k:
                            #print(u, t, h-k)
                            rcc[k,:] = cc[u,:]
                t=0
        
            #print()
            #print(h)
            #print(np.asarray(rcc))

            # Polygon Area
            aa = PolygonArea(h, rcc)
            px_area[i,1] = aa
            #print(3, N, 'area ', aa, aa/(px[0]*px[1]),'\n \n')
            #input('enter')
            #plt.show()
            #plt.close()

            cc[:,:]   = -999.
            rcc[:,:]  = -999.
            acos0[:]  = -999.



    return np.asarray(px_area) #, np.asarray(intersects)
