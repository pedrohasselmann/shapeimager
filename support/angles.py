#!/usr/bin/python
# -*- coding: utf-8 -*-

__version__="1.0.0"

__doc__='''
  ############################ ILLUMINATION ANGLES FUNCTIONS ##############################
  ### Filename: angles.py
  ### Author: Pedro H. A. Hasselmann
  ###
  ### Compute the illumination angles from incidence, emergence, phase and azimuth angles.
  ### Spherical coordinates and local solar time are included.
  ###
  #########################################################################################
'''
### FUNCTIONS ####

def phase(i, e, phi):
  ''' Phase Angle.
  
     Parameters
     ==========
     i: incidence, e: emergence, phi: azimuth 
  '''
  from numpy import sign, arccos, cos, sin
  return sign(e)*arccos( cos(i)*cos(e) + sin(i)*sin(e)*cos(phi) )


def azimuth(i, e, phase):
  ''' Azimuth.
  
     Parameters
     ==========
     i: incidence, e: emergence, phase: phase angle 
  '''
  from numpy import cos, arccos, sin
  return arccos( (cos(phase) -cos(i)*cos(e))/(sin(i)*sin(e)) )


def cos_illum_lat(i, e, phi):
  ''' 
     Luminance Latitude.
     
     Parameters
     ==========
     i: incidence, e: emergence, phi: azimuth
  '''
  from numpy import cos, arccos, sin, tan, sqrt 
  sin2e_sin2i = sin(2e0*e) * sin(2e0*i)
  sinie2 = ( sin(i + e) )**2
  cosphi2 = ( cos(phi/2e0) )**2
    
  term1 = sinie2 -cosphi2*sin2e_sin2i
    
  return  sqrt( term1/(term1 + ( sin(i)*sin(e)*sin(phi) )**2 ) )


def cos_illum_lon(e, coslat):
  ''' 
     Luminance Longitude. 
  '''
  from numpy import cos
  # Change sign if not contained by the photometric equator
  #sign = ones(ph.shape)
  #boolean = (ph < meridian_phase).values
  #sign[boolean] = -1e0
  #if isnan(meridian_phase): sign = -1e0 
  return cos(e)/coslat#*sign


def spherical_coord(X, Y, Z):
  ''' 
  Spherical coordinates: latitude, longitude, radius
  
  Paramters
  =========
  X, Y, Z 

  https://stackoverflow.com/questions/4116658/faster-numpy-cartesian-to-spherical-coordinate-conversion
  '''
  from numpy import degrees, sqrt, arctan2, arccos
  xy = X**2 + Y**2
  r = sqrt(xy + Z**2)
  lon = arctan2(Y, X)
  lat = arctan2(Z, sqrt(xy))
  return lat, lon, r


def orthorectification(l, alt0, alt1, res, c=(1014,1014)):
   '''
      Orthorectification.
      Correct a length l by distortion caused by viewing angle and distance.
      
      l : array of (X0,Y0,X1,Y1)
      alt0, alt1 : altitude in respect to points (X0, Y0) and (X1, Y1)
      res : pixel angular resolution (x_res, y_res)
      c : central reference (generally the middle of the image)
   '''
   from numpy import array, int32, sqrt, cos, sin, fabs, zeros, where

   ax1 = res[0]*(l[:,[0,2]]-c[0])    # X-column image
   ax2 = res[1]*(l[:,[1,3]]-c[1])    # Y-column image
   ax11, ax22 = ax1.copy(), ax2.copy()

   # Re-order the vectors
   # azimuth vector length must be larger than theta vector length
   c = where(ax1[:,0]<ax1[:,1])
   ax11[:,0][c]=ax1[:,1][c]
   ax11[:,1][c]=ax1[:,0][c]

   c = where(ax2[:,0]>ax2[:,1])
   ax22[:,0][c]=ax2[:,1][c]
   ax22[:,1][c]=ax2[:,0][c]
   az, theta = ax11.copy(), ax22.copy()
    
   c=fabs(az[:,0]-az[:,1])>fabs(theta[:,0]-theta[:,1])
   az[c,:]=ax22[c,:]
   theta[c,:]=ax11[c,:] 
    
   # Orthorectified length
   r2 = alt0**2 +alt1**2 -2e0*alt0*alt1*(sin(theta[:,0])*sin(theta[:,1])*cos(az[:,0]-az[:,1]) +cos(theta[:,0])*cos(theta[:,1]))

   return sqrt(r2), theta, az


def solid_angle_tri(v1_v2, d):
    '''
       Compute projected triangular Solid Angle.
       Unit: stereoradians.

       Parameters
       ==========
       v1_v2 : 2D frame-projected coordinates.
       d     : 1D observer distance at facet center.
    '''
    from numpy import float32, sum, sqrt
    
    s1 = sqrt(sum((v1_v2[:,0,:]-v1_v2[:,1,:])**2, axis=1))
    s2 = sqrt(sum((v1_v2[:,0,:]-v1_v2[:,2,:])**2, axis=1))
    s3 = sqrt(sum((v1_v2[:,1,:]-v1_v2[:,2,:])**2, axis=1))

    p = (s1+s2+s3)/2e0

    omega = sqrt(p*(p-s1)*(p-s2)*(p-s3))/(d**2)
    return omega.astype(float32)


# END
