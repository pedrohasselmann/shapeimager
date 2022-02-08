#!/usr/bin/python
# -*- coding: utf-8 -*-

__version__="0.2"

doc = '''
  #########################################################################################
  #  author: Pedro H. A. Hasselmann, Rio de Janeiro (LESIA-OBSPM/ON-MCT)
  #  script: shapeimager.py
  #
  #  2D Shape Model Imager.
  #
  # Requirements:
  # cython 0.29, pandas 1.2.4, numpy 1.20, scipy 1.6, astropy 4.2, matplotlib 3.3, PIL
  #
  ##########################################################################################
'''
from .. import *
from ..main import *
#from imager_utils import raster_c, raytracing_loop_c, raytracing_c, weighted_average_c

###############
## 2D IMAGER ## 
###############
class Imager(object):
  '''
     2D Imager Class.
     Reconstruction of surfaces onto Field-of-Views.
     
     From Instrument specifications, positioning vectors, and shape models/DTMs the target image is reconstructed.
  '''
  def __init__(self, 
               S, 
               cam_matrix,
               boresight,
               source_pos, 
               obs_pos,
               visible=False, 
               illuminated=False,
               raytrace=False,
               shaded=0,
               occ=0,
               ):
    '''
       Remove non-illuminated facets and hemispherically-hidden facets.
       
       Parameters
       ==========
       S          : shape model class
       cam_matrix : Camera Rotation Matrix (numpy.array, 3x3)
       boresight  : Camera boresight vector (numpy.array, 1x2)
       ccd        : ccd pixel format (integer)
       FOV        : FOV quaternion
       source_pos : light source position vector
       obs_pos    : observer position vector
       
       -Boolean Flags-
       visible    : hemispherical visibility
       illuminated: direct illumination
       shaded     : shadow cast

       -Float Flags-
       precision : Degree of precision in the hidden_by function.
       It divides the average projected separation of facets. 
       As it get higher, more precise the limb and shaded profiles.
    '''
    from numpy import uint32, int32, arange
    self.cam_matrix = cam_matrix
    self.boresight = boresight
    self.source_pos = source_pos
    self.obs_pos = obs_pos

    # SHAPE MODEL & GEOMETRY
    assert hasattr(S, 'n_')
    #S.normal_vector(True)

    self._camera_framing = S.camera_framing
    self._plane          = S.plane
    self._irradiance     = S.irradiance

    self.inc_ = S.project(source_pos)[0]
    self.emi_ = S.project(-obs_pos)[0]
    self.pha_ = S.phase(obs_pos, -source_pos)[0]
    #self.distance_ = S.distance(-obs_pos)

    # FACET ID
    self.facetid_ = arange(1,S.facets_.shape[0]+1).copy(order='C').astype(uint32)

    # FIRST-ORDER VISIBILITY & ILLUMINATION
    c1 = self.inc_<0.49999*3.14159265359 #(pi/2)
    c2 = self.emi_<0.49999*3.14159265359 #(pi/2)

    # HIDDEN_BY
    self._set_raytrace(occ, shaded, raytrace, c1, c2)

    # FLAGS
    self._set_flags(visible, illuminated)

    # CAMERA FRAME
    self._framing(cam_matrix, boresight, -obs_pos)

  ########################
  ### FLAGS & SETTINGS ###
  ########################
  def _set_flags(self, visible, illuminated):
    from numpy import ones, bool_
    # FIRST-ORDER VISIBILITY & ILLUMINATION
    c1 = self.inc_<0.49999*3.14159265359 #(pi/2)
    c2 = self.emi_<0.49999*3.14159265359 #(pi/2)

    # FLAGS
    if   visible == True  and  illuminated == False:  
      self.active_ = c2
      
    elif visible == False and  illuminated == True:   
      self.active_ = c1
      
    elif visible == True  and  illuminated == True:
      self.active_=(c2&c1)
       
    elif visible == False and  illuminated == False:
      self.active_ = ones(self.inc_.shape[0]).astype(bool_)



  def _set_raytrace(self, occ, shaded, raytrace, c1, c2):
    # HIDDEN_BY (OCCULTED/OVERLAYED FACETS)
    if occ!=0:
      print('occultation')
      self.occult_ = self._hidden_by(-self.obs_pos, occ, c2, raytrace=True)
    else:
      self.occult_ = True
    
    # HIDDEN_BY (SHADED FACETS)
    if shaded !=0:
        print('shaded by')
        self.shaded_ = self._hidden_by(1e-4*self.source_pos, shaded, c1, raytrace)
    else:
        self.shaded_ = True

  # Resetting flags
  def reset(visible=True, illuminated=True, occ=0, shaded=0):
    self.set_flags(visible, illuminated)
    self.set_raytrace(occ, shaded)
    self._framing(self.cam_matrix, self.boresight, -self.obs_pos)
    
  ####################
  ### CAMERA SCENE ###
  ####################
  def _framing(self, cam_matrix, boresight, pos):
    '''
       Shape Model in the Camera frame.
    '''
    assert hasattr(self, 'active')
    assert hasattr(self, 'occult')
    assert hasattr(self, 'shaded')

    self.v1_v2_, self.d_ = self._camera_framing(cam_matrix, boresight, pos)

    # Imaged facets
    in_frame = self.shaded_&self.occult_&self.active_
    self.facetid_ = self.facetid_[in_frame]
    #self.v1_v2_ = v1_v2[self.facetid_-1,:,:]
    #self.d_ = d[self.facetid_-1,:]
    print('Active #facets ', self.facetid_.shape)



  ######################
  ### CORE FUNCTIONS ###
  ######################
  @staticmethod
  def _raster(vert, facetid, frame, grid, mode=u'pixelation', pools=4):
    '''
      Image rasterization with multiprocessing.
      
      -Input-
      =======
      vert
      facetid: ndarray, facet index
      frame  : 2-tuple, frame pixel size
      grid   : 2-tuple of ndarrays
      mode   : string, pixelation (facet number) or occultation (count)
      pools : int, number of active pools (check your processor)
      
    '''
    from .imager_utils import raster_c
    import multiprocessing as multiproc
    from scipy import stats
    from collections import deque
    from numpy import int32, uint32, float32, array, arange, split, array_split, vstack
    from functools import partial
    from time import time
    
    start = time()
    # Partitioning of vertices into pixels/cases (front)
    vert = vert.transpose(2,0,1).reshape(2,-1).copy(order='C') 
    image_vertex = stats.binned_statistic_2d(vert[1,:], vert[0,:], None, 'count', bins=grid, expand_binnumbers=True)
    vert_pix = array(split(image_vertex.binnumber.T, image_vertex.binnumber.shape[1]/3), order='F', dtype=float32).T
    image_vertex = 0
    print('_raster overhead time ', time()-start)

    # Rasterization
    start = time()
    image, facet_pix  = raster_c((0,vert_pix.shape[2]),
                                vert_pix, 
                                facetid.copy(order='C'), 
                                frame, 
                                mode)
    print('_raster computation time ', time()-start)

    if facet_pix.shape[0] == 0:
        raise Exception("No facets are visible or illuminated.")

    return image, facet_pix[facet_pix.all(axis=1)] 
    #######################################################


  def _hidden_by(self, o, prec, active, raytrace=False):
    '''
       Identify and flag facets hidden by other facets in respect to a given viewing vector.
       
       -Input-
       =======
       o : ndarray(1,3), viewing vector
       prec : int, precision
       active : ndarray(boolean), flags of active facets
       raytrace : bool
       
       -Output-
       ========
       ndarray(boolean) : refined active facet flags after ray-tracing
    '''
    from numpy import bool_, uint8, uint32, int32, float32, nanmin, nanmax, \
     nanmean, diff, fabs, hstack, vstack, arange, flatnonzero, zeros, in1d, unique, arange
    from time import time
    import matplotlib.pyplot as plt
    #print('active facets ',active[active==1].shape)
    
    start = time()
    v1_v2, z, base = self._plane(o, o)     # project to the viewing vector frame.
    fc = vstack((v1_v2.T.mean(1), z.mean(1).flatten())).copy(order='C')
    z = z[active,:].mean(axis=1).flatten().argsort(kind='mergesort')[::-1].copy(order='C') # Distance per facet, ascending order
    vert = v1_v2[active,...][z,...]   # Vertices: Convert 3D array to 2D array
    print('_hidden_by overhead time ',time() -start)

    # facet sizes
    delta = diff(vert, axis=1).transpose(2,0,1).reshape(2,-1)

    # Grid resolution based on the sizes of projected facets
    #v = v1_v2[self.facetid_-1,:,:].transpose(2,0,1).reshape(2,-1)
    #plt.plot(v[0,:], v[1,:], '.', markersize=0.3, color='black')
    #plt.show()
    i1 = nanmin(v1_v2, axis=0).min(axis=0) # X, Y inferior values
    i2 = nanmax(v1_v2, axis=0).max(axis=0) # X, Y upper values
    
    dX = fabs(delta[0,...])[delta[0,...]!=0.0]
    #print(dX.min(), dX.mean(), prec*(i2[0]-i1[0])/dX.mean())
    dY = fabs(delta[1,...])[delta[1,...]!=0.0]
    #print(dY.min(), dY.mean(), prec*(i2[1]-i1[1])/dY.mean())

    grid_x = arange(i1[0],i2[0], nanmean(dX)/prec)
    grid_y = arange(i1[1],i2[1], nanmean(dY)/prec)

    # Rasterization
    frame = (int(prec*(i2[1]-i1[1])/nanmean(dY)), int(prec*(i2[0]-i1[0])/nanmean(dX)))
    front_image, facet_pix = self._raster(vert, 
                                          self.facetid_[active][z], 
                                          frame, 
                                          (grid_y,grid_x), 
                                          mode=u'pixelation')
    z = 0
    vert = 0
    
    # Separate front and back facets
    front_facets = unique(front_image.ravel()[flatnonzero(front_image)]).copy(order='C')
    #print('first-order front facets ',front_facets.shape)

    ############### BEGIN RAY-TRACING #################################
    if raytrace == True:
      from .raytracing import raytracing_loop1_c, raytracing_loop2_c

      interc_bool = in1d(facet_pix[...,2], front_facets, assume_unique=False, invert=True) # Select only pixel containing intercepted facets
      interc_pix = unique(facet_pix[interc_bool,:2], axis=0).copy(order='C')  # list of intercepted pixels

      # Get all intercepted facets (front and back) per pixel
      interc_bool = inNd(facet_pix[...,(0,1)].copy(order='C'), interc_pix, assume_unique=False, invert=False)
      interc_pix = unique(facet_pix[interc_bool,:], axis=0).copy(order='C')   # list of intercepted pixels
      interc_facets = unique(interc_pix[...,2]).copy(order='C') # list of intercepted facets
      unique_bool = in1d(front_facets, interc_facets, assume_unique=True, invert=True)
      #print('intercepted pixels ',interc_pix.shape)
    
      # Build image only with the intercepted pixels
      #interc_image = self._raster(v1_v2[interc_facets-1,...], 
      #                            interc_facets, 
      #                            frame, 
      #                            (grid_y,grid_x), 
      #                            mode=u'pixelation')[0]
      # test
      #f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, sharex=True)
      #ax1.imshow(fc[2,front_image-1]-fc[2,...].min(), interpolation='none')
      #ax2.imshow(fc[2,interc_image-1]-fc[2,...].min(), interpolation='none')
      #plt.show()


      # Ray-tracing (refine, again, and again)
      # Pandas required
      # https://stackoverflow.com/questions/48875601/group-by-apply-with-pandas-and-multiprocessing
      # https://stackoverflow.com/questions/34967759/pandas-dataframe-groupby-apply-function-returning-an-array-and-map-back-the-re
      start = time()
      #interc_df = pd.DataFrame(interc_pix[...,(0,1)], index=interc_pix[...,2], columns=('y','x'), dtype=uint32)
      #interc_df['vis'] = interc_df.index
      #interc_df = interc_df.groupby(by=['y','x'], as_index=False, sort=False, dropna=False)

      interc_df = pd.DataFrame(interc_pix, columns=('y','x','f'))
      interc_df.sort_values(['y','x'], ascending=[True, True], inplace=True, kind='mergesort')
      memsize = int(interc_df.groupby(by=['y','x']).count().max())
      #print(interc_df.head(30))
      #print(interc_df.reset_index(level=0).head(10))
      #print(interc_df.agg(lambda x: print(x)))

      #interc = hstack()
      #interc = hstack(interc_df['vis'].agg(raytracing_loop_c,
      #                                       fc  = fc,                        # fc
      #                                       fxy = v1_v2.T.copy(order='C'),   #fxy
      #                                       engine='cython'                     
      #                                      )['vis'].to_numpy())

      ifacets = interc_df['f'].to_numpy()-1
      interc = raytracing_loop2_c(interc_df.T.to_numpy(dtype=uint32, copy=False).copy(order='C'), # fp
                                  fc[:,ifacets].copy(order='C'),                                  # fc
                                  v1_v2[ifacets,...].T.copy(order='C'),                           # fxy
                                  buff_size = memsize
                                  )

      # GroupBy
      interc_s = pd.Series(interc[1,...], index=interc[0,...], dtype=float32).groupby(level=0).mean().round(0)
      print('ray-tracing computation time ',time() -start)
      #print(interc.head(10))
      front_facets = hstack((front_facets[unique_bool], interc_s[interc_s>0].index)).astype(uint32)
      #print('second-order front facets ', front_facets.shape)
    ############### END RAY-TRACING #################################

    # Output: Boolean array: refined active facets 
    active = zeros(active.shape[0], dtype=bool_, order='C')
    active[front_facets-1] = True
    return active


  ##################
  ### MATPLOTLIB ###
  ##################
  def plot_v(self, FOV, ccd, name, factor=1, save=False):
    '''
      Plot frame-projected shape model vertices.

      -Input-
      factor : int, dpi resolution multiplying factor
    '''
    import matplotlib.pyplot as plt
    from matplotlib.transforms import Bbox
    from matplotlib.ticker import NullLocator
    from matplotlib import patches
    from numpy import array
    assert hasattr(self, 'v1_v2_')

    self.ccd = ccd
    qua = FOV[4]
    dpi = factor*96

    # PLOT -- BUILD EQUIVALENT IMAGE
    v = self.v1_v2_[self.facetid_-1,:,:].transpose(2,0,1).reshape(2,-1) # Convert 3D array to 2D array
    
    # https://stackoverflow.com/questions/13714454/specifying-and-saving-a-figure-with-exact-size-in-pixels
    fig = plt.figure(figsize=(ccd[0]/dpi, ccd[1]/dpi), dpi=dpi)
    frame = fig.add_subplot(111)
    frame.plot(v[0,:], v[1,:], '.', markersize=0.3, color='black')
    
    if FOV[0] == "RECTANGLE":
      frame.set_ylim(qua[:,1].min(),qua[:,1].max())
      frame.set_xlim(qua[:,0].min(),qua[:,0].max())

    elif FOV[0] == "CIRCLE":
      frame.set_ylim(-qua[0][0],qua[0][0])
      frame.set_xlim(-qua[0][0],qua[0][0])

      # Make a circle
      circ = patches.Circle((0, 0), qua[0][0], facecolor='none', edgecolor='red', fill=False, linestyle='-')
      frame.add_patch(circ) # Plot the outline

    if save:
      # https://stackoverflow.com/questions/11837979/removing-white-space-around-a-saved-image-in-matplotlib
      plt.gca().set_axis_off()
      plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
      plt.margins(0,0)
      plt.gca().xaxis.set_major_locator(NullLocator())
      plt.gca().yaxis.set_major_locator(NullLocator())
      plt.savefig(path.join(home,prod,name), dpi=dpi, pad_inches=0, frameon=False)
      plt.close()
    else:
      plt.draw()
      plt.show()
      plt.close()

  ###############
  ### IMAGING ###
  ###############
  def imaging(self, FOV, ccd):
    '''
       Image rasterization on the viewing instrument plane.
       
       -Input-
       =======
       FOV : spice.getfov object (spiceypy output)
       ccd : 2-tuple, ccd pixel size
       
       -Output-
       =======
       self:
          facet_image_
          facet_pix_
          solid_angle_ (per pixel size)
    '''
    from numpy import bool_, uint32, int32, float64, float32, arange, unique, ogrid, vstack, hstack, meshgrid
    import matplotlib.pyplot as plt
    assert hasattr(self, 'v1_v2_')

    v1_v2 = self.v1_v2_[self.facetid_-1,:,:]
    d = self.d_[self.facetid_-1,:].mean(axis=1).flatten().argsort(kind='mergesort')[::-1].copy(order='C') # Distance per facet, descending order

    self.ccd_ = ccd
    qua = FOV[4]

    if FOV[0] == "RECTANGLE":
      print('FOV window (rad)', (qua[:,1].max()-qua[:,1].min()), 
                                (qua[:,0].max()-qua[:,0].min()))
    
      self.pix_ = (
                   (qua[:,0].max()-qua[:,0].min())/ccd[0], 
                   (qua[:,1].max()-qua[:,1].min())/ccd[1]
                  )
      print('pixel angular size (rad)', self.pix_)
      print()
    
      grid_y = arange(qua[:,1].min(),qua[:,1].max(), self.pix_[1], dtype=float32) #[qua[0,0] + pix*j for j in range(ccd)]
      grid_x = arange(qua[:,0].min(),qua[:,0].max(), self.pix_[0], dtype=float32) #[qua[0,0] + pix*j for j in range(ccd)]

    elif FOV[0] == "CIRCLE":
      print('FOV window (rad)', 2*qua[0][0])
      print()
    
      self.pix_ = (
                   2*qua[0][0]/ccd[0], 
                   2*qua[0][0]/ccd[1]
                  )
      print('pixel angular size (rad)', self.pix_)
    
      grid_y = arange(-qua[0][0],qua[0][0], self.pix_[1], dtype=float32) #[qua[0,0] + pix*j for j in range(ccd)]
      grid_x = arange(-qua[0][0],qua[0][0], self.pix_[0], dtype=float32) #[qua[0,0] + pix*j for j in range(ccd)]

    # rasterization
    #input('imaging _raster enter')
    image, facet_pix = self._raster(v1_v2[d,...],
                                   self.facetid_[d],
                                   ccd,
                                   (grid_y, grid_x),
                                   mode=u'pixelation')
    facet_pix[:,0] = ccd[1] -facet_pix[:,0] -1
    
    if FOV[0] == "CIRCLE":
      a=ccd[1]/2
      b=ccd[0]/2
      y,x = ogrid[-a:ccd[1]-a,-b:ccd[0]-b]
      footprint = y*y/(a*a) + x*x/(b*b) >1e0
      image[footprint] = 0e0
      footprint, yx = image_to_table(footprint)
      yx = yx[footprint.astype(bool_)==False]
      facet_pix = facet_pix[inNd(facet_pix[...,(0,1)].copy(order='C'), 
                                yx.copy(order='C'), 
                                assume_unique=False, 
                                invert=False)]

    # output
    self.solid_angle_emi_ = solid_angle_tri(self.v1_v2_, 1e0)
    self.facet_image_ = image.astype(int32)
    # facet per pixel is in descending order in respect to Observer distance.
    self.facet_pix_ = unique(facet_pix, axis=0)  # sort: Y X facet
    self.facet_pix_.view('u4,u4,u4').sort(axis=0, kind='mergesort', order='f2') 

    # Line clipping : Get the facet ratios intersected by every pixel
    from .lineclipping import clipping_loop_c
    from time import time

    #print('#facets on the frame: ', self.facet_pix_[:,2].shape)
  
    start = time()
    px_area = clipping_loop_c(self.facet_pix_, 
                              self.v1_v2_[self.facet_pix_[:,2]-1,:,:].copy(order='C'),
                              self.pix_, 
                              grid_x.copy(order='C'), 
                              grid_y[::-1].copy(order='C')
                              )
    print('clipping computation time ', time() -start)

    # Remove empty pixels from facet-pixel array
    self.facet_area_pix_ = px_area[px_area[:,1]!=0.].astype(float64)
    self.facet_pix_ = self.facet_pix_[px_area[:,1]!=0.]
    
    print('#facets on the frame (empty removed): ', self.facet_pix_[:,2].shape)





  ###############################
  ### PIXEL ONTO OBJECT FRAME ###
  ###############################
  def onto_target_frame(self):
     '''
        Given an array of pixels (x, y) in radians in the camera frame, obtain their coordinates (X, Y, Z) in the target frame.
        
        Output
        ======
        XYZ_ : ndarray, images with pixel-wise cartesian coordinates
     '''
     from numpy import c_, uint32, float32, asarray, matmul, zeros
     from numpy.linalg import inv
     
     assert hasattr(self, 'd_')
     assert hasattr(self, 'facet_pix_')

     # Facet per pixel in camera variables (X_rad, Y_rad, d)
     fp_d = self.d_[self.facet_pix_[...,2]-1].mean(axis=1).flatten()
     fp_rad  = ((self.facet_pix_[...,(0,1)]+0.5).T*fp_d).T*asarray(self.pix_)[::-1]
     #fp_d = self.d_
     #fp_rad = self.v1_v2_*self.d_
     fp = c_[fp_rad, fp_d].T#.reshape(3,-1)
     #print(fp.reshape(3,-1).T)

     # Intrinsic & Extrinsic Matrix inversion
    
     # Intrinsic Matrix
     f =  self.boresight[2]
     y0 = self.boresight[1]
     x0 = self.boresight[0]
     K = array([[ f,   0,   x0],
                [ 0,   f,   y0],
                [ 0,   0,   1]], dtype=float32)
    
     I = inv(matmul(K,self.cam_matrix))
     
     XYZ = I.dot(fp).T -self.obs_pos
     #from mpl_toolkits.mplot3d import Axes3D
     #import matplotlib.pyplot as plt
     #fig = plt.figure()
     #ax = fig.add_subplot(111, projection='3d')
     #ax.scatter(XYZ[...,0], XYZ[...,1], XYZ[...,2], c='k', marker='.')
     #plt.show()
     
     # To image format:
     XYZ_ = zeros((3, self.ccd_[0], self.ccd_[1]), dtype=float32)
     yx = pd.MultiIndex.from_arrays(self.facet_pix_[...,(0,1)].T.astype(uint32), names=('y','x'))
     table = pd.DataFrame(XYZ, columns=('X','Y','Z'), index=yx, dtype=float32)
     for n, c in enumerate(table.columns):
       XYZ_[n,...] = table_to_image(table[c], self.ccd_[::-1], b=0e0)[::-1]
     #plt.imshow(XYZ_[...,0], interpolation='none')
     #plt.show()
     return XYZ_.astype(float32)



  ###########################################
  ### OVERLAY PROPERTY ON THE IMAGE PIXEL ###
  ###########################################
  def broadcast1(self, X, plot=False):
    '''
       Broadcast: Overlay property, version 1.
    
       Get facet-wise property (phase, emi, inc, band_depth, ...)
       and display it as an image into the instrument frame.

       This function is applied directly to the self.facet_image_
       and is to be considered a rough property approximation for pixels with several facets inside.

       Parameters
       ==========
       X : ndarray, facet-wise, property to be overlayed
    '''
    assert hasattr(self, 'facet_image_')
    from numpy import nonzero, isnan, arange, where, maximum

    # replace NaNs
    mask = isnan(X)
    idx = where(~mask,arange(mask.shape[0]),0)
    maximum.accumulate(idx,axis=0, out=idx)
    X[mask] = X[idx[mask]]

    # Ximage
    ximage = X[self.facet_image_-1]
    ximage[self.facet_image_==0] = 0

    # plot - test
    if plot:
      import matplotlib.pyplot as plt
      plt.gca().set_axis_off()
      plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
      plt.margins(0,0)
      plt.rcParams['figure.facecolor'] = 'white'
      plt.imshow(ximage, interpolation='none', cmap='Greys', alpha=0.99)
      plt.gca().invert_yaxis()
      plt.show()
      #plt.savefig(path.join(home,prod,'{}_property.png'.format(obj)), dpi=2*96, pad_inches=0, frameon=False)
      plt.clf()

    return ximage




  def broadcast2(self, X, plot=False, thresh=0.99):
    '''
       Broadcast: Overlay property, version 2.

       Get facet-wise property (phase, emi, inc, band_depth, ...)
       and display it as an image into the instrument frame.

       Using the precise facet-pixel ratios, an image is composed with higher precision.
       This function uses self.facet_pix_ and self.facet_area_pix_ .

       Parameters
       ==========
       X : ndarray, facet-wise, property to be overlayed

    '''
    from .imager_utils import weighted_average_c
    from numpy.core import defchararray
    from numpy import uint32, float32, unicode_, arange, stack, meshgrid, unique
    from time import time
    import matplotlib.pyplot as plt
    assert hasattr(self, 'facet_pix_')
    assert hasattr(self, 'ccd_')
    assert hasattr(self, 'facet_area_pix_')
    
    start = time()
    # Data Frame: Convert pixel position to string ndarray to preserve unicity
    df = pd.DataFrame(self.facet_pix_, columns=('y','x','f'))
    df['a'] = self.facet_area_pix_[:,1]/(self.pix_[0]*self.pix_[1])
    #print(df.query('a>{}'.format(thresh),inplace=False).head(50))
    df.sort_values(['y','x'], ascending=[True, True], inplace=True, kind='mergesort')
    #print(df.head(50))

    # Precise value of X per pixel
    X_pix = weighted_average_c(df.to_numpy(dtype=float32).copy(order='C'), 
                                           X[df['f']-1].astype(float32).copy(order='C')
                                           )
    print('broadcast computation time ',time() -start)
    # To multi-index DataFrame
    idx = pd.MultiIndex.from_arrays(X_pix[1:,(0,1)].astype(uint32).T, names=['y', 'x'])
    X_pix = pd.Series(X_pix[1:,2], index=idx)
    #print(X_pix)

    
    # Build image back from pixel position
    ximage = table_to_image(X_pix, self.ccd_[::-1], b=0e0)[::-1]
    
    
    # plot - test
    if plot:
      import matplotlib.pyplot as plt
      plt.gca().set_axis_off()
      plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
      plt.margins(0,0)
      plt.rcParams['figure.facecolor'] = 'white'
      plt.imshow(ximage, interpolation='none', cmap='Greys', alpha=0.99)
      plt.gca().invert_yaxis()
      plt.show()
      #plt.savefig(path.join(home,prod,'{}_property.png'.format(obj)), dpi=2*96, pad_inches=0, frameon=False)
      plt.clf()
   
    return ximage




  ##################
  ### IRRADIANCE ###
  ##################
  def irradiance(self, Sw):
     '''
       Irradiance at the instrument distance
       
       Input
       =====
       Sw : Source spectral irradiance at given wavelength at source position
     '''
     assert hasattr(self, 'pix_')
     self.solid_angle_inc_ =  self._irradiance(Sw, self.source_pos, self.obs_pos)/(self.pix_[0]*self.pix_[1])


# END
