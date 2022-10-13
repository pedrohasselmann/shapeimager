# -*- coding: utf-8 -*-
#!/usr/bin/python

__version__="0.2"

doc = '''
  #########################################################################################
  #  author: Pedro H. Hasselmann, Rio de Janeiro (LESIA-OBSPM/ON-MCT)
  #  script: shapemodel.py
  #
  # Shape Model Manipulation
  #
  # Load .OBJ, .PLY, .VER shape model.
  #
  # Requirements:
  # numpy, scipy, astropy, pandas, skimage, matplotlib
  #
  ##########################################################################################
'''
from .. import *
from ..main import *

###############################################################################
############################ SHAPE MODEL TOOL #################################
class ShapeModel(object):

  def __init__(self,shapefile, comments=0, pickle=True):
    import itertools
    from glob import glob
    from numpy import float32, int32, uint32, fabs, array, split, savez, load, unique, where, delete
    
    npz = path.join(home,aux,shapefile.split('.')[0]+'.npz')
    shapefile = path.join(home,aux,shapefile)
    print(shapefile)
    self.shapefile_ = shapefile
    print(npz)

    if len(glob(npz))==1:
      with load(npz, allow_pickle=True) as f:
          self.vertices_ = f["vertices"].copy(order='C')
          self.facets_   = f["facets"].copy(order='C')
          self.array_    = f["array3D"]
          self.facet_index_ = f["facet_index"].astype(uint32).copy(order='C')

    else:
      if shapefile[-3:].lower() == 'obj':
        self.load_obj(shapefile, comments)
    
      if shapefile[-3:].lower() == 'ply':
        self.load_ply(shapefile)
    
      if shapefile[-5:].lower() == 'shape':
        self.load_sha(shapefile)

      if shapefile[-3:].lower() == 'ver':
        self.load_ver(shapefile)

    # LICIACube - Cone re-orientation
    #self.vertices_[:,0] = -self.vertices_[:,0]
    #self.vertices_[:,1] = -self.vertices_[:,1]


    if len(glob(npz))==0:
      # Final database format : 3D-array packing all triangle vertices
      vertices = pd.DataFrame(self.vertices_, index=range(self.vertices_.shape[0]), columns=['X','Y','Z'])
      df = vertices.loc[self.facets_.ravel()]
      df['facet_index'] = list(itertools.chain.from_iterable(itertools.repeat(x, 3) for x in range(self.facets_.shape[0])))
      df.set_index('facet_index', inplace=True)
      #print(df.loc[0])
      self.facet_index_ = unique(df.index).astype(uint32)+1
      self.array_ = array(split(df.values, len(self.facet_index_)), order='C', dtype=float32)
    
      if pickle == True:
        savez(npz,
          vertices=self.vertices_.copy(order='C'),
          facets=self.facets_.copy(order='C'),
          array3D=self.array_,
          facet_index=self.facet_index_)

    print('shape model is loaded......',self.vertices_.shape, self.facets_.shape)


################################
## LOAD SHAPE MODEL FUNCTIONS ##
################################

  def load_obj(self, shapefile, comments=11):
    ''' 
       Load Shape Model OBJ.

       Uses pandas.read_csv
    '''
    from pandas import read_csv
    from numpy import float32, uint32, split
  
    #{'names':('X','Y','Z'), 'formats':('f4','f4','f4')}
    self.vertices_ = read_csv(shapefile, dtype=float32, 
                      comment='f', usecols=(1,2,3), skiprows=comments, header=None, delim_whitespace=1).values.copy(order='C')
    #{'names':('V1','V2','V3'), 'formats':('i4','i4','i4')}
    self.facets_ = read_csv(shapefile, dtype=uint32,
                    comment='v', usecols=(1,2,3), skiprows=comments, header=None, delim_whitespace=1).values.copy(order='C') -1


  def load_ply(self, shapefile):
    '''
      Load Shape Model PLY.

      Uses plyfile package.
    '''
    from plyfile import PlyData, PlyElement
    from numpy import  array, asarray, matrix, float32, uint32, split

    mesh = PlyData.read(open(shapefile, 'rb'))
    print(shapefile)
    print(mesh)

    self.vertices_ = array(map(lambda x: list(x), mesh['vertex'].data), order='C', dtype=float32)
    self.facets_   = array(map(lambda x: list(x[0]), mesh['face'].data), order='C', dtype=uint32) #-1

  def load_sha(self,shapefile):
    '''
       Load Shape Model .shape

       Uses pandas.read_csv
    '''
    from pandas import read_csv
    from numpy import float32, uint32, split

    f = open(shapefile, 'r+')
    delimiter = int(f.readline().split(' ')[0]) -1
    mesh = read_csv(f, dtype=float32, skiprows=1, usecols=(0,1,2)).values
    
    self.vertices_ = mesh[:delimiter,:]
    self.facets_ = mesh[delimiter:,:].astype(uint32)

  def load_ver(self,shapefile):
    '''
       Load Shape Model .ver

       Uses collections.deque
    '''
    from collections import deque
    from numpy import uint32, float32, array
    
    with open(shapefile, 'r+') as f:
       N = f.readline().split("   ")
       N.remove(u'')
       N_vertices = int(N[0])
       N_facets = int(N[1])
       print(N[0], N[1])
       # Read Vertices
       V = deque()
       l=0
       while l<N_vertices:
           V.append(f.readline().split(" "))
           l+=1
    
       self.vertices_ = array(V, dtype=float32)
    
       # Read facets
       F = deque()
       k,l=0,1
       while k<N_facets:
        if l%2==0:  
          line = f.readline().split(" ")
          line=filter(lambda a: a!='', line)
          F.append(line)
          k+=1
        else:
          f.readline()
        l+=1

       self.facets_ = array(F, order='C', dtype=uint32) -1

######################
## ARRAY FORMATTING ##
######################
  def rebase(self, v):
    '''
       DEPRECATED -- Use imager_utils.rebase (cythonized).
       Re-format projected vertices to facet-oriented index.
       v : coordinates per vertice - (M, N)

       Final database format : 3D-array packing all triangle vertices
    '''
    import itertools
    from numpy import float64, float32, uint32, array, split, empty

    vertices = pd.DataFrame(v, index=range(self.vertices_.shape[0]))
    df = vertices.loc[self.facets_.ravel()]
    df['facet_index'] = list(itertools.chain.from_iterable(itertools.repeat(x, 3) for x in range(self.facets_.shape[0])))
    df.set_index('facet_index', inplace=True)
    #print(df.loc[0])
    return array(split(df.values, self.facet_index_.shape[0]), dtype=float64, order='C')
      

##################
## SAVE FORMATS ##
##################

  def save_obj(self):
    '''
       Save as Shape Model .OBJ

       Uses numpy only.
    '''
    from numpy import around, savetxt

    obj = open(path.join(aux,self.shapefile_[:-4]+'.obj'),'w+')

    obj.write(''.join(['##############################################\n',
            '# 3D Shape model                               \n',
            '# --------------                               \n',
            '# This mesh has been reformatted by the \n',
            '# ShapeModel.save_obj (P. H. Hasselmann, LESIA)   \n',
            '#                                             \n',
            '#                                             \n',
            '# Vertices: '+str(self.vertices_.shape[1])+' \n',
            '# Faces: '+str(self.facets_.shape[0])+'  \n',
            '#                                              \n',
            '##############################################\n']))

    #print(facets)
    savetxt(obj, around(self.vertices_, 8), fmt=str('v %2.8e %2.8e %2.8e'))
    savetxt(obj, self.facets_, fmt=str('f %8i %8i %8i'))
  
    obj.close()
    print('.obj ....OK.')

###############
## SELECTION ##
###############
  def select(self, f):
    '''
      Select facets by their index.
    '''
    self.facet_index = tuple(f)
    self.array  = self.array_[f,...]
    self.facets = self.facets_[f,...]
    print('new shape ',self.array_.shape, self.facets_.shape)

#################
## COORDINATES ##
#################

  def spherical_coords(self, pickle=False):
    '''
       Cartesian to Spherical coordinates.
    '''
    import itertools
    from numpy import float32, int32, split
        
    # Spherical coordinates
    v = self.vertices_
    lat, lon, r = spherical_coord(v[:,0], v[:,1], v[:,2])

    coords = pd.DataFrame({'lat':lat.T, 'lon':lon.T}, index=range(lat.size))
    df = coords.loc[self.facets_.ravel()]
    df['facet_index'] = list(itertools.chain.from_iterable(itertools.repeat(x, 3) for x in range(self.facets_.shape[0])))
    df.set_index('facet_index', inplace=True)
    coords = array(split(df.values, self.facet_index_.shape[0]),dtype=float32)

    self.lat_ = coords[:,:,0]
    self.lon_ = coords[:,:,1] 
    
    if pickle == True:
      self.lat_.to_pickle(self.shapefile.split('.')[-1]+'_lat.pkl')
      self.lon_.to_pickle(self.shapefile.split('.')[-1]+'_lon.pkl')

##################
## NORMAL VECTOR ##
###################

  def normal_vector(self, pickle=False):
    '''
       Compute facet centroid and normal vector.

       Save calculations into a pickle file for speed up.
    '''
    from numpy import array, isnan, transpose, sqrt, cross, arccos, clip, sum,  gradient, float32, int32, split, load

    shapepath = self.shapefile_.split('.')[0]

    if path.exists(shapepath+'_normalvector.pkl'):
      self.n_ = load(shapepath+'_normalvector.pkl', allow_pickle=True)
      self.tilt_ = load(shapepath+'_tilt.pkl', allow_pickle=True)
      return

    # Normal Vector in respect to the center of the facet
    der = array(gradient(self.array_, axis=1)) # Derivatives
    # Colunms X Y Z
    n = cross(der[:,0,:],der[:,2,:]) # Normal Vector
    
    # Tilt between normal vector and a vertice
    v = self.array_[:,0,:] # Vertices 1
    v = transpose(v.T/sqrt(sum(v*v, axis=1)))
    n = transpose(n.T/sqrt(sum(n*n, axis=1)))
    tilt = arccos(clip(sum(n*v, axis=1), -1e0, 1e0)) # Tilt
    tilt[tilt >= 3.1415926/2e0] = 3.1415926 -tilt[tilt >= 3.1415926/2e0]
    tilt[isnan(tilt)] = 0e0
      
    # Save into self
    self.n_ = n
    self.tilt_ = tilt

    if pickle == True:
      self.n_.dump(shapepath+'_normalvector.pkl')
      self.tilt_.dump(shapepath+'_tilt.pkl')


  def normal_vector_ellipsoid(self, axes, pickle=False):
    '''
      Compute the normal vector by replacing the local surface element tilt (facet)
      for a smooth tangent plane represented by a equivalent ellipsoid of same radius.
      
      Tangent Plane at point (x0, y0, z0):
      
      x*x0/a**2 + y*y0/b**2 + z*z0/c**2 = 1
      
      Normal vector to such plane:
      
      n = [x0/a**2, y0/b**2, z0/c**2]
    '''
    from numpy import array, isnan, transpose, sqrt, cross, arccos, clip, sum,  gradient, float32, int32, split, load

    shapepath = self.shapefile_.split('.')[0]

    if path.exists(shapepath+'_normalvectorellip.pkl'):
      self.n_ = load(shapepath+'_normalvectorellip.pkl', allow_pickle=True)
      return

    n = self.array_/axes**2
    self.n_ = transpose(n.T/sqrt(sum(n*n, axis=1)))

    if pickle == True:
      self.n_.dump(shapepath+'_normalvectorellip.pkl')


##########################
## GEOMETRIC PROPERTIES ##
##########################

  def project(self, v):
    '''
       Compute the direction vector and angle of a shape model vertices towards a given source.
       The vector-direction must be referenced to the shape model cartesian origin.
       Thus, incidence and emergence angles may be now calculated with precision.
       
       Angles in radians.
    '''
    from numpy import array, isnan, transpose, sqrt, cross, arccos, clip, sum, float32, float64
    assert hasattr(self, 'n_')
    
    n = self.n_
    d = self.array_[:,0,:]
    v = array(v, dtype=float64) +d
    v = transpose(v.T/sqrt(sum(v*v, axis=1)))
    #n = transpose(n.T/sqrt(sum(n*n, axis=1)))
    
    v_n = sum(n*v, axis=1).astype(float32)
    angle = arccos(clip(v_n, -1e0, 1e0)).astype(float32) # Tilt
    #angle[angle >= pi/2e0] = pi - angle[angle >= pi/2e0]
    angle[isnan(angle)] = -3.1415926/2e0    
    return angle, v_n


  def phase(self, v, o):
    '''
       Compute phase angle between two re-based vectors 'v' and 'o'.
       
       Angles are given in radians.
    '''
    from numpy import array, isnan, transpose, degrees, sqrt, cross, arccos, clip, sum, \
    float32, float64
    assert hasattr(self, 'array_')

    phase_at_center = arccos(sum(transpose(v.T/sqrt(sum(v*v)))*transpose(o.T/sqrt(sum(o*o)))
                                                            )
                                                        )
    print('phase angle at body center is {: 3.3f}'.format(degrees(phase_at_center)))
    
    d = self.array_[:,0,:]
    v = array(v, dtype=float64) +d
    o = array(o, dtype=float64) +d
    
    v = transpose(v.T/sqrt(sum(v*v, axis=1)))
    o = transpose(o.T/sqrt(sum(o*o, axis=1)))
 
    v_o = sum(o*v, axis=1)
    angle = arccos(clip(v_o, -1e0, 1e0)) # Tilt
    #angle[angle >= pi/2e0] = pi - angle[angle >= pi/2e0]
    angle[isnan(angle)] = -3.1415926/2e0    
    return angle, v_o



  def distance(self, o):
     from numpy import float32, sqrt, sum
     return sqrt(sum((self.array_-o)**2, axis=2))



  def plane(self, o, n, e1=[1,0,0], e2=[0,1,0]):
    '''
       Reproject all shapemodel vertices onto a given plane.
      
       Projected vertice satisfies the relation:
       P = o + v1*e1 + v2*e2 + s*n

       s = Dot(n, v-O)
       v1 = Dot(e1, v-O)    
       v2 = Dot(e2, v-O)
       
       Parameters
       ==========
       n : vector that defines the plane
       o : origin of the plane
       
       v1, v2 : (X, Y), Cartesian coordinates of the shape model projected onto the (e1,e2) plane
       https://stackoverflow.com/questions/33658620/generating-two-orthogonal-vectors-that-are-orthogonal-to-a-particular-direction

       REBASE optimized.
    '''
    from numpy import float64, float32, asarray, transpose, cross, sum, vstack
    from math import sqrt
    from .imager_utils import rebase_c
    from time import time
    
    n = asarray(n)
    o = asarray(o)

    start = time()
    # Define the base of the plane orthogonal to n (e1,e2)
    e1 = asarray(e1, dtype=float64)  # take an arbitrary vector
    #e2 = array(e2, dtype=float32)  # take an arbitrary vector

    # Normalize vectors
    n = transpose(n.T/sqrt(sum(n*n, axis=0)))
    #o = transpose(o.T/sqrt(sum(o*o, axis=0)))
    
    e1 -= e1.dot(n)*n                   # make it orthogonal to n
    e1 /= sqrt(sum(e1*e1, axis=0))      # normalize it
    e2 = cross(n, e1)

    #print('orthogonality e1',e1.dot(n))
    #print('orthogonality e2',e2.dot(e1), e2.dot(n))

    # Project onto Plane
    d = self.vertices_ -o  # Columns X Y Z (Vertices/rows 1 2 3)
    s=sum(-d*n, axis=1)
    v1=sum(d*e1, axis=1)
    v2=sum(d*e2, axis=1)
    print('shapemodel.plane calculation time', time() - start)
    
    # Re-format the plane-projected shape model to facet-oriented index.
    start = time()
    v1_v2 = rebase_c(vstack([v1, v2]).T.copy(order='C'), self.facets_)
    sr = rebase_c(s.reshape(-1, 1), self.facets_)

    print('shapemodel.rebase overhead time', time() - start)

    return v1_v2, sr, (-n, e1, e2)

  def camera_framing(self, R, n, o):
    '''
       Re-project the 3D shape model onto the camera perspective.
    
       http://ksimek.github.io/2012/08/22/extrinsic/
       
       Parameters
       ==========
       R : Camera Rotation matrix (numpy.array, 3x3)
       n : Boresight vector (numpy.array, 1x3)
       o : camera origin vector (numpy.array, 1x3)
    '''
    from numpy import int32, float32, transpose, sqrt, median, argmax, fabs, sum
    from .imager_utils import rebase_c
    from math import fabs

    # boresight
    #f = sqrt(sum(n*n, axis=0)) #focal length
    #nn = transpose(n.T/f)
    
    # Intrinsic Matrix
    f =  n[2]
    y0 = n[1]
    x0 = n[0]
    K = array([[ f,   0,   x0],
               [ 0,   f,   y0],
               [ 0,   0,   1]], dtype=float32)

    # Extrinsic Matrix
    d = (self.vertices_ -o).copy(order='C')
    X = K.dot(R.dot(d.T))
    
    # Re-format the plane-projected shape model to facet-oriented index.
    i=(1,0)
    j=2
    v1_v2 = rebase_c((X[i,:]/X[j,:]).T.copy(order='C'), self.facets_)
    d = rebase_c(X[j,:].reshape(-1, 1).copy(order='C'), self.facets_)
    
    return v1_v2, d


  def irradiance(self, Sw, s, o):
     '''
       Irradiance at the instrument distance
       
       Input
       =====
       Sw : Source spectral irradiance at given wavelength at source position
     '''
     from numpy import float32, transpose, sum, sqrt
     source_at_obs = sqrt(sum(o**2))*transpose(s.T/sqrt(sum(s**2)))
  
     # SHAPE MODEL
     # Area projection & solid angle
     v1_v2, d, base = self.plane(source_at_obs, s)
     solid_angle_inc_ =  solid_angle_tri(v1_v2/d, 1e0)
     return float32(Sw)*solid_angle_inc_



################
### PATCHING ###
################
  def patches(self, n, semiaxis, seed=121189, plot=False, save=True):
        '''
        '''
        from scipy import stats
        from scipy.spatial import transform
        from numpy import random, linspace, arange, ones, dstack, newaxis, \
        fabs, sqrt, argmin, diag, cross, zeros, log, radians, degrees, \
        arccos, cos, sin, savez, isnan, isinf

        # load shape model
        S = load_shapemodel(shafile)
        vecs = self.array_.mean(1)

        # generate random points in the space : https://tutorial.math.lamar.edu/Classes/CalcIII/EqnsOfLines.aspx
        rnd = random.RandomState(seed)
        v = 2*rnd.rand(n,3) -1.
        print(v)

        # trace line connecting the points to the shape model origin
        t = linspace(0,2,100)
        x = v[:,0]*t[:,None]
        y = v[:,1]*t[:,None]
        z = v[:,2]*t[:,None]
        lines = dstack((x,y,z))
        print(lines.shape)

        # from traced line, selected the closest facets
        samp = list()
        probs = zeros(S.array_.shape[0])
        z = array((0.,0.,1.))
        semiaxis = diag(semiaxis)
        for j in range(n):
                l = lines[:,j,:]
                subt  = sqrt((vecs[:, newaxis] -l)**2).sum(2).min(1)
                print(subt.shape)
                idx = argmin(subt)
                print(idx)
                vc = vecs[idx,:]
                nv = self.n_[idx,:]
                print('vc', vc, v[j,:])
                print('nv', nv)

                # from distribution around points, project to shape model
                # build covariance: 
                #https://math.stackexchange.com/questions/1956699/getting-a-transformation-matrix-from-a-normal-vector
                #https://robotics.stackexchange.com/questions/2556/how-to-rotate-covariance
                rotaxis = cross(nv, z)
                rotrad = arccos(nv.dot(z))
                print('rotaxis ',rotaxis)
                print('rotrad ', degrees(rotrad),' deg')

                quat = cos(0.5*rotrad)*array([1., 0., 0., 0.])
                quat[1:] = rotaxis[::-1]*sin(0.5*rotrad+1e-8)
                R = transform.Rotation.from_quat(quat)
                print('R ', R.as_matrix())
                cov = R.as_matrix() @ semiaxis @ R.as_matrix().T
                print('cov ',cov)

                # Probility distribution of the patch
                distr = stats.multivariate_normal(vc, cov, random_state=seed, allow_singular=True)
                #print(distr.pdf(vc))
                prob = distr.pdf(vecs)
                #print('prob j = {} : '.format(j), prob)
                #print('probs max min : ', prob.max(), prob.min())
                samp.extend(distr.rvs(size=100))

                probs = probs + fabs(prob)
                # end loop

        prob_weight = fabs(1e0/log(probs))
        prob_weight[isinf(prob_weight)] = 0.0

        samp = array(samp)
        print('sampling : ', samp.shape)
        print('prob total : ', prob_weight)
        print('probs max mean min : ', prob_weight.max(), prob_weight.mean(), probs.min())

        if plot == True:
              import plotly.graph_objects as go
              import trimesh

              fig = go.Figure(data=[go.Scatter3d(
              x=samp[:,0],
              y=samp[:,1],
              z=samp[:,2],
              mode='markers',
              marker=dict(
              size=8,
              colorscale='Viridis',   # choose a colorscale
              opacity=0.8
              )
              )])
              fig.add_trace(go.Mesh3d(
                        x=S.vertices_[:, 0], 
                        y=S.vertices_[:, 1], 
                        z=S.vertices_[:, 2],
                        i=S.facets_[:, 0],
                        j=S.facets_[:, 1],
                        k=S.facets_[:, 2],
                        intensity=fabs(1e0/probs), 
                        colorscale='Spectral', 
                        opacity=.5,
                        flatshading=True,
                        #lighting=dict(specular=0.2),
                        alphahull=1))
              fig.show()


        if save == True: savez('patches_gaussians_n{}_{}'.format(n, shafile[:-4]), values=prob_weight)
        return prob_weight



####################
### PARTITIONING ###
####################
  def partitioning(self, dlat=5, dlon=5, dxy=0.020, dz=0.100, pole=(-90,90), save=False, mode='yield'):
    '''
       Partitioning the shape model in equal block.
       The blocks may intercept and have variable local surface sizes.
       https://math.stackexchange.com/questions/225614/tangent-plane-to-sphere
       https://stackoverflow.com/questions/37670658/python-dot-product-of-each-vector-in-two-lists-of-vectors
        
        P4____________P5
     P0 /|           /|
       x-|----------x P1
       | |          | |
       | |          | |
       | |   x      | |
       | |____N ____|_|
       |/ P6        |/ P7
       x------------x
     P2              P3

       Tangent Plane: 2(x-xp) + 2(y-yp) + 2(z-zp) = 0
       Normal vector (x, y, z)
       https://math.stackexchange.com/questions/2679926/unit-normal-of-sphere-in-cartesian-coordinates
       http://citadel.sjfc.edu/faculty/kgreen/vector/Block3/flux/node6.html


       Two modes:
       yield : return generator. Output: (corners, box_content)
       return : charge (base, cartesian, corners, box_content) into self.patches

    '''
    from numpy import int32, float32, isnan, arange, meshgrid, where, sin, cos, array, logical_and, \
    sqrt, cross, transpose, sum, zeros, ones, einsum, radians, stack, vstack, fabs, dot, all, savez
    from numpy import linalg
    from scipy.interpolate import griddata
    from collections import deque
    #print(self.vertices_.max(0))
    #print(self.vertices_.min(0))

    #cart_array = self.vertices_
    lat, lon, r = array(spherical_coord(self.vertices_[:,0], self.vertices_[:,1], self.vertices_[:,2]), dtype=float32)
    #print(lat.max(), lat.min())
    #print(lon.max(), lon.min())

    # Grid -- Correct
    glat = radians(arange(pole[0],pole[1],dlat))#[1:]
    glon = radians(arange(-180,180,dlon))#[1:]
    glat, glon = meshgrid(glat, glon)
    gr = griddata(array((lat,lon)).T,r,(glat,glon), method='nearest')

    # Spherical to Cartesian -- Correct
    x = gr*cos(glat)*cos(glon)
    y = gr*cos(glat)*sin(glon)
    z = gr*sin(glat)
    #print(x.max(), x.min())
    #print(y.max(), y.min())
    #print(z.max(), z.min())
    # Normal Vectors -- Correct
    xyz = vstack((x,y,z)).reshape(3, -1)
    n = xyz/sqrt(sum(xyz*xyz, axis=0))

    # Orthogonal base -- Correct
    dim = n.shape[1:]
    e1 = vstack((zeros(dim), zeros(dim), ones(dim))) # vector parallel to longitude line
    e1 -= einsum('ij, ij->i', e1.T, n.T)*n   # make it orthogonal to n (einsum: multidimensional dot operation)
    e1 /= sqrt(sum(e1*e1, axis=0))            # normalize it
    e2 = cross(n.T, e1.T).T                   # vector parallel to latitude line
    e2 /= sqrt(sum(e2*e2, axis=0))            # normalized it

    # Continue....
    # Partitioning shape model vertices into boxes
    #Box ---  dX dY dZ
    corners  = array((( dxy, dxy, dz),
                      (-dxy, dxy, dz),
                      ( dxy,-dxy, dz),
                      (-dxy,-dxy, dz),
                      ( dxy, dxy, -dz),
                      (-dxy, dxy, -dz),
                      ( dxy,-dxy, -dz),
                      (-dxy,-dxy, -dz)), dtype=float32)

    # Rotation Matrix from Global frame to "surface frame"
    B0 = array(((1,0,0),
                (0,1,0),
                (0,0,1)), dtype=float32)
    #print(corners)
    corners_global = deque()
    matrix2local = deque()
    for j in range(n.shape[1]):
       B1 = vstack((e2[:,j],e1[:,j],n[:,j]))
       # Global --> Local
       M1 = linalg.solve(B1,B0)
       # Local --> Global
       M2 = linalg.solve(B0,B1)

       corners_global.append(dot(M1, corners.T).T + xyz[:,j])
       matrix2local.append(M2)

    corners_global = array(corners_global, dtype=float32)
    matrix2local = array(matrix2local, dtype=float32)

    # Binning the shape model:
    N_boxes = matrix2local.shape[0]
    print(N_boxes, ' patches')
    patch_content = deque()
    #del self.vertices_, self.facets_
    for m in range(N_boxes):
      
      M = matrix2local[m,...]

      v = self.array_ -xyz[:,m]
      p = dot(M.T, v.T).T

      boxyz = all((p[...,0]<dxy)&(p[...,0]>-dxy)&(p[...,1]<dxy)&(p[...,1]>-dxy)&(p[...,2]<dz)&(p[...,2]>-dz),axis=1)

      if mode == 'return':
        patch_content.append(self.facet_index_[boxyz])

      yield M, self.facet_index_[boxyz]

    if (save == True)&(mode=='return'):
      savez(self.shapefile_.split('.')[0]+'_patch{}_{}_{}_{}.npz'.format(dlon,dlat,dxy,dz), corners=corners_global, patches=patch_content)
      
    if mode == 'return':
      self.patches = ((e1, e2, n), xyz, corners_global, patch_content)
    


# END
