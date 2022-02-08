# -*- coding: utf-8 -*-
#!/usr/bin/python

__version__="0.2"

doc = '''
  #########################################################################################
  #  author: Pedro H. Hasselmann, Rio de Janeiro (LESIA-OBSPM/ON-MCT)
  #  script: shapeviewer.py
  #
  # Shape Model 3D Visualization.
  #
  # Requirements:
  # numpy, scipy, astropy, pandas, mayavi, TVTK, skimage, matplotlib
  #
  ##########################################################################################
'''
from .. import *
from ..main import *

########################
## 3D GRAPHICAL TOOLS ## 
########################
class Viewer:
  '''
       3D Shape Viewer
       
       Parameters
       ==========
       S : loaded Shape Model Class
  '''
  from os import name, environ
# First, and before importing any Enthought packages, set the ETS_TOOLKIT
# environment variable to qt4, to tell Traits that we will use Qt.
  print(name)
  if name != 'nt':
    environ['ETS_TOOLKIT'] = 'qt4'
# By default, the PySide binding will be used. If you want the PyQt bindings
# to be used, you need to set the QT_API environment variable to 'pyqt'
  if name == 'nt':
    environ['QT_API'] = 'pyside'
  
  def __init__(self, 
               S, 
               size=(600,600)
               ): # Import loaded shape model class
    from tvtk.api import tvtk
    from mayavi.mlab import figure
    self.S=S
    self.fig = figure(size=size)
    self.mesh= tvtk.PolyData(points=S.vertices_, polys=S.facets_) 

  def cell_data(self,params,label):
    '''
        Load values to be printed into the shape model.
        
        Parameters
        ==========
        params : numpy.array or pandas.Series
        label  : str
    '''
    if isinstance(params, pd.Series):
      params = params.reindex(range(1,self.S.facets_.shape[0]), fill_value=-999, copy=True).values
      
    self.mesh.cell_data.scalars= params
    self.mesh.cell_data.scalars.name= label  

  def bodyrepr(self,repres='surface',scale=(None,None),colormap=None,nlabel=5, alpha=0.99):
    from mayavi.mlab import pipeline, colorbar, scalarbar
    mesh=pipeline.surface(self.mesh, vmax=scale[1], vmin=scale[0], colormap=colormap, representation=repres, opacity=alpha)
    mesh.module_manager.scalar_lut_manager.reverse_lut = True
    #scalarbar(title=self.mesh.cell_data.scalars.name, nb_labels=nlabel)
    cb = colorbar(title=self.mesh.cell_data.scalars.name, nb_labels=nlabel)
    cb.label_text_property.font_family = 'times'
    cb.label_text_property.bold = 0
    cb.label_text_property.font_size=16

  def vector(self,v,o,**args):
    from mayavi.mlab import quiver3d
    from numpy import vstack, hstack
    v=list(hstack((o,v)).T) 
    quiver3d(*v,**args)

  def latlon(self, dlat, dlon):
    from mayavi.mlab import text3d
    from scipy.interpolate import griddata
    from numpy import float32, array, radians, degrees, sin, cos, meshgrid, arange, vstack, array2string
    
    lat, lon, r = array(spherical_coord(self.S.vertices_[:,0], self.S.vertices_[:,1], self.S.vertices_[:,2]), dtype=float32)
    #print(lat.max(), lat.min())
    #print(lon.max(), lon.min())

    # Grid -- Correct
    glat = radians(arange(-90,90,dlat))#[1:]
    glon = radians(arange(-180,180,dlon))#[1:]
    glat, glon = meshgrid(glat, glon)
    gr = griddata(array((lat,lon)).T,r,(glat,glon), method='nearest')+0.030
    coord_txt = vstack((degrees(glat),degrees(glon))).reshape(2, -1).T.astype(int).astype(str)
    coord_txt = [', '.join(row) for row in coord_txt]

    # Spherical to Cartesian -- Correct
    x = gr*cos(glat)*cos(glon)
    y = gr*cos(glat)*sin(glon)
    z = gr*sin(glat)
    # Normal Vectors -- Correct
    cart = vstack((x,y,z)).reshape(3, -1)

    print(cart.shape)
    for j in range(cart.shape[1]):
      text3d(cart[0,j], cart[1,j], cart[2,j], coord_txt[j], color=(0,0,0), scale=0.01, orient_to_camera=True)
    
    

  def camera(self,az,evel,d,roll=None):
    from numpy import array
    from mayavi.mlab import view, orientation_axes
    orientation_axes(line_width=6.0)
    self.v = view(az,90-evel,d,roll=roll,focalpoint='auto')

  def view(self):
     from mayavi.mlab import draw, show, orientation_axes
     orientation_axes(line_width=6.0)
     draw()
     show()
  
  def clear(self):
     from mayavi.mlab import clf, close
     clf()

  def close(self):
     from mayavi.mlab import clf, close
     close(all=True)
        
  def screenshot(self,figname,label):
     import wx
     from mayavi.mlab import savefig, draw, screenshot, title
     if name != 'nt': 
       locale = wx.Locale.GetSystemLanguage()
            
     title(label,size=0.8,height=0.92,line_width=20.0, color=(0.9,0,0))
     savefig(path.join(home,prod,figname+'.png'),size=(600,600))
     image = screenshot()
    
     if name != 'nt':
       app = wx.App(redirect=False)
       app.locale = wx.Locale(locale)
     return image

# END
