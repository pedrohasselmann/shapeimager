# -*- coding: utf-8 -*-
#!/usr/bin/python

# Escrevendo em python3 e usando python2.6:
from __future__ import print_function, unicode_literals, absolute_import, division

__version__="1.0.0"

__doc__='''
  ############################ HAPKE PLOTTING TOOL ################################
  ### Author: Pedro Henrique A. Hasselmann
  ### Last Modified: March, 2016
  ###
  ### Tool for displaying astronomical images of resolved bodies
  ###
  ### requirements: matplotlib, matplotlib-scalebar
  ###
  #################################################################################
'''
######################## GLOBAL IMPORT ##########################################

from source import *
from numpy import degrees, radians, tan

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

##
## BEGIN
##
from ..output.plot import Plot, font
font['size']=20

class ImagePlot(Plot):

   __doc__ = ''' Python Class dealing with 2D image display. 
                 Matplotlib interface.                  
              '''

   Plot.matplotlib.rcParams.update({'font.size': font['size']})
   imgplot = []

   def get_image(self, imagename, channel, treat=False, **args):
      '''
         Load .PDS & .FITS formats.
      '''
      if not treat:
         hdu = read_image(imagename, channel, **args)
         self.image = hdu.data #, hdu.header
      else:
         self.image = image_treated(imagename, channel, **args)

   def get_array(self, array):
      '''
        Load from numpy.array
      '''
      self.image = array # (Y, X)
            
   def get_figure(self, imagename, s=None):
      '''
         Load .tiff, .jpg & .png format
         
         imagename   --> figure
         imageheader --> image containing header
         s           --> cropping
      '''
      print(imagename)
      self.image = plt.imread(imagename)
      if s:
        self.image = self.image[2048-s[1]:2048-s[0], s[2]:s[3]] # (Y, X)
      #self.header = read_image(imageheader, 1).header

   def imshow(self, n, scale=['quar',0,100], factor=1, **args):
      from numpy import isnan, isinf, nanpercentile, dstack

      if scale[0] == 'quar':
        inf = factor*nanpercentile(self.image, scale[1])
        sup = factor*nanpercentile(self.image, scale[2])
        
      elif scale[0] == 'abs':
        inf = scale[1]
        sup = scale[2]
        
      self.current = self.frame[n].imshow(factor*self.image, **args)
      self.imgplot.append(self.current)
      self.imgplot[-1].set_clim(inf,sup)

   def threshold(self, value, boolean):
     from numpy import ma, greater_equal, less_equal, equal
     if boolean == 'ge': boo = greater_equal
     if boolean == 'le': boo = less_equal
     if boolean == 'eq': boo = equal
     self.image = ma.masked_where(boo(self.image, value), self.image)

   def interpolate(self, sigma=(2,2)):
      import scipy.ndimage as ndimage
      self.image = ndimage.gaussian_filter(self.image, sigma=sigma, order=0)

   def framework(self, n, colorbar_title, alt, pix, ticks=False, color='white', format='%.2f', side="right"):
      '''
         Build image framework: scalebar, colorbar, ticks, dimension,...
      '''
      from mpl_toolkits.axes_grid1 import make_axes_locatable                        
      from matplotlib_scalebar.scalebar import ScaleBar
      
      #pix = float(self.header['ROSETTA:VERTICAL_RESOLUTION'].split(' ')[-2])
      #alt = float(self.header['SPACECRAFT_ALTITUDE'].split(' ')[-2])
      
      if alt!=None and pix!=None:
        pix_scale = tan(pix)*alt
      
        bar = ScaleBar(pix_scale, frameon=False, color=color, pad=0.8, location='upper right',
                     length_fraction=0.125, height_fraction=0.02,
                     font_properties={'size':36})
      
        self.frame[n].add_artist(bar)
      
      if colorbar_title!=None:
        # colorbar
        # Create divider for existing axes instance
        divider = make_axes_locatable(self.frame[n])
        cax = divider.append_axes(side, size="3%", pad=0.12)

        # Create colorbar in the appended axes
        cbar = plt.colorbar(self.current, cax=cax, format=format)#, extend='min')
        
        tick_locator = ticker.MaxNLocator(nbins=8)
        cbar.locator = tick_locator
        
        side_args={}
        if side == "left": 
            side_args["labelleft"]=True
            side_args["left"]=True 
            side_args["labelright"]=False
            side_args["right"]=False
            labelpad=-68
        else:
            side_args["labelleft"]=False
            side_args["left"]=False
            side_args["labelright"]=True
            side_args["right"]=True
            labelpad=26
            
        cbar.ax.tick_params(labelsize=int(0.8*font['size']), **side_args)
        cbar.update_ticks()
      
        cbar.set_label(colorbar_title, 
                         fontsize=font['size'], 
                         rotation=270, 
                         labelpad=labelpad)

      if ticks==False:
        self.frame[n].set_xticklabels([])
        self.frame[n].set_yticklabels([])

      # Ticks
      #loc = Plot.ticker.MultipleLocator(base=100) # this locator puts ticks at regular intervals
      #self.frame[n].xaxis.set_major_locator(loc)
      #self.frame[n].yaxis.set_major_locator(loc)
      
      #self.fig.set_frameon = False
      
      #fr.set_xticklabels(map(lambda x: int(x),fr.get_xticks()), fontsize=16)
      #fr.set_yticklabels(map(lambda x: int(x),fr.get_yticks()), fontsize=16)
   
   def set_anchoredtext(self, n, text, loc=3, fontsize=18):
      from matplotlib.offsetbox import AnchoredText
      # TEXT
      anchored_text = AnchoredText(text,
                       loc=loc, prop=dict(size=fontsize), frameon=False,
                       bbox_to_anchor=(0., 1.),
                       bbox_transform=self.frame[n].transAxes
                       )
      anchored_text.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
      self.frame[n].add_artist(anchored_text)
   
   def set_coords(self, n, coords, lat, lon, marker, tol=0.3):
      '''
         Posit a given sequence of latitude & longitude coordinate. 
      '''
      from numpy import int32, where, logical_and, isclose, median
      
      X, Y = [], []
      for i, c in enumerate(coords):
         print(c)
         XY = where(logical_and(isclose(lat,c[0],atol=tol), isclose(lon,c[1],atol=tol)))
         Y = lat.shape[0]-XY[0]
         X = XY[1]
    
         pos = (1700, 200+i*35)
         try:
           prec=int(max(X)-min(X)+max(Y)-min(Y))/2
           self.frame[n].scatter(X, Y, s=(3*72./self.fig.dpi)**2, marker=marker, facecolors='none', edgecolors='r')
           self.text(n, str(c), pos, fontsize=12, horiz_align='left', vertical_align='top')
         except ValueError:
           pass


   def draw_lines(self, n, coords, labels, x_lbl, y_lbl, **args):
      '''
         Draw lines and print a corresponding text.
      '''
      from numpy import array, around
      fontsize=12
      
      for i, c in enumerate(coords):
        X = c[[0,2]]
        Y = self.image.shape[0]-array(c[[1,3]])
        
        self.frame[n].plot(X, Y, **args)
        self.text(n, str(around(labels[i], 1)), (x_lbl, self.image.shape[0]-y_lbl+35*i), fontsize=fontsize, horiz_align='right', vertical_align='top')


   def contour(self, n, intervals=5, **args):
      CS = self.frame[n].contour(self.image, intervals, linestyle='--', corner_mask=True, **args)
      self.frame[n].clabel(CS, inline=True, inline_spacing=font['size'], fontsize=font['size'], linestyle='--', fmt='%.0f')

   def overplot_image(self, original, filtered, filter_name, **kargs):

      fig, (ax1, ax2) = plt.subplots(ncols=2, sharex=True, sharey=True, **kargs)
      ax1.imshow(original, cmap='gray')
      ax1.set_title('original')
      ax1.axis('off')
      ax1.set_adjustable('box-forced')
      ax2.imshow(filtered, cmap='gray')
      ax2.imshow(original, cmap='Blues', alpha=0.5)
      ax2.set_title(filter_name)
      ax2.axis('off')
      ax2.set_adjustable('box-forced')

   def compare(self, *images, **kwargs):
      """
        Utility function to display images side by side.
    
        Parameters
        ----------
        image0, image1, image2, ... : ndarrray
           Images to display.
        labels : list
           Labels for the different images.
      """
      from numpy import array
    
      f, axes = plt.subplots(1, len(images), **kwargs)
      axes = narray(axes, ndmin=1)
    
      labels = kwargs.pop('labels', None)
      if labels is None:
        labels = [''] * len(images)
    
      for n, (image, label) in enumerate(zip(images, labels)):
         axes[n].imshow(image, interpolation='nearest', cmap='gray')
         axes[n].set_title(label)
         axes[n].axis('off')


# END
