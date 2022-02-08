# -*- coding: utf-8 -*-
#!/usr/bin/python

__version__="0.1.0"

__doc__='''
  ############################ MAIN PLOT TOOL ################################
  ### Author: Pedro Henrique A. Hasselmann
  ### Last Modified: August, 2019
  ###
  ### Tool for graphical plotting using matplotlib.
  ###
  ###
  ### Class: Plot
  ###
  ### Attributes:
  ###          
  ###  
  ###   
  ###     
  ###             
  ###
  #############################################################################
'''
####################################################################################
from .. import *
from ..main import *


from itertools import cycle

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D

##########################################################################################

#http://www.randalolson.com/wp-content/uploads/percent-bachelors-degrees-women-usa.png
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
def tableau(start):
  for i in tableau20[start:]:    
    r, g, b = i    
    yield (r / 255., g / 255., b / 255.)

color_cycle = cycle(('#FE2E2E','#B40404','#08088A','#2E9AFE','green','#FF8000','#086A87','purple','#8A0829','magenta','#FF0040','#0B3B0B','black'))
color_space = lambda x: cycle(cm.Paired(linspace(0, 1, x)))
symbols = cycle(('^','o','d','s','*','p','8','+','1','x'))
linestl = cycle(('solid', 'dashed', 'dashdot', 'dotted'))


##############################################################
################### GRAPHIC PLOTS ############################
##############################################################
class Plot:
  __doc__ = ''' Python Class dealing with plot display. 
                Matplotlib interface. 
            '''
  import matplotlib
  from numpy import linspace
  global matplotlib, linspace

  def __init__(self, ncols=1, nlines=1,
               size=(10,8),
               titles=['',],
               suptitle='', 
               fcolor='k',
               fontsize=16,
               cycle_color = cm.Spectral_r(linspace(0,1,15)),
               background='white', 
               shared=[False,False],
               nosep=False,
               projection=None):

     font = {'family' : 'serif',
             'weight' : 'bold',
             'size'   : fontsize}
     self.fontsize = fontsize
     self.cmap = cycle_color
     matplotlib.rc('font', **font)
     #matplotlib.rc('text', usetex=True)
     matplotlib.rcParams['agg.path.chunksize'] = 1000
     from cycler import cycler
     matplotlib.rcParams['axes.prop_cycle'] = cycler(color=cycle_color)

     self.fig, frame = plt.subplots(ncols, nlines, figsize=size, sharey=shared[1], sharex=shared[0], 
                                    subplot_kw={'adjustable':'box','projection':projection})
     self.fcolor =fcolor
     self.background = background
     self.titles = titles
     
     if suptitle is not None:
       plt.suptitle(suptitle, fontsize=int(font['size'])+4, fontweight='bold')
      
     if ncols == 1 and nlines == 1: 
             self.frame = [frame, ]
     else:
             self.frame = frame.flatten()
         
     for n, fr in enumerate(self.frame):
       fr.set_facecolor(background)
       fr.set_alpha(0.8)
       fr.minorticks_on()
       fr.tick_params('both', length=9, width=1.8, which='major', direction='inout', color=fcolor, top=True, right=True)
       fr.tick_params('both', length=4, width=0.9, which='minor', direction='inout', color=fcolor, top=True, right=True)

       fr.set_title(titles[n], fontweight='bold', fontsize=int(font['size']))

       if (nosep==True):#&(shared[0]==True):
         plt.setp(fr.get_xticklabels(), visible=False)
         plt.subplots_adjust(hspace=.0)
       #if (nosep==True)&(shared[1]==True):
         plt.setp(fr.get_yticklabels(), visible=False)
         plt.subplots_adjust(wspace=.0)


     plt.xticks(fontsize = font['size'])
     plt.yticks(fontsize = font['size'])

     #if nosep==True:
     #   plt.subplots_adjust(hspace=.0, wspace=.0)

  def hide(self,n):
     self.frame[n].axis('off')

  def canvas(self, n, on_off, which):
    if on_off=="on": flag=True
    if on_off=="off": flag=False
    self.frame[n].spines[which].set_visible(flag)
    args=dict(); args[which] = on_off
    self.frame[n].tick_params(which='major',**args)
    self.frame[n].tick_params(which='minor',**args)

  def yaxis_right(self, n):
    self.frame[n].yaxis.set_label_position("right")
    self.frame[n].yaxis.tick_right()
    plt.setp(self.frame[n].get_yticklabels(), visible=True)
    
  def show(self, tight=False):    
     if tight: plt.tight_layout()
     plt.draw()
     plt.show()
     plt.clf()
  
  def clear(self):
    self.fig.clear()

  def close(self):
    plt.cla()
    plt.clf()
    plt.close()

  def save(self, pathname, dpi=300, tight=False, **args):    
     if tight: plt.tight_layout()
     plt.savefig(pathname, dpi=dpi, bbox_inches='tight', **args)

  def ylim(self,minimum,maximum):
     for fr in self.frame:
        fr.set_ylim(minimum, maximum)
        fr.set_ybound(minimum, maximum)

  def xlim(self,minimum,maximum):
     for fr in self.frame:
        fr.set_xlim(minimum, maximum)
        fr.set_xbound(minimum, maximum)

  def xlabel(self,n,x, **args):
    if 'fontsize' not in args: args['fontsize'] = self.fontsize+2
    self.frame[n].set_xlabel(x, fontweight='bold', **args)

  def ylabel(self,n,y, **args):
    if 'fontsize' not in args: args['fontsize'] = self.fontsize+2
    self.frame[n].set_ylabel(y, fontweight='bold', **args)

  def label_axes(self,x,y):
     for fr in self.frame:
        fr.set_xlabel(x, fontsize=self.fontsize, fontweight='bold')
        fr.set_ylabel(y, fontsize=self.fontsize, fontweight='bold')

  def colorbar(self, figure, title, ticks=15, vmin=None, vmax=None, **args):
     # Color Bar

     self.fig.subplots_adjust(right=0.75)
     ax = self.fig.add_axes([0.78, 0.15, 0.025, 0.7])
     cm = plt.cm.ScalarMappable(cmap=self.cmap)
     cm.set_array(figure)
     cm.set_clim(vmin, vmax)
     cbar = self.fig.colorbar(cm, cax=ax, orientation='vertical', pad=.04, aspect=12, fraction=0.03)
     if title != None: cbar.ax.set_ylabel(title, fontsize=28, fontweight=font['weight'], labelpad=18.0)

     tick_locator = ticker.MaxNLocator(nbins=ticks)
     cbar.locator = tick_locator
     cbar.ax.tick_params(labelsize=int(1.5*self.fontsize))
     cbar.update_ticks()

  def legend(self, n, color='k', **args):
     l = self.frame[n].legend(loc=0, frameon=False, numpoints=1, fontsize=self.fontsize-2, **args)
     plt.setp(l.get_texts(), color=color)

  def text(self, n, write, pos, color='k', **args):
              self.frame[n].text(
                        pos[0],
                        pos[1], 
                        write, 
                        color=color, 
                        fontsize=font['size']-2,
                        fontweight='bold',
                        horizontalalignment=args['horiz_align'],
                        verticalalignment=args['vertical_align']
                        #bbox=dict(facecolor='white', edgecolor='white', alpha=0.7)
                        )

  def points3d(self, nnn, x, y, z, **args):
     frame = self.fig.add_subplot(nnn, projection='3d')
     frame.plot_surface(x, y, z, **args)
     return frame

  def scatter(self, n, x, y, c, cmap='Spectral', **args):
    draw = self.frame[n].scatter(x, y, c=c, marker='o', edgecolors='none', cmap=cmap, **args)
    return draw

  def plotcurve(self, n, x, y, linewidth=2.5, **args):
    draw = self.frame[n].plot(x, y, linewidth=linewidth, **args)
    return draw

  def hist(self, n, x, kde=True, bandwidth=0.2, **args):
    from scipy.stats import gaussian_kde
    from numpy import linspace

    if kde==True:
      k = gaussian_kde(x, bandwidth)
      points = linspace(x.min(), x.max(), 100)
      Z = k.evaluate(points)
      self.frame[n].plot(points, Z, linewidth=2.5, color='k')
    else:
      self.frame[n].hist(x, linewidth=1.5, edgecolor='black', **args)


  def imshow(self, n, image, **args):
    im = self.frame[n].imshow(image[::-1], interpolation='none', **args)
    return im




  def plot_kde(self, x, wt, bandwidth, xlabels, N=3000, **args):
    from numpy import tile, nonzero, zeros, sort, argsort, sqrt, diagonal, mean, linspace, ones, percentile
    from scipy.stats import gaussian_kde
    
    #k = gaussian_kde(x, bandwidth, wt)
    #tot = k.integrate_box(x.min(1),x.max(1))
    #std = sqrt(diagonal(k.covariance))

    #grid = k.resample(N)
    #print(grid.shape)
    #Z = k.evaluate(grid)
    #indx = nonzero(Z == Z.max())[0][0]

    iqr = {}
    for n, lbl in enumerate(xlabels):
      k = gaussian_kde(x[n,...], bandwidth, wt)
      tot = k.integrate_box_1d(x[n,...].min(), x[n,...].max())
      points = linspace(x[n,...].min(), x[n,...].max(), 500)
      Z = k.evaluate(points)/tot
    
      std = sqrt(k.covariance)[0][0]
      mode = points[nonzero(Z==Z.max())][0]
      sm = k.resample(N)
      p5 = percentile(sm, 5)
      p25 = percentile(sm, 25)
      p50 = percentile(sm, 50)
      p75 = percentile(sm, 75)
      p95 = percentile(sm, 95)

      #sorting = argsort(grid[n,...])
      #indx = nonzero(wt==wt.max())
      self.frame[n].scatter(x[n,...],zeros(x.shape[1]), s=3, color='k')
      self.frame[n].plot(points, Z/Z.max(), linewidth=2.5, color='k', **args)
      #self.frame[n].plot(x[n,indx][0]*ones(3),linspace(Z.min(),Z.max(),3), linestyle='-.',linewidth=2.5, color='red',)
      self.frame[n].plot(mode*ones(3),linspace(0.0,1.0,3), linestyle='-.',linewidth=2.5, color='k',)
      #self.frame[n].hist(x[n,...], bins=30, histtype='stepfilled',linewidth=3, **args)
      self.xlabel(n, lbl, fontsize=18)
      
      #self.legend(n, 'mode: {:.3f} std: {:.4f}'.format(mode, std))
      iqr[lbl] = {'mode':mode, 'std':std, 'p5':p5, 'p25':p25, 'p50':p50, 'p75':p75, 'p95':p95}

    return iqr
    
# END
