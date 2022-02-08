# -*- coding: utf-8 -*-
#!/usr/bin/python
# Author: Pedro H. A. Hasselmann

__version__="1.0.0"

__doc__='''
  ############## Structuralize Data for Facet manipulation #######################
  ### Author: Pedro H. A. Hasselmann
  ### Last modified: March, 2016
  ### 
  ### Main module for structuralizing data to pandas.DataFrame & HDF5
  ###
  ### 
  ### Retrieve : Radiance & irradiance calculation at facet basis.
  ###            For data on pandas.DataFrame.
  ###
  #################################################################################
'''    

# GLOBAL IMPORT
from .. import *
from ..main import *

# Retrieve Radiance/Flux at facet basis
class Retrieve(object):
   '''
      From image converted to pandas.DataFrame, 
      calculate the radiance/irradiance/flux at facet basis.
      
      Adapted for Store class output.
   '''
   def __init__(self, partition, image_df=None, radiance_label):
     '''
        Input
        =====
        image_df : pandas.DataFrame, geometric pixel-wise information [pixels, facets, incidence, emergence, ...]
     '''
     from numpy import float32
     self.radiance_label = radiance_label
     self.image_df_ = image_df
     # Initialize new dataframe for adding irradiance and radiance arrays
     #self.radiance_df_ = image_df[["utctime","y","x","facet"]].set_index("facet", verify_integrity=False)

     # memory collect
     self.c = collect(1)
   
   @property
   def image_df_(self):
     return self._image_df_
   @image_df_.setter
   def image_df_(self, df):
     df["id"] = df.index
     self._image_df_ = df

   def irradiance(self, Sw, Sha, spc, ccd=(1024,1024)):
     '''
         From all stored geometric & measured information,
         we estimate the irradiance per facet.
         
         Formulation:
         dFj = Sw * dSij/Spix
         
         RADF : dIj/dFj
         
         Si : solid angle at source
         Sw : Solar flux at given wavelength
         Spix : Solid angle of the pixel/boresight
         
         Input
         =====
         sha : shape model filename
         spc : geo.position.from_spice class with mate-kernel loaded.
         ccd : ccd pixel grid
     '''
     from numpy import float32, sum, sqrt, transpose
  
     spc.load_time(self.image_df_["utctime"].loc[0])
     
     # Instrument boresight & pixel resolution
     qua = spc.instrument_frame()[0][4]
     pix_size = (
             (qua[:,0].max()-qua[:,0].min())/ccd[0], 
             (qua[:,1].max()-qua[:,1].min())/ccd[1]
             )
  
     sun = spc.solar_coord(spc.body_frame)[0]
     sc = spc.sc_coord(spc.body_frame)[0]
     
     facet_index = self.image_df_["facet"].drop_duplicates(keep='first').values
     s_source =  Sha.irradiance(Sw, sun, sc)/(pix_size[0]*pix_size[1])
     s_source = pd.Series(s_source[facet_index-1],
                          index = facet_index,
                          name ='s_source',
                          dtype=float32)
     
     # plot
     #plt.hist(s_source, bins=30, alpha=0.7, label='incidence')
     
     self.image_df_["irradiance"] = float32(Sw)*s_source
     self.s_source_eff_ = self.image_df.groupby(("y","x"))["irradiance"].sum()/float32(Sw)


   def isotropic_radiance(self, radiance_label, light_diffusion=1e-10):
     '''
         From all stored geometric & measured information,
         we estimate the radiance factor (I/F) per facet.
         
         The radiance enrgy is isotropic, equally redistributed into all facets.
         
         Formulation:
         dIj = Ipix * Sepix - Ifacet * sum(Sej-1)
         dFj = Sw * dSij
         
         Ifacet = (Ipix * Sepix)/sum(Sej-1)
         
         RADF : dIj/dFj
         
         Se : solid angle at observer
         Si : solid angle at source
         Sw : Solar flux at given wavelength
         Ifacet : Average intensity of the facets
         Ipix : Intensity in the pixel/boresight
         Spix : Solid angle of the pixel/boresight
     '''
     from diskfunc import isoradiance_c
     from numpy import uint32, float32, hstack, any
     from numpy.core import defchararray

     image_df = self.image_df_.set_index("facet", verify_integrity=False)
     facet_pix_count = self.image_df.index.value_counts()
     # facets' solid angle sharing with more than one pixel are divided by the number of duplicates (Approximation!)
     image_df["s"] = image_df.groupby(level=0)["s"].mean()/facet_count

     #plot
     #plt.hist(image_df["s"]*facet_count, bins=30, alpha=0.4, label='emergence')
     #plt.ylabel("N", fontsize=14)
     #plt.xlabel("Solid Angle (str)", fontsize=14)
     #plt.legend(loc=0, fontsize=14)
     #plt.xticks(fontsize=14)
     #plt.show()
     #plt.clf()
     
     # Pixel basis
     pixels =  image_df.groupby(("y","x"))
     # Compute isotropic radiance  from pixel to facet basis
     radiance =  hstack(pixels[[radiance_label,"s"]].apply(isoradiance_c,
                                                    diff=light_diffusion,
                                                    label=self.radiance_label).values)
     print(radiance.shape)
     self.image_df_["isoradiance"] = pd.Series(radiance[1,...], index=radiance[0,...].astype(uint32), dtype=float32)
     
     #plot
     #plt.hist(self.s_eff_, bins=30, alpha=0.4, label='emergence')
     #plt.ylabel("N", fontsize=14)
     #plt.xlabel("Effective Solid Angle (str)", fontsize=14)
     #plt.legend(loc=0, fontsize=14)
     #plt.xticks(fontsize=14)
     #plt.show()
     #plt.clf()

   def point_spread_function(self, psf):
     '''
        TO DO:
        Apply the instrument Point-spread-function to obtain the exact energy distribution 
        contribution of every facet to the pixels during the image/boresight integration.
        
        Compute the proportion of energy contribution for every facet.
     '''
     pass

# END
