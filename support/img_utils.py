# -*- coding: utf-8 -*-

__version__ = '0.1.0'

__doc__ = '''
##############################################################################
#  author: Pedro H. A. Hasselmann
#  script: img_utils.py
#
#  A script for dealing with images of different formats.
#  IMG (PDS format binary)
#  SAV (IDL)
#  FIT (IAU)
#  TABLE
#
#  dependencies:
#  numpy, matplotlib, astropy, gdal, PyPDS
#
#
#
#
###############################################################################
'''

###########################################
## Treat and intensity renormalize image ##
###########################################

def treatimage(img, 
               s=([0,None],[0,None]),
               scale    = (0,99.99),
               contrast = 0.5,
               sigma    = 2.0,
               erode    = False,
               r = 2,
               **args):
  '''
     Rescales intensity and smooths down.
     
     Parameters
     ----------
     s : 2-tuple - 2-list ([int,int],[int,int])
      Crop limits
      
     contrast : float
      intensity contrast 
  
     erode : boolean
      Active the function to degrading resolution.
     
     r : int
      Pixel radius of the convolutive disk
  '''                                    
  from skimage.exposure import rescale_intensity, equalize_adapthist
  from numpy import isnan, isinf, nanpercentile, dstack

  # Cropping
  img = img[::-1][s[0][0]:s[0][1], s[1][0]:s[1][1]] # (Y, X)
  img[isnan(img)] = 1.0
  img[isinf(img)] = 1.0
  img[img<0e0] = 0e0


  inf = nanpercentile(img, scale[0])
  sup = nanpercentile(img,scale[1])
  print('lower: ', inf, 'upper: ', sup)
  
  # Equalization
  #img_denoised = denoise_bilateral(img, sigma_spatial=sigma, bins=100, multichannel=False)
  img_rescaled = rescale_intensity(img, in_range=(inf,sup))
  #img_equalized = equalize_adapthist(img, clip_limit=contrast)

  if erode:
      from skimage.filters import rank
      from skimage.morphology import square

      img_rescaled = rank.median(img_rescaled/img_rescaled.max(), square(r))
      img_rescaled = img_rescaled/img_rescaled.max()

  #img_equalized[img_equalized > 1.0] = 1.0
  img_rescaled[img_rescaled > 1.0] = 1.0
  print('Rescaled', 'lower: ', img_rescaled.min(), 'upper: ', img_rescaled.max())
  
  return dstack((img, img_rescaled))#, img_equalized))


def rebin(arr, new_shape):
    """Rebin 2D array arr to shape new_shape by averaging."""
    from numpy import nanmean
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return nanmean(arr.reshape(shape), axis=(-1,1))

##########################
## READ IMAGES UTILITY ###
##########################

def image_treated(imagename, channel, **args):
   from numpy import float32
   import argparse

   # Image Reference
   ref = read_image(imagename, channel=channel)
   cube = treatimage(ref.data, **args).astype(float32, order='C',copy=False)
   return argparse.Namespace(**{'data':cube, 'header':ref.header})


def read_image(imagename, channel, **args): 

   if imagename[-4:].lower() =='.fts' or imagename[-4:].lower() =='.fit' or imagename[-5:] =='.fits':
     return read_fits(imagename, channel, **args)
   
   elif imagename[-4:] =='.IMG' or imagename[-4:] =='.img':   
     return read_PDS(imagename, channel, **args)


################
## FITS Image ##
################
def read_fits(imagename, channel, treat=False, **treat_args):
   ''' 
      Decode fits-image or cubes using astropy.fits. 
   '''
   from astropy.io import fits
   from numpy import float32, dstack
   image = fits.open(imagename, ignore_missing_end=True)
   #print(image.info())
   
   if channel != None:
     return image[channel]

   # only for cube-images.   
   else:
     channel = 1
     cube = list()
     header = list()
     
     while True:
       try:

         if treat:
            img = treatimage( image[channel].data.astype(float32), **treat_args)[1]
         else:
            img = image[channel].data.astype(float32)
         
         cube.append(img)
         try:
           header.append(image[channel].header)
         except NameError:
           header.append(None)

         channel += 1
      
       except IndexError:
         import argparse
         return argparse.Namespace(**{'data':dstack(cube), 'header':header})

# end class

##############################################################
# PDS-format image -- HEADER, IMAGE, SIGMA_MAP, QUALITY_MAP ##
##############################################################
class read_PDS:
 
   from numpy import array, float16
   from collections import deque
   global array, deque, float16
 
   def __init__(self, imagename, channel=1, sigma_map=False):
     '''
      Read the PDS-format (.IMG) image using PyPDS.
      GDAL - Geospatial Data Abstraction Library --> http://www.gdal.org/index.html
     '''
     from .pds.core.common import open_pds
     from .pds.core.parser import Parser
     import struct
     from osgeo import gdal
     from osgeo.gdalconst import GA_ReadOnly, GDT_Float32

     self.imagename = imagename

     # Mode
     if sigma_map:
       return self.sigma_map()
     else:
       pass


     # Data
     data  = gdal.Open(imagename, GA_ReadOnly)
     print(data)
     drive = data.GetDriver()
     band = data.GetRasterBand(channel)
     #print('Driver: ', drive.ShortName,'/', drive.LongName)
     #print('Size is',data.RasterXSize,'x',data.RasterYSize,'x',data.RasterCount,
     #'type is',gdal.GetDataTypeName(band.DataType))
     
     image = deque()
     for Y in xrange(band.YSize):
       try:
          scan  = band.ReadRaster( 0, Y, band.XSize, 1, band.XSize, 1, GDT_Float32 )
          image.append( array( struct.unpack(str('f' * band.XSize), scan) ) )
        
       except struct.error:
          break

     self.data = array(image)
   
     # Header
     parser = Parser()
     self.header = parser.parse(open_pds(imagename))


   def sigma_map(self):
     '''
      Decode the sigma map into readable Numpy array.   
     '''
     import struct
     
     bytelist = deque()
     image = deque()
     header = dict()
     obj, flag = '', False

     # Begin read .IMG file
     with open(self.imagename, 'r') as data:
       for n, line in enumerate(data):
          
          # Header
          try:
            if line.split('=')[0].strip() == 'OBJECT':
              obj = line.split('=')[1].strip()+':'
     
            elif line.split('=')[0].strip() == 'END_OBJECT':
              obj = ''
            
            
            header[obj+line.split('=')[0].strip()] = line.split('=')[1].strip()

          except IndexError:
            pass

          except UnicodeDecodeError:
            flag = True

          # Collect bytes 
          if flag: bytelist.append(line)

       bytemap = str("").join(bytelist)

       #image_offset = int(header['RECORD_BYTES']) * (int(header['^SIGMA_MAP_IMAGE']) - 1)
       #print('offset: ', image_offset)

       w = int(header['SIGMA_MAP_IMAGE:LINE_SAMPLES'])
       #l = int(header['SIGMA_MAP_IMAGE:LINES'])
       #band = int(header['SIGMA_MAP_IMAGE:BANDS'])
       #byt = int(header['SIGMA_MAP_IMAGE:SAMPLE_BITS'])
       #typ = str(header['SIGMA_MAP_IMAGE:SAMPLE_TYPE'])
       
       # Decode bytes data
       raster = map(str("").join, zip(*[iter(bytemap)]*struct.calcsize(str('f' * w))))

       image = array([array( struct.unpack(str('f' * w), segment) ) for segment in raster])
             
       self.data = array(image)[:2048]
       self.sigma = array(image)[2048:2*2048]
       self.header = header


   def __str__(self):
     import json
     return json.dumps(self.header, indent=2)
     #return reduce(lambda x,y: str(x)+'\n'+str(y)+'\n', self.header.items())
# End class

def to_fits(filename, image):
   from astropy.io import fits
   from numpy import float32

   hdu = fits.PrimaryHDU(image.astype(float32))
   #hdulist = fits.HDUList([hdu])
   hdu.writeto(filename)#, output_verify='warn')
   #hdulist.close()

###############################
### TABLE - IMAGE FUNCTIONS ###
###############################
def image_to_table(image):
   '''
      2D image array to pandas.DataFrame table [Y, X, value].
   '''
   import pandas as pd
   from numpy import float32, uint32, arange, stack, meshgrid
   x, y = image.shape
   # Build image back from pixel position
   # https://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays      
   index = (arange(y, dtype=uint32), arange(x, dtype=uint32))
   yx = stack(meshgrid(*index), -1).reshape(-1, 2)[:,::-1]
   data = pd.Series(image.flatten(), index=pd.MultiIndex.from_arrays(yx.T, names=('y','x')), dtype=float32)
   return data, yx

def table_to_image(table, dpi, b=0e0):
   '''
      Transform pandas.DataFrame table of pixel-index & values into an image.
      
      Parameters
      ==========
      table     : Pandas.Series, Multi-index ('y','x') & value
      dpi       : 2-item tuple, image dimension
   '''
   import pandas as pd
   from numpy import any, float32, uint32, arange, stack, meshgrid
   x, y = dpi
   if any(table.index.duplicated(keep=False), axis=0) == True:
     table = table.groupby(level=(0,1)).mean()
    # Build image back from pixel position
    # https://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays
   index = (arange(y, dtype=uint32), arange(x, dtype=uint32))
   yx = stack(meshgrid(*index), -1).reshape(-1, 2)[:,::-1]
   background= pd.Series(b, index=pd.MultiIndex.from_arrays(yx.T, names=['y','x']), dtype=float32)
   return background.add(table, fill_value=b).values.reshape(dpi)#.values[:dpi[0]*dpi[1]].reshape(dpi)

# SAV Cube
def sav_cube(filename, mode):
   '''
      Read IDL SAV format.
   '''
   from scipy.io import readsav
   data = readsav(filename)
   return data[mode]


# END
