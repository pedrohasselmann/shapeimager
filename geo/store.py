# -*- coding: utf-8 -*-
#!/usr/bin/python
# Author: Pedro H. A. Hasselmann

################################################
#  Structuralize Data for Facets manipulation. #
################################################

__version__="0.2"

__doc__='''
  ############## Structuralize Data for Facet manipulation #######################
  ### Author: Pedro H. A. Hasselmann
  ### Last modified: April, 2021
  ### 
  ### Main module for structuralizing data to pandas.DataFrame & HDF5
  ###
  ### Store : Structuralization
  ### 
  ###        
  ###
  #################################################################################
'''    
from .. import *
from ..main import *

## BINNING
def binning(data, col_data, cols_bin,  nbins, kernel="mean"):
  from scipy.stats import binned_statistic_dd, scoreatpercentile
  from numpy import float32, stack, concatenate, meshgrid, newaxis, degrees

  data = data.reset_index()
  cols_bin = cols_bin.copy()

  upper = scoreatpercentile(data[cols_bin], per= 99, axis=0)
  lower = scoreatpercentile(data[cols_bin], per= 0, axis=0)
  print(upper)
  print(lower)
  print('angle step (degrees): ', (upper-lower)/array(nbins))

  binned_data, edges, locator = binned_statistic_dd(data[cols_bin].values, 
                                                    data[col_data].values, 
                                                    bins=nbins,
                                                    range=list(zip(lower.tolist(),upper.tolist())), 
                                                    statistic=kernel)
  central_point = [array([sum(edges[n][i:i+2])/2.0 for i in range(0, len(edges[n]))])[:-1] for n in range(0,len(nbins))]
  mesh = stack(meshgrid(*central_point, indexing='ij'), axis=len(nbins))
  binned = concatenate((mesh, binned_data[...,newaxis]), axis=len(nbins)).astype(float32)
  
  cols_bin.append(col_data)
  return pd.DataFrame(binned.reshape(-1, binned.shape[-1]), columns=cols_bin).dropna(how='any')





# NPZ
def to_npz(label, utctime, Imager, **items):
  from numpy import savez, vstack

  #print(Imager.facet_area_pix_)
  items.update({'dist':Imager.d_.mean(1).flatten(),
                'inc' :Imager.inc_, 
                'emi' :Imager.emi_, 
                'pha' :Imager.pha_,
                'facetid' :Imager.facetid_})

  if hasattr(Imager, "facet_image_"):
    try:
      to_fits('facet_'+label+'.fit', Imager.facet_image_)
    except IOError:
      pass

    items.update({'s'              :Imager.solid_angle_emi_,
                  #'s_source' :Imager.solid_angle_inc_,
                  'facet_area_pix_' :Imager.facet_area_pix_,
                  'facet_pix'      :Imager.facet_pix_})


  print('geo_'+label+'_'+utctime.replace(':','.')+'.npz')
  savez('geo_'+label+'_'+utctime.replace(':','.')+'.npz',
         utctime=utctime,
         vectors=vstack((Imager.source_pos, 
                         Imager.obs_pos, 
                         Imager.cam_matrix, 
                         Imager.boresight)),
         **items)






# SHAPEMODEL PATCHING DATA
def __(utctime, store, cols, facets):
      with pd.HDFStore(path.join(home,core,store), mode='r') as hdf:
        image_data = hdf.select(utctime)[cols]
        facet_rows = image_data.merge(facets, on=['facet'])
        #facet_rows = image_data[image_data.facet.isin(facets)]
        hdf.close()
      return facet_rows


def store_shape_patches(storage, columns, patch_params, pools=6):
   '''
     Store into the necessary format for the local phase modeling.
     https://stackoverflow.com/questions/40684896/multiprocessing-hdf5-file-read-what-is-better-multiple-connections-or-repeated
     http://matthewrocklin.com/blog/work/2015/03/16/Fast-Serialization
     https://stackoverflow.com/questions/23945493/a-faster-alternative-to-pandas-isin-function
     
     Shape model is break into patches and data for every patch is gathered into pickle.
   
   '''
   import cPickle as cpkl
   from functools import partial
   # Output Store filename
   #pklname = path.join(home,prod,'patches_{}_{}_{}.pkl'.format(obj, fi, patch_params['sha'][-8:-4]))
   #storeout = open(pklname, mode='ab')
   directory = make_dir(path.join(home,prod,'patches_{}_{}_{}'.format(obj, fi, patch_params['sha'][-8:-4])))
   columns.extend(['facet','y','x','utctime'])

   # Load HDF5
   S = ShapeModel(patch_params['sha'], comments=71)
   patch_params.pop('sha', None)
   patches = S.partitioning(mode='yield',**patch_params)

   c = collect(100)
   
   for N, p in enumerate(patches):

     print('patch ',N+1)
     print('#facets ',p[1].shape[0])
     print(directory)
     pklname = path.join(directory,'patch_{}.pkl'.format(N+1))
     if path.isfile(pklname)==False:
       storeout = open(pklname, mode='ab')
       facets_df = pd.DataFrame(p[1], columns=['facet'])
     
       facet_rows = list()
       for sto in storage:
          print(sto)
          with pd.HDFStore(path.join(home,core,sto), mode='r') as hdf:
            utctimes = hdf.keys()#pd.read_csv(path.join(home,core,sto.replace('h5','txt')), delim_whitespace=1, usecols=(0,)).values.flatten().tolist()
            hdf.close()
          # Seach for facets at every image (multi-processing)
          with poolcontext(pools, init_worker) as pool:
              try:
                facet_rows.extend(pool.map(partial(__, store=sto, cols=columns, facets=facets_df), utctimes))
              except KeyboardInterrupt:
                   print("Caught KeyboardInterrupt, terminating workers")
                   pool.terminate()
                   pool.join()
       #raw_input()
       facet_data = pd.concat(facet_rows, ignore_index=True)#.groupby(['y','x','utctime']).mean()
       facet_data['patch'] = N+1
       #facet_data.rename(columns={'residue':'weight'}, inplace=True)
       facet_data.to_pickle(storeout, compression=None, protocol=2)
       storeout.close()
       #print(facet_data.head(5))
       print('#pixels',facet_data.shape[0])
       c.next()
        
    #hdf.close()
    #print('done...')




# STORING
class Store(object):
   '''
     Store into HDF5/Pickle.
     Adapted for output from source.geo.shapeimager
   '''
   def __init__(self, columns, label='test'):  
     
     # Open HDF5 file
     self.partition = 'geo_data'
     columns.extend(('y','x','facet','utctime'))
     self.columns = set(columns)
     self.store = pd.HDFStore('{}.h5'.format(label), mode='a')
     #self.store = open(path.join(home,prod,obj+'_'+fi+'_{}.pkl'.format(label)), 'ab')
     pd.set_option('io.hdf.default_format','table')
     # memory collect
     self.c = collect(1)

   def image_dataframe(self, value_name, image, geo, offset=(0,0), threshold=5e-5):
      from numpy import uint32, int32, float32, asarray
      
      # 1. Multi-index: Match data to the facets for every pixel
      # data on pandas.DataFrame
      data = image_to_table(image)[0].to_frame(value_name)
      data.query('{}>{}'.format(value_name,threshold), inplace=True)
      
      # facet_pix on pandas.DataFrame
      yx = geo["facet_pix"][...,(0,1)]-asarray(offset, dtype=uint32)
      yx = pd.MultiIndex.from_arrays(yx.T, names=('y','x'))
      facet_pix = pd.DataFrame(geo["facet_pix"][...,2], columns=('facet',), index=yx)
      facet_pix['area'] = geo['facet_area_pix_'][:,1].astype(float32)

      # plot - test
      #import matplotlib.pyplot as plt
      #plt.imshow(table_to_image(data["IF"], (1024,1024), b=0e0), cmap='Greys', interpolation='none', alpha=0.99)
      #plt.imshow(table_to_image(facet_pix["facet"], (1024,1024), b=0e0), interpolation='none', cmap='Blues', alpha=0.5)
      #plt.title('Pixel Correspondence')
      #plt.show()

      # Merge both table and conserve duplicity
      # https://stackoverflow.com/questions/49786769/need-to-handle-a-concatenated-dataframe-with-non-unique-multi-index
      facet_pix = facet_pix.merge(data, left_index=True, right_on=('y','x'), how='outer')
      facet_pix.dropna(how='any', axis=0, inplace=True)
      
      # 2. Add geo columns
      dictio = dict()
      dictio['utctime'] = asarray([geo['utctime'],]*facet_pix.index.size, dtype=str)
      dictio['y']       = asarray(facet_pix.index.get_level_values('y'), dtype=uint32)
      dictio['x']       = asarray(facet_pix.index.get_level_values('x'), dtype=uint32)
      dictio['facet']   = facet_pix['facet'].values.astype(uint32)
      dictio['area']     = facet_pix['area'].values.astype(float32)
      dictio[value_name] = facet_pix[value_name].values.astype(float32)
      
      for c in self.columns:
         try:
           #print(c)
           dictio[c] = geo[c][dictio['facet']-1].astype(float32)
         except (IndexError, KeyError):
           #print('fail...',c)
           pass
       
      self.image_df_ = pd.DataFrame.from_dict(dictio)
      print(self.image_df_.head(10))
      print(self.image_df_.shape)


   def storing_block(self, compression=3):
      self.image_df_.to_hdf(self.store, self.partition, format='table', data_columns=True, append=True, complib='blosc', complevel=compression)
      #raw_input(self.store[self.name].head(10))
      next(self.c)

   def storing_image(self, compression=3):
      self.image_df_.to_hdf(self.store, self.image_df_["utctime"].iloc[0], format='fixed', data_columns=True, append=False, complib='blosc', complevel=compression)
      #raw_input(self.store[self.name].head(10))
      next(self.c)

   def storing_image2(self, protocol=3):
     self.image_df_.to_pickle(self.store, None, protocol)

   def close(self):
      self.store.close()



# END
