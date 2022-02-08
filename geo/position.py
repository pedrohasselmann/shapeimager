#!/usr/bin/python
# -*- coding: utf-8 -*-
# Filename: position.py
# Retrieve target, observer and camera positioning from SPICE kernels or image headers
# Author: Pedro H. A. Hasselmann

version = "0.2"

doc = '''
  #########################################################################################
  #  author: Pedro H. Hasselmann, Rio de Janeiro (LESIA-OBSPM/ON-MCT)
  #  script: position.py
  #
  # Position and pointing vectors retrieval for OASIS
  #
  #
  # dependencies:
  #
  # Numpy, Scipy, Astropy, spiceypy, PyPDS
  #
  ##########################################################################################

   Essential vectors for OASIS:
   
   OBJ:
   X-coordinate in EQU J2000 (AU)
   Y-coordinate in EQU J2000 (AU)
   Z-coordinate in EQU J2000 (AU)
   Euler angle #1 (deg)
   Euler angle #2 (deg)
   Euler angle #3 (deg)
   rotational rate (deg/day)
   date at which angle #2 = VALUE_E2 (JD)

   OBS:
   # Spacecraft
   bodycentric-X S/C position in EQU J2000 (AU)
   bodycentric-Y S/C position in EQU J2000 (AU)
   bodycentric-Z S/C position in EQU J2000 (AU)
   
   # Camera
   OBS_RA
   OBS_DEC
   OBS_ROLL

'''
from .. import *
from ..main import *

def from_header(pds_image,order='AU'):
   '''
      Obtain positioning from PDS image header.
      
      time = START_TIME
      
      Julian Day (JD2000)
      
   '''
   from numpy import array
   from astropy.time import Time
   from pds.core.common import open_pds
   from pds.core.parser import Parser

   
   # load headers
   parser = Parser()
   header = parser.parse(open_pds(pds_image))
   
   FRAME = path.basename(pds_image)
   
   # Set time in JD2000
   TIME = header['START_TIME']
   
   t = Time(TIME, format='isot', scale='utc')
   print(t, t.jd)
   
   # Point to the necessary header keys
   SC_COORD = header['SC_COORDINATE_SYSTEM']  # Cartesian <km> 
   CAM_COORD = header['CAMERA_COORDINATE_SYSTEM']  # Cartesian <km>
   #IMAGE_COORD = header['IMAGE_POI']    # Georeferencing
   
   #for i in SC_COORD.items(): print(i)
   #for i in CAM_COORD.items(): print(i)
   #for i in IMAGE_COORD.items(): print(i)

   CAM =header['INSTRUMENT_ID']
    
   SC_TARGET_POS = header['SC_TARGET_POSITION_VECTOR']  # <km>
   SC_TARGET_VEL = header['SC_TARGET_VELOCITY_VECTOR']  # <km/s>
   SC_SUN_POS = header['SC_SUN_POSITION_VECTOR']        # <km>
   PHASE_ANGLE = header['PHASE_ANGLE']                  # <deg>
   LATITUDE = header['SUB_SPACECRAFT_LATITUDE']
   LONGITUDE = header['SUB_SPACECRAFT_LONGITUDE']
   SOLAR_DISTANCE = header['SPACECRAFT_SOLAR_DISTANCE']
   ALTITUDE = header['SPACECRAFT_ALTITUDE']
   TARGET_RA = header['RIGHT_ASCENSION']
   TARGET_DEC = header['DECLINATION']

   # Convert from string to float64
   SC_TARGET_POS_2 = array(map(float, [s.translate(None, " (),") for s in SC_COORD['ORIGIN_OFFSET_VECTOR'].split('<km>')][:-1]))
  
   SC_TARGET_POS = array(map(float, [s.translate(None, " (),") for s in SC_TARGET_POS.split('<km>')][:-1]))
   SC_SUN_POS = array(map(float, [s.translate(None, " (),") for s in SC_SUN_POS.split('<km>')][:-1]))
   #SC_TARGET_VEL = map(float, [s.translate(None, " (),") for s in SC_TARGET_VEL.split('<km/s>')][:-1])

   #print(SC_TARGET_POS_2) # SC_TARGET 2   <-- it seems the most correct one
   #print(SC_TARGET_POS)   # SC_TARGET ORIG
   #print(SC_SUN_POS)  

   if order == 'AU':
     import astropy.units as u
     au_km = 149597870.659999996424 #u.au.to(u.km)
     
     TARGET_SUN = (SC_TARGET_POS + SC_TARGET_POS_2)/au_km
     SC_TARGET_POS_2 = SC_TARGET_POS_2/au_km
     
     print('Target heliocentric vector (UA).....',TARGET_SUN)
     print('S/C heliocentric vector (UA).......',SC_TARGET_POS_2)
   
   elif order == 'km':
     TARGET_SUN = SC_TARGET_POS + SC_TARGET_POS_2
     print('Target heliocentric vector (km).....',TARGET_SUN)
     print('S/C heliocentric vector (km)........',SC_TARGET_POS_2)

   # Camera Spherical Coordinates
   CAM_POS = array(map(float, [s.translate(None, " (),") for s in CAM_COORD['ORIGIN_OFFSET_VECTOR'].split('<km>')][:-1]))
 
   return t, SC_TARGET_POS_2, TARGET_SUN, CAM_POS

######################################
## POSITION FROM NAIF/SPICE KERNELS ##
######################################
class from_spice:

  import spiceypy as spice
  global spice

  def __init__(self,
               body, # Body code
               obs,  # Spacecraft code
               ins,  # Instrument code
               body_frame = None, # Body Frame name, str
               ins_frame = None   # Instrument frame name, str
               ):
    '''
       __init__(self, body, obs, ins, body_frame = None, ins_frame=None)
       
       Parameters
       ==========
       body        : int, Body/target SPICE code
       obs         : int, Spacecraft SPICE code
       ins         : int, Instrument SPICE code
       body_frame  : str, Body/target frame name
       ins_frame   : str, instrument frame name
    '''
    print(spice.tkvrsn("toolkit"))
   
    #if obs == -37:      self.obs_frame='HAYABUSA2_SC_BUS_PRIME'
    #if obs == -64:      self.obs_frame='ORX_SPACECRAFT'
    #if obs == -226:     self.obs_frame='ROS_SPACECRAFT'
   
    #if body == 1000012: self.body_frame='67P/C-G_CK'
    #if body == 2162173: self.body_frame='RYUGU_FIXED'
    #if body == 2101955: self.body_frame='IAU_BENNU'
    
    #if ins == -226111: self.ins_frame='ROS_OSIRIS_NAC'
    #if ins == -226112: self.ins_frame='ROS_OSIRIS_WAC'
    #if ins == -64360:  self.ins_frame='ORX_OCAMS_POLYCAM'
    #if ins == -37100:  self.ins_frame='HAYABUSA2_ONC-T'

    self.body = str(body)
    self.obs = str(obs)
    self.ins = str(ins)

    self.body_frame = body_frame
    self.ins_frame = ins_frame

  def furnish(self,tm):
    spice.furnsh(path.join(home,kern,tm))

    if self.body_frame is None:
        self.body_frame = spice.bodc2n(int(self.body))
    else:
        self.body_frame = str(self.body_frame)

    if self.ins_frame is None:
        self.ins_frame = spice.bodc2n(int(self.ins))
    else:
        self.ins_frame = str(self.ins_frame)

    print('#KERNELS.......',spice.ktotal('ALL'))
    print('BODY...........',self.body, spice.bodc2n(int(self.body)))
    print('OBSERVER.......',self.obs, spice.bodc2n(int(self.obs)))
    print('INSTRUMENT.....',self.ins, spice.bodc2n(int(self.ins)))
    print(spice.frmnam(int(self.ins)))   
    print("spice.furnsh()...........SPICE Kernels...............loaded.\n")

  @staticmethod
  def kclear():
    spice.kclear()

  def load_time(self,utctime):#pds_image=None,time=False):
    '''
       Load time stamp.
       
       time in ISOT format.
    '''
    #from source.support.pds.core.common import open_pds
    #from source.support.pds.core.parser import Parser
    from astropy.time import Time
    
    # load header
    #if pds_image:
    #  parser = Parser()
    #  header = parser.parse(open_pds(pds_image))
      
    # Set time
    print()
    self.utctime = Time(utctime, format='isot', scale='utc')
    print('UTC (PDS Image)........',self.utctime)
    print('JD2000 (astropy).......',self.utctime.jd)


    jd2000 = 2451545e0
    #et = (self.utctime.jd-jd2000)*86400e0
    et = spice.str2et(self.utctime.isot)
    print('ET (SPICE)......',et)
    print('JD2000 (SPICE)......',jd2000+et/86400e0)

    slck = spice.sce2s(int(self.obs), et)
    print('S/C clock ticks.....',slck)

    #calet = spice.etcal(et)
    #print(calet)
    calet = spice.timout(et, 'YYYY-MON-DDTHR:MN:SC ::TDB')
    print(calet)

    print('VISIBILITY..........',spice.fovtrg(self.ins,self.body,'POINT',' ','LT+S',self.obs,et),'\n')

    self.et = et
    self.slck = slck

  def load_clock(self, slck):
    
    et = spice.scs2e(int(self.obs), str(slck))
    print('SLCK clock.....',slck)
    print('ET (SPICE)......',et)

    calet = spice.timout(et, 'YYYY-MON-DDTHR:MN:SC ::TDB')
    print(calet)
    #print('VISIBILITY..........',spice.fovtrg(self.ins,self.body,'POINT',' ','LT+S',self.obs,et),'\n')

    self.slck = slck
    self.et = et

  def local_time_ellipsoid(self):
    '''
      1. Calculate subsolar latitude & longitude
      2. Local time in respect to longitude=0
      
      Output
      ======
      subsolar latitude & longitude, greenwich local time, rotational direction
    '''
    from numpy import degrees, radians, sign
    assert hasattr(self, 'et')
    
    subsolar_pos0, et0, subsc_pos0 = spice.subslr('INTERCEPT/ELLIPSOID',self.body,self.et,self.body_frame,'LT+S',self.obs)
    subsolar_pos1, et1, subsc_pos1 = spice.subslr('INTERCEPT/ELLIPSOID',self.body,self.et+1e0,self.body_frame,'LT+S',self.obs)

    ss_r, ss_lat, ss_lon = spice.recsph(subsolar_pos1)
    localtime1 = -degrees(ss_lon)/15e0
    ss_r, ss_lat, ss_lon = spice.recsph(subsolar_pos1)
    localtime0 = -degrees(ss_lon)/15e0

    if sign(localtime0) == -1.0: localtime0=24e0+localtime0
    
    return (ss_r, 90-degrees(ss_lat), -degrees(ss_lon)), localtime0, sign(localtime1-localtime0)

  def solar_coord(self, ref_frame='J2000'):
    '''
       1. Solar XYZ (in KM)
       2. Solar RA and DEC in respect to the resting frame of the body. 
    '''
    from numpy import degrees, radians
    # positions [X, Y, Z]
    sun_target, ltime1 = spice.spkpos( 'SUN', self.et, ref_frame, 'LT+S', self.body)
    sun_ra_dec = spice.recsph(sun_target)
    #print('Solar RA and DEC.........',degrees(sun_ra_dec[1:]))
    return sun_target, sun_ra_dec

  def sc_coord(self, ref_frame='J2000'):
    '''
       1. S/C XYZ (in KM)
       2. S/C RA and DEC in respect to the resting frame of the body. 
    '''
    from numpy import degrees, radians
    # positions [X, Y, Z]
    sc_target, ltime1 = spice.spkpos(self.body, self.et, ref_frame, 'LT+S', self.obs)
    sc_ra_dec = spice.recsph(sc_target)[1:]
    #print('S/C DEC and RA.........',degrees(sc_ra_dec))
    return sc_target, sc_ra_dec

  def geo_ellipsoid(self, surf=[0,0,0]):
    '''
       As a rough ellipsoid, we compute the phase, emergence and incidence angles.
    '''
    from numpy import degrees, radians
    angles = spice.ilumin('ELLIPSOID', self.body, self.et, self.body_frame, 'LT+S', self.obs, surf)
    print('BODY-CENTER PHASE ANGLE ', degrees(angles[2]))
    
    return angles[2], angles[3], angles[4]

  def instrument_frame(self):
    '''
        Reconstructing the boresight frame from instrument kernels.
        
        1. Get Field-of-View;
        2. Compute boresight in spacecraft-frame;
        3. Compute the Instrument Tranformation Matrix in body-frame;
        
        Output
        ======
        Field-of-View Matrix, 
        Camera pointing vector,
        Camera X-axis,
        Instrument Transformation Matrix
    '''
    from numpy import float32, int32, array, degrees, radians, fabs, empty

    FOV = spice.getfov(int(self.ins), 4, shapelen=256, framelen=256)
    qua = FOV[4] #quaternion

    distance, ltime1 = spice.spkpos(self.ins_frame, self.et, self.body_frame, 'LT+S', self.obs)
    #print('SC on target body-reference',sc_target)
    if FOV[0] == "RECTANGLE":
      boresight = qua.sum(0)/4e0 #FOV[2]
    elif FOV[0] == "CIRCLE":
      boresight = qua[0]

    #R = spice.pxform(self.ins_frame,self.obs_frame, self.et) # Instrument to S/C
    #boresight_in_sc = spice.mxv(R, boresight) # Camera pointing vector in S/C frame

    #Instrument Transformation Matrix with the correction for the light time.
    M_ins = spice.pxform(self.body_frame, self.ins_frame, self.et-ltime1) # target to camera
    #M_targ = spice.pxform(self.ins_frame, self.body_frame, self.et-ltime1) # camera to target      
    #boresight_in_targ = spice.mxv(M_targ, boresight) # Camera pointing vector in target frame
    
    # Camera Euler Angles in body-frame
    #ins_euler = spice.m2eul(M, 3, 1, 3)
    #ra   = ins_euler[2]  
    #dec  = ins_euler[1]
    #roll = ins_euler[0]
    #print('RA',degrees(ra),'DEC',degrees(dec),'ROLL',degrees(roll))

    print('Camera Instrument Rotation Matrix')
    print(M_ins)
    print('Boresight in Instrument')
    print(boresight)
    return FOV, M_ins, boresight

  def to_oasis(self):
    from numpy import float64, array, radians, degrees
    from math import pi, sqrt

    # positions [X, Y, Z]
    sc_pos, ltime1 = spice.spkpos( 'SUN', self.et, 'J2000', 'LT+S', self.obs)
    sc_target, ltime2 = spice.spkpos(self.body, self.et, 'J2000', 'LT+S', self.obs)

    # Light-aberration correction to the ET timestamp
    calet = spice.timout(self.et+ltime1, 'YYYY-MON-DDTHR:MN:SC ::TDB')
    #print(calet)


    # Instrument
    m_ins = spice.pxform('J2000',self.ins_frame,self.et-ltime2)
    #print("Instrument Rotation Matrix")
    #print(m_ins) 
    ins_euler = spice.m2eul(m_ins, 3, 1, 3)
    #print("Instrument Euler Angles")
    #print(ins_euler)

    ra   = ins_euler[2]  -pi/2e0  +2e0*pi
    dec  = -ins_euler[1] +pi/2e0 
    roll = -ins_euler[0] +pi/2e0 #+ drol
    #ra  = ra-pi
    #dec = -dec

    #print('RA',degrees(ra),'DEC',degrees(dec),'ROLL',degrees(roll))

    # Body Euler Angles
    #spice.ckgp(int(self.body), self.ticks, 1e3, self.body_frame)
    m_body = spice.pxform('J2000',self.body_frame,self.et-ltime2)
    #print("Body Rotation Matrix")
    #print(m_body) 
    body_euler = degrees(spice.m2eul(m_body, 3, 1, 3))
    #print("Body Euler Angles")
    #print(body_euler)

    # To Astronomical Unity
    sc_pos_au     = array([-spice.convrt(x, 'KM', 'AU') for x in sc_pos]    ,dtype=float64)
    sc_target_au  = array([-spice.convrt(x, 'KM', 'AU') for x in sc_target] ,dtype=float64)

    # Target_to_Sun = -Sc_to_Sun +Sc_to_Target
    # the correct light time calculation is:
    # TARGET ---> S/C (light time 1) & S/C ---> SUN (light time 2)   
    target_pos_au = -sc_target_au +sc_pos_au

    #print('S/C position vector (AU): ')
    #print(sc_pos_au)
    #print('target position vector (AU): ')
    #print(target_pos_au)
    #print('\n')

    #dist_au = spice.convrt( sqrt(sum(target_pos*target_pos)), 'KM', 'AU' )
    #print('Helio distance (UA)',dist_au)
   
    return sc_pos_au, target_pos_au, body_euler, degrees([ra,dec,roll])

# END
