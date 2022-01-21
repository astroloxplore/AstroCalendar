#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import unicode_literals

'''
PROGRESSIONS
This program (calendar) calculates on the sky transits for a period of 
given time (default for next 2 months). 

The natal of a person is used to mark the repeated interactions.

Transiting Moon's aspects are calculated and plotted first. 
Other transiting planets' aspects are first stored, and then plotted.


This is a working variant that has following functions:
+ moon signs and aspects
+ aspects of transit-transit
+ transiting aspects that are present in natal are marked with thicker lines
+ Lilith is implemented
+ natal and transit orbs implemented

USE:
python ../my_astro_program/transit_na_nebe.py -i me1.zbs 

TO DO: 
1) dubliaz dni mesiaca i polnyh let cheloveka +
2) smena znaka planetoj +
3) smena doma planetoj + 
4) smena napravlenija (ot ME i starshe) + needs checking
5) progressivnyje ASC i MC: these two make a conjunctions with every planet in the chart every year, so no need to plot

'''


from optparse import OptionParser, OptionGroup
import math
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

import datetime

import swisseph as swe
import util

import matplotlib


def parse_commandline():
    parser = OptionParser()

    # in file
    parser.add_option("-i", "--input",  dest="input_file", metavar="<file>",\
    help="input file: records in zbs format")
    parser.set_defaults(input_file='record.zbs')
    
    # flags to select what aspects to calculate
    parser.add_option("-a", "--asp",  dest="asp_set", type="int",\
    help="Aspects: 0 for default set [0, 60, 90, 120, 180], 1 for major [0, 30, 45, 60, 90, 120, 135, 150, 180], 2 for all [0.0, 30.0, 36.0, 40.0, 45.0, 60.0, 72.0, 80.0, 90.0, 100.0, 108.0, 120.0, 135.0, 144.0, 150.0, 180.0], 3 for my default set [0, 30, 60, 90, 120, 150, 180]")
    parser.set_defaults(asp_set=0)
    
    parser.add_option("-t", "--time",  dest="time_interval", type="int",\
    help="A period of time (in days)")
    parser.set_defaults(time_interval=60) 
    
    # flags to select orbs
    parser.add_option("-b", "--orb",  dest="orbis", type="float",\
    help="default is 1 (100%), to use wider orbs: 1.1 - 10% wider, 1.5 - 50% wider and so on")
    parser.set_defaults(orbis=1)
    
    parser.add_option("-e", "--expand",  dest="transorbis", type="float",\
    help="default is 1 (100%), to use wider orbs: 1.1 - 10% wider, 1.5 - 50% wider and so on")
    parser.set_defaults(transorbis=1)
    
    # house system
    parser.add_option("-s", "--housesystem",  dest="house_sys", type="str",\
    help="Select house system: P for Placidus, K for Koch")
    parser.set_defaults(house_sys="P")
        
    
    # choose ruling system
    parser.add_option("-r", "--rule",  dest="ruling", type="str",\
    help="ruling systems: f - full, s - septener, d - day, n - night, m - my style")
    parser.set_defaults(ruling="m")
    
    parser.add_option("-o", "--output",  dest="output_file", metavar="<file>",\
    help="name for output graph")
    parser.set_defaults(output_file='progressions')
    
    # get the options:
    (options, args) = parser.parse_args()

    # check for any leftover command line arguments:
    if len(args):
        warning("Ignoring additional arguments "+str(args))
    
    # clean up:
    del(parser)
    return options



# ---------------- EPHEMERIS PARAMS ---------------------

SE_JUL_CAL = 0
SE_GREG_CAL = 1

#planet numbers for the ipl parameter in swe_calc()
SE_ECL_NUT = -1      

SE_SUN = 0       
SE_MOON = 1       
SE_MERCURY = 2       
SE_VENUS = 3       
SE_MARS = 4       
SE_JUPITER = 5       
SE_SATURN = 6       
SE_URANUS = 7       
SE_NEPTUNE = 8       
SE_PLUTO = 9       
SE_MEAN_NODE = 10      
SE_TRUE_NODE = 11
SE_MEAN_APOG = 12      
SE_OSCU_APOG = 13    
SE_EARTH = 14      
SE_CHIRON = 15      
SE_PHOLUS = 16      
SE_CERES = 17      
SE_PALLAS = 18      
SE_JUNO = 19      
SE_VESTA = 20      
SE_INTP_APOG = 21      
SE_INTP_PERG = 22    

SE_NPLANETS = 23      

SE_AST_OFFSET = 10000
SE_VARUNA = (SE_AST_OFFSET + 20000)

SE_FICT_OFFSET = 40
SE_FICT_OFFSET_1 = 39
SE_FICT_MAX = 999 
SE_NFICT_ELEM = 15

SE_COMET_OFFSET = 1000

SE_NALL_NAT_POINTS = (SE_NPLANETS + SE_NFICT_ELEM)

SE_ASC = 0
SE_MC = 1
SE_ARMC = 2
SE_VERTEX = 3
SE_EQUASC = 4    # "equatorial ascendant"
SE_COASC1 = 5    # "co-ascendant" (W. Koch)
SE_COASC2 = 6    # "co-ascendant" (M. Munkasey)
SE_POLASC = 7    # "polar ascendant" (M. Munkasey)
SE_NASCMC = 8


# flag bits for parameter iflag in function swe_calc()
# The flag bits are defined in such a way that iflag = 0 delivers what one
# usually wants:
#    - the default ephemeris (SWISS EPHEMERIS) is used,
#    - apparent geocentric positions referring to the true equinox of date
#      are returned.
# If not only coordinates, but also speed values are required, use 
# flag = SEFLG_SPEED.
#
# The 'L' behind the number indicates that 32-bit integers (Long) are used.

SEFLG_JPLEPH = 1       # use JPL ephemeris
SEFLG_SWIEPH = 2       # use SWISSEPH ephemeris
SEFLG_MOSEPH = 4       # use Moshier ephemeris

SEFLG_HELCTR = 8          # return heliocentric position
SEFLG_TRUEPOS = 16         # return true positions, not apparent
SEFLG_J2000 = 32             # no precession, i.e. give J2000 equinox
SEFLG_NONUT = 64             # no nutation, i.e. mean equinox of date
SEFLG_SPEED3 = 128     # speed from 3 positions (do not use it, SEFLG_SPEED is faster and more precise.)
SEFLG_SPEED = 256         # high precision speed
SEFLG_NOGDEFL = 512     # turn off gravitational deflection
SEFLG_NOABERR =1024    # turn off 'annual' aberration of light
SEFLG_EQUATORIAL = (2*1024)    # equatorial positions are wanted
SEFLG_XYZ = (4*1024)    # cartesian, not polar, coordinates
SEFLG_RADIANS = (8*1024)    # coordinates in radians, not degrees
SEFLG_BARYCTR = (16*1024)   # barycentric positions
SEFLG_TOPOCTR = (32*1024)   # topocentric positions
SEFLG_SIDEREAL = (64*1024)   # sidereal positions
SEFLG_ICRS = (128*1024)  # ICRS (DE406 reference frame)

SE_SIDBITS = 256
# for projection onto ecliptic of t0
SE_SIDBIT_ECL_T0 = 256
# for projection onto solar system plane 
SE_SIDBIT_SSY_PLANE = 512

# default ephemeris used when no ephemeris flagbit is set
SEFLG_DEFAULTEPH = SEFLG_SWIEPH




# ----------- houses ------------
# default house system is Placidus
hsystems = ('P', 'K', 'R', 'C', 'E', 'W', 'X', 'M', 'H', 'T', 'B', 'O')
hsysnames = {'P': 'Placidus', 'K': 'Koch', 'R': 'Regiomontanus', \
'C':'Campanus', 'E': 'Equal', 'W': 'WholeSign', 'X': 'Meriadian(Axial)',\
'M': 'Morinus', 'H': 'Horizontal', 'T': 'Topocentric(Page-Polich)',\
'B': 'Alchabitius', 'O': 'Porphyry'}
    
def_hsys = 'P'


# for houses
ASC, MC, ARMC, VERTEX, EQUASC, COASC, COASC2, POLARASC = range(0, 8)


# aspects
# major: 0, 30, 60, 90, 120, 150, 180
# carmic: 20, 40, 80, 100
# artistic: 36, 72, 108, 144
# additional: 45, 135
#Aspects = [0.0, 20.0, 30.0, 36.0, 40.0, 45.0, 60.0, 72.0, 80.0, 90.0, \
#100.0, 108.0, 120.0, 135.0, 144.0, 150.0, 180.0]
#ASPECT_NUM = 17    # we need all of them, but for now aspmatrix
                        # has only orbs for 11
    
#Aspects = [0.0, 30.0, 45.0, 60.0, 72.0, 90.0, 120.0, 135.0, 144.0, 150.0, 180.0]
ASPECT_NUM = 11    
Aspects = [0, 60, 90, 120, 180]    


# Natal orbs
orbs = {'SU': {0.0:10.0,30.0:1.5,60.0:6.0,90.0:8.0,120.0:8.0,150.0:1.0,180.0:10.0}, \
'MO': {0.0:10.0,30.0:1.5,60.0:6.0,90.0:8.0,120.0:8.0,150.0:1.0,180.0:10.0}, \
'ME': {0.0:8.0,30.0:1.0,60.0:5.0,90.0:6.0,120.0:6.0,150.0:1.0,180.0:8.0}, \
'VE': {0.0:8.0,30.0:1.0,60.0:5.0,90.0:6.0,120.0:6.0,150.0:1.0,180.0:8.0}, \
'MA': {0.0:8.0,30.0:1.0,60.0:5.0,90.0:6.0,120.0:6.0,150.0:1.0,180.0:8.0}, \
'JU': {0.0:6.0,30.0:1.0,60.0:4.0,90.0:5.0,120.0:5.0,150.0:1.0,180.0:6.0}, \
'SA': {0.0:6.0,30.0:1.0,60.0:4.0,90.0:5.0,120.0:5.0,150.0:1.0,180.0:6.0}, \
'UR': {0.0:5.0,30.0:1.0,60.0:4.0,90.0:4.0,120.0:4.0,150.0:1.0,180.0:5.0}, \
'NE': {0.0:5.0,30.0:1.0,60.0:4.0,90.0:4.0,120.0:4.0,150.0:1.0,180.0:5.0}, \
'PL': {0.0:5.0,30.0:1.0,60.0:4.0,90.0:4.0,120.0:4.0,150.0:1.0,180.0:5.0} }


orbs2 = {0: {0.0:10.0,30.0:1.5,60.0:6.0,90.0:8.0,120.0:8.0,150.0:1.0,180.0:10.0}, \
1: {0.0:10.0,30.0:1.5,60.0:6.0,90.0:8.0,120.0:8.0,150.0:1.0,180.0:10.0}, \
2: {0.0:8.0,30.0:1.0,60.0:5.0,90.0:6.0,120.0:6.0,150.0:1.0,180.0:8.0}, \
3: {0.0:8.0,30.0:1.0,60.0:5.0,90.0:6.0,120.0:6.0,150.0:1.0,180.0:8.0}, \
4: {0.0:8.0,30.0:1.0,60.0:5.0,90.0:6.0,120.0:6.0,150.0:1.0,180.0:8.0}, \
5: {0.0:6.0,30.0:1.0,60.0:4.0,90.0:5.0,120.0:5.0,150.0:1.0,180.0:6.0}, \
6: {0.0:6.0,30.0:1.0,60.0:4.0,90.0:5.0,120.0:5.0,150.0:1.0,180.0:6.0}, \
7: {0.0:5.0,30.0:1.0,60.0:4.0,90.0:4.0,120.0:4.0,150.0:1.0,180.0:5.0}, \
8: {0.0:5.0,30.0:1.0,60.0:4.0,90.0:4.0,120.0:4.0,150.0:1.0,180.0:5.0}, \
9: {0.0:5.0,30.0:1.0,60.0:4.0,90.0:4.0,120.0:4.0,150.0:1.0,180.0:5.0} }

# Transit orb = 1 for all (default)

progress_orb = {'SU':2,'MO':1.5}
        
# -------- planets rule signs --------
# full        
planet_rules_signs = {'Sun': [120], 'Moon': [90], \
'Mars': [0, 210], 'Venus': [30, 180], \
'Mercury': [60, 150], 'Jupiter': [240, 330],\
'Saturn': [270, 300], 'Uranus': [300, 270], \
'Neptune': [330, 240], 'Pluto': [210, 0], 'true Node': [], 'mean Node': [], 'mean Apogee': [], 'osc. Apogee': []}

# septener only
planet_rules_signs_sep = {'Sun': [120], 'Moon': [90], \
'Mars': [0, 210], 'Venus': [30, 180], \
'Mercury': [60, 150], 'Jupiter': [240, 330],\
'Saturn': [270, 300], 'Uranus': [], \
'Neptune': [], 'Pluto': [], 'true Node': [], 'mean Node': [], 'mean Apogee': [], 'osc. Apogee': []}

# one sign strictly starting from Jupiter # day birth
planet_rules_one_sign_day = {'Sun': [120], 'Moon': [90], \
'Mars': [0], 'Venus': [30, 180], \
'Mercury': [60, 150], 'Jupiter': [240],\
'Saturn': [270], 'Uranus': [300], \
'Neptune': [330], 'Pluto': [210], 'true Node': [], 'mean Node': [], 'mean Apogee': [], 'osc. Apogee': []}

# one sign strictly starting from Jupiter # night birth
planet_rules_one_sign_night = {'Sun': [120], 'Moon': [90], \
'Mars': [210], 'Venus': [30, 180], \
'Mercury': [60, 150], 'Jupiter': [330],\
'Saturn': [300], 'Uranus': [270], \
'Neptune': [240], 'Pluto': [0], 'true Node': [], 'mean Node': [], 'mean Apogee': [], 'osc. Apogee': []}

# minors for both signs, vyshki dlia odnogo znaka
planet_rules_my_style = {'Sun': [120], 'Moon': [90], \
'Mars': [0,210], 'Venus': [30, 180], \
'Mercury': [60, 150], 'Jupiter': [240,330],\
'Saturn': [270,300], 'Uranus': [300], \
'Neptune': [330], 'Pluto': [210], 'true Node': [], 'mean Node': [], 'mean Apogee': [], 'osc. Apogee': []}







# -------------- shifts and encodings ---------------------
uni_signs = {'Ari': u'\u2648', 'Tau': u'\u2649', 'Gem': u'\u264A', 'Can': u'\u264B', \
'Leo': u'\u264C', 'Vir': u'\u264D', 'Lib': u'\u264E', 'Sco': u'\u264F',\
'Sag': u'\u2650', 'Cap': u'\u2651', 'Aqu': u'\u2652', 'Pis': u'\u2653'}

uni_pls = {'Sun': u'\u2609', 'Moon': u'\u263D', 'Mercury': u'\u263F', 'Venus': u'\u2640', 'Mars': u'\u2642', \
'Jupiter': u'\u2643', 'Saturn': u'\u2644', 'Uranus': u'\u2645', 'Neptune': u'\u2646', 'Pluto': u'\u2647', \
'true Node': u'\u260A', 'mean Node': u'\u260A', 'mean Apogee':u'\u26B8', 'Chiron':u'\u26B7'}

uni_pls_2 = {'SU': u'\u2609', 'MO': u'\u263D', 'ME': u'\u263F', 'VE': u'\u2640', 'MA': u'\u2642', \
'JU': u'\u2643', 'SA': u'\u2644', 'UR': u'\u2645', 'NE': u'\u2646', 'PL': u'\u2647', \
'NN': u'\u260A', 'Li': u'\u26B8', 'Ch': u'\u26B7'}

#uni_asps = {0: u'\u260C', 60: u"\u26B9", 90: u'\u25A1', 120: u'\u25B3', 180: u'\u260D'}

uni_asps = {0: u'\u260C', 60: u'\u2731', 90: u'\u25A1', 120: u'\u25B3', 180: u'\u260D',\
20: '20'+u'\u00B0', 30: u'\u26BA', 36: '36'+u'\u00B0', 40: '40'+u'\u00B0', 45: '45'+u'\u00B0', 72: '72'+u'\u00B0', 80: '80'+u'\u00B0', \
100: '100'+u'\u00B0', 108: '108'+u'\u00B0', 135: '135'+u'\u00B0', 144: '144'+u'\u00B0', 150: u'\u26BB'}

    
    
# --------additional--------------
    
#calendars
GREGORIAN, JULIAN = 0, 1    # Gregorian is default

#times
ZONE = 0
GREENWICH = 1
LOCALMEAN = 2
LOCALAPPARENT = 3
    
HOURSPERDAY = 24.0

swisseph_pla_to_int = {'Sun':0,'Moon':1,'Mercury':2,'Venus':3,'Mars':4,\
    'Jupiter':5,'Saturn':6,'Uranus':7,'Neptune':8,'Pluto':9, 'Node':11, 'Lilith':12}

#list_of_houses = [1,2,3,4,5,6,7,8,9,10,11,12]

list_of_houses = ['1','2','3','4','5','6','7','8','9','10','11','12']

pl_to_pl2 = {'Sun':'SU','Moon':'MO','Mercury':'ME','Venus':'VE','Mars':'MA',\
    'Jupiter':'JU','Saturn':'SA','Uranus':'UR','Neptune':'NE','Pluto':'PL', 'Node':'NN'}




# make a flag what to return
def get_houses(tjd_ut, flag, geolat, geolon, hsys, obl):
    """Calculates the cusps of the Houses"""

    HOUSE_NUM = 12
    ASC, MC, ARMC, VERTEX, EQUASC, COASC, COASC2, POLARASC = range(0, 8)
    LON, LAT, RA, DECL = 0, 1, 2, 3
    obl = obl

    cusps, ascmc = swe.houses_ex(tjd_ut, geolat, geolon, bytes(hsys, 'utf-8'), flag)
    #print self.cusps
    ascra, ascdecl, dist = swe.cotrans(ascmc[ASC], 0.0, 1.0, -obl)             
    mcra, mcdecl, dist = swe.cotrans(ascmc[MC], 0.0, 1.0, -obl)
    ascmc2 = ((ascmc[ASC], 0.0, ascra, ascdecl), (ascmc[MC], 0.0, mcra, mcdecl))
       
    #zdAsc=90.0, zdMC=0.0
    #poleAsc=lat, poleMC=0.0
    qasc = math.degrees(math.asin(math.tan(math.radians(ascdecl))*math.tan(math.radians(geolat))))
    regioMPAsc = ascra-qasc
    regioMPMC = mcra

    cuspstmp = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
    for i in range(HOUSE_NUM):
       cuspstmp[i][0], cuspstmp[i][1], dist = swe.cotrans(cusps[i], 0.0, dist, -obl)
          
    cusps2 = ((cuspstmp[0][0], cuspstmp[0][1]), (cuspstmp[1][0], cuspstmp[1][1]), (cuspstmp[2][0], cuspstmp[2][1]), (cuspstmp[3][0], cuspstmp[3][1]), (cuspstmp[4][0], cuspstmp[4][1]), (cuspstmp[5][0], cuspstmp[5][1]), (cuspstmp[6][0], cuspstmp[6][1]), (cuspstmp[7][0], cuspstmp[7][1]), (cuspstmp[8][0], cuspstmp[8][1]), (cuspstmp[9][0], cuspstmp[9][1]), (cuspstmp[10][0], cuspstmp[10][1]), (cuspstmp[11][0], cuspstmp[11][1]))

    return cusps, ascmc, ascmc2





    
   # ------------ GET PLANET POSITION -----------    
    
def get_planet(tjd_ut, pId, flag, lat, ascmc2, raequasc, nolat, obl):
    """Data of a Planet"""

    #Speculum #Common
    LONG, LAT, RA, DECL = 0, 1, 2, 3

    #Placidus
    ADLAT, SA, MD, HD, TH, HOD, PMP, ADPH, POH, AODO = 4, 5, 6, 7, 8, 9, 10, 11, 12, 13

    #Regiomontanian/Campanian - for Primary Directions 
    RMD, RHD, ZD, POLE, Q, W, CMP, RMP, AZM, ELV  = 4, 5, 6, 7, 8, 9, 10, 11, 12, 13

    #data[x]
    DIST, SPLON, SPLAT, SPDIST = 2, 3, 4, 5
    
    #dataEqu
    RAEQU, DECLEQU, DISTEQU, SPRAEQU, SPDECLEQU, SPDISTEQU = 0, 1, 2, 3, 4, 5

    speculums = None
    computePlacidianSpeculum = False

    #if (ecl == None):
    data, rflag = swe.calc_ut(tjd_ut, pId, flag)
    dataEqu, rflag = swe.calc_ut(tjd_ut, pId, flag+SEFLG_EQUATORIAL)
    # data : [longitude, latitude, distance, speed in long, speed in lat, speed in dist]
    name = swe.get_planet_name(pId)
       
    '''else:
       data = tuple(ecl)
       dataEqu = tuple(equ)
       name = 'DescNode'''

    if nolat:
       #print 'no lat given'
       # data : [longitude, latitude, distance, speed in long, speed in lat, speed in dist]
       data[1] = 0.0
       ra, decl, dist = swe.cotrans(data[LONG], 0.0, 1.0, -obl)
       dataEqu = (ra, decl, dataEqu[DISTEQU], dataEqu[SPRAEQU], dataEqu[SPDECLEQU], dataEqu[SPDISTEQU])

    if lat != None:
       #print 'use some other files'
       #placspec.py and regiospec should be used instead, remove these!
       speculums = []
       computePlacidianSpeculum = True


    if computePlacidianSpeculum:
       #print "compute placidian speculum"
       ramc = ascmc2[MC][RA]
       raic = ramc+180.0
       if raic > 360.0:
          raic -= 360.0

       eastern = True
       if ramc > raic:
          if dataEqu[RAEQU] > raic and dataEqu[RAEQU] < ramc:
             eastern = False
       else:
          if (dataEqu[RAEQU] > raic and dataEqu[RAEQU] < 360.0) or (dataEqu[RAEQU] < ramc and dataEqu[RAEQU] > 0.0):
             eastern = False

       #adlat
       adlat = 0.0
       val = math.tan(math.radians(lat))*math.tan(math.radians(dataEqu[DECLEQU]))
       if math.fabs(val) <= 1.0:
          adlat = math.degrees(math.asin(val))

       #md
       med = math.fabs(ramc-dataEqu[RAEQU])

       if med > 180.0:
          med = 360.0-med
       icd = math.fabs(raic-dataEqu[RAEQU])
       if icd > 180.0:
          icd = 360.0-icd

       md = med

       #hd
       aoasc = ramc+90.0
       if aoasc >= 360.0:
          aoasc -= 360.0

       dodesc = raic+90.0
       if dodesc >= 360.0:
          dodesc -= 360.0

       aohd = dataEqu[RAEQU]-adlat
       hdasc = aohd-aoasc
       if hdasc < 0.0:
          hdasc *= -1
       if hdasc > 180.0:
          hdasc = 360.0-hdasc 

       dohd = dataEqu[RAEQU]+adlat
       hddesc = dohd-dodesc
       if hddesc < 0.0:
          hddesc *= -1
       if hddesc > 180.0:
          hddesc = 360.0-hddesc 

       hd = hdasc
       if hddesc < hdasc:
          hd = hddesc
          hd *= -1

       #sa (southern hemisphere!?)
       dsa = 90.0+adlat
       nsa = 90.0-adlat

       abovehorizon = True
       if med > dsa:
          abovehorizon = False

       sa = dsa
       if not abovehorizon:
          sa = -nsa #nocturnal if negative
          md = icd
          md *= -1

       #TH(TemporalHour)
       th = sa/6.0

       #HOD(HourlyDistance)
       hod = 0.0
       if th != 0.0:
          hod = md/math.fabs(th)

       #pmp
       pmp = 0.0
       tmd = md
       if tmd < 0.0:
          tmd *= -1

       pmpsa = sa
       if pmpsa < 0.0:
          pmpsa *= -1

       if not abovehorizon and eastern:
          pmp = 90.0-90.0*(tmd/pmpsa)
       elif not abovehorizon and not eastern:
          pmp = 90.0+90.0*(tmd/pmpsa)
       elif abovehorizon and not eastern:
          pmp = 270.0-90.0*(tmd/pmpsa)
       elif abovehorizon and eastern:
          pmp = 270.0+90.0*(tmd/pmpsa)

       #adphi
       tval = math.fabs(sa)
       adphi = 0.0
       if tval != 0.0:
          adphi = math.fabs(tmd)*adlat/tval

       #phi
       tval = math.tan(math.radians(dataEqu[DECLEQU]))
       phi = 0.0
       if tval != 0.0:
          phi = math.degrees(math.atan(math.sin(math.radians(adphi))/tval))

       #ao/do (southern hemisphere!?)
       if eastern:
          ao = dataEqu[RAEQU]-adphi
       else:
          ao = dataEqu[RAEQU]+adphi
          ao *= -1 #do if negative

       speculums.append((data[LONG], data[LAT], dataEqu[RAEQU], dataEqu[DECLEQU], adlat, sa, md, hd, th, hod, pmp, adphi, phi, ao))
       #speculums.append((dataEqu[LONG], dataEqu[LAT], dataEqu[RAEQU], dataEqu[DECLEQU], adlat, sa, md, hd, th, hod, pmp, adphi, phi, ao))
       #print data[LONG]
       #print dataEqu[LONG]

       return speculums

    

def calc_planet_houses(pl_poss, cusps, ruling_mode):
    planets_in_houses = {}
    
    #planet_list = ['Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto', 'NNode', 'Lilith']
    
    # add first house at the end
    #cusps = cusps + (cusps[0],) # tuple
    cusps.append(cusps[0])
    
    for pl in pl_poss:
        if isinstance(pl, int):
            pos = pl_poss[pl]
            planet_house = []
            for i in range(0,12):
                            cusp1, cusp2 = cusps[i], cusps[i+1]
                        # check for circle beginning
                            if cusp1 > cusp2: # cusp1 in pisces
                                    if pos < cusp2: # planet in Aries
                                            planet_house.append(i+1)
                                    elif pos > cusp1: # planet in Pisces
                                            planet_house.append(i+1)
                            else:
                                    if cusp1 <= pos and pos < cusp2: 
                                            #print pos, i+1
                                            planet_house.append(i+1)
                            
                                                                            
            # do check for house 13th aka 1st
            for h in planet_house:
                    if h==13:
                            planet_house[planet_house.index(h)]=1
            if len(planet_house)!=0:
                    planets_in_houses[swe.get_planet_name(pl)] = [planet_house[0]]
      
    
    # ruling system 
    if ruling_mode=='f':
            planet_rules = planet_rules_signs
    elif ruling_mode=='s':
            planet_rules = planet_rules_signs_sep
    elif ruling_mode=='d':
            planet_rules = planet_rules_one_sign_day
    elif ruling_mode=='n':
            planet_rules = planet_rules_one_sign_night
    elif ruling_mode=='m':
            planet_rules = planet_rules_my_style
    else:
            planet_rules = planet_rules_signs
        
    #print planet_rules    
            
    # ruling
    for p in planets_in_houses:
        
        if ruling_mode=='n': # check where SUN is
            if planets_in_houses['Sun'] in [7,8,9,10,11,12]:
                planet_rules = planet_rules_one_sign_day
        
        try:
            rul_signs = planet_rules[p] # 1 or 2 signs
        except KeyError:
            rul_signs = []
            
        rul_houses = []
        for s in rul_signs:
            # example Leo: 120, so s=120 for Sun
        
            for i in range (0,12):
                cusp1, cusp2 = cusps[i], cusps[i+1]
                
                
                # if cuspid gets into this sign
                if s <= cusp1 and cusp1 <= s+30:
                    rul_houses.append(cusps.index(cusp1)+1)
                # check included signs
                if cusp1 <= s and s+30 <= cusp2: 
                    #print 'vkliucenyj znak mezdu kuspidami', cusp1, cusp2
                    rul_houses.append(cusps.index(cusp1)+1)
                    
                    
                if cusp2 < cusp1: # the edge 360-0
                    # variant 1: cusp1 is in Aqu/Pis and cusp2 is in Ari
                    if cusp1 <= s and util.normalize(s+30) <= cusp2:
                        rul_houses.append(cusps.index(cusp1)+1) # Pisces included
                    # variant 2: cusp1 is in Pis and cusp2 is in Tau
                    if s <= cusp1 and util.normalize(s+30) <= cusp2:
                        rul_houses.append(cusps.index(cusp1)+1) # Aries included
                    
                    
        # remove duplicates 
        rul_houses = list(set(rul_houses))
        planets_in_houses[p].append(rul_houses)
    return planets_in_houses



def check_overflow(year, month, day, time):

    #check over/underflow
    if time >= HOURSPERDAY:
        time -= HOURSPERDAY
        year, month, day = util.incrDay(year, month, day)
        overflow = True
    elif time < 0.0:
        time += HOURSPERDAY
        year, month, day = util.decrDay(year, month, day)
        overflow = True
    else: 
        overflow = False

    return year, month, day, time

    

def find_asp(pl1, asp_exact):
            
    is_asp = False
    try:
        pl1_orbs = orbs2[pl1] # in orbs2 planets incoded as integers
    except KeyError:
        pl1_orbs = {0.0:1.0,30.0:1.0,60.0:1.0,90.0:1.0,120.0:1.0,150.0:1.0,180.0:1.0}
    asp_type = None
    diff_with_exact = None        
            
    for a in pl1_orbs: # for every type of aspect get ranges
        ranges = []
        try:
            orbis = pl1_orbs[a]
        except KeyError:
            orbis = 1.0
        if a==0:
            ranges = [a, a+orbis]
        elif a==180:
            ranges = [a-orbis, a]
        else:
            ranges = [a-orbis, a+orbis]
                    
        if ranges[0]<=asp_exact and asp_exact<=ranges[1]: 
            #print 'this is the aspect!'
            asp_type = a
            #how_exact = 1-abs(a-asp_exact)/orbis
            diff_with_exact = abs(a-asp_exact)
            is_asp = True
            break    
    
    return is_asp, asp_type, diff_with_exact


# function to find transit_asp
def find_tr_asp(asp_exact, trans_orbis, Aspects):
            
    is_asp = False
    #pl1_orbs = orbs[pl1]
    asp_type = None
    diff_with_exact = None        
            
    #for a in pl1_orbs: # for every type of aspect get ranges
    for a in Aspects:
        ranges = []
        #orbis = pl1_orbs[a]
                
        if a==0:
            ranges = [a, a+trans_orbis]
        elif a==180:
            ranges = [a-trans_orbis, a]
        else:
            ranges = [a-trans_orbis, a+trans_orbis]
                    
        if ranges[0]<=asp_exact and asp_exact<=ranges[1]: 
            #print 'this is the aspect!'
            asp_type = a
            #how_exact = 1-abs(a-asp_exact)/orbis
            diff_with_exact = abs(a-asp_exact)
            is_asp = True
            break    
    
    return is_asp, asp_type, diff_with_exact
    


def find_available_track(track_order,track_occupancy,from_when,till_when, otstup):
   
    ok_track = 6
    # ---------- find available track ----------                
    for f in track_order:
            coords = track_occupancy[f]
            
            if len(coords)==0:
                    ok_track = f
                    track_occupancy[f] = [from_when, till_when]
                    break
            else: 
                    no_overlap_on_the_track = True
                    # coords can have [x1,y1,x2,y2]
                    for x in range(0,len(coords),2):
                            min_coord, max_coord = coords[x], coords[x+1]
                            if from_when > max_coord+otstup or till_when < min_coord-otstup:
                                    pass
                            else:
                                    no_overlap_on_the_track = False
                            
                    if no_overlap_on_the_track:
                            ok_track = f
                            coords.append(from_when)
                            coords.append(till_when)
                            track_occupancy[ok_track] = coords
                            break
    return ok_track
    
    
    
    
def plot_aspect(plt, natal_asps, asp_coords, ok_track, min_days, max_days, otstup, k):
    
    from_when = min(asp_coords)
    till_when = max(asp_coords)
    
    x = [ok_track, ok_track]
    y = [from_when, till_when]
         
    tr_pl_name = k[0]
    n_pl_name = k[1]
    asp_type = k[2]
    try:
        text = k[3]
    except IndexError:
        text=''
    #print text
    
    text_position = from_when+(till_when-from_when)/2 # in the middle of an aspect line
    if text_position>max_days:
        text_position = max_days-otstup
    elif text_position<min_days:
        text_position = min_days+otstup
    try:
        plt.text(ok_track, text_position-otstup, uni_pls[tr_pl_name],  size=10)
    except KeyError:
        plt.text(ok_track, text_position-otstup, tr_pl_name[0],  size=10)
    
    
    if n_pl_name=='1':
        if asp_type==180:
            n_pl_name = '7'
            asp_type = 0
        plt.text(ok_track, text_position+otstup, n_pl_name+text,  size=10)
    elif n_pl_name=='10':
        if asp_type==180:
            n_pl_name = '4'
            asp_type = 0
        plt.text(ok_track, text_position+otstup, n_pl_name+text,  size=10)
    else:
        try:
            plt.text(ok_track, text_position+otstup, uni_pls[n_pl_name]+text,  size=10)
        except KeyError: #in case natal obj is a cusp
            plt.text(ok_track, text_position+otstup, n_pl_name+text,  size=10)
    
        
            
    if asp_type in [30,60,120]:
        col = 'g'
    elif asp_type in [90,150,180]:
        col = 'r'
    elif asp_type in [0]:
        col = 'b'
    else:
        col = 'k'    
    plt.text(ok_track, text_position, uni_asps[asp_type],  size=10)
    
    
    
    lwidth = 1
    if (tr_pl_name, n_pl_name) in natal_asps or (n_pl_name,tr_pl_name) in natal_asps:
        lwidth = 2
                                    
    plt.plot(x, y, color=col, linewidth=lwidth, alpha=0.9)
    
    
    
# -----------------------------------------------------------------------------
# ----------------------------------- MAIN PROGRAM ----------------------------
# -----------------------------------------------------------------------------
        
            
def main():
    
    matplotlib.rc('font', family='FreeSerif')
    
    # get command line options
    options = parse_commandline()

    # open the zbs file
    try:
        records = open(options.input_file).read().splitlines() 
        data = records[0].split(';')
        #s_transit = records[1].split(';')
        s_transits = records[1:]
    except IOError:
        print ("ERROR: cannot open or read input file:", input_file)
        exit(-1)
    
    
    if options.orbis:
        # natal orbs
        for planet in orbs:
            asp_orbs = orbs[planet]
            for asp in asp_orbs:
                new_orb = asp_orbs[asp]*options.orbis
                orbs[planet][asp] = new_orb
    
    trans_orbis = options.transorbis
    
    
    # ---------------- set ASPECTS set -----------------
    # default set Aspects = [0, 60, 90, 120, 180]

    if options.asp_set==1:
        Aspects = [0, 30, 45, 60, 90, 120, 135, 150, 180] # major
    elif options.asp_set==2:
        Aspects = [0.0, 30.0, 36.0, 40.0, 45.0, 60.0, 72.0, 80.0, 90.0, 100.0, 108.0, 120.0, 135.0, 144.0, 150.0, 180.0]
    elif options.asp_set==3:
        Aspects = [0, 30, 60, 90, 120, 150, 180]
    else:
        Aspects = [0, 60, 90, 120, 180]
    
    
    # ---------------- ruling system ---------------
    if options.ruling=='f':
        planet_rules = planet_rules_signs
    elif options.ruling=='s':
        planet_rules = planet_rules_signs_sep
    elif options.ruling=='d':
        planet_rules = planet_rules_one_sign_day
    elif options.ruling=='n':
        planet_rules = planet_rules_one_sign_night
    elif options.ruling=='m':
        planet_rules = planet_rules_my_style
    else:
        planet_rules = planet_rules_signs
      

    

    # --------------- Start: getting NATAL data --------------------
    name = data[0]
    
    # ------- date of birth --------
    date = data[1].split('.')
    year, month, day = int(date[2]), int(date[1]), int(date[0])

    # ------- time of birth --------
    time = data[2].split(':')
    hour = float(time[0])
    if len(time)>=2:
        minute = float(time[1])
    else:
        minute = 0.0
    if len(time)==3:
        second = float(time[2])
    else:
        second = 0.0
    
    # ------- place of birth --------
    place = data[4]
    zbslat, zbslon = data[5], data[6]
    
    east, north = False, False
    z_lon = re.findall('[a-zA-Z]+', zbslon)
    z_lat = re.findall('[a-zA-Z]+', zbslat)
                
    if "e" in z_lon or "E" in z_lon:
        east = True      
    if "n" in z_lat or "N" in z_lat:
        north = True            
        
    lat_list = re.findall('\d+', zbslat)
    lon_list = re.findall('\d+', zbslon)
    
    deglat, minlat = int(lat_list[0]), int(lat_list[1])
    deglon, minlon = int(lon_list[0]), int(lon_list[1])
    seclat, seclon = 0, 0

    #altitude = 100
        
    lat = deglat+minlat/60.0+seclat/3600.0
    lon = deglon+minlon/60.0+seclon/3600.0
        
    if not north:
        lat *= -1.0

    if not east:
        lon *= -1.0

    # 'cal' means calendar (gregorian=0 or julian=1)
    # zt is zonetime (ZONE = 0, GREENWICH = 1, LOCALMEAN = 2, LOCALAPPARENT = 3)
    # zh is zonehour, zm is zoneminute
    
    origyear = year
    origmonth = month
    origday = day
    bc = False # before christ
    cal = GREGORIAN # calendar
    zt = 0
    
    zone = data[3].strip()
    
    if re.findall('\:', zone):
        zone_list = zone.split(':')
        sign = zone_list[0][0]
        zh = int(zone_list[0][1:])
        zm = int(zone_list[1])
        zs = int(zone_list[2])
    else:
        sign = zone[0] # +/-
        zh = int(zone[1:])
        zm, zs = 0, 0
    
    if sign=='+':
        plus = True
    else:
        plus = False
    
    
    daylightsaving = False # default is False

    time = hour+minute/60.0+second/3600.0
    
    if daylightsaving:
        time -= 1.0
        dhour -= 1
    if time < 0.0:
        time += HOURSPERDAY
        year, month, day = util.decrDay(year, month, day)
        dhour += int(HOURSPERDAY)
        dyear, dmonth, dday = year, month, day
            
    if zt == 0:
        ztime = zh+zm/60.0
        if plus:
            time-=ztime
        else:
            time+=ztime
        
    if bc:
        year = 1-self.year
        
    #check over/underflow
    if time >= HOURSPERDAY:
        time -= HOURSPERDAY
        year, month, day = util.incrDay(year, month, day)
    elif time < 0.0:
            time += HOURSPERDAY
            year, month, day = util.decrDay(year, month, day)
    
    dyear, dmonth, dday, dhour, dmin, dsec = year, month, day, hour, minute, second    
    
    calflag = SE_GREG_CAL
    if cal == JULIAN:
        calflag = SE_JUL_CAL

    hsys = options.house_sys
    
    if hsys in hsystems:
        hsys = hsys
    else:
        hsys = hsystems[0]
    
    # additional stuff
    htype = 1 #natal
    notes = 'notes'
    nolat = False
    #swe_set_topo(place.lon, place.lat, place.altitude)
    #swe_set_topo(lon, lat, altitude)
    hflag = 0
    fsflag = 0
    pflag = SEFLG_SWIEPH+SEFLG_SPEED
    astflag = SEFLG_SWIEPH
        
    # --------------- End: getting NATAL data -----
    count_years = 0
    
    
    
    
    n_signs = {1: 'Ari', 2: 'Tau', 3: 'Gem', 4: 'Can', 5: 'Leo',\
6: 'Vir', 7: 'Lib', 8: 'Sco', 9: 'Sag', 10: 'Cap',\
11: 'Aqu', 12: 'Pis'}
    planets = [0,1,2,3,4,5,6,7,8,9,11,12,15]    # 11 is mean node, 12 lilith
    n_pls = {0: 'Sun', 1:'Moon', 2:'Mercury', 3:'Venus', 4:'Mars', 5:'Jupiter', 6:'Saturn', 7:'Uranus', 8:'Neptune', 9:'Pluto', 11:'Mean Node', 12:'Lilith', 15:'Chiron'}
        
    natal_asps = []



    # ------------ calculate NATAL ----------------------
    
    year, month, day, time = check_overflow(dyear, dmonth, dday, time)    
    #print 'year, month, day, stime', year, month, day, stime
                                    
    # -------------horoscope--------------
    jd = swe.julday(year, month, day, time, calflag) # julian date
    d = swe.deltat(jd)
    obl, rflag = swe.calc(jd+d, SE_ECL_NUT, 0)

    # -------------houses--------------
    cusps, ascmc, ascmc2 = get_houses(jd, hflag, lat, lon, hsys, obl[0])
    raequasc, declequasc, dist = swe.cotrans(ascmc[SE_EQUASC], 0.0, 1.0, -obl[0])
    options_meannode = True 
    #print('Natal cusps:',cusps, ascmc)       
    
    # -------------planets--------------       
    for pId in n_pls:
          pl_info = get_planet(jd, pId, pflag, lat, ascmc2, raequasc, nolat, obl[0])
          n_pls[pId] = pl_info
          #print pId, pl_info, '\n'
    
    # ---------- check degree overflow ----------
    for pId in n_pls:
       pl_info = n_pls[pId]
       pl_norm_pos = util.normalize(pl_info[0][0])
       n_pls[pId] = pl_norm_pos
    cusps = [util.normalize(c) for c in cusps]
    
    #print n_pls
    
    # ---------calculate natal aspects--------------
    for i in range(0,len(planets)-2):
            p1 = planets[i]
            
            for j in range(i+1,len(planets)):
                    p2 = planets[j]
                            
                    if p1!=p2:
                    
                            p1_pos = n_pls[p1]
                            p2_pos = n_pls[p2]
                    
                            if p1_pos>=p2_pos:
                                    diff = p1_pos - p2_pos
                                    if diff>=180: # further apart or at the boundary of pisces-aries
                                            diff = abs(p1_pos-360) + p2_pos
                            else:
                                    diff = p2_pos - p1_pos
                                    if diff>=180:
                                            diff = abs(p2_pos - 360) + p1_pos    
                    
                            is_asp, asp_type, diff_with_exact = find_asp(p1, diff)
                                    
                            if is_asp:
                                    if asp_type in Aspects:
                                        natal_asps.append((swe.get_planet_name(p1),swe.get_planet_name(p2)))

    
    
    
    
    
    # for every month that we put in
    for s in s_transits:
            s_transit = s.split(';')
            # --------------- Start: getting TRANSIT data --------------------
            #name = s_transit[0]
            
            # ------- date of starting transit --------
            tr_date = s_transit[1].split('.')
            tr_year, tr_month, tr_day = int(tr_date[2]), int(tr_date[1]), int(tr_date[0])

            # ------- time of first transit --------
            '''
            tr_time = s_transit[2].split(':')
            tr_hour = float(tr_time[0])
            if len(tr_time)>=2:
                    tr_minute = float(tr_time[1])
            else:
                    tr_minute = 0.0
            if len(tr_time)==3:
                    tr_second = float(tr_time[2])
            else:
                    tr_second = 0.0
            '''
            #tr_hour, tr_minute, tr_second = 0.0, 0.0, 0.0 # transit starts at midnight
            # transit of progression starts at the same time as birth time
            tr_hour, tr_minute, tr_second = hour, minute, second
            
            # ------- place of transit --------
            tr_place = s_transit[4]
            tr_zbslat, tr_zbslon = s_transit[5], s_transit[6]
            
            tr_east, tr_north = False, False
            tr_z_lon = re.findall('[a-zA-Z]+', tr_zbslon)
            tr_z_lat = re.findall('[a-zA-Z]+', tr_zbslat)
                                    
            if "e" in tr_z_lon or "E" in tr_z_lon:
                    tr_east = True      
            if "n" in tr_z_lat or "N" in tr_z_lat:
                    tr_north = True            
                    
            tr_lat_list = re.findall('\d+', tr_zbslat)
            tr_lon_list = re.findall('\d+', tr_zbslon)
            
            tr_deglat, tr_minlat = int(tr_lat_list[0]), int(tr_lat_list[1])
            tr_deglon, tr_minlon = int(tr_lon_list[0]), int(tr_lon_list[1])
            tr_seclat, tr_seclon = 0, 0

            #altitude = 100
                    
            tr_lat = tr_deglat+tr_minlat/60.0+tr_seclat/3600.0
            tr_lon = tr_deglon+tr_minlon/60.0+tr_seclon/3600.0
                    
            if not tr_north:
                    tr_lat *= -1.0

            if not tr_east:
                    tr_lon *= -1.0
            
            tr_origyear = tr_year
            tr_origmonth = tr_month
            tr_origday = tr_day
            bc = False # before christ
            cal = GREGORIAN # calendar
            tr_zt = 0
            
            tr_zone = s_transit[3].strip()
            
            if re.findall('\:', tr_zone):
                    tr_zone_list = tr_zone.split(':')
                    tr_sign = tr_zone_list[0][0]
                    tr_zh = int(tr_zone_list[0][1:])
                    tr_zm = int(tr_zone_list[1])
                    tr_zs = int(tr_zone_list[2])
            else:
                    tr_sign = tr_zone[0] # +/-
                    tr_zh = int(tr_zone[1:])
                    tr_zm, tr_zs = 0, 0
            
            if tr_sign=='+':
                    tr_plus = True
            else:
                    tr_plus = False    
            
            daylightsaving = False # default is False

            tr_time = tr_hour+tr_minute/60.0+tr_second/3600.0
            
                    
            if daylightsaving:
                    tr_time -= 1.0
                    tr_dhour -= 1
            if tr_time < 0.0:
                    tr_time += HOURSPERDAY
                    tr_year, tr_month, tr_day = util.decrDay(tr_year, tr_month, tr_day)
                    tr_dhour += int(HOURSPERDAY)
                    tr_dyear, tr_dmonth, tr_dday = tr_year, tr_month, tr_day
                            
            if tr_zt == 0:
                    tr_ztime = tr_zh+tr_zm/60.0
                    if tr_plus:
                            tr_time-=tr_ztime
                    else:
                            tr_time+=tr_ztime
                    
                    
            #check over/underflow
            if tr_time >= HOURSPERDAY:
                    tr_time -= HOURSPERDAY
                    tr_year, tr_month, tr_day = util.incrDay(tr_year, tr_month, tr_day)
            elif tr_time < 0.0:
                            tr_time += HOURSPERDAY
                            tr_year, tr_month, tr_day = util.decrDay(tr_year, tr_month, tr_day)
            
            tr_dyear, tr_dmonth, tr_dday, tr_dhour, tr_dmin, tr_dsec = tr_year, tr_month, tr_day, tr_hour, tr_minute, tr_second
                    
            # --------------- End: getting TRANSIT data -----

            

            

            
            # -------------- PREPARING THE PLOT ----------------
            fig = plt.figure()    
            ax = fig.add_subplot(111)
                                    
            #ax = plt.axes()
            plt.xticks([])
            
            #plt.xticks(timeticks, xticklabels, rotation=60, fontsize=7)
            #diffe = (timelim2 - timelim1) / 7
            
            interval = options.time_interval # in days
            min_days, max_days = 0, options.time_interval
            otstup = 0.012*max_days
                                    
            ax.set_xlim(0,38)
            ax.set_ylim(interval, 0)
            
            day_label = [0]
            
            # add lines to separate days
            for i in range(0,max_days+1):
                    plt.axhline(y=i, xmin=0, xmax=1, alpha=0.2) # xmin and xmax take only relative val 0 and 1 (1=whole range, 100%)
                    #day_label.append(i-0.5)
                    day_label.append(i)
            
                                                            
            #track_order = [6,8,10,12,14,16,18,20,22,24,26,28]
            # tracks 2 and 4 for Moon
            #track_occupancy = {2:[],4:[],6:[],8:[],10:[],12:[],14:[],16:[],18:[],20:[],22:[],24:[],26:[],28:[]}
            free_track = 8
            
            track_order = [8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38]
            track_occupancy = {2:[],4:[],6:[],7:[],8:[],9:[],10:[],11:[],12:[],13:[],14:[],15:[],16:[],17:[],18:[],19:[],20:[],21:[],22:[],23:[],24:[],25:[],26:[],27:[],28:[],29:[],30:[],31:[],32:[],33:[],34:[],35:[],36:[],37:[],38:[]}
            
            moon_track = 6 # moon track =2, but aspects may overlap, so region [1.5,4.5] is assign to them
            # ----------- END preparing to plot ------------
            
            
            
            
            # -------- START ----- calculate TRANSIT 1 ------------
            tr1_year, tr1_month, tr1_day, tr1_time = tr_year, tr_month, tr_day, tr_time          

            # -------------horoscope--------------
            tr1_jd = swe.julday(tr1_year, tr1_month, tr1_day, tr1_time, calflag) # julian date
            tr1_d = swe.deltat(tr1_jd)
            tr1_obl, tr1_rflag = swe.calc(tr1_jd+tr1_d, SE_ECL_NUT, 0)    

            # -------------houses--------------
            tr1_cusps, tr1_ascmc, tr1_ascmc2 = get_houses(tr1_jd, hflag, tr_lat, tr_lon, hsys, tr1_obl[0])
            tr1_raequasc, tr1_declequasc, tr1_dist = swe.cotrans(tr1_ascmc[SE_EQUASC], 0.0, 1.0, -tr1_obl[0])
            options_meannode = True        

            # -------------planets--------------       
            #tr1_pls = []    
            tr1_pls = {}
            for pId in planets:
                pl = get_planet(tr1_jd, pId, pflag, tr_lat, tr1_ascmc2, tr1_raequasc, nolat, tr1_obl[0])
                tr1_pls[pId] = pl
                
            # ---------- check degree overflow ----------
            for pId in tr1_pls:
                pl_info = tr1_pls[pId]
                #print pl_info
                pl_norm_pos = util.normalize(pl_info[0][0])
                tr1_pls[pId] = pl_norm_pos
            
            # ------------ set time labels
            #ylabels = ['', str(tr_origday)+'.'+str(tr_origmonth)]
            ylabels = ['', str(tr_origday)+'('+str(count_years)+')']
            count_years+=1
            
            weekdays = {0:'Mon',1:'Tue',2:'Wed',3:'Thu',4:'Fri',5:'Sat',6:'Sun'}
            months = {1:'January', 2:'February',3:'March',4:'April',5:'May',6:'June',\
            7:'July', 8:'August',9:'September',10:'October',11:'November',12:'December'}
            
            # -------- MAKE LABELS USING TRANSIT 2 (END OF THE PERIOD) ------------
            tr2_year, tr2_month, tr2_day, tr2_time = tr_origyear, tr_origmonth, tr_origday, 0.0
            
            prog_year=tr2_year
            for i in range(1,interval+1): # +1 because we need not 00:00, but 23:59:59
                    tr2_year, tr2_month, tr2_day = util.incrDay(tr2_year, tr2_month, tr2_day)
                    prog_year+=1
                    ylabels.append(str(tr2_day)+'('+str(count_years)+','+str(prog_year)+')')
                    count_years+=1
            
            plt.yticks(day_label, ylabels)
            
            pls_signs_title = ''
            transit_pls = [9,8,7,6,5,4,3,2,0] # exclude Moon=1 for now, excluding NN,Lilith,Chiron 
            
            for t in transit_pls:
                tr1_pl_pos = util.normalize(tr1_pls[t])
                tr1_pl_sign = n_signs[tr1_pl_pos//30+1]
                
                pls_signs_title = pls_signs_title + uni_pls[swe.get_planet_name(t)]+uni_signs[tr1_pl_sign]+'  '
            
            #plt.title(months[tr_origmonth]+' '+str(tr_origyear)+'     '+pls_signs_title+' '+tr_zbslat+';'+tr_zbslon, loc='left', fontsize=12)
            #plt.title(months[tr_origmonth]+' '+str(tr_origyear)+'     '+pls_signs_title, loc='left', fontsize=12)
            
            plt.text(0, -otstup*2.4, months[tr_origmonth]+' '+str(tr_origyear)+' Progressions',  size=10)
            plt.text(0, -otstup*1.2, pls_signs_title,  size=10)
            
            
            
            
            
            #transit_pls = [9,8,7,6,5,4,3,2,1,0] #[0,1,2,3,4,5,6,7,8,9,11]    # 11 is mean node
            transit_pls = [9,8,7,6,5,4,3,2,0] # exclude Moon=1 for now, excluding NN,Lilith,Chiron 
                                            
            transits_in = []    
            transits_in_track = {}
            moon_transits_in1 = [] # for 1 and 2 days
            moon_transits_in2 = [] # for 2 and 3 days (one day overlap)
            moon_transits_in = []
            step = 1
            
            #natal_obj = [0,1,2,3,4,5,6,7,8,9,11,12, 'Asc', 'MC']
            natal_obj = [0,1,2,3,4,5,6,7,8,9, 'Asc', 'MC'] # NN, Lilith, Chiron
            

            tr_planet_sign_change = [] # [ME, Ari, when] # mercury entered Aries
            tr_planet_house_change = [] # [ME, 1, when] # mercury entered 1
            tr_planet_direction_change = [] 

            # ------------- DAY BY DAY TRANSIT CALCULATION FOR MOON ------------------
            moon_sign_before = n_signs[tr1_pls[1]//30+1]
            moon_sign_ymin = 0

           
            moon_house_ymin = 0
            moon_first_h_in_plotted=False

            tr2_year, tr2_month, tr2_day, tr2_time = tr1_year, tr1_month, tr1_day, tr1_time
            from_when, till_when = min_days, max_days

            moon_asp_shift = 1

            # one MOON sign at the top
            plt.text(-0.75, moon_sign_ymin-otstup*0.3, uni_pls['Moon'],  size=18)

            for i in range(1,max_days+1,step): # i = every day

                    moon_transits_in = moon_transits_in1
                    if i not in [1,2] and i%2==1:
                           #print i, moon_transits_in
                           moon_transits_in1 = [] # renew Moon aspects that have happend every 2 days
                           moon_transits_in = moon_transits_in2
                    if i not in [1,2] and i%2==0:
                           moon_transits_in2 = [] # renew Moon aspects that have happend every 2 days
                           moon_transits_in = moon_transits_in1

                    #print 'Day', i, '\n'
                    # -------------calculate each day
                    tr2_year, tr2_month, tr2_day = util.incrDay(tr2_year, tr2_month, tr2_day)

                    # -------------horoscope--------------
                    tr2_jd = swe.julday(tr2_year, tr2_month, tr2_day, tr2_time, calflag)
                    tr2_d = swe.deltat(tr2_jd)
                    tr2_obl, tr2_rflag = swe.calc(tr2_jd+tr2_d, SE_ECL_NUT, 0)

                    # -------------houses--------------
                    tr2_cusps, tr2_ascmc, tr2_ascmc2 = get_houses(tr2_jd, hflag, tr_lat, tr_lon, hsys, tr2_obl[0])
                    tr2_raequasc, tr2_declequasc, tr2_dist = swe.cotrans(tr2_ascmc[SE_EQUASC], 0.0, 1.0, -tr2_obl[0])
                    options_meannode = True 

                    # -------------planets--------------    
                    tr2_pls = {}
                    for pId in planets:
                           pl = get_planet(tr2_jd, pId, pflag, tr_lat, tr2_ascmc2, tr2_raequasc, nolat, tr2_obl[0])
                           tr2_pls[pId] = pl
                           
                    # ---------- check degree overflow ----------
                    for pId in tr2_pls:
                        pl_info = tr2_pls[pId]
                        pl_norm_pos = util.normalize(pl_info[0][0])
                        tr2_pls[pId] = pl_norm_pos
            
                    
                    
                    # --------- Calculate Moon (every half an hour) ---------
                    # make a list of Moon positions for every hour
                    modiff = tr2_pls[1] - tr1_pls[1]
                    
                    if modiff<0: # further apart or at the boundary of pisces-aries
                        modiff = tr2_pls[1]+360 - tr1_pls[1]
                    
                    moon_path_per_day = modiff
                    #print 'Moon path per day:', moon_path_per_day
                    
                    moon_per_hour = float(moon_path_per_day) / 48 # every half an hour
                    
                    moon_pos1 = tr1_pls[1]
                    
                    moon_poss = [moon_pos1]
                    
                    for m in range(1,47): # 00:00, 00:30, ... , 23:30, 24:00 = 00:00 of another day
                        moon_pos1+=moon_per_hour
                        if moon_pos1>=360:
                                moon_pos1 = 0
                        moon_poss.append(moon_pos1)
                    
                    #print moon_poss
                    
                    
                    
                    
                    for m_pos in moon_poss: # pozicii Luny v techenii odnih sutok
                            
                            moon_sign_now = n_signs[m_pos//30+1]
                            if moon_sign_now!=moon_sign_before: # Moon went to new sign, plot the line and sign
                                    moon_sign_ymax = i-1 + float((moon_poss.index(m_pos)+1))/48
                                    
                                    rect_col = ''
                                    
                                    if moon_sign_before in ['Ari', 'Leo', 'Sag']:
                                            rect_col = 'r'
                                    elif moon_sign_before in ['Can', 'Sco', 'Pis']:
                                            rect_col = 'b'
                                    elif moon_sign_before in ['Tau', 'Vir', 'Cap']:
                                            rect_col = 'g'
                                    elif moon_sign_before in ['Gem', 'Lib', 'Aqu']:
                                            rect_col = 'y'
                                    else:
                                            rect_col = 'k'
                                    
                                    
                                    #plt.axvspan(ymin=moon_sign_ymin, ymax=moon_sign_ymax, xmin=0, xmax=4, alpha=0.05, color=rect_col)
                                    
                                    
                                    #ymin and ymax are relative to 0 and 1
                                    moon_sign_ymin_rel = moon_sign_ymin/max_days
                                    moon_sign_ymax_rel = moon_sign_ymax/max_days
                                    # and it starts from the bottom 0,0 => so I had to reverse
                                    plt.axvspan(ymin=1-moon_sign_ymin_rel, ymax=1-moon_sign_ymax_rel, xmin=0, xmax=7, alpha=0.05, color=rect_col)
                                    plt.text(0.1, moon_sign_ymin+otstup*1.5, uni_pls['Moon']+uni_signs[moon_sign_before],  size=12)
                                    
                                    moon_sign_ymin = moon_sign_ymax
                                    moon_sign_before = moon_sign_now
                                    
                            
                            # check Moon aspects to natal planets
                            for n in [0,1,2,3,4,5,6,7,8,9]:
                                    try:
                                        n_pl_name = swe.get_planet_name(n)
                                    except TypeError:
                                        n_pl_name=n
                                    if n==11:
                                        n_pl_pos = util.normalize(n_pls[10]) # NN
                                    elif n==12:
                                        n_pl_pos = util.normalize(n_pls[11]) # Lilith
                                    elif n=='Asc':
                                        n_pl_pos = util.normalize(n_pls[12]) # Asc 
                                    elif n=='MC':
                                        n_pl_pos = util.normalize(n_pls[13]) # MC
                                    else:
                                        n_pl_pos = util.normalize(n_pls[n])
                                    
                                    if m_pos>=n_pl_pos:
                                        diffm = m_pos - n_pl_pos
                                        if diffm>=180: # further apart or at the boundary of pisces-aries
                                            diffm = abs(m_pos-360) + n_pl_pos
                                    else:
                                        diffm = n_pl_pos - m_pos
                                        if diffm>=180:
                                            diffm = abs(n_pl_pos - 360) + m_pos                                        
                                            
                                    is_asp_mo, asp_type_mo, diff_with_exact_mo = find_tr_asp(diffm, trans_orbis, Aspects)
                                    
                                    if is_asp_mo:
                                            if ('Moon', n_pl_name,asp_type_mo) not in moon_transits_in:
                                                    moon_transits_in.append(('Moon', n_pl_name,asp_type_mo))
                                                    
                                                    if i==1:
                                                            moon_transits_in1.append(('Moon', n_pl_name,asp_type_mo))
                                                    else:
                                                            moon_transits_in1.append(('Moon', n_pl_name,asp_type_mo))
                                                            moon_transits_in2.append(('Moon', n_pl_name,asp_type_mo))
                                                    
                                                    where_moon_tr = i-1 + float((moon_poss.index(m_pos)+1))/48
                                                    
                                                    # moon text size
                                                    mo_text_size = 5
                                                    if ('Moon', n_pl_name) in natal_asps or (n_pl_name,'Moon') in natal_asps:
                                                            mo_text_size = 7
                                                    
                                                    if n_pl_name=='Asc' or n_pl_name=='MC':
                                                            plt.text(moon_track-moon_asp_shift, where_moon_tr, uni_asps[asp_type_mo]+n_pl_name[0],  size=mo_text_size)
                                                    else:
                                                            plt.text(moon_track-moon_asp_shift, where_moon_tr, uni_asps[asp_type_mo]+uni_pls[n_pl_name],  size=mo_text_size)
                                                    
                                                    if moon_track-moon_asp_shift<=2.4:
                                                            moon_asp_shift=0
                                                    else:
                                                            moon_asp_shift+=0.8
                                                    
                                                    #print 'MOON aspect with', n_pl_name, where_moon_tr, moon_asp_shift
                    
                    # second transit becomes first
                    tr1_pls = tr2_pls        
            

            # --------------- End: calculating the chart for every step -----


            # plot last moon sign
            
            if moon_sign_now in ['Ari', 'Leo', 'Sag']:
                    rect_col = 'r'
            elif moon_sign_now in ['Can', 'Sco', 'Pis']:
                    rect_col = 'b'
            elif moon_sign_now in ['Tau', 'Vir', 'Cap']:
                    rect_col = 'g'
            elif moon_sign_now in ['Gem', 'Lib', 'Aqu']:
                    rect_col = 'y'
            else:
                    rect_col = 'k'
            moon_sign_ymin_rel = moon_sign_ymin/max_days
            moon_sign_ymax_rel = moon_sign_ymax/max_days
            plt.axvspan(ymin=1-moon_sign_ymin_rel, ymax=0, xmin=0, xmax=7, alpha=0.05, color=rect_col)
            plt.text(0.1, moon_sign_ymin+otstup*1.5, uni_pls['Moon']+uni_signs[moon_sign_before],  size=12)





            
            

            # ------------- DAY BY DAY TRANSIT CALCULATION FOR ALL TRANSITS ------------------


            # -------- TRANSIT 1 DECREASED BY 1 hour (to check planets direction)------------
            tr1_year, tr1_month, tr1_day, tr1_time = tr_year, tr_month, tr_day, tr_time
            tr1_year, tr1_month, tr1_day, tr1_time = util.subtractHour(tr1_year, tr1_month, tr1_day, tr1_time)


            # -------------horoscope--------------
            tr1_jd = swe.julday(tr1_year, tr1_month, tr1_day, tr1_time, calflag) # julian date
            tr1_d = swe.deltat(tr1_jd)
            tr1_obl, tr1_rflag = swe.calc(tr1_jd+tr1_d, SE_ECL_NUT, 0)    

            # -------------houses--------------
            tr1_cusps, tr1_ascmc, tr1_ascmc2 = get_houses(tr1_jd, hflag, tr_lat, tr_lon, hsys, tr1_obl[0])
            tr1_raequasc, tr1_declequasc, tr1_dist = swe.cotrans(tr1_ascmc[SE_EQUASC], 0.0, 1.0, -tr1_obl[0])
            options_meannode = True        

            # -------------planets--------------       
            #tr1_pls = []    
            tr1_pls = {}
            for pId in planets:
                pl = get_planet(tr1_jd, pId, pflag, tr_lat, tr1_ascmc2, tr1_raequasc, nolat, tr1_obl[0])
                tr1_pls[pId] = pl
                
            # ---------- check degree overflow ----------
            for pId in tr1_pls:
             pl_info = tr1_pls[pId]
             #print pl_info
             pl_norm_pos = util.normalize(pl_info[0][0])
             tr1_pls[pId] = pl_norm_pos

            #print tr1_pls, tr1_cusps
            tr1_pls_in_houses = calc_planet_houses(tr1_pls, list(tr1_cusps), options.ruling)

            
            
            transit_pls = [9,8,7,6,5,4,3,2,0] # exclude Moon=1 for now, exclude NN,LIlith,Ch
            pl_direction = {0:'D',1:'D',2:'D',3:'D',4:'D',5:'D',6:'D',7:'D',8:'D',9:'D',11:'D',12:'D'}
            #natal_obj = [0,1,2,3,4,5,6,7,8,9,11,12, 'Asc', 'MC']            
            transit_pls_rev = [0,2,3,4,5,6,7,8,9] # exclude NN,LIlith,Ch
            
            transits_in = {}    
            transits_in2 = {}
            transits_in_track = {}
            
            
            tr2_year, tr2_month, tr2_day, tr2_time = tr1_year, tr1_month, tr1_day, tr1_time
            from_when, till_when = min_days, max_days

            max_hours = max_days*24
            step = 1
            
            #print '\nAll transits\n'
            
            for i in range(0,max_hours+1,step): # i = every 1 hour should be enough
                    #print 'Day', i/24 ,'Hour', i, '\n'

                    # -------------calculate each hour
                    tr2_year, tr2_month, tr2_day, tr2_time = util.addHour(tr2_year, tr2_month, tr2_day, tr2_time)

                    # -------------horoscope--------------
                    tr2_jd = swe.julday(tr2_year, tr2_month, tr2_day, tr2_time, calflag)
                    tr2_d = swe.deltat(tr2_jd)
                    tr2_obl, tr2_rflag = swe.calc(tr2_jd+tr2_d, SE_ECL_NUT, 0)

                    # -------------houses--------------
                    tr2_cusps, tr2_ascmc, tr2_ascmc2 = get_houses(tr2_jd, hflag, tr_lat, tr_lon, hsys, tr2_obl[0])
                    tr2_raequasc, tr2_declequasc, tr2_dist = swe.cotrans(tr2_ascmc[SE_EQUASC], 0.0, 1.0, -tr2_obl[0])
                    options_meannode = True 

                    # -------------planets--------------    
                    tr2_pls = {}
                    for pId in planets:
                           pl = get_planet(tr2_jd, pId, pflag, tr_lat, tr2_ascmc2, tr2_raequasc, nolat, tr2_obl[0])
                           tr2_pls[pId] = pl
                           
                    # ---------- check degree overflow ----------
                    for pId in tr2_pls:
                        pl_info = tr2_pls[pId]
                        pl_norm_pos = util.normalize(pl_info[0][0])
                        tr2_pls[pId] = pl_norm_pos
            
            
            
                    # ---------- check directions of planets ---------        
                    for p in [0,1,2,3,4,5,6,7,8,9]:#n_pls:
                        
                        if isinstance(p, int): # if it is a planet incoded with int
                            pos1 = tr1_pls[p]
                            pos2 = tr2_pls[p]
                            
                            current_pl_direction = 'D' # default
                            
                            if pos2 > pos1 and abs(pos2-pos1)<180: # D
                                current_pl_direction='D'
                            elif pos2 < pos1 and abs(pos2-pos1)>180: # pisces-aries step
                                current_pl_direction='D'
                            else:
                                current_pl_direction='R'
                            
                            
                            try:
                                if pl_direction[p]!=current_pl_direction:
                                    where_change = float(i)/24
                                    tr_planet_direction_change.append([swe.get_planet_name(p),current_pl_direction,where_change])
                            except KeyError:
                                pass
                                
                            
                            pl_direction[p] = current_pl_direction
                        
                        #print 'pl status:', pl_direction, '\n'
                    
                    
                        # -------- check if transiting planet changed a sign and house ---------- START ---------
                    
                        tr1_pl_pos = util.normalize(tr1_pls[p])
                        tr2_pl_pos = util.normalize(tr2_pls[p])
                        
                        tr1_pl_sign = n_signs[tr1_pl_pos//30+1]
                        tr2_pl_sign = n_signs[tr2_pl_pos//30+1]
                        
                        #print p, tr1_pl_sign, tr2_pl_sign
                        if tr1_pl_sign!=tr2_pl_sign and p!=1: # Moon is separate
                            where_change = float(i)/24
                            tr_planet_sign_change.append([swe.get_planet_name(p),tr2_pl_sign,where_change])
                            
                        
                        
                        tr2_pls_in_houses = calc_planet_houses(tr2_pls, list(cusps), options.ruling) # using cusps of natal
                        tr1_pl_h = tr1_pls_in_houses[swe.get_planet_name(p)][0] # dict {ME:[1, [2,3]]}
                        tr2_pl_h = tr2_pls_in_houses[swe.get_planet_name(p)][0]
                        
                        
                        if tr1_pl_h!=tr2_pl_h and p!=1:
                            where_change = float(i)/24
                            tr_planet_house_change.append([swe.get_planet_name(p),tr2_pl_h,where_change])
                        
                        #tr_planet_house_change
                        
                        # -------- check if transiting planet changed a sign and house -------- END ---------
            
            
            
            
            
            
                    # any aspects in between transiting planets?
                    for k in range(0,len(transit_pls_rev)-1):
                    #for t in transit_pls:
                            t = transit_pls_rev[k]
                            tr2_pl_name = swe.get_planet_name(t)
                            tr2_pl_pos = util.normalize(tr2_pls[t])
                            #tr2_pl_sign = n_signs[tr2_pl_pos//30+1]
                    
                    
                            for j in range(k+1,len(transit_pls_rev)):
                            #for n in transit_pls:
                                n = transit_pls_rev[j]
                        
                                update = False
                                
                                if t!=n: # e.g., no ME-ME in transits
                                
                                    try:
                                        n_pl_name = swe.get_planet_name(n)
                                    except TypeError:
                                        n_pl_name=n
                                    n_pl_pos = util.normalize(tr2_pls[n])
                            
            
                                    if tr2_pl_pos>=n_pl_pos:
                                        diff2 = tr2_pl_pos - n_pl_pos
                                        if diff2>=180: # further apart or at the boundary of pisces-aries
                                            diff2 = abs(tr2_pl_pos-360) + n_pl_pos
                                    else:
                                        diff2 = n_pl_pos - tr2_pl_pos
                                        if diff2>=180:
                                            diff2 = abs(n_pl_pos - 360) + tr2_pl_pos                                        
                                    
                                    
                                    # orbs for transits are the same = 1 degree for example
                                    is_asp2, asp_type2, diff_with_exact2 = find_tr_asp(diff2, trans_orbis, Aspects)    
                                    
                                    if is_asp2:
                                        #print t,n, tr2_pl_name, n_pl_name, is_asp2, asp_type2, diff_with_exact2
                                        #print tr2_pl_pos, n_pl_pos
                                        
                                        for_hours = i
                                        if for_hours==24:
                                            for_hours=0
                                        where_asp = float(i)/24
                                
                                        if (pl_direction[t], tr2_pl_name, n_pl_name, asp_type2) not in transits_in:
                                            transits_in[(pl_direction[t], tr2_pl_name, n_pl_name, asp_type2)] = [where_asp]
                                        else: # update
                                            asps_coords = transits_in[(pl_direction[t], tr2_pl_name, n_pl_name, asp_type2)]
                                            asps_coords.append(where_asp)
                                            transits_in[(pl_direction[t], tr2_pl_name, n_pl_name, asp_type2)] = asps_coords
                                        
                                                    
                                                    
                                                    
                    # any aspects in between transiting and natal planets?
                    for k in range(0,len(transit_pls_rev)-2): # from transiting sun to mars, others move 
                        
                    #for t in transit_pls:
                        t = transit_pls_rev[k]
                        tr2_pl_name = swe.get_planet_name(t)
                        #print tr2_pl_name
                        tr2_pl_pos = util.normalize(tr2_pls[t])
                        #tr2_pl_sign = n_signs[tr2_pl_pos//30+1]
                
                
                        for j in range(0,len(planets)-4): # no need for NN,Li,Ch
                            n = planets[j]
                        
                            update = False
                            
                            #if t!=n: # e.g., no ME-ME in transits
                            # in progression can be SU to natal SU aspects
                            
                            try:
                                n_pl_name = swe.get_planet_name(n)
                            except TypeError:
                                n_pl_name=n
                            n_pl_pos = util.normalize(n_pls[n])
                    
    
                            if tr2_pl_pos>=n_pl_pos:
                                diff2 = tr2_pl_pos - n_pl_pos
                                if diff2>=180: # further apart or at the boundary of pisces-aries
                                    diff2 = abs(tr2_pl_pos-360) + n_pl_pos
                            else:
                                diff2 = n_pl_pos - tr2_pl_pos
                                if diff2>=180:
                                    diff2 = abs(n_pl_pos - 360) + tr2_pl_pos                                        
                            
                            
                            # orbs for transits are the same = 1 degree for example
                            is_asp2, asp_type2, diff_with_exact2 = find_tr_asp(diff2, trans_orbis, Aspects)    
                            
                            if is_asp2:
                                #print t,n, tr2_pl_name, n_pl_name, is_asp2, asp_type2, diff_with_exact2
                                #print tr2_pl_pos, n_pl_pos
                                
                                for_hours = i
                                if for_hours==24:
                                    for_hours=0
                                where_asp = float(i)/24
                        
                                if (pl_direction[t], tr2_pl_name, n_pl_name, asp_type2) not in transits_in2:
                                    transits_in2[(pl_direction[t], tr2_pl_name, n_pl_name, asp_type2)] = [where_asp]
                                else: # update
                                    asps_coords = transits_in2[(pl_direction[t], tr2_pl_name, n_pl_name, asp_type2)]
                                    asps_coords.append(where_asp)
                                    transits_in2[(pl_direction[t], tr2_pl_name, n_pl_name, asp_type2)] = asps_coords

                    '''                        
                    # any aspects of progressing ASC / MC to natal planets?
                    for k in ['ASC', 'MC']: # from transiting sun to mars, others move 
                        
                        tr2_pl_name = k
                        
                        if k=='ASC':
                            tr2_pl_pos = util.normalize(tr2_cusps[0])
                        elif k=='MC':
                            tr2_pl_pos = util.normalize(tr2_cusps[9])
                        else:
                            print 'Not ASC, not MC? Something is wrong!'
                
                        for j in range(0,len(planets)-2):
                            n = planets[j]
                        
                            update = False
                            
                            #if t!=n: # e.g., no ME-ME in transits
                            # in progression can be SU to natal SU aspects
                            
                            try:
                                n_pl_name = swe_get_planet_name(n)
                            except TypeError:
                                n_pl_name=n
                            n_pl_pos = util.normalize(n_pls[n])
                    
    
                            if tr2_pl_pos>=n_pl_pos:
                                diff2 = tr2_pl_pos - n_pl_pos
                                if diff2>=180: # further apart or at the boundary of pisces-aries
                                    diff2 = abs(tr2_pl_pos-360) + n_pl_pos
                            else:
                                diff2 = n_pl_pos - tr2_pl_pos
                                if diff2>=180:
                                    diff2 = abs(n_pl_pos - 360) + tr2_pl_pos                                        
                            
                            
                            # orbs for transits are the same = 1 degree for example
                            is_asp2, asp_type2, diff_with_exact2 = find_tr_asp(diff2, trans_orbis, Aspects)    
                            
                            if is_asp2 and asp_type2==0.0:
                                
                                for_hours = i
                                if for_hours==24:
                                    for_hours=0
                                where_asp = float(i)/24
                        
                                if ('D', tr2_pl_name, n_pl_name, asp_type2) not in transits_in2:
                                    transits_in[('D', tr2_pl_name, n_pl_name, asp_type2)] = [where_asp]
                                    
                                else: # update
                                    asps_coords = transits_in2[('D', tr2_pl_name, n_pl_name, asp_type2)]
                                    asps_coords.append(where_asp)
                                    transits_in[('D', tr2_pl_name, n_pl_name, asp_type2)] = asps_coords
                    '''
                    
                    
                    
                    
                    
                    
                    
            
                    tr1_pls = tr2_pls    
                    tr1_pls_in_houses = tr2_pls_in_houses
                    
            
            #print transits_in
            #print transits_in2
            
            
            
            # ----------------- plot all other transits (except Moon) ------------------------
            
            for t in transit_pls:
                    tr_pl_name = swe.get_planet_name(t)
                    
                    tr_pl_asps = {}
                    tr_pl_asps2 = {}
                    tr_pl_asps_plotted = {}
                    tr_pl_asps_plotted2 = {}
                    single_asp_type = [] # (tr_pl,n_pl,asp_type)
                    single_asp_type2 = []
                    
                    for k in transits_in:
                            
                        if tr_pl_name==k[1]:
                            tr_pl_asps[k] = transits_in[k]
                            
                            tr_pl_asps_plotted[(k[1],k[2],k[3])] = 0
                            
                            if (k[1],k[2],k[3]) not in single_asp_type:
                                single_asp_type.append((k[1],k[2],k[3]))
                                
                    for k in transits_in2:
                            
                        if tr_pl_name==k[1]:
                            tr_pl_asps2[k] = transits_in2[k]
                            
                            tr_pl_asps_plotted2[(k[1],k[2],k[3])] = 0
                            
                            if (k[1],k[2],k[3]) not in single_asp_type2:
                                single_asp_type2.append((k[1],k[2],k[3]))
                                                    
                    
                    # ------------ TRANSIT to TRANSIT PLANETS (SKY TRANSIT) ------------
                    for k in single_asp_type:
                            
                        if tr_pl_asps_plotted[k]==0: # that aspect has not been plotted yet
                                
                                
                            tr_pl_D = ('D',k[0],k[1],k[2])
                            tr_pl_R = ('R',k[0],k[1],k[2])

                            
                            if tr_pl_D in tr_pl_asps and tr_pl_R in tr_pl_asps: # then we can try to plot them on one track
                                
                                asp_coords_merged = tr_pl_asps[tr_pl_D] + tr_pl_asps[tr_pl_R]
                            
                                from_when = min(asp_coords_merged)
                                till_when = max(asp_coords_merged)
                                
                                ok_track = find_available_track(track_order, track_occupancy,from_when,till_when, otstup) 
                                
                                asp_coords = sorted(asp_coords_merged)
                                gap_found = False
                                # check if there is a gap due to change in direction
                                i0=0
                                for i1 in range(0,len(asp_coords)-1):
                                    if (asp_coords[i1+1] - asp_coords[i1]) > 2: # if change direction once
                                        #print 'Gap', k, asp_coords[i1], asp_coords[i1+1]
                                        
                                        asp_coords1 = asp_coords[i0:i1+1]
                                        #asp_coords2 = asp_coords[i1+1:]
                                        i0=i1+1
                                        #print asp_coords1
                                        # -------------- plot aspect ---------------
                                        #print 'Plotting', tr_pl_D, tr_pl_R, min(asp_coords1), max(asp_coords1)
                                        plot_aspect(plt, natal_asps, asp_coords1, ok_track, min_days, max_days, otstup, k)
                                        #plot_aspect(plt, natal_asps, asp_coords2, ok_track, min_days, max_days, otstup, k)
                                        
                                        gap_found = True
                                        #break
                                if gap_found:
                                    asp_coords2 = asp_coords[i0:]
                                    #print asp_coords2, 'LAST'
                                    plot_aspect(plt, natal_asps, asp_coords2, ok_track, min_days, max_days, otstup, k)    
                                    #print 'Plotting', tr_pl_D, tr_pl_R, min(asp_coords2), max(asp_coords2)
                                else:
                                    # -------------- plot aspect ---------------
                                    asp_coords1 = tr_pl_asps[tr_pl_D]
                                    asp_coords2 = tr_pl_asps[tr_pl_R]
                                    plot_aspect(plt, natal_asps, asp_coords1, ok_track, min_days, max_days, otstup, k)
                                    plot_aspect(plt, natal_asps, asp_coords2, ok_track, min_days, max_days, otstup, k)
                                    #print 'Plotting both', tr_pl_D, tr_pl_R, min(asp_coords1), max(asp_coords1), 'and', min(asp_coords2), max(asp_coords2)
                                    
                        
                                # flag them as used    
                                tr_pl_asps_plotted[k]==1
                                    
                            else:
                                try:
                                    asp_coords = tr_pl_asps[tr_pl_D]
                                except KeyError:
                                    asp_coords = tr_pl_asps[tr_pl_R]
                                #print k, asp_coords
                                from_when = min(asp_coords)
                                till_when = max(asp_coords)
                                
                                #print asp_coords
                                ok_track = find_available_track(track_order, track_occupancy,from_when,till_when, otstup) 
                                
                                gap_found = False
                                # check if there is a gap due to change in direction
                                i0 = 0
                                for i1 in range(0,len(asp_coords)-1):
                                    if (asp_coords[i1+1] - asp_coords[i1]) > 2: # if change direction once
                                        #print 'Gap', k, asp_coords[i1], asp_coords[i1+1]
                                        
                                        asp_coords1 = asp_coords[i0:i1+1]
                                        #print asp_coords1
                                        #asp_coords2 = asp_coords[i1+1:]
                                        i0 = i1+1
                                        # -------------- plot aspect ---------------
                                        #print 'Plotting', tr_pl_D, tr_pl_R, min(asp_coords1), max(asp_coords1)
                                        plot_aspect(plt, natal_asps, asp_coords1, ok_track, min_days, max_days, otstup, k)
                                        #plot_aspect(plt, natal_asps, asp_coords2, ok_track, min_days, max_days, otstup, k)
                                        
                                        gap_found = True
                                        #break
                                if gap_found:
                                    asp_coords2 = asp_coords[i0:]
                                    plot_aspect(plt, natal_asps, asp_coords2, ok_track, min_days, max_days, otstup, k)
                                    #print 'Plotting', tr_pl_D, tr_pl_R, min(asp_coords2), max(asp_coords2)
                                else:
                                    plot_aspect(plt, natal_asps, asp_coords, ok_track, min_days, max_days, otstup, k)
                                    #print 'Plotting full', tr_pl_D, tr_pl_R, min(asp_coords), max(asp_coords)
                                            
                                            
                    # ------------ TRANSIT to NATAL PLANETS ------------
                    for k in single_asp_type2:
                            
                        if tr_pl_asps_plotted2[k]==0: # that aspect has not been plotted yet
                                
                                
                            tr_pl_D = ('D',k[0],k[1],k[2])
                            tr_pl_R = ('R',k[0],k[1],k[2])

                            
                            if tr_pl_D in tr_pl_asps2 and tr_pl_R in tr_pl_asps2: # then we can try to plot them on one track
                                    
                                asp_coords_merged = tr_pl_asps2[tr_pl_D] + tr_pl_asps2[tr_pl_R]
                            
                                from_when = min(asp_coords_merged)
                                till_when = max(asp_coords_merged)
                                
                                ok_track = find_available_track(track_order, track_occupancy,from_when,till_when, otstup) 
                                
                                asp_coords = sorted(asp_coords_merged)
                                gap_found = False
                                # check if there is a gap due to change in direction
                                i0=0
                                for i1 in range(0,len(asp_coords)-1):
                                    if (asp_coords[i1+1] - asp_coords[i1]) > 2: # if change direction once
                                        #print 'Gap', k, asp_coords[i1], asp_coords[i1+1]
                                        
                                        asp_coords1 = asp_coords[i0:i1+1]
                                        #asp_coords2 = asp_coords[i1+1:]
                                        i0=i1+1
                                        #print asp_coords1
                                        # -------------- plot aspect ---------------
                                        #print 'Plotting', tr_pl_D, tr_pl_R, min(asp_coords1), max(asp_coords1)
                                        plot_aspect(plt, natal_asps, asp_coords1, ok_track, min_days, max_days, otstup, k+('N',))
                                        #plot_aspect(plt, natal_asps, asp_coords2, ok_track, min_days, max_days, otstup, k)
                                        
                                        gap_found = True
                                        #break
                                if gap_found:
                                    asp_coords2 = asp_coords[i0:]
                                    #print asp_coords2, 'LAST'
                                    plot_aspect(plt, natal_asps, asp_coords2, ok_track, min_days, max_days, otstup, k+('N',))    
                                    #print 'Plotting', tr_pl_D, tr_pl_R, min(asp_coords2), max(asp_coords2)
                                else:
                                    # -------------- plot aspect ---------------
                                    asp_coords1 = tr_pl_asps2[tr_pl_D]
                                    asp_coords2 = tr_pl_asps2[tr_pl_R]
                                    plot_aspect(plt, natal_asps, asp_coords1, ok_track, min_days, max_days, otstup, k+('N',))
                                    plot_aspect(plt, natal_asps, asp_coords2, ok_track, min_days, max_days, otstup, k+('N',))
                                    #print 'Plotting both', tr_pl_D, tr_pl_R, min(asp_coords1), max(asp_coords1), 'and', min(asp_coords2), max(asp_coords2)
                                    
                        
                                # flag them as used    
                                tr_pl_asps_plotted2[k]==1
                                    
                            else:
                                try:
                                    asp_coords = tr_pl_asps2[tr_pl_D]
                                except KeyError:
                                    asp_coords = tr_pl_asps2[tr_pl_R]
                                #print k, asp_coords
                                from_when = min(asp_coords)
                                till_when = max(asp_coords)
                                
                                #print asp_coords
                                ok_track = find_available_track(track_order, track_occupancy,from_when,till_when, otstup) 
                                
                                gap_found = False
                                # check if there is a gap due to change in direction
                                i0 = 0
                                for i1 in range(0,len(asp_coords)-1):
                                    if (asp_coords[i1+1] - asp_coords[i1]) > 2: # if change direction once
                                        #print 'Gap', k, asp_coords[i1], asp_coords[i1+1]
                                        
                                        asp_coords1 = asp_coords[i0:i1+1]
                                        #print asp_coords1
                                        #asp_coords2 = asp_coords[i1+1:]
                                        i0 = i1+1
                                        # -------------- plot aspect ---------------
                                        #print 'Plotting', tr_pl_D, tr_pl_R, min(asp_coords1), max(asp_coords1)
                                        plot_aspect(plt, natal_asps, asp_coords1, ok_track, min_days, max_days, otstup, k+('N',))
                                        #plot_aspect(plt, natal_asps, asp_coords2, ok_track, min_days, max_days, otstup, k)
                                        
                                        gap_found = True
                                        #break
                                if gap_found:
                                    asp_coords2 = asp_coords[i0:]
                                    plot_aspect(plt, natal_asps, asp_coords2, ok_track, min_days, max_days, otstup, k+('N',))
                                    #print 'Plotting', tr_pl_D, tr_pl_R, min(asp_coords2), max(asp_coords2)
                                else:
                                    plot_aspect(plt, natal_asps, asp_coords, ok_track, min_days, max_days, otstup, k+('N',))
                                    #print 'Plotting full', tr_pl_D, tr_pl_R, min(asp_coords), max(asp_coords)
                                #plot_aspect(plt, natal_asps, asp_coords, ok_track, min_days, max_days, otstup, k+('N',))                    
                                            
                                            
                                            
            #print tr_planet_sign_change
            for c in tr_planet_sign_change:
                pl_name, pl_sign, when_changed = c[0], c[1], c[2]
                #print pl_name, pl_sign, when_changed
                ok_track = find_available_track(track_order, track_occupancy,when_changed-otstup*2,when_changed+otstup*2, otstup) 
                
                plt.text(ok_track, when_changed, uni_pls[pl_name],  size=12)
                plt.text(ok_track, when_changed+otstup, uni_signs[pl_sign], color="black", size=12)    
                plt.plot([ok_track,ok_track], [when_changed, when_changed+otstup], color="darkorange", linewidth=2)
                
            #print tr_planet_house_change
            for c in tr_planet_house_change:
                pl_name, pl_house, when_changed = c[0], c[1], c[2]
                
                ok_track = find_available_track(track_order, track_occupancy,when_changed-otstup*2,when_changed+otstup*2, otstup) 
                
                plt.text(ok_track, when_changed, uni_pls[pl_name],  size=10)
                plt.text(ok_track, when_changed+otstup, str(pl_house), color="black", size=10)    
                plt.plot([ok_track,ok_track], [when_changed, when_changed+otstup], color="darkorange", linewidth=1)
                    
            for c in tr_planet_direction_change:
                pl_name, pl_direction, when_changed = c[0], c[1], c[2]
                
                ok_track = find_available_track(track_order, track_occupancy,when_changed-otstup*2,when_changed+otstup*2, otstup) 
                
                plt.text(ok_track, when_changed, uni_pls[pl_name],  size=12)
                plt.text(ok_track, when_changed+otstup, pl_direction, color="black", size=12)    
                plt.plot([ok_track,ok_track], [when_changed, when_changed+otstup], color="darkorange", linewidth=2)
                                    
            '''                        
            # plot progressing ASC and MC making conjunctions
            for c in ['ASC', 'MC']:
                
                    tr_pl_name = c
                    
                    tr_pl_asps = {}
                    tr_pl_asps_plotted = {}
                    single_asp_type = [] # (tr_pl,n_pl,asp_type)
                    
                    for k in transits_in:
                        print k
                        if tr_pl_name==k[1]:
                            print tr_pl_name
                            tr_pl_asps[k] = transits_in[k]
                            
                            tr_pl_asps_plotted[(k[1],k[2],k[3])] = 0
                            
                            if (k[1],k[2],k[3]) not in single_asp_type:
                                single_asp_type.append((k[1],k[2],k[3]))


                    for k in single_asp_type:
                            
                        if tr_pl_asps_plotted[k]==0: # that aspect has not been plotted yet
                                
                                
                            tr_pl_D = ('D',k[0],k[1],k[2]) # cusps move only directly
                            #tr_pl_R = ('R',k[0],k[1],k[2])

                                
                            asp_coords = tr_pl_asps[tr_pl_D]
                        
                            from_when = min(asp_coords)
                            till_when = max(asp_coords)
                            
                            ok_track = find_available_track(track_order, track_occupancy,from_when,till_when, otstup) 
                                    
                            print 'plos'        
                            # -------------- plot ---------------
                            plot_aspect(plt, natal_asps, asp_coords, ok_track, min_days, max_days, otstup, k)

            '''
            
            
            

            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            
            #plt.show()
            
            fig.set_size_inches(8, 12)
                            
            if len(str(tr_origmonth))==1:
                month_n = '0'+ str(tr_origmonth)
            else:
                month_n = str(tr_origmonth)
            if len(str(tr_origday))==1:
                day_n = '0'+ str(tr_origday)
            else:
                day_n = str(tr_origday)
            
            #fig.savefig(options.output_file+'_'+str(tr_origyear)+'_'+months[tr_origmonth][:3]+str(tr_origday)+'.png', dpi=300, bbox_inches='tight')
            fig.savefig(str(tr_origyear)+'_'+month_n+'_'+day_n+'_'+options.output_file+'.png', dpi=300, bbox_inches='tight')
            plt.close(fig)

    
    
    
        

if __name__ == "__main__":
    main()        
