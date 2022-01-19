#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import unicode_literals

'''
This script calculates natal positions and transits, and makes a plot - a calendar for a time period in years.

Input text file example (ZBS format / ZET astroprocessor, first line - natal, 
second line - start point to calculate transits from, better to set to 1st of Jan):

t; 01.01.1992; 12:20:00; +3; Vilnius, Lithuania; 54n41; 25e19; -; ;
t; 1.1.2000; 00:00:00; +3; Vilnius, Lithuania; 54n41; 25e19; -; ;


Calendar includes:

+ aspects of transit planets to natal planets/cusps, which are:
++ yearly transits of higher and social planets (Pluto, Neptune, Uranus, Saturn, Jupiter)
++ yearly transits of Node, Lilith and Chiron if selected with parameter -f

+ Retrograde planets transits are implemented (direct and retro are marked with joined line)

Default parameters are:
+ transits to all planets and ASC/MC 
+ orb for transits 1 degree (can be expanded)
+ time period 10 years


Running the script:

python3 transit_calendar_years.py -i 'test.zbs'


Running the script with additional parameters:

python3 transit_calendar_years.py -i 'test.zbs' -t 5 -e 3 -c 0
# -t -> time period, 5 years
# -e 3 -> expand transit orb to 3 degrees
# -c 0 -> transits to all cusps

LIMITATIONS:

- transits of Node are limited to only 0 and 180 aspects, and transits of
Lilith and Chiron only to 0 aspect

- no transits of NN, Li, Ch to NN, Li, Ch


'''



from optparse import OptionParser, OptionGroup
import math
import re
import numpy as np
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
    
    usage = "python3 transit_calendar.py -i input.zbs [options]"
    version = "1.0"
    description = "%prog makes a transit calendar years"
    epilog = "Written by Svetlana Stoma / @astroloxplore (2022)"
    parser = OptionParser(usage=usage, description=description,
                          version="%prog "+version, epilog=epilog)
    # in file
    parser.add_option("-i", "--input",  dest="input_file", metavar="<file>",\
    help="Input file: record in ZBS format (ZET astroprocessor). Example: \nt; 01.01.1992; 12:20:00; +3; Vilnius, Lithuania; 54n41; 25e19; -; ;")
    parser.set_defaults(input_file='input.zbs')
    
    # flags to select what aspects to calculate
    parser.add_option("-a", "--asp",  dest="asp_set", type="int",\
    help="Aspects: 0 for short set [0, 60, 90, 120, 180], 1 for major [0, 30, 45, 60, 90, 120, 135, 150, 180], 2 for all [0.0, 30.0, 36.0, 40.0, 45.0, 60.0, 72.0, 80.0, 90.0, 100.0, 108.0, 120.0, 135.0, 144.0, 150.0, 180.0], 3 for my default set [0, 30, 60, 90, 120, 150, 180]")
    parser.set_defaults(asp_set=0)
    
    parser.add_option("-t", "--time",  dest="time_interval", type="int",\
    help="A period of time (in years)")
    parser.set_defaults(time_interval=10) 
    
    parser.add_option("-e", "--expand",  dest="transorbis", type="float",\
    help="default is 1 (100%), to use wider orbs: 1.1 - 10% wider, 1.5 - 50% wider and so on")
    parser.set_defaults(transorbis=1)
    
    # house system
    parser.add_option("-s", "--housesystem",  dest="house_sys", type="str",\
    help="Select house system: P for Placidus, K for Koch")
    parser.set_defaults(house_sys="P")
    
    # analyse all or only Asc/Dsc and MC/IC
    parser.add_option("-c", "--cusps",  dest="cusps_selection", type="int",\
    help="Default is 0 = calculate aspects to all cusps, 1 = only to ASC/Dsc and MC/IC, 2 = only cusps, given in comments (e.g. 5 if event is child birth, 8 if event is related to death), and 3 = to ASC/Dsc and MC/IC and cusps in comments")
    parser.set_defaults(cusps_selection=1)

    # choose ruling system
    parser.add_option("-r", "--rule",  dest="ruling", type="str",\
    help="ruling systems: f - full, s - septener, d - day, n - night, m - my style")
    parser.set_defaults(ruling="m")
    
    parser.add_option("-p", "--personal_info",  dest="print_pers_info", type="int",\
    help="print date of birth at the bottom of the calendar. Default is 0 = do Not print birt data, just print natal positions")
    parser.set_defaults(print_pers_info=0)
    
    parser.add_option("-f", "--ficti",  dest="fictif", type="int",\
    help="default: not to calculate transits of Node, Lilith and Chiron")
    parser.set_defaults(fictif=0)
    
    parser.add_option("-o", "--output",  dest="output_file", metavar="<file>",\
    help="name for output graph")
    parser.set_defaults(output_file='transitkalendar_years')
    
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
SEFLG_SWIEPH = 2      # use SWISSEPH ephemeris  
SEFLG_SPEED = 256        # high precision speed

SE_ASC = 0
SE_MC = 1
SE_ARMC = 2
SE_VERTEX = 3
SE_EQUASC = 4    # "equatorial ascendant"
SE_COASC1 = 5    # "co-ascendant" (W. Koch)
SE_COASC2 = 6    # "co-ascendant" (M. Munkasey)
SE_POLASC = 7    # "polar ascendant" (M. Munkasey)
SE_NASCMC = 8

SEFLG_EQUATORIAL = (2*1024)    # equatorial positions are wanted







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
major_aspects = [0, 60, 90, 120, 180]



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
'true Node': u'\u260A', 'mean Node': u'\u260A', 'mean Apogee':u'\u26B8', 'Chiron':u'\u26B7', 'Rahu': u'\u260A', 'Ketu': u'\u260B'}

uni_pls_2 = {'SU': u'\u2609', 'MO': u'\u263D', 'ME': u'\u263F', 'VE': u'\u2640', 'MA': u'\u2642', \
'JU': u'\u2643', 'SA': u'\u2644', 'UR': u'\u2645', 'NE': u'\u2646', 'PL': u'\u2647', \
'NN': u'\u260A', 'Li': u'\u26B8', 'Ch': u'\u26B7', 'SN': u'\u260B'}

uni_asps = {0: u'\u260C', 60: u'\u2731', 90: u'\u25A1', 120: u'\u25B3', 180: u'\u260D',\
20: '20'+u'\u00B0', 30: u'\u26BA', 36: '36'+u'\u00B0', 40: '40'+u'\u00B0', 45: '45'+u'\u00B0', 72: '72'+u'\u00B0', 80: '80'+u'\u00B0', \
100: '100'+u'\u00B0', 108: '108'+u'\u00B0', 135: '135'+u'\u00B0', 144: '144'+u'\u00B0', 150: u'\u26BB'}
# sextile replace by heavy asterisk       
 
    
    
    
    
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
    


def find_tr_asp(asp_exact, trans_orbis):
            
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
    ok_track = 1
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
                            if from_when > max_coord+2*otstup or till_when < min_coord-2*otstup:
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
    
    if tr_pl_name=='true Node' and asp_type==180:
        tr_pl_name='Ketu'
        asp_type=0
        
    if n_pl_name=='true Node' and asp_type==180:
        n_pl_name='Ketu'
        asp_type=0
    
    text_position = from_when+(till_when-from_when)/2 # in the middle of an aspect line
    if text_position>max_days:
          text_position = max_days-otstup
    elif text_position<min_days:
          text_position = min_days+otstup
    plt.text(ok_track, text_position-otstup, uni_pls[tr_pl_name],  size=10)
    
    
    lwidth = 1
    
    opp_cusp = {'1':'7','2':'8','3':'9','4':'10','5':'11','6':'12','7':'1','8':'2','9':'3','10':'4','11':'5','12':'6'}
    
    # if transit is to cusp
    if n_pl_name in ['1','2','3','4','5','6','7','8','9','10','11','12']:
       if asp_type==180:
          n_pl_name = opp_cusp[n_pl_name]
          asp_type = 0
       if asp_type==60:
          n_pl_name = opp_cusp[n_pl_name]
          asp_type = 120
       plt.text(ok_track, text_position+otstup, n_pl_name,  size=10)
       lwidth = 0.75
       linestyle_asp='--'
    else:
       try:
          plt.text(ok_track, text_position+otstup, uni_pls[n_pl_name],  size=10)
       except KeyError: #in case natal obj is a cusp
          plt.text(ok_track, text_position+otstup, n_pl_name,  size=10)
          
          
    if asp_type in [60,120]:
          col = 'g'
    elif asp_type in [90,180]:
          col = 'r'
    elif asp_type in [0]:
          col = 'b'
    elif asp_type in [30]:
          col = 'lightgreen'
    elif asp_type in [150]:
          col = 'orange'
    else:
          col = 'k'    
    plt.text(ok_track, text_position, uni_asps[asp_type],  size=10)
    
    
    
    if (tr_pl_name, n_pl_name) in natal_asps or (n_pl_name,tr_pl_name) in natal_asps:
          lwidth = 2
          
    if asp_type in major_aspects:
       linestyle_asp='-'
    else:
       linestyle_asp='--'
                             
    plt.plot(x, y, color=col, linewidth=lwidth, alpha=0.9, linestyle=linestyle_asp)
    
    
    
    
    
    
    
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
        s_transit = records[1].split(';')
    except IOError:
        print ("ERROR: cannot open or read input file:", input_file)
        exit(-1)
    
    
    
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
    
    
    # ---------------- set SIGN RULING system ---------------
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




    e_list = data    
        
    # check the cusps' selection: all (0) or just Asc/Dsc and MC/IC (1)
    if options.cusps_selection==0:
        e_houses = [1,2,3,4,5,6]
                
    elif options.cusps_selection==1:
        e_houses = [1,10]
                
    elif options.cusps_selection==2: # only those that are in the comment
        e_houses = e_list[8].strip().split(',')
        e_houses = [x for x in e_houses if x]
        #if len(e_houses):
            #e_houses = [int(i) for i in e_houses]
                
    elif options.cusps_selection==3:
        e_houses = e_list[8].strip().split(',')
        e_houses = [int(x) for x in e_houses if x]
                
        #if len(e_houses):
        #    e_houses = [int(i) for i in e_houses]
        e_houses = [1,10] + e_houses
    e_houses = list(set(e_houses))
    
    print (e_houses)    


    
    
    # --------------- Start: getting TRANSIT data (first given time point) --------------------
    #name = s_transit[0]
    
    # ------- date of starting transit --------
    tr_date = s_transit[1].split('.')
    tr_year, tr_month, tr_day = int(tr_date[2]), int(tr_date[1]), int(tr_date[0])

    # ------- time of first transit --------
    tr_time = s_transit[2].split(':')
    tr_hour = float(tr_time[0])
    
    tr_hour, tr_minute, tr_second = 0.0, 0.0, 0.0 # transit starts at midnight
    
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

    

    
    n_signs = {1: 'Ari', 2: 'Tau', 3: 'Gem', 4: 'Can', 5: 'Leo',\
6: 'Vir', 7: 'Lib', 8: 'Sco', 9: 'Sag', 10: 'Cap',\
11: 'Aqu', 12: 'Pis'}
        
    planets = [0,1,2,3,4,5,6,7,8,9,11,12,15]    # 11 is mean node, 12 lilith=mean apogee, 15=chiron
    
    n_pls = {0: 'Sun', 1:'Moon', 2:'Mercury', 3:'Venus', 4:'Mars', 5:'Jupiter', 6:'Saturn', 7:'Uranus', 8:'Neptune', 9:'Pluto', 11:'Mean Node', 12:'Lilith', 15:'Chiron'}
    #planets = [0,1,2,3,4,5,6,7,8,9,10]    # 10 is true node
    
    
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
            
                #p = swe_get_planet_name(p1)
                is_asp, asp_type, diff_with_exact = find_asp(p1, diff)
                    
                if is_asp:
                    if asp_type in Aspects:
                        natal_asps.append((p1,swe.get_planet_name(p2)))
                    
    #print natal_asps, '\n'
    
  
    
    
    # -------------- SETTING PLOT PARAMETERS ----------------
    fig = plt.figure()    
    ax = fig.add_subplot(111)
                
    #ax = plt.axes()
    plt.xticks([])
    
    #plt.xticks(timeticks, xticklabels, rotation=60, fontsize=7)
    #diffe = (timelim2 - timelim1) / 7
    
    interval = options.time_interval # in days
    min_years, max_years = 0, options.time_interval
    otstup = 0.02*max_years
                
    ax.set_xlim(0,41)
    ax.set_ylim(interval, 0)
    
    year_label = [0]
    
    # add lines to separate days
    for i in range(1,max_years+1):
        plt.axhline(y=i, xmin=0, xmax=1, alpha=0.2) # xmin and xmax take only relative val 0 and 1 (1=whole range, 100%)
        #year_label.append(i-0.5)
        year_label.append(i)
    
                        
    ok_track = 1
    
    track_order = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40]
    track_occupancy = {1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],10:[],11:[],12:[],13:[],14:[],15:[],16:[],17:[],18:[],19:[],20:[],21:[],22:[],23:[],24:[],25:[],26:[],27:[],28:[],29:[],30:[],31:[],32:[],33:[],34:[],35:[],36:[],37:[],38:[],39:[],40:[]}
    
    # ----------- END preparing to plot ------------
    
    
    
    
    #weekdays = {0:'Mon',1:'Tue',2:'Wed',3:'Thu',4:'Fri',5:'Sat',6:'Sun'}
    months = {1:'January', 2:'February',3:'March',4:'April',5:'May',6:'June',\
    7:'July', 8:'August',9:'September',10:'October',11:'November',12:'December'}
    
    
    
    
    

    # -------- START ----- calculate TRANSIT 1 ------------
    tr1_year, tr1_month, tr1_day, tr1_time = tr_year, tr_month, tr_day, tr_time          

    # -------------horoscope--------------
    tr1_jd = swe.julday(tr1_year, tr1_month, tr1_day, tr1_time, calflag) # julian date
    tr1_d = swe.deltat(tr1_jd)
    tr1_obl, tr1_rflag = swe.calc(tr1_jd+tr1_d, SE_ECL_NUT, 0)    

          
    
    ylabels = [months[tr_origmonth][:3]+' '+str(tr_origyear)]
    
    
    otstup = 0.015*max_years
    
    # -------- MAKE LABELS USING END OF THE PERIOD ------------
    tr2_year, tr2_month, tr2_day, tr2_time = tr_origyear, tr_origmonth, tr_origday, 0.0
    
    year_list = [tr2_year] # for joining repeating transits later
    
    for i in range(1,interval+1): # +1 because we need not 00:00, but 23:59:59
        tr2_year = tr2_year+1
        ylabels.append(months[tr2_month][:3]+' '+str(tr2_year))
        year_list.append(tr2_year)
    
    plt.yticks(year_label, ylabels)
    #plt.title('Directions: '+months[tr_origmonth][:3]+' '+str(tr_origyear)+'-'+str(tr2_year), loc='left')
    
    plt.text(0, -otstup*2, 'Natal: '+ data[1] + ', ' +  data[2],  size=12)
    plt.text(0, -otstup, 'Transits: '+months[tr_origmonth][:3]+' '+str(tr_origyear)+'-'+str(tr2_year),  size=12)

    

    

    # ------------- DAY BY DAY TRANSIT CALCULATION FOR ALL TRANSITS ------------------
    transits_in = {}    
    transits_in_track = {}
    step = 1
    
    natal_obj = [0,1,2,3,4,5,6,7,8,9,11,12,15]
    #transit_pls = [9,8,7,12,11,15,6,5,4] # exclude Sun,Moon,Merc,Venus,Mars 
    
    if options.fictif!=0:
        transit_pls = [9,8,7,6,5,15,11,12] # slow moving pls, and social pls, and node, lilith, chiron
    else:
        transit_pls = [9,8,7,6,5] # exclude all fict objects
        
       
        
    natal_cusps = e_houses
    
    
    
    # -------- TRANSIT 1 DECREASED BY 1 hours (to check planets direction)------------
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
    
    
    
    planet_list = ['Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']
    
    # --------------- START: CALCULATING TRANSITS OF SLOW PLANETS AND PLOTTING ----------------------------
    # --------------- checking direction of the slow planets ------------------------
    pl_direction = {0:'D',1:'D',2:'D',3:'D',4:'D',5:'D',6:'D',7:'D',8:'D',9:'D',11:'D',12:'D',15:'D'}

    tr2_year, tr2_month, tr2_day, tr2_time = tr1_year, tr1_month, tr1_day, tr1_time
    from_when, till_when = min_years, max_years

    max_months = max_years*12
    max_days = max_years*365
    step = 1
    
    #print '\nAll transits\n'
    
    for i in range(0,max_days+1,step): # i = every 1 month should be enough
    
        
        tr2_year, tr2_month, tr2_day = util.incrDay(tr2_year, tr2_month, tr2_day)
        
        tr2_year, tr2_month, tr2_day, tr2_time = check_overflow(tr2_year, tr2_month, tr2_day, tr2_time)    
        
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
        for p in n_pls:
            if isinstance(p, int): # if it is a planet incoded with int
                pos1 = tr1_pls[p]
                pos2 = tr2_pls[p]
                
                if pos2 > pos1 and abs(pos2-pos1)<180: # D
                        pl_direction[p]='D'
                elif pos2 < pos1 and abs(pos2-pos1)>180: # pisces-aries step
                        pl_direction[p]='D'
                else:
                        pl_direction[p]='R'
            
        #print 'pl status:', pl_direction, '\n'
    
    

    
        # any aspects with natal?
    
        for t in transit_pls:
            tr2_pl_name = swe.get_planet_name(t)
            tr2_pl_pos = util.normalize(tr2_pls[t])
            
            #tr2_pl_sign = n_signs[tr2_pl_pos//30+1]
        
        
            for n in natal_obj:
                update = False
                try:
                    n_pl_name = swe.get_planet_name(n)
                except TypeError:
                    n_pl_name = n
                n_pl_pos = n_pls[n]
        
                if tr2_pl_pos>=n_pl_pos:
                    diff2 = tr2_pl_pos - n_pl_pos
                    if diff2>=180: # further apart or at the boundary of pisces-aries
                        diff2 = abs(tr2_pl_pos-360) + n_pl_pos
                else:
                    diff2 = n_pl_pos - tr2_pl_pos
                    if diff2>=180:
                        diff2 = abs(n_pl_pos - 360) + tr2_pl_pos
                    
                    
                # orbs for transits are the same
                is_asp2, asp_type2, diff_with_exact2 = find_tr_asp(diff2, trans_orbis)    
    
                if is_asp2:
                    
                    skip_aspect = False
                    
                    if tr2_pl_name in planet_list: # transit pl -> check natal pls
                        if n_pl_name in ['true Node', 'mean Apogee', 'Chiron']:
                            if n_pl_name in ['true Node'] and asp_type2 in [0,180]:
                                skip_aspect = False
                            elif n_pl_name in ['mean Apogee', 'Chiron'] and asp_type2 in [0]:
                                skip_aspect = False
                            else:
                                skip_aspect = True
                            
                    else: # NN, Li, Ch transiting
                        if tr2_pl_name in ['true Node'] and asp_type2 in [0,180] and n_pl_name in planet_list:
                            skip_aspect = False
                        elif tr2_pl_name in ['Chiron', 'mean Apogee'] and asp_type2 in [0] and n_pl_name in planet_list:
                            skip_aspect = False
                        else:
                            skip_aspect = True
                    
                
                    
                    
                    if skip_aspect==False:
                        #for_months = i
                        #if for_months==12:
                        #    for_months=0
                        
                        where_asp = float(i)/365
                        #print pl_direction[t], tr2_pl_name, n_pl_name, asp_type2, where_asp, i
                        
                        if (pl_direction[t], tr2_pl_name, n_pl_name, asp_type2) not in transits_in:
                            transits_in[(pl_direction[t], tr2_pl_name, n_pl_name, asp_type2)] = [where_asp]
                                
                        else: # update
                            asps_coords = transits_in[(pl_direction[t], tr2_pl_name, n_pl_name, asp_type2)]
                            asps_coords.append(where_asp)
                            transits_in[(pl_direction[t], tr2_pl_name, n_pl_name, asp_type2)] = asps_coords
                                            
                        
                    
            # --------- transits to natal cusps
            
            for c in natal_cusps:
                update = False
                
                n_cusp_pos = cusps[c-1]
        
                if tr2_pl_pos>=n_cusp_pos:
                    diff2 = tr2_pl_pos - n_cusp_pos
                    if diff2>=180: # further apart or at the boundary of pisces-aries
                        diff2 = abs(tr2_pl_pos-360) + n_cusp_pos
                else:
                    diff2 = n_cusp_pos - tr2_pl_pos
                    if diff2>=180:
                        diff2 = abs(n_cusp_pos - 360) + tr2_pl_pos                                        
                    
                    
                # orbs for transits are the same
                is_asp2, asp_type2, diff_with_exact2 = find_tr_asp(diff2, trans_orbis)    
    
                if is_asp2:
                                    
                    skip_aspect = False
                    
                    if tr2_pl_name in planet_list: # transit pl -> check natal pls
                        if c in [1, 10]: # any asp to ASC or MC
                            skip_aspect = False
                        else:
                            if asp_type2 in [0]:
                                skip_aspect = False
                            else:
                                skip_aspect = True
                            
                    else: # NN, Li, Ch transiting
                        
                        if asp_type2 in [0]:
                            skip_aspect = False
                        else:
                            skip_aspect = True
                    
                    
                    if skip_aspect==False:
                        #for_months = i
                        #if for_months==12:
                        #    for_months=0
                        where_asp = float(i)/365
                        
                        if (pl_direction[t], tr2_pl_name, str(c), asp_type2) not in transits_in:
                            transits_in[(pl_direction[t], tr2_pl_name, str(c), asp_type2)] = [where_asp]
                        else: # update
                            asps_coords = transits_in[(pl_direction[t], tr2_pl_name,  str(c), asp_type2)]
                            asps_coords.append(where_asp)
                            transits_in[(pl_direction[t], tr2_pl_name,  str(c), asp_type2)] = asps_coords        
                    
    
        tr1_pls = tr2_pls    
    
       

    
    # ----------------- plot slow transits ------------------------
    
    for t in transit_pls:
        tr_pl_name = swe.get_planet_name(t)
        
        tr_pl_asps = {}
        tr_pl_asps_plotted = {}
        single_asp_type = [] # (tr_pl,n_pl,asp_type)
        
        
        # ------- looking for aspect duplicates -------------
        
        for k in transits_in:
            
            if tr_pl_name==k[1]:
                tr_pl_asps[k] = transits_in[k]
                
                tr_pl_asps_plotted[(k[1],k[2],k[3])] = 0
                
                if (k[1],k[2],k[3]) not in single_asp_type:
                    single_asp_type.append((k[1],k[2],k[3]))
                
                #print t, single_asp_type
                
        for k in single_asp_type:

            if tr_pl_asps_plotted[k]==0: # that aspect has not been plotted yet
                
                tr_pl_D = ('D',k[0],k[1],k[2])
                tr_pl_R = ('R',k[0],k[1],k[2])
                                    
                                
                
                if tr_pl_D in tr_pl_asps and tr_pl_R in tr_pl_asps: 
                    
                    # combine direct and retro periods
                    asp_coords = tr_pl_asps[tr_pl_D] + tr_pl_asps[tr_pl_R]
                    #print k, asp_coords
                    sorted_asp_coords = asp_coords
                    sorted_asp_coords.sort()
                    #print sorted_asp_coords
                    min_from_when = sorted_asp_coords[0]
                    
                    asp_repeats = []
                    one_asp_repeat_coords = [min_from_when]
                    
                    for a in range(1, len(sorted_asp_coords)):
                        if sorted_asp_coords[a]-min_from_when > 1.5: # then it is another period of that aspect
                            # add to asp_repeats
                            asp_repeats.append(one_asp_repeat_coords)
                            one_asp_repeat_coords = [sorted_asp_coords[a]]
                        else:
                            #print 'else'
                            one_asp_repeat_coords.append(sorted_asp_coords[a])
                        min_from_when = sorted_asp_coords[a]
                    
                    asp_repeats.append(one_asp_repeat_coords) # the last one manually   
                        
                    # print 'asp_repeats',asp_repeats
                    
                    
                    for asp_coords in asp_repeats: # consists of coords
                        from_when = min(asp_coords)
                        till_when = max(asp_coords)
                        
                        ok_track = find_available_track(track_order, track_occupancy,from_when,till_when, otstup) 
                        
                        plot_aspect(plt, natal_asps, asp_coords, ok_track, min_years, max_years, otstup, k)
                
                                               
                                        
                    tr_pl_asps_plotted[k]==1
                    
                else:
                    try:
                        asp_coords = tr_pl_asps[tr_pl_D]
                    except KeyError:
                        asp_coords = tr_pl_asps[tr_pl_R]
                
                    #print k, asp_coords
                    sorted_asp_coords = asp_coords
                    sorted_asp_coords.sort()
                    #print sorted_asp_coords
                    min_from_when = sorted_asp_coords[0]
                    
                    asp_repeats = []
                    one_asp_repeat_coords = [min_from_when]
                    
                    for a in range(1, len(sorted_asp_coords)):
                        if sorted_asp_coords[a]-min_from_when > 1.5: # then it is another period of that aspect
                            # add to asp_repeats
                            asp_repeats.append(one_asp_repeat_coords)
                            one_asp_repeat_coords = [sorted_asp_coords[a]]
                        else:
                            #print 'else'
                            one_asp_repeat_coords.append(sorted_asp_coords[a])
                        min_from_when = sorted_asp_coords[a]
                    
                    asp_repeats.append(one_asp_repeat_coords) # the last one manually   
                        
                    #print 'asp_repeats',asp_repeats
                    
                    
                    for asp_coords in asp_repeats: # consists of coords
                        from_when = min(asp_coords)
                        till_when = max(asp_coords)
                        
                        ok_track = find_available_track(track_order, track_occupancy,from_when,till_when, otstup) 
                        
                        plot_aspect(plt, natal_asps, asp_coords, ok_track, min_years, max_years, otstup, k)
                         
    # --------------- END: CALCULATING TRANSITS OF SLOW PLANETS AND PLOTTING ----------------------------    
        
        
        
     
        
        

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    # birth information
    if options.print_pers_info!=0:
        plt.text(0, max_years+otstup*1.2, 'Natal: '+ data[1] + ', ' +  data[2] + ', ' +data[4],  size=10)
    else:
        plt.text(0, max_years+otstup*1.2, 'Natal ',  size=10)
    
    
    # natal positions
    nat_pls_houses = calc_planet_houses(n_pls, cusps, options.ruling)
    
    nat_pls_title = ''
    nat_pls_upr = ''
    
    for t in planets:
        t_name = swe.get_planet_name(t)
        nat_pl_pos = util.normalize(n_pls[t])
        nat_pl_sign = n_signs[nat_pl_pos//30+1]

                
        t_house = nat_pls_houses[t_name][0]
        u_houses = nat_pls_houses[t_name][1]
        u_houses2 = [str(u) for u in u_houses]
        u_houses = ",".join(u_houses2)
        
        nat_pls_title = nat_pls_title + uni_pls[swe.get_planet_name(t)]+uni_signs[nat_pl_sign]+','+str(t_house)+'  '
        if t < 10:
            nat_pls_upr = nat_pls_upr + uni_pls[swe.get_planet_name(t)]+'u'+str(u_houses)+'  '
    
    
    plt.text(0, max_years+otstup*2.4, nat_pls_title, fontsize=10)
    plt.text(0, max_years+otstup*1.5*2.4, nat_pls_upr, fontsize=10)

    
    fig.set_size_inches(8, 12)
            
    fig.savefig(options.output_file+'_'+str(tr_origyear)+'_'+months[tr_origmonth][:3]+str(tr_origday)+'.png', dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close(fig)
    

    
    
    
        

if __name__ == "__main__":
    main()        

