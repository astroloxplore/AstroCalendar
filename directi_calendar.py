#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import unicode_literals


'''
This script calculates natal positions and directions, and makes a plot - a calendar.

Input text file example (ZBS format / ZET astroprocessor, first line - natal, 
second line - start point to calculate directions from, better to set to 1st of Jan):

t; 01.01.1992; 12:20:00; +3; Vilnius, Lithuania; 54n41; 25e19; -; ;
t; 1.1.2015; 00:00:00; +3; Vilnius, Lithuania; 54n41; 25e19; -; ;



Running the script:

python3 directi_calendar.py -i 'test_transit_calendar_years.zbs'


Running the script with additional parameters:

python3 directi_calendar.py -i 'test_transit_calendar_years.zbs' -d '0,1,2' -c 0
# -d -> directions of planets to cusps, of cusps to planets and of planets to planets
# -c 0 -> directions to all cusps

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
import sys

import matplotlib

def parse_commandline():
    parser = OptionParser()
    
    usage = "python3 directi_calendar.py -i input.zbs [options]"
    version = "1.0"
    description = "%prog makes a calendar of directions (years)"
    epilog = "Written by Svetlana Stoma / @astroloxplore (2022)"
    parser = OptionParser(usage=usage, description=description,
                          version="%prog "+version, epilog=epilog)
    # in file
    parser.add_option("-i", "--input",  dest="input_file", metavar="<file>",\
    help="Input file: record in ZBS format (ZET astroprocessor)")
    parser.set_defaults(input_file='input.zbs')
    
    
    # flags to select what aspects to calculate
    parser.add_option("-a", "--asp",  dest="asp_set", type="int",\
    help="Aspects: 0 for default set [0, 60, 90, 120, 180], 1 for major [0, 30, 45, 60, 90, 120, 135, 150, 180], 2 for all [0.0, 30.0, 36.0, 40.0, 45.0, 60.0, 72.0, 80.0, 90.0, 100.0, 108.0, 120.0, 135.0, 144.0, 150.0, 180.0], 3 for my default set [0, 60, 90, 120, 150, 180]")
    parser.set_defaults(asp_set=0)
    
    # INTERVAL and STEP, example 30 30 1 = 30 min back and 30 min forward, calc every 1 min
    parser.add_option("-t", "--time",  dest="time_interval", type="int",\
    help="A period of time (in years)")
    parser.set_defaults(time_interval=10) 
    
    # flags to select orbs
    parser.add_option("-b", "--transorb",  dest="orbis", type="float",\
    help="default is 1 (100%), to use wider orbs: 1.1 - 10% wider, 1.5 - 50% wider and so on")
    parser.set_defaults(orbis=1)
    
    # flags to select orbs
    parser.add_option("-y", "--dirorb",  dest="dirorbis", type="float",\
    help="default is 1 (+/-1 degree/year before exact aspect), to use stricter orbs: 0.5 (+/-0.5 degree/year)")
    parser.set_defaults(dirorbis=1.0)
    
    # house system
    parser.add_option("-s", "--housesystem",  dest="house_sys", type="str",\
    help="Select house system: P for Placidus, K for Koch")
    parser.set_defaults(house_sys="P")
    

    # analyse all or only Asc/Dsc and MC/IC
    parser.add_option("-c", "--cusps",  dest="cusps_selection", type="int",\
    help="Default is 0 = calculate aspects to all cusps, 1 = only to ASC/Dsc and MC/IC, 2 = only cusps, given in comments (e.g. 5 if event is child birth, 8 if event is related to death), and 3 = to ASC/Dsc and MC/IC and cusps in comments")
    parser.set_defaults(cusps_selection=1)

    
    # analyse D PL to N cusps (0), D cusps to N pl (1), PL to PL (2) or any combination of these
    parser.add_option("-d", "--ds",  dest="directions", type="str",\
    help="0 = calculate D planets to N cusps, 1 = D cusps to N pl, 2 = D planets to N planets, 3 = ingressions, 0,1,2,3 for all")
    parser.set_defaults(directions="0,1") # full "0 1 2", only pl to pl "2"
    
    # directions type
    parser.add_option("-e", "--dirtype",  dest="dir_type", type="int",\
    help="0 = symbolical directions (zodiac), 1 - Solar Arc (zodiac), 2 - Primary Directions ")
    parser.set_defaults(dir_type=1) # 
    
    # choose ruling system
    parser.add_option("-r", "--rule",  dest="ruling", type="str",\
    help="ruling systems: f - full, s - septener, d - day, n - night, m - my style")
    parser.set_defaults(ruling="m")
    
    parser.add_option("-o", "--output",  dest="output_file", metavar="<file>",\
    help="name for output graph")
    parser.set_defaults(output_file='directi_kalendar')
    
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
'Neptune': [330, 240], 'Pluto': [210, 0], 'true Node': [], 'mean Node': []}

# septener only
planet_rules_signs_sep = {'Sun': [120], 'Moon': [90], \
'Mars': [0, 210], 'Venus': [30, 180], \
'Mercury': [60, 150], 'Jupiter': [240, 330],\
'Saturn': [270, 300], 'Uranus': [], \
'Neptune': [], 'Pluto': [], 'true Node': [], 'mean Node': []}

# one sign strictly starting from Jupiter # day birth
planet_rules_one_sign_day = {'Sun': [120], 'Moon': [90], \
'Mars': [0], 'Venus': [30, 180], \
'Mercury': [60, 150], 'Jupiter': [240],\
'Saturn': [270], 'Uranus': [300], \
'Neptune': [330], 'Pluto': [210], 'true Node': [], 'mean Node': []}

# one sign strictly starting from Jupiter # night birth
planet_rules_one_sign_night = {'Sun': [120], 'Moon': [90], \
'Mars': [210], 'Venus': [30, 180], \
'Mercury': [60, 150], 'Jupiter': [330],\
'Saturn': [300], 'Uranus': [270], \
'Neptune': [240], 'Pluto': [0], 'true Node': [], 'mean Node': []}

# minors for both signs, vyshki dlia odnogo znaka
planet_rules_my_style = {'Sun': [120], 'Moon': [90], \
'Mars': [0,210], 'Venus': [30, 180], \
'Mercury': [60, 150], 'Jupiter': [240,330],\
'Saturn': [270], 'Uranus': [300], \
'Neptune': [330], 'Pluto': [210], 'true Node': [], 'mean Node': []}





# -------------- shifts and encodings ---------------------
uni_signs = {'Ari': u'\u2648', 'Tau': u'\u2649', 'Gem': u'\u264A', 'Can': u'\u264B', \
'Leo': u'\u264C', 'Vir': u'\u264D', 'Lib': u'\u264E', 'Sco': u'\u264F',\
'Sag': u'\u2650', 'Cap': u'\u2651', 'Aqu': u'\u2652', 'Pis': u'\u2653'}

uni_pls = {'Sun': u'\u2609', 'Moon': u'\u263D', 'Mercury': u'\u263F', 'Venus': u'\u2640', 'Mars': u'\u2642', \
'Jupiter': u'\u2643', 'Saturn': u'\u2644', 'Uranus': u'\u2645', 'Neptune': u'\u2646', 'Pluto': u'\u2647', \
'true Node': u'\u260A', 'mean Node': u'\u260A', 'mean Apogee':u'\u26B8'}

uni_pls_2 = {'SU': u'\u2609', 'MO': u'\u263D', 'ME': u'\u263F', 'VE': u'\u2640', 'MA': u'\u2642', \
'JU': u'\u2643', 'SA': u'\u2644', 'UR': u'\u2645', 'NE': u'\u2646', 'PL': u'\u2647', \
'NN': u'\u260A', 'Li': u'\u26B8'}

#uni_asps = {0: u'\u260C', 60: u"\u26B9", 90: u'\u25A1', 120: u'\u25B3', 180: u'\u260D'}

uni_asps = {0: u'\u260C', 60: u'\u26B9', 90: u'\u25A1', 120: u'\u25B3', 180: u'\u260D',\
20: '20'+u'\u00B0', 30: u'\u26ba', 36: '36'+u'\u00B0', 40: '40'+u'\u00B0', 45: '45'+u'\u2220', 72: u'\u0051', 80: '80'+u'\u00B0', \
100: '100'+u'\u00B0', 108: '108'+u'\u00B0', 135: u'\u26bc', 144: u'\u0062'+u'\u0051', 150: u'\u26BB'}
# sextile can be replace by heavy asterisk    60: u'\u2731'
    
    
    
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

osi = {1:7, 2:8, 3:9, 4:10, 5:11, 6:12, 7:1, 8:2, 9:3, 10:4, 11:5, 12:6}



# make a flag what to return
def get_houses(tjd_ut, flag, geolat, geolon, hsys, obl):
    """Calculates the cusps of the Houses"""

    HOUSE_NUM = 12
    ASC, MC, ARMC, VERTEX, EQUASC, COASC, COASC2, POLARASC = range(0, 8)
    LON, LAT, RA, DECL = 0, 1, 2, 3
    obl = obl
    #res, cusps, ascmc = swe.houses_ex(tjd_ut, flag, geolat, geolon, ord(hsys))
    
    cusps, ascmc = swe.houses_ex(tjd_ut, geolat, geolon, bytes(hsys, 'utf-8'), flag)
    #cusps, ascmc = swe.houses_ex(tjd_ut, geolat, geolon, hsys.encode("ascii"), flag)
    #cusps, ascmc = swe.houses(tjd_ut, geolat, geolon, bytes(hsys, 'utf-8'))

    ascra, ascdecl, dist = swe.cotrans(ascmc[ASC], 0.0, 1.0, -obl)                
    mcra, mcdecl, dist = swe.cotrans(ascmc[MC], 0.0, 1.0, -obl)
    ascmc2 = ((ascmc[ASC], 0.0, ascra, ascdecl), (ascmc[MC], 0.0, mcra, mcdecl))
        
    #zdAsc=90.0, zdMC=0.0
    #poleAsc=lat, poleMC=0.0
    qasc = math.degrees(math.asin(math.tan(math.radians(ascdecl))*math.tan(math.radians(geolat))))
    regioMPAsc = ascra-qasc
    regioMPMC = mcra

    cuspstmp = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
    for i in range(0,HOUSE_NUM):
        #cuspstmp[i][0], cuspstmp[i][1], dist = swe.cotrans(cusps[i+1], 0.0, dist, -obl)
        cuspstmp[i][0], cuspstmp[i][1], dist = swe.cotrans(cusps[i], 0.0, dist, -obl)
            
    cusps2 = ((cuspstmp[0][0], cuspstmp[0][1]), (cuspstmp[1][0], cuspstmp[1][1]), (cuspstmp[2][0], cuspstmp[2][1]), (cuspstmp[3][0], cuspstmp[3][1]), (cuspstmp[4][0], cuspstmp[4][1]), (cuspstmp[5][0], cuspstmp[5][1]), (cuspstmp[6][0], cuspstmp[6][1]), (cuspstmp[7][0], cuspstmp[7][1]), (cuspstmp[8][0], cuspstmp[8][1]), (cuspstmp[9][0], cuspstmp[9][1]), (cuspstmp[10][0], cuspstmp[10][1]), (cuspstmp[11][0], cuspstmp[11][1]))

    #return cusps[1:], ascmc, ascmc2
    #print(cusps)
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
    #print name, data, dataEqu
        

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

        #speculums.append((data[LONG], data[LAT], dataEqu[RAEQU], dataEqu[DECLEQU], adlat, sa, md, hd, th, hod, pmp, adphi, phi, ao))
        #speculums.append((dataEqu[LONG], dataEqu[LAT], dataEqu[RAEQU], dataEqu[DECLEQU], adlat, sa, md, hd, th, hod, pmp, adphi, phi, ao))
        speculums.append((data[LONG], data[LAT], dataEqu[LONG], dataEqu[LAT]))
        #print data[LONG]
        #print dataEqu[LONG]

        return speculums

    


def calc_planet_houses(planets, pl_poss, cusps_pl_h, ruling_mode):
    planets_in_houses = {}
    
    # add first house at the end
    cusps_pl_h.append(cusps_pl_h[0])
    #print(len(cusps_pl_h))
    
    
    for pos in pl_poss:
        planet_house = []
        for i in range(0,12):
            cusp1, cusp2 = cusps_pl_h[i], cusps_pl_h[i+1]
            # check for circle beginning
            if cusp1 > cusp2: # cusp1 in pisces
                if pos < cusp2: # planet in Ari
                    planet_house.append(i+1)
                elif pos > cusp1: # planet in Pis
                    planet_house.append(i+1)
            else:
                if cusp1 <= pos and pos < cusp2: 
                    #print pos, i+1
                    planet_house.append(i+1)
                
        #print planet_house                                     
        # do check for house 13th aka 1st
        for h in planet_house:
            if h==13:
                planet_house[planet_house.index(h)]=1
        i = pl_poss.index(pos)
        planets_in_houses[swe.get_planet_name(planets[i])] = [planet_house[0]]
      
    
    # ruling system 
    if ruling_mode=='f':
            planet_rules = planet_rules_signs
    elif ruling_mode=='s':
            planet_rules = planet_rules_signs_sep
    elif ruling_mode=='d':
            planet_rules = planet_rules_one_sign_day
    elif ruling_mode=='n':
            planet_rules = planet_rules_one_sign_night
    else:
            planet_rules = planet_rules_signs
        
            
    # ruling
    for p in planets_in_houses:
        
        if ruling_mode=='n': # check where SUN is
            if planets_in_houses['Sun'] in [7,8,9,10,11,12]:
                planet_rules = planet_rules_one_sign_day
        
        rul_signs = planet_rules[p] # 1 or 2 signs
        rul_houses = []
        for s in rul_signs:
            # example Leo: 120, so s=120 for Sun
        
            for i in range (0,11):
                cusp1, cusp2 = cusps_pl_h[i], cusps_pl_h[i+1]
                
                
                # if cuspid gets into this sign
                if s <= cusp1 and cusp1 <= s+30:
                    rul_houses.append(cusps_pl_h.index(cusp1)+1)
                # check included signs
                if cusp1 <= s and s+30 <= cusp2: 
                    #print 'vkliucenyj znak mezdu kuspidami', cusp1, cusp2
                    rul_houses.append(cusps_pl_h.index(cusp1)+1)
                    
                    
                if cusp2 < cusp1: # the edge 360-0
                    # variant 1: cusp1 is in Aqu/Pis and cusp2 is in Ari
                    if cusp1 <= s and util.normalize(s+30) <= cusp2:
                        rul_houses.append(cusps_pl_h.index(cusp1)+1) # Pis included
                    # variant 2: cusp1 is in Pis and cusp2 is in Tau
                    if s <= cusp1 and util.normalize(s+30) <= cusp2:
                        rul_houses.append(cusps_pl_h.index(cusp1)+1) # Ari included
                    
                    
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
    pl1_orbs = orbs[pl1]
    asp_type = None
    diff_with_exact = None        
            
    for a in pl1_orbs: # for every type of aspect get ranges
        ranges = []
        orbis = pl1_orbs[a]
                
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
    


def find_tr_asp(asp_exact, orbis):
            
    is_asp = False
    #pl1_orbs = orbs[pl1]
    asp_type = None
    diff_with_exact = None        
            
    #for a in pl1_orbs: # for every type of aspect get ranges
    for a in Aspects:
        ranges = []
        #orbis = pl1_orbs[a]
                
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
    



def calc_pd_diff(jd, pd_jd, calflag):
        
        diff = pd_jd - jd # diff in days
        d_full_years = int(diff//365.2421904) # full years
        
        ra = swe.revjul(jd, calflag)
        ra = list(ra)
        ra[0] = ra[0]+d_full_years
        jd_with_full_yy =  swe.julday(ra[0], ra[1], ra[2], ra[3], calflag)
        
        diff = pd_jd - jd_with_full_yy # diff in days
        days_in_last_year = swe.julday(ra[0]+1, ra[1], ra[2], ra[3], calflag) - jd_with_full_yy # days from last bday to the pd date
        
        d_part_year = float(diff%365.2421904)/days_in_last_year # rest of the year

        #print d_full_years,d_part_year
        return d_full_years+d_part_year
            
        



def get_directions(jd, pd_jd, pls, cusps, calflag, dir_type, solar_arc):    

        diff = calc_pd_diff(jd, pd_jd, calflag)
        
        
        if dir_type==0 or dir_type==1: # if symbolical directions or Solar Arc (ecliptic)
            diff = diff * solar_arc
        elif dir_type==2: 
            diff = diff 
        else:
            pass
        
        pd_planets, pd_cusps = [], []
            
        if diff>=0:
            # planet directions (PD positions on given date)
            for pl in pls:
                pd_planets.append(pl[0][0] + diff) # planet position at birth + pd
            for cu in cusps:
                pd_cusps.append(cu + diff)
        return pd_planets, pd_cusps    
            



    
    
    
def calc_pd_n_aspects(direct, natal, rflag, pls_in_h, Aspects, dir_orbis):
        
    n_signs = {1: 'Ari', 2: 'Tau', 3: 'Gem', 4: 'Can', 5: 'Leo',\
6: 'Vir', 7: 'Lib', 8: 'Sco', 9: 'Sag', 10: 'Cap',\
11: 'Aqu', 12: 'Pis'}
    default_asps = [0, 60, 90, 120, 180]    
    
    d_direct = {}
    
    for d in direct:
        for n in natal:
            #calc diff
            d_norm = util.normalize(d)
            n_norm = util.normalize(n)
            if d>=n:
                diff = d-n
                if diff>=180: # further apart or at the boumdary of pisces-aries
                    diff = abs(d-360) + n
            else:
                diff = n-d
                if diff>=180:
                    diff = abs(n-360) + d                                        
            
                    
            # check asp
            for a in Aspects:
                
                orbis = 1.0
                
                if a in default_asps: # if SUN or Moon then orbis changes
                    if rflag==0: # pl to pl
                        d_planet = direct.index(d) # SUN is at 0, Moon at 1
                        n_planet = natal.index(n)
                        if d_planet == 0 or n_planet == 0:
                            orbis = 2.0
                        elif d_planet == 1 or n_planet == 1:
                            orbis = 1.5
                        else:
                            orbis = 1.0
                        
                    elif rflag==1: # planet to cusp    
                        d_planet = direct.index(d) # SUN is at 0, Moon at 1
                        if d_planet == 0:
                            orbis = 2.0
                        elif d_planet == 1:
                            orbis = 1.5
                        else:
                            orbis = 1.0
                    else: # cusp to planet    
                        n_planet = natal.index(n) # SUN is at 0, Moon at 1
                        if n_planet == 0:
                            orbis = 2.0
                        elif n_planet == 1:
                            orbis = 1.5
                        else:
                            orbis = 1.0
                        
                        
                orbis = orbis * dir_orbis # elsi nuzhny bolee tochnye aspekty    
            
                if a-orbis <= diff and diff<= a+orbis:
                    
                    d_sign = n_signs[d_norm//30+1]
                    n_sign = n_signs[n_norm//30+1]
                
                    if rflag==0: # planet to planet
                        d_pl_name = swe.get_planet_name(direct.index(d))
                        d_pl_info = pls_in_h[d_pl_name]
                        n_pl_name = swe.get_planet_name(natal.index(n))
                        n_pl_info = pls_in_h[n_pl_name]
                        #print 'D', d_pl_name, d_sign, ' ', round(diff,2), ' ', n_pl_name, n_sign, n_pl_info
                        # dict[key is (D pl, N pl, aspect)] = [(sign for D pl), (sign, house, and ruled houses for N pl), (float aspect)]
                        d_direct[(d_pl_name, n_pl_name, a)] = [[(d_sign, d_pl_info), (n_sign, n_pl_info), (diff)]]
                        
                    elif rflag==1: # planet to cusp    
                        d_pl_name = swe.get_planet_name(direct.index(d))
                        d_pl_info = pls_in_h[d_pl_name]
                        ncuspid = natal.index(n)+1
                        #print 'D', d_pl_name, d_sign, ' ', round(diff,2), ' ', 'K'+str(ncuspid)
                        # dict[key is (D pl, cusp, aspect)] = [(sign for D pl), (n_sign), (float aspect)]
                        d_direct[(d_pl_name, ncuspid, a)] = [[(d_sign, d_pl_info), (n_sign), (diff)]]
                    
                    else: # cusp to planet
                        dcuspid = direct.index(d)+1
                        n_pl_name = swe.get_planet_name(natal.index(n))
                        n_pl_info = pls_in_h[n_pl_name]
                        #print 'D', 'K'+str(dcuspid), ' ', round(diff,2), ' ', n_pl_name, n_sign, n_pl_info
                        # dict[key is (D cusp, N pl, aspect)] = [(d_sign), (sign, house, and ruled houses for N pl), (float aspect)]
                        d_direct[(dcuspid, n_pl_name, a)] = [[(d_sign), (n_sign, n_pl_info), (diff)]]
                    break
    #print d_direct
    return d_direct









    
    
    
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
    
    
    
    trans_orbis = options.orbis
    dir_orbis = options.dirorbis
    
    
    # ---------------- set ASPECTS set -----------------
    # default set Aspects = [0, 60, 90, 120, 180]

    if options.asp_set==1:
        Aspects = [0, 30, 45, 60, 90, 120, 135, 150, 180] # major
    elif options.asp_set==2:
        Aspects = [0.0, 30.0, 36.0, 40.0, 45.0, 60.0, 72.0, 80.0, 90.0, 100.0, 108.0, 120.0, 135.0, 144.0, 150.0, 180.0]
    elif options.asp_set==3:
        Aspects = [0, 60, 90, 120, 150, 180]
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

    
    pl_to_cu, cu_to_pl, pl_to_pl, ingressions = False, False, False, False
    d_types = options.directions.split(',')
    if "0" in d_types:
        pl_to_cu = True
    if "1" in d_types:
        cu_to_pl = True
    if "2" in d_types:
        pl_to_pl = True
    if "3" in d_types:
        ingressions = True
        
    
    
    
    e_list = data    
        
    # check the cusps' selection: all (0) or just Asc/Dsc and MC/IC (1)
    if options.cusps_selection==0:
        e_houses = [1,2,3,4,5,6]
                
    elif options.cusps_selection==1:
        e_houses = [1,10]
                
    elif options.cusps_selection==2: # only those that are in the comment
        e_houses = e_list[8].strip().split(',')
        e_houses = [x for x in e_houses if x]
        if len(e_houses):
            e_houses = [int(i) for i in e_houses]
                
    elif options.cusps_selection==3:
        e_houses = e_list[8].strip().split(',')
        e_houses = [x for x in e_houses if x]
                
        if len(e_houses):
            e_houses = [int(i) for i in e_houses]
        e_houses = [1,10] + e_houses
    e_houses = list(set(e_houses))
    
    print ('Natal:',data)
    print ('Selected cusps', e_houses)
    print ('Direction from:', s_transit)
    
    # --------------- Start: getting DIRECTION start data --------------------
    #name = s_transit[0]
    
    # ------- date of starting transit --------
    tr_date = s_transit[1].split('.')
    tr_year, tr_month, tr_day = int(tr_date[2]), int(tr_date[1]), int(tr_date[0])

    # ------- time of first transit --------
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
    planets = [0,1,2,3,4,5,6,7,8,9]    # 11 is mean node, 12 lilith
    pls_names = ['Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn','Uranus', 'Neptune', 'Pluto']
    pls_names_short = ['SU', 'MO', 'ME', 'VE', 'MA', 'JU', 'SA', 'UR', 'NE', 'PL', 'NN', 'Li']

    
    natal_asps = []



    # ------------ calculate NATAL ----------------------
    
    year, month, day, time = check_overflow(dyear, dmonth, dday, time)    
    #print 'year, month, day, time', year, month, day, time
                        
 
    # --------------- horoscope ----------------
    jd = swe.julday(year, month, day, time, calflag) # julian date
    d = swe.deltat(jd)
    obl, rflag = swe.calc(jd+d, SE_ECL_NUT, 0)

    # ------------- houses --------------
    cusps, ascmc, ascmc2 = get_houses(jd, hflag, lat, lon, hsys, obl[0])
    #print(cusps)
    raequasc, declequasc, dist = swe.cotrans(ascmc[SE_EQUASC], 0.0, 1.0, -obl[0])
    options_meannode = True 
        
    
    # -------------planets--------------        
    pls = []    
    for pId in planets:
        pl = get_planet(jd, pId, pflag, lat, ascmc2, raequasc, nolat, obl[0])
        pls.append(pl)
        
    # ---------- check degree overflow ----------
    n_pls = [util.normalize(pl[0][0]) for pl in pls]
    cusps = [util.normalize(c) for c in cusps]
    
    # ---------planet's house and houses it rules----------
    pls_in_h = calc_planet_houses(planets, n_pls, list(cusps), options.ruling)
   
   

    l = len(n_pls)
    i=0
    for i in range(0,10):
        d = n_pls[i]
        for k in range(i+1,l):
            n = n_pls[k]
            if d>=n:
                diff = d-n
                if diff>=180: # further apart or at the boundary of pisces-aries
                    diff = abs(d-360) + n
            else:
                diff = n-d
                if diff>=180:
                    diff = abs(n-360) + d    

            is_asp, asp_type, diff_with_exact = find_asp(pls_names_short[i], diff)
            
            if is_asp:
                natal_asps.append((pls_names_short[i],pls_names_short[k]))
            
        i+=1

    
    
    
    
    # -----------------SOLAR ARCS--------------------
    solar_arc = 1.0 # one degree of Solar Arc    
    if options.dir_type==1:
        # calculate solar arc
        # decrease natal by one day
        ddyear, ddmonth, ddday = util.decrDay(dyear, dmonth, dday)
        pyear, pmonth, pday, ptime = check_overflow(ddyear, ddmonth, ddday, time)    
                    
        # -------------horoscope--------------
        solar_jd = swe.julday(pyear, pmonth, pday, ptime, calflag) # julian date
        solar_d = swe.deltat(solar_jd)
        solar_obl, solar_rflag = swe.calc(solar_jd+solar_d, SE_ECL_NUT, 0)    

        # -------------houses--------------
        solar_cusps, solar_ascmc, solar_ascmc2 = get_houses(solar_jd, hflag, lat, lon, hsys, solar_obl[0])
        solar_raequasc, solar_declequasc, solar_dist = swe.cotrans(solar_ascmc[SE_EQUASC], 0.0, 1.0, -solar_obl[0])
        options_meannode = True 

    
        # -------------planets--------------        
        solar_pls = []    
        for pId in planets:
            pl = get_planet(solar_jd, pId, pflag, lat, solar_ascmc2, solar_raequasc, nolat, solar_obl[0])
            solar_pls.append(pl)
        
        # ---------- check degree overflow ----------
        solar_n_pls = [util.normalize(pl[0][0]) for pl in solar_pls]    
        solar_arc = abs(solar_n_pls[0] - n_pls[0])
        print ('Solar Arc', solar_arc, solar_n_pls[0], n_pls[0])
        
    elif options.dir_type==2:
        pass
        # calculate primary directions # to do
    else:
        solar_arc = 1.0
    
    
    
    

    
    # -------------- PREPARING THE PLOT ----------------
    fig = plt.figure()    
    ax = fig.add_subplot(111)
                
    plt.xticks([])
    
    interval = options.time_interval # in years
    min_years, max_years = 0, options.time_interval
                
    ax.set_xlim(0,46)
    ax.set_ylim(interval, 0)
    
    year_label = [0]
    
    # add lines to separate days
    for i in range(1,max_years+1):
        plt.axhline(y=i, xmin=0, xmax=1, alpha=0.2) # xmin and xmax take only relative val 0 and 1 (1=whole range, 100%)
        #year_label.append(i-0.5)
        year_label.append(i)
    
                        
    free_track = 1
    
    track_order = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45]
    track_occupancy = {1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],10:[],11:[],12:[],13:[],14:[],15:[],16:[],17:[],18:[],19:[],20:[],21:[],22:[],23:[],24:[],25:[],26:[],27:[],28:[],29:[],30:[],31:[],32:[],33:[],34:[],35:[],36:[],37:[],38:[],39:[],40:[],41:[],42:[],43:[],44:[],45:[]}
    
    # ----------- END preparing to plot ------------
    
    
    
    #weekdays = {0:'Mon',1:'Tue',2:'Wed',3:'Thu',4:'Fri',5:'Sat',6:'Sun'}
    months = {1:'January', 2:'February',3:'March',4:'April',5:'May',6:'June',\
    7:'July', 8:'August',9:'September',10:'October',11:'November',12:'December'}
    
    
    
    # -------- START ----- DIRECTION 1 ------------
    tr1_year, tr1_month, tr1_day, tr1_time = tr_year, tr_month, tr_day, tr_time            
    
    # -------------horoscope--------------
    tr1_jd = swe.julday(tr1_year, tr1_month, tr1_day, tr1_time, calflag) # julian date
    tr1_d = swe.deltat(tr1_jd)
    tr1_obl, tr1_rflag = swe.calc(tr1_jd+tr1_d, SE_ECL_NUT, 0)    

  
    
    ylabels = [months[tr_origmonth][:3]+' '+str(tr_origyear)]
    
    
    otstup = 0.015*max_years
    
    # -------- MAKE LABELS USING END OF THE PERIOD ------------
    tr2_year, tr2_month, tr2_day, tr2_time = tr_origyear, tr_origmonth, tr_origday, 0.0
    
    for i in range(1,interval+1): # +1 because we need not 00:00, but 23:59:59
        tr2_year = tr2_year+1
        ylabels.append(months[tr2_month][:3]+' '+str(tr2_year))
    
    plt.yticks(year_label, ylabels)
    #plt.title('Directions: '+months[tr_origmonth][:3]+' '+str(tr_origyear)+'-'+str(tr2_year), loc='left')
    
    plt.text(0, -otstup*2, 'Natal: '+ data[1] + ', ' +  data[2],  size=12)
    plt.text(0, -otstup, 'Directions: '+months[tr_origmonth][:3]+' '+str(tr_origyear)+'-'+str(tr2_year),  size=12)

    
    

    transit_pls = planets
                    
    transits_in = []    
    transits_in_track = {}
    moon_transits_in = []
    step = 1
    
    natal_obj = [0,1,2,3,4,5,6,7,8,9,11,12, 'Asc', 'MC']
    

    
    

    

    # ------------- DAY BY DAY DIRECTION CALCULATION FOR ALL DIR ------------------
    # all commented as 'TRANSITs', cause it is based on transit calendar script
    
    
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
    tr1_pls = []    
    for pId in planets:
        pl = get_planet(tr1_jd, pId, pflag, tr_lat, tr1_ascmc2, tr1_raequasc, nolat, tr1_obl[0])
        tr1_pls.append(pl)
        
    # ---------- check degree overflow ----------
    tr1_pls = [util.normalize(pl[0][0]) for pl in tr1_pls]
    
    
    
    
    #transit_pls = [9,8,7,12,6,5,4,3,2,0] # exclude Moon=1 for now 
    pl_direction = {0:'D',1:'D',2:'D',3:'D',4:'D',5:'D',6:'D',7:'D',8:'D',9:'D',11:'D',12:'D'}
    #natal_obj = [0,1,2,3,4,5,6,7,8,9,11,12, 'Asc', 'MC']            
    
    transits_in = {}    
    transits_in_track = {}
    
    
    tr2_year, tr2_month, tr2_day, tr2_time = tr1_year, tr1_month, tr1_day, tr1_time
    from_when, till_when = min_years, max_years

    max_months = max_years*12
    step = 1
    
    #print '\nAll directions\n'
    
    for i in range(0,max_months+1,step): # i = every 1 month should be enough
    
        tr2_year, tr2_month = util.incrMonth(tr2_year, tr2_month)
        
        
        tr2_year, tr2_month, tr2_day, tr2_time = check_overflow(tr2_year, tr2_month, tr2_day, tr2_time)    
        
        pd_jd = swe.julday(tr2_year, tr2_month, tr2_day, 12, calflag) # default time 12 pm
        
        # -----get direction time and approx PL or CUSP directed position---
        d_pls, d_cusps = get_directions(jd, pd_jd, pls, cusps, calflag, options.dir_type, solar_arc)
                        
        d_pls = [util.normalize(pl) for pl in d_pls]
        d_cusps = [util.normalize(c) for c in d_cusps]
        #print  tr2_year, tr2_month, d_pls        
                
        if pl_to_pl:    
            # -----------D planets to N planets---------
            ds = calc_pd_n_aspects(d_pls, n_pls, 0, pls_in_h, Aspects, dir_orbis)    
            
            # return ds
            # d_direct[(d_pl_name, n_pl_name, a)] = [[(d_sign, d_pl_info), (n_sign, n_pl_info), (diff)]]
                    
                
            for d in ds:
                        
                for_hours = i # which month it is at
                if for_hours==12: # next year
                    for_hours=0
                where_asp = float(i)/12
                
                d_pl_name = d[0]
                n_pl_name = d[1]
                asp_type = d[2]
                    
                if ('D', d_pl_name, n_pl_name, asp_type) not in transits_in:
                    transits_in[('D', d_pl_name, n_pl_name, asp_type)] = [where_asp]
                else: # update
                    asps_coords = transits_in[('D', d_pl_name, n_pl_name, asp_type)]
                    asps_coords.append(where_asp)
                    transits_in[('D', d_pl_name, n_pl_name, asp_type)] = asps_coords        
            
        
        
        if pl_to_cu:
            # -------------D planets to N cusps-------
            ds = calc_pd_n_aspects(d_pls, cusps, 1, pls_in_h, Aspects, dir_orbis)
            # d_direct[(d_pl_name, ncuspid, a)] = [[(d_sign, d_pl_info), (n_sign), (diff)]]
            
            for d in ds:
                        
                for_hours = i # which month it is at
                if for_hours==12: # next year
                    for_hours=0
                where_asp = float(i)/12
                
                d_pl_name = d[0]
                ncuspid = d[1]
                asp_type = d[2]
                    
                if ncuspid in e_houses:    
                    
                    if ('D', d_pl_name, ncuspid, asp_type) not in transits_in:
                        transits_in[('D', d_pl_name, ncuspid, asp_type)] = [where_asp]
                    else: # update
                        asps_coords = transits_in[('D', d_pl_name, ncuspid, asp_type)]
                        asps_coords.append(where_asp)
                        transits_in[('D', d_pl_name, ncuspid, asp_type)] = asps_coords
    
    
            
        
        if cu_to_pl:    
            # -----------D cusps to N planets---------
            ds = calc_pd_n_aspects(d_cusps, n_pls, 2, pls_in_h, Aspects, dir_orbis)
            
            for d in ds:
                        
                for_hours = i # which month it is at
                if for_hours==12: # next year
                    for_hours=0
                where_asp = float(i)/12
                
                dcuspid = d[0]
                n_pl_name = d[1]
                asp_type = d[2]
                    
                if dcuspid in e_houses:    
                    
                    if ('D', dcuspid, n_pl_name, asp_type) not in transits_in:
                        transits_in[('D', dcuspid, n_pl_name, asp_type)] = [where_asp]
                    else: # update
                        asps_coords = transits_in[('D', dcuspid, n_pl_name, asp_type)]
                        asps_coords.append(where_asp)
                        transits_in[('D', dcuspid, n_pl_name, asp_type)] = asps_coords
                        
                        
        if ingressions:    
            # -----------D cusps changes sign (both directions)---------
            # ds = calc_pd_n_aspects(d_cusps, n_pls, 2, pls_in_h, Aspects)
            
            ds = {}
            # take every natal cusp, and calc when it changes signs
            for d in d_cusps:
                cuspid_index = d_cusps.index(d)
                
                if d%30 < 1: # forward ingression
                    d_cuspid = cuspid_index +1
                    entered_sign = n_signs[util.normalize(d)//30+1]
                    
                    ds[(d_cuspid, entered_sign, 0)] = [[(d_cuspid), (entered_sign), (d%30)]]
                
                natal_cusp = cusps[cuspid_index]
                if (natal_cusp - (d-natal_cusp))%30 < 1:  #backward ingression
                    d_cuspid = cuspid_index +1
                    entered_sign_index = util.normalize((natal_cusp - (d-natal_cusp)))//30
                    if entered_sign_index==0:
                        entered_sign_index=12
                    entered_sign = n_signs[entered_sign_index]
                    
                    ds[(d_cuspid, entered_sign, 0)] = [[(d_cuspid), (entered_sign), (d%30)]]
            
            for d in ds:
                        
                for_hours = i # which month it is at
                if for_hours==12: # next year
                    for_hours=0
                where_asp = float(i)/12
                
                dcuspid = d[0]
                n_sign = d[1]
                asp_type = d[2]
                    
                if dcuspid in e_houses:    
                    
                    if ('D', dcuspid, n_sign, asp_type) not in transits_in:
                        transits_in[('D', dcuspid, n_sign, asp_type)] = [where_asp]
                    else: # update
                        asps_coords = transits_in[('D', dcuspid, n_sign, asp_type)]
                        asps_coords.append(where_asp)
                        transits_in[('D', dcuspid, n_sign, asp_type)] = asps_coords
                        
                    
                        
    
    # ----------------- plot all directions ------------------------
    the_planets = ['Sun','Moon','Mercury','Venus','Mars','Jupiter','Saturn','Uranus','Neptune','Pluto', 'mean Node'] + e_houses
    
    for tr_pl_name in the_planets: # for planets and selected houses
        #tr_pl_name = swe_get_planet_name(t)
        
        tr_pl_asps = {}
        tr_pl_asps_plotted = {}
        single_asp_type = [] # (tr_pl,n_pl,asp_type)
        
        for k in transits_in:
            
            if tr_pl_name==k[1]:
                tr_pl_asps[k] = transits_in[k]
                
                tr_pl_asps_plotted[(k[1],k[2],k[3])] = 0
                
                if (k[1],k[2],k[3]) not in single_asp_type:
                    single_asp_type.append((k[1],k[2],k[3]))
                
            
        for k in single_asp_type:
            
            if tr_pl_asps_plotted[k]==0: # that aspect has not been plotted yet
                
                
                tr_pl_D = ('D',k[0],k[1],k[2])
                tr_pl_R = ('R',k[0],k[1],k[2])

                
                if tr_pl_D in tr_pl_asps and tr_pl_R in tr_pl_asps: # then we need to merge those and plot on one track
                    
                    asp_coords = tr_pl_asps[tr_pl_D] + tr_pl_asps[tr_pl_R]
                
                    # flag them as used    
                    tr_pl_asps_plotted[k]==1
                else:
                    try:
                        asp_coords = tr_pl_asps[tr_pl_D]
                    except KeyError:
                        asp_coords = tr_pl_asps[tr_pl_R]
                
                from_when = min(asp_coords)
                till_when = max(asp_coords)
                
                    
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
                            if from_when > max_coord+1 or till_when < min_coord-1:
                                pass
                            else:
                                no_overlap_on_the_track = False
                            
                        if no_overlap_on_the_track:
                            ok_track = f
                            coords.append(from_when)
                            coords.append(till_when)
                            track_occupancy[ok_track] = coords
                            break
                        
                    
                    
                # -------------- plot ---------------
                
                x = [ok_track, ok_track]
                y = [from_when, till_when]
                                
                
                text_position = from_when+(till_when-from_when)/2 # in the middle of an aspect line
                if text_position>max_years:
                    text_position = max_years-otstup
                elif text_position<min_years:
                    text_position = min_years+otstup
                
                asp_type = k[2]
                
                if tr_pl_name not in uni_pls: # then it is a cusp
                    
                    tr_pl_name_cuspid = osi[tr_pl_name]
                    change = False
                    if asp_type==180:
                        asp_type=0
                        change = True
                    if asp_type==60:
                        asp_type=120
                        change = True
                    if change:
                        plt.text(ok_track, text_position-otstup, tr_pl_name_cuspid,  size=12)
                    else:
                        plt.text(ok_track, text_position-otstup, tr_pl_name,  size=12)
                else:
                    plt.text(ok_track, text_position-otstup, uni_pls[tr_pl_name],  size=12)
                
                
                
                
                
                
                
                n_pl_name = k[1]
                
               
                if n_pl_name not in uni_pls: # then it is a cuspid
                    try:
                        plt.text(ok_track, text_position+otstup, uni_signs[n_pl_name],  size=12) # for ingressions?
                    except KeyError: # for cusps
                        n_pl_name_cuspid = osi[n_pl_name]
                        change = False
                        if asp_type==180:
                            asp_type=0
                            change = True
                        if asp_type==60:
                            asp_type=120
                            change = True
                        if change:
                            plt.text(ok_track, text_position+otstup, n_pl_name_cuspid,  size=12)
                        else:
                            plt.text(ok_track, text_position+otstup, n_pl_name,  size=12)
                else:
                    plt.text(ok_track, text_position+otstup, uni_pls[n_pl_name],  size=12)
                
                
                
                if asp_type in [30,60,120]:
                    col = 'g'
                elif asp_type in [90,180,150]:
                    col = 'r'
                elif asp_type in [0]:
                    col = 'b'
                else:
                    col = 'k'
                
                plt.text(ok_track, text_position, uni_asps[asp_type],  size=12)
                
                
                
                
                
                
                lwidth = 1
                if (tr_pl_name, n_pl_name) in natal_asps or (n_pl_name,tr_pl_name) in natal_asps:
                    lwidth = 2
                if (tr_pl_name, n_pl_name) in natal_asps or (n_pl_name,tr_pl_name) in natal_asps:
                    lwidth = 2    
                                
                plt.plot(x, y, color=col, linewidth=lwidth, alpha=0.9)
                
                


    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    
    fig.set_size_inches(8, 12)
            
    fig.savefig(options.output_file+'_'+str(tr_origyear)+'_'+months[tr_origmonth][:3]+str(tr_origday)+'_'+str(int(hour))+'_'+str(int(minute))+'_'+str(int(second))+'.png', dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close(fig)
    

    
    
    
        

if __name__ == "__main__":
    main()        
