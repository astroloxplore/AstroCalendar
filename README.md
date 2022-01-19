# AstroCalendar

## AstroCalendar of transits (months)

**transit_calendar.py** is a script that generates individual astrological calendar of transits for every given month.

Example of input file format is in test_transit_calendar_2021.zbs

Running the script (default parameters):
**python3 transit_calendar.py -i test_transit_calendar_2021.zbs**

Running the script (with additional parameters):
**python3 transit_calendar.py -i test_transit_calendar_2021.zbs -t 31 -e 3 -c 0 -b 1.2 -a 3**
* -t -> time period in days (default is 31 days)
* -e 3 -> expand transit orb to 3 degrees (default is 1)
* -c 0 -> transits to all cusps (default - only transits to Asc/MC)
* -b 1.2 -> expand orb of natal aspects by 20%
* -a 3 -> use a different set of aspects (default is [0,60,90,120,180])


## AstroCalendar of transits (years)

**transit_calendar_years.py** is a script that generates individual astrological calendar of transits for a period of time in years.

Running the script (default parameters):
**python3 transit_calendar_years.py -i test_transit_calendar_years.zbs**

Running the script (with additional parameters):
**python3 transit_calendar_years.py -i test_transit_calendar_years.zbs -t 5 -e 3 -c 0**
* -t -> time period, 5 years (default is 10 years)
* -e 3 -> expand transit orb to 3 degrees
* -c 0 -> transits to all cusps

## DEPENDENCIES
**util.py** (provided) from Morinus
