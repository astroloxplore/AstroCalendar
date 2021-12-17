# AstroCalendar of transits

transit_calendar.py is a script that generates individual astrological calendar of transits.
It works as follows:
- calculates natal positions of planets, cuspids for a given date, time, and place of birth (provided in the input file)
- calculates transits for given period of time in months (also in the input file)
- plots each month's transits

Input file format (ZBS format / ZET astroprocessor, first line - natal, other - each month):
- t; 01.01.1992; 12:20:00; +3; Vilnius, Lithuania; 54n41; 25e19; -; ;
- t; 1.1.2022; 00:00:00; +3; Vilnius, Lithuania; 54n41; 25e19; -; ;
- t; 1.2.2022; 00:00:00; +3; Vilnius, Lithuania; 54n41; 25e19; -; ;


Running the script (default parameters):
- python3 transit_calendar.py -i 'test_transit_calendar_2021.zbs'

Running the script (with additional parameters):
- python3 transit_calendar.py -i 'test_transit_calendar_2021.zbs' -t 31 -e 3 -c 0 -b 1.2 -a 3
- -t -> time period, one month is recommended (default)
- -e 3 -> expand transit orb to 3 degrees (default is 1)
- -c 0 -> transits to all cusps (default - only transits to Asc/MC)
- -b 1.2 -> expand orb of natal aspects by 20%
- -a 3 -> use a different set of aspects (default is [0,60,90,120,180])

DEPENDENCIES:
util.py (provided) from Morinus
