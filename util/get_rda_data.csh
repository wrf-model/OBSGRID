#!/bin/csh

### This is ONLY for the NCAR super computer. 
### This script will get all the data needed for OBSGRID and place it in one 
### big file called 'rda_obs'. After this step, run the get_rda_data.f Fortran 
### code to re-order the data into 3hourly input files
 

set start_YY = 2010
set start_MM = 3
set start_DD = 6
set start_HH = 0
set end_YY = 2010
set end_MM = 3
set end_DD = 7
set end_HH = 0

set YY=`printf %04d $end_YY`
set MM=`printf %02d $end_MM`
set DD=`printf %02d $end_DD`
set HH=`printf %02d $end_HH`
set end_date = ${YY}${MM}${DD}${HH}

set YY=`printf %04d $start_YY`
set MM=`printf %02d $start_MM`
set DD=`printf %02d $start_DD`
set HH=`printf %02d $start_HH`
set get_date = ${YY}${MM}${DD}${HH}


while ( $get_date <= $end_date ) 

  echo "Now getting data for " ${get_date}
  cat /glade/p/rda/data/ds351.0/little_r/OBS:${get_date} >> rda_obs
  cat /glade/p/rda/data/ds461.0/little_r/SURFACE_OBS:${get_date} >> rda_obs

  @ start_HH = $start_HH + 6
  if ( ${start_HH} == 24 ) then
    set start_HH = 0
    @ start_DD = $start_DD + 1
    set DD=`printf %02d $start_DD`
  endif
  set HH=`printf %02d $start_HH`
  set get_date = ${YY}${MM}${DD}${HH}

end
