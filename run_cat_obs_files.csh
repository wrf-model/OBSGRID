#!/bin/csh
#set echo

## Run this script to create the single OBS_DOMAIN101 file needed for WRF Obs nudging

unalias ls

foreach fil (` ls -1 OBS_DOMAIN1?? `)
  echo "appending " $fil
  mv $fil ${fil}_sav
  cat ${fil}_sav >> OBS_DOMAIN101
end
