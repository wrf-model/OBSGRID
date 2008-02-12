PROGRAM gts_cleaner

   CHARACTER(LEN=600) :: first_line
   CHARACTER(LEN=200) :: second_line
   CHARACTER(LEN=21)  :: third_line

   INTEGER , PARAMETER :: input = 9 , output = 19 , temp = 29
   INTEGER :: ok_1 , ok_2 , ok_3

   LOGICAL :: trouble

   OPEN ( UNIT= input,FILE='gts_data'      ,STATUS='old'    ,FORM='FORMATTED')
   OPEN ( UNIT=output,FILE='gts_data_clean',STATUS='unknown',FORM='FORMATTED')

   read_all_stations : DO

      !  Open a temporary file for each station.

      OPEN ( UNIT=  temp,FILE='gts_data_temp' ,STATUS='unknown',FORM='FORMATTED')

      !  Initialize each station's data to OK.

      trouble = .FALSE.

      !  Read first line.

      READ(UNIT=input,IOSTAT=ok_1,FMT='(A)') first_line

      !  If this was not a successful read, this is the end of the file, probably.

      IF(ok_1 .NE. 0 ) THEN
         PRINT '(A)','Found end of data for all stations.'
         EXIT read_all_stations
      END IF

      !  If this was successful read, are there any '*' in the string?

      ok_1 = INDEX ( first_line , '*' )
      IF ( ok_1 .EQ. 0 ) THEN
         WRITE(UNIT=temp,FMT='(A)') first_line
      ELSE
         trouble = .TRUE.
         PRINT '(A,A)','Found trouble in first line with ',first_line(41:80)
      END IF

      read_this_station : DO 

         !  Read second line, well, multiple lines in the second format.
   
         READ(UNIT=input,IOSTAT=ok_2,FMT='(A)') second_line
   
         !  If this was not a successful read, this is the end of the station
         !  info.  Backup and read it with the last format.
   
         IF(ok_2 .NE. 0 ) THEN
            PRINT '(A)','Found end of this station in an unusual way.'
            BACKSPACE(UNIT=input)
            EXIT read_this_station
         END IF
   
         !  If this was successful read, are there any '*' in the string?
   
         ok_2 = INDEX ( second_line , '*' )
         IF ( ok_2 .EQ. 0 ) THEN
            WRITE(UNIT=temp,FMT='(A)') second_line
         ELSE
            trouble = .TRUE.
            PRINT '(A,A)','Found trouble in middle line with ',first_line(41:80)
         END IF

         !  Was this the flag value of "no more levels"?

         IF ( second_line(1:2) .EQ. '-7' ) THEN
            EXIT read_this_station
         END IF

      END DO read_this_station
      
      !  Read third line.

      READ(UNIT=input,IOSTAT=ok_3,FMT='(A)') third_line

      !  If this was not a successful read, this is just plain confusing, so quit.

      IF(ok_3 .NE. 0 ) THEN
         PRINT '(A)','Found an error in the third line, unexpected.'
         EXIT read_all_stations
      END IF

      !  If this was successful read, are there any "*" in the string?

      ok_3 = INDEX ( third_line , '*' )
      IF ( ok_3 .EQ. 0 ) THEN
         WRITE(UNIT=temp,FMT='(A)') third_line
      ELSE
         trouble = .TRUE.
         PRINT '(A,A)','Found trouble in last line with ',first_line(41:80)
      END IF

      !  If there was no trouble with this ob, then we can put it in the
      !  good output file, otherwise, we discard it.

      IF ( .NOT. trouble ) THEN
   
         REWIND(UNIT=temp)

         READ(UNIT=temp,FMT='(A)') first_line
         WRITE(UNIT=output,FMT='(A)') first_line
         PRINT '(A,A)','Processing ',first_line(41:80)
         
         again : DO 
            READ(UNIT=temp,FMT='(A)') second_line
            WRITE(UNIT=output,FMT='(A)') second_line
            IF ( second_line(1:2) .EQ. '-7' ) THEN
               EXIT again
            END IF
         END DO again

         READ(UNIT=temp,FMT='(A)') third_line
         WRITE(UNIT=output,FMT='(A)') third_line

         !  Close the temporary file

         CLOSE (UNIT=  temp)

      ELSE

         !  Close the temporary file, throw this suspect data away, attack
         !  the next sequence of obs.

         CLOSE (UNIT=  temp)
         PRINT '(A)','Found errors in ob, discarding it.'

      END IF

   END DO read_all_stations

   CLOSE (UNIT= input)
   CLOSE (UNIT=output)
   CLOSE (UNIT=  temp)

END PROGRAM gts_cleaner
