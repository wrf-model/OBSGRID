MODULE ob_access

USE observation

CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE query_ob ( obs , date , time , &
request_variable , request_level , request_qc_max , request_p_diff , &
value , qc )

   USE date_pack

   IMPLICIT NONE

   TYPE (report)           :: obs
   INTEGER                 :: date , time
   CHARACTER ( LEN =   8 ) :: request_variable
   REAL                    :: request_level
   INTEGER                 :: request_qc_max   , &
                              request_p_diff
   REAL                    :: value
   INTEGER                 :: qc

   TYPE ( measurement ) , POINTER :: next_one

   LOGICAL                 :: still_more

   REAL                    :: r
   LOGICAL                 :: not_missing
   INTEGER                 :: ds , time_error
   CHARACTER (LEN=19)      :: obs_date , analysis_date , analysis_p1_date , analysis_m1_date
   LOGICAL                 :: close

   INCLUDE 'error.inc'
   INCLUDE 'missing.inc'
   INTERFACE
      INCLUDE 'error.int'
   END INTERFACE

   not_missing ( r )  = ( ABS ( r - missing_r ) .GT. 1. )

   !  Compute the difference in times bewteen the requested
   !  data and time, and the observation's data and time.

   obs_date( 1: 5) = obs%valid_time%date_char( 1: 4) // '-'
   obs_date( 6: 8) = obs%valid_time%date_char( 5: 6) // '-'
   obs_date( 9:11) = obs%valid_time%date_char( 7: 8) // '_'
   obs_date(12:14) = obs%valid_time%date_char( 9:10) // ':'
   obs_date(15:17) = obs%valid_time%date_char(11:12) // ':'
   obs_date(18:19) = obs%valid_time%date_char(13:14)

   WRITE ( analysIs_date , '(I4.4,"-",I2.2,"-",I2.2,"_",I2.2,":",I2.2,":",I2.2)' ) &
   date / 10000 ,  &
   ( date - (date / 10000 ) * 10000 ) / 100 , &
   date - ( date / 100 ) * 100 , &
   time / 10000 , &
   ( time - ( time / 10000 ) * 10000 ) / 100, &
   time - ( time / 100 ) * 100

   IF ( obs%info%platform(36:39) .EQ. '    ' ) THEN
      IF ( ( obs%info%platform( 1:11) .EQ. 'FM-97 AIREP'     ) .OR. &
           ( obs%info%platform( 1:15) .EQ. 'FM-36 TEMP SHIP' ) .OR. &
           ( obs%info%platform( 1:10) .EQ. 'FM-35 TEMP'      ) .OR. &
           ( obs%info%platform( 1:11) .EQ. 'FM-88 SATOB'     ) ) THEN
         ds = 3600
      ELSE
         ds = 1800
      END IF
   ELSE
      READ (obs%info%platform(36:39),IOSTAT=time_error,FMT='(I4)') ds
      IF ( time_error .NE. 0 ) THEN
         IF ( ( obs%info%platform( 1:11) .EQ. 'FM-97 AIREP'     ) .OR. &
              ( obs%info%platform( 1:15) .EQ. 'FM-36 TEMP SHIP' ) .OR. &
              ( obs%info%platform( 1:10) .EQ. 'FM-35 TEMP'      ) .OR. &
              ( obs%info%platform( 1:11) .EQ. 'FM-88 SATOB'     ) ) THEN
            ds = 3600
         ELSE
            ds = 1800
         END IF
      END IF
   END IF

   CALL geth_newdate ( analysis_p1_date , analysis_date ,    ds ) 
   CALL geth_newdate ( analysis_m1_date , analysis_date , -1*ds )

   IF ( ( obs_date .LE. analysis_p1_date ) .AND. &
        ( obs_date .GE. analysis_m1_date ) ) THEN
      close = .TRUE.
   ELSE
      close = .FALSE.
   END IF

   !  Is the the correct time, or not.

   right_time : IF ( close ) THEN

      !  Data is different for the single values (surface, aircraft, etc) and the 
      !  vertically stacked data (soundings, satellite winds, etc), so we need
      !  to differentiate with respect to those.
   
      slp_vs_others : SELECT CASE ( request_variable ) 
   
         CASE ( 'PMSL    ' ) slp_vs_others
            IF ( ( not_missing ( obs%ground%slp%data ) ) .AND. &
                 ( obs%ground%slp%qc .LT. request_qc_max ) ) THEN
               value = obs%ground%slp%data
               qc    = obs%ground%slp%qc
            END IF 
   
         CASE DEFAULT slp_vs_others
   
            !  All of the data that is not in the "ground" section is stored in
            !  a linked list.  We can not assume that the first level is the
            !  surface, nor can we assume that the first level exists. Initialize 
            !  the linked list trace pointer.  We will follow the linked list as 
            !  long as there is still_more to do.
   
            next_one => obs%surface
            still_more = ASSOCIATED ( next_one )
    
            not_slp_which_level : SELECT CASE ( NINT ( request_level ) ) 
   
               CASE ( 1001 ) not_slp_which_level
   
                  !  Does the level exist?  Loop through the linked list of the data
                  !  to find a correct pressure.  Surface data will be handled separately,
                  !  as it is defined when the elevation is equal to the geopotential
                  !  height at a level.
   
                  search_for_sfc : DO WHILE ( still_more ) 
                     found_sfc : IF ( (       eps_equal ( next_one%meas%height%data , obs%info%elevation , 1. ) ) .AND. &
                                      ( .NOT. eps_equal ( next_one%meas%height%data , missing_r          , 1. ) ) .AND. &
                                      ( .NOT. eps_equal ( obs%info%elevation        , missing_r          , 1. ) ) ) THEN
   
                        !  Does the data exist on this level?  We found the right
                        !  pressure level, now check the existence of the 
                        !  variable, and if the qc info is correct.  Check has 
                        !  to be done for each variable.
   
                        which_sfc_var : SELECT CASE ( request_variable ) 
   
                           CASE ( 'UU      ' ) which_sfc_var
                              IF ( ( not_missing ( next_one%meas%u%data ) ) .AND. &
                                                 ( next_one%meas%u%qc .LT. request_qc_max ) ) THEN
                                 value = next_one%meas%u%data
                                 qc    = next_one%meas%u%qc
                              END IF
   
                           CASE ( 'VV      ' ) which_sfc_var
                              IF ( ( not_missing ( next_one%meas%v%data ) ) .AND. &
                                                 ( next_one%meas%v%qc .LT. request_qc_max ) ) THEN
                                 value = next_one%meas%v%data
                                 qc    = next_one%meas%v%qc
                              END IF
   
                           CASE ( 'TT      ' ) which_sfc_var
                              IF ( ( not_missing ( next_one%meas%temperature%data ) ) .AND. &
                                                 ( next_one%meas%temperature%qc .LT. request_qc_max ) ) THEN
                                 value = next_one%meas%temperature%data
                                 qc    = next_one%meas%temperature%qc
                              END IF
   
                           CASE ( 'RH      ' ) which_sfc_var
                              IF ( ( not_missing ( next_one%meas%rh%data ) ) .AND. &
                                                 ( next_one%meas%rh%qc .LT. request_qc_max ) ) THEN
                                 value = next_one%meas%rh%data
                                 qc    = next_one%meas%rh%qc
                              END IF
               
                        END SELECT which_sfc_var      
   
                        !  We found the surface data, and we either have the data or not.  
                        !  Either way we need to exit this loop which will try to keep
                        !  finding the surface layer.
   
                        EXIT search_for_sfc
   
                     ELSE IF ( ( next_one%meas%height%data - obs%info%elevation .GT. 10. ) .AND. &
                               ( .NOT. eps_equal ( next_one%meas%height%data , missing_r , 1. ) ) .AND. &
                               ( .NOT. eps_equal ( obs%info%elevation        , missing_r , 1. ) ) ) THEN found_sfc
   
                        !  Once the geopotential height is higher than the elevation, since the data 
                        !  is vertically ordered (surface to top of atmosphere), we will not find
                        !  the surface above here.
   
                        EXIT search_for_sfc
   
                     ELSE found_sfc
   
                        !  We did not find the surface level, so go back and keep looking. 
   
                        next_one => next_one%next
                        still_more = ASSOCIATED ( next_one )
   
                     END IF found_sfc
   
                  END DO search_for_sfc
   
               CASE DEFAULT not_slp_which_level
   
                  !  Initialize the linked list trace pointer.  We will
                  !  follow the linked list as long as there is 
                  !  still_more to do.
   
                  next_one => obs%surface
                  still_more = ASSOCIATED ( next_one )
   
                  !  Does the level exist?  Loop through the linked list of the data
                  !  to find a correct pressure.
   
                  search_for_level : DO WHILE ( still_more ) 
                     found_level : IF ( ABS ( next_one%meas%pressure%data - request_level * 100 ) .LT. request_p_diff ) THEN
   
                        !  Does the data exist on this level?  We found the right
                        !  pressure level, now check the existence of the 
                        !  variable, and if the qc info is correct.  Check has 
                        !  to be done for each variable.
   
                        which_var : SELECT CASE ( request_variable ) 
   
                           CASE ( 'UU      ' ) which_var
                              IF ( ( not_missing ( next_one%meas%u%data ) ) .AND. &
                                                 ( next_one%meas%u%qc .LT. request_qc_max ) ) THEN
                                 value = next_one%meas%u%data
                                 qc    = next_one%meas%u%qc
                              END IF
   
                           CASE ( 'VV      ' ) which_var
                              IF ( ( not_missing ( next_one%meas%v%data ) ) .AND. &
                                                 ( next_one%meas%v%qc .LT. request_qc_max ) ) THEN
                                 value = next_one%meas%v%data
                                 qc    = next_one%meas%v%qc
                              END IF
   
                           CASE ( 'TT      ' ) which_var
                              IF ( ( not_missing ( next_one%meas%temperature%data ) ) .AND. &
                                                 ( next_one%meas%temperature%qc .LT. request_qc_max ) ) THEN
                                 value = next_one%meas%temperature%data
                                 qc    = next_one%meas%temperature%qc
                              END IF
   
                           CASE ( 'RH      ' ) which_var
                              IF ( ( not_missing ( next_one%meas%rh%data ) ) .AND. &
                                                 ( next_one%meas%rh%qc .LT. request_qc_max ) ) THEN
                                 value = next_one%meas%rh%data
                                 qc    = next_one%meas%rh%qc
                              END IF
   
                        END SELECT which_var
   
                        !  Inside this IF, we found the level.  Either we have the data or
                        !  not, but either way we must exit the level searching loop (this
                        !  effectively exits this routine).
   
                        EXIT search_for_level
   
                     ELSE found_level
     
                        !  We didn't find the correct level yet, so go to the next location
                        !  in the linked list.
     
                        next_one => next_one%next
                        still_more = ASSOCIATED ( next_one ) 
   
                     END IF found_level
                     
                  END DO search_for_level
   
            END SELECT not_slp_which_level
               
      END SELECT slp_vs_others

   ELSE right_time
      !obs%info%discard = .TRUE.                          !!! 2012-12-17 - cB - to allow data - not qc'ed to enter OBS_DOMAIN101 file
      obs%surface%meas%pressure%qc =  no_qc_possible
      obs%surface%meas%height%qc =  no_qc_possible
      obs%surface%meas%temperature%qc =  no_qc_possible
      obs%surface%meas%u%qc =  no_qc_possible
      obs%surface%meas%v%qc =  no_qc_possible
      obs%surface%meas%rh%qc =  no_qc_possible
      error_number = 00351001
      error_message(1:31) = 'query_ob                       '
      error_message(32:)  = ' Wrong time for observation ' // &
      TRIM ( obs%location%id ) // ', at time = ' // obs%valid_time%date_char(1:12) // &
      ' (ccyymmddhhmm).'  
      fatal = .false.
      listing = .false.
! foo
!    CALL error_handler ( error_number , error_message ,  &
!    fatal , listing )
!    print*,obs%surface%meas%temperature%data, obs%surface%meas%temperature%qc
   END IF right_time
   
END SUBROUTINE query_ob

!------------------------------------------------------------------------------

SUBROUTINE store_ob ( obs , &
request_variable , request_level , request_p_diff , &
value , qc )

   IMPLICIT NONE

   TYPE (report)           :: obs
   CHARACTER ( LEN =   8 ) :: request_variable
   REAL                    :: request_level
   INTEGER                 :: request_p_diff
   REAL                    :: value
   INTEGER                 :: qc

   TYPE ( measurement ) , POINTER :: next_one

   LOGICAL                 :: still_more

   CHARACTER ( LEN = 4 )   :: plevel_char

   INCLUDE 'error.inc'
   INTERFACE
      INCLUDE 'error.int'
   END INTERFACE

   !  Data is different for the single values (surface, aircraft, etc) and the 
   !  vertically stacked data (soundings, satellite winds, etc), so we need
   !  to differentiate with respect to those.

   slp_vs_others : SELECT CASE ( request_variable ) 

      CASE ( 'PMSL    ' ) slp_vs_others
         obs%ground%slp%data = value
         obs%ground%slp%qc   = qc

      CASE DEFAULT slp_vs_others

         !  All of the data that is not in the "ground" section is stored in
         !  a linked list.  We can not assume that the first level is the
         !  surface. Initialize the linked list trace pointer.  We will 
         !  follow the linked list as long as there is still_more to do.
         !  If we run out of data, that is a FATAL error, since we should
         !  be able to find where this observation came from.

         next_one => obs%surface
         still_more = ASSOCIATED ( next_one )
 
         not_slp_which_level : SELECT CASE ( NINT ( request_level ) ) 

            CASE ( 1001 ) not_slp_which_level

               !  Loop through the linked list of the data
               !  to find the correct pressure.  Surface data will be handled separately,
               !  as it is defined when the elevation is equal to the geopotential
               !  height at a level.

               search_for_sfc : DO WHILE ( still_more ) 
                  found_sfc : IF ( (       eps_equal ( next_one%meas%height%data , obs%info%elevation , 1. ) ) .AND. &
                                   ( .NOT. eps_equal ( next_one%meas%height%data , missing_r          , 1. ) ) .AND. &
                                   ( .NOT. eps_equal ( obs%info%elevation        , missing_r          , 1. ) ) ) THEN

                     !  Does the data exist on this level?  We found the right
                     !  pressure level, now check the existence of the 
                     !  variable, and if the qc info is correct.  Check has 
                     !  to be done for each variable.

                     which_sfc_var : SELECT CASE ( request_variable ) 

                        CASE ( 'UU      ' ) which_sfc_var
                           next_one%meas%u%data = value
                           next_one%meas%u%qc   = qc

                        CASE ( 'VV      ' ) which_sfc_var
                           next_one%meas%v%data = value
                           next_one%meas%v%qc   = qc

                        CASE ( 'TT      ' ) which_sfc_var
                           next_one%meas%temperature%data = value
                           next_one%meas%temperature%qc   = qc

                        CASE ( 'RH      ' ) which_sfc_var
                           next_one%meas%rh%data = value
                           next_one%meas%rh%qc   =qc
            
                     END SELECT which_sfc_var      

                     !  We found the surface data layer, and stored the value and QC 
                     !  information. We need to exit this loop which will try to keep
                     !  finding the surface layer.

                     EXIT search_for_sfc

                  ELSE IF ( ( next_one%meas%height%data - obs%info%elevation .GT. 10. ) .AND. &
                            ( .NOT. eps_equal ( next_one%meas%height%data , missing_r , 1. ) ) .AND. &
                            ( .NOT. eps_equal ( obs%info%elevation        , missing_r , 1. ) ) ) THEN found_sfc

                     !  Once the geopotential height is higher than the elevation, since the data 
                     !  is vertically ordered (surface to top of atmosphere), we will not find
                     !  the surface above here.

                     EXIT search_for_sfc

                  ELSE found_sfc

                     !  We did not find the surface level, so go back and keep looking. 

                     next_one => next_one%next
                     still_more = ASSOCIATED ( next_one )

                     !  Here is where we panic if we ran out of data to search.
                     !  We must find the correct level, since it used to exist.

                     bad_news_1 : IF ( .NOT. still_more ) THEN
                        error_number = 00352001
                        error_message(1:31) = 'store_ob                       '
                        error_message(32:)  = ' Could not find surface level for ' // &
                        TRIM ( obs%location%id ) // ', for variable ' // request_variable // '.'
                        fatal = .TRUE.
                        listing = .false.
                        CALL error_handler ( error_number , error_message ,  &
                        fatal , listing )
                     END IF bad_news_1

                  END IF found_sfc

               END DO search_for_sfc

            CASE DEFAULT not_slp_which_level

               !  Initialize the linked list trace pointer.  We will follow the 
               !  linked list as long as there is still_more to do.  We should find
               !  the correct level prior to running out of data.

               next_one => obs%surface
               still_more = ASSOCIATED ( next_one )

               !  Loop through the linked list of the data to find a correct pressure.

               search_for_level : DO WHILE ( still_more ) 
                  found_level : IF ( ABS ( next_one%meas%pressure%data - request_level * 100 ) .LT. request_p_diff ) THEN

                     !  We found the right pressure level, so just shove the data in.

                     which_var : SELECT CASE ( request_variable ) 

                        CASE ( 'UU      ' ) which_var
                           next_one%meas%u%data = value
                           next_one%meas%u%qc   = qc

                        CASE ( 'VV      ' ) which_var
                           next_one%meas%v%data = value
                           next_one%meas%v%qc   = qc

                        CASE ( 'TT      ' ) which_var
                           next_one%meas%temperature%data = value
                           next_one%meas%temperature%qc   = qc

                        CASE ( 'RH      ' ) which_var
                           next_one%meas%rh%data = value
                           next_one%meas%rh%qc = qc

                     END SELECT which_var

                     !  Inside this IF, we found the level.  We have stored the data, 
                     !  so we must exit the level searching loop (this effectively exits 
                     !  this routine).

                     EXIT search_for_level

                  ELSE found_level
  
                     !  We didn't find the correct level yet, so go to the next location
                     !  in the linked list.
  
                     next_one => next_one%next
                     still_more = ASSOCIATED ( next_one ) 

                     !  Here is where we panic if we ran out of data to search.
                     !  We must find the correct level, since it used to exist.

                     bad_news_2 : IF ( .NOT. still_more ) THEN
                        WRITE ( plevel_char , FMT='(i4)' ) NINT ( request_level ) 
                        error_number = 00352002
                        error_message(1:31) = 'store_ob                       '
                        error_message(32:)  = ' Could not find level ' // plevel_char // &
                        ' hPa for ' // &
                        TRIM ( obs%location%id ) // ', for variable ' // request_variable // '.'
                        fatal = .TRUE.
                        listing = .false.
                        CALL error_handler ( error_number , error_message ,  &
                        fatal , listing )
                     END IF bad_news_2

                  END IF found_level
                  
               END DO search_for_level

         END SELECT not_slp_which_level
            
   END SELECT slp_vs_others

END SUBROUTINE store_ob

!------------------------------------------------------------------------------
   
END MODULE ob_access
