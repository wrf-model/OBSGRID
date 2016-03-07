MODULE ob_access

USE observation
!BPR BEGIN
USE get_fg_at_point
!BPR END

CONTAINS

!------------------------------------------------------------------------------
SUBROUTINE query_ob ( obs , date , time , &
request_variable , request_level , request_qc_max , request_p_diff , &
!BPR BEGIN
!value , qc )
value , qc , fg_3d_t , fg_3d_h )
!BPR END

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
!BPR BEGIN
   REAL , INTENT(IN) , OPTIONAL , DIMENSION(:,:,:) :: fg_3d_t, &
                                                      fg_3d_h
   LOGICAL                 :: mark_as_false_extend
!BPR END

   TYPE ( measurement ) , POINTER :: next_one

   LOGICAL                 :: still_more

   REAL                    :: r
   LOGICAL                 :: not_missing
   INTEGER                 :: ds , time_error
   CHARACTER (LEN=19)      :: obs_date , analysis_date , analysis_p1_date , analysis_m1_date
   LOGICAL                 :: close

!BPR BEGIN
   REAL                    :: curpdiff, minpdiff
   INTEGER                 :: curlev, uselev, maxlev
   LOGICAL                 :: isin

!BPR END
!BPR BEGIN
   REAL                    :: rh_value, t_value, rh_qc, t_qc
   REAL                    :: t_layer_average, t_value_at_sea_level
   REAL                    :: psfc_value, psfc_qc
   REAL                    :: elevation_value
   INCLUDE 'constants.inc'
!BPR END

   INCLUDE 'error.inc'
   INCLUDE 'missing.inc'
   INTERFACE
      INCLUDE 'error.int'
   END INTERFACE

   not_missing ( r )  = ( ABS ( r - missing_r ) .GT. 1. )

!BPR BEGIN
   mark_as_false_extend = .false.
!BPR END

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

                           !BPR BEGIN
                           !Added so that we can perform QC on dewpoint to improve QC of RH
                           CASE ( 'DEWPOINT' ) which_sfc_var
                              IF ( ( ( not_missing ( next_one%meas%rh%data ) ) .AND. &
                                                   ( next_one%meas%rh%qc .LT. request_qc_max ) ) .AND. &
                                   ( ( not_missing ( next_one%meas%temperature%data ) ) .AND. &
                                                   ( next_one%meas%temperature%qc .LT. request_qc_max ) ) ) THEN
                                 rh_value = next_one%meas%rh%data
                                 rh_qc    = next_one%meas%rh%qc
                                 t_value = next_one%meas%temperature%data
                                 t_qc    = next_one%meas%temperature%qc
                                 !Calculate dewpoint using formula in
                                 !obs_sort_module.F90
                                 !Since log(0) is undefined, if RH is near zero then set it to
                                 !some small value
                                 IF( rh_value .LT. very_small_rh ) THEN
                                  value = 1. / ( 1./t_value - Rv_over_L * LOG ( very_small_rh / 100. ) )
                                 ELSE
                                  value = 1. / ( 1./t_value - Rv_over_L * LOG ( rh_value / 100. ) )
                                 END IF

                                 !Since purpose is to improve QC of RH, pull RH QC values
                                 qc = rh_qc
                              END IF

                           CASE ( 'PSFC    ' ) which_sfc_var
                              !Although ob header stores surface pressure in obs%ground%psfc%data
                              !Obsgrid seems to use the pressure in the surface level of the ob
                              IF ( ( not_missing ( next_one%meas%pressure%data ) ) .AND. &
                                   ( next_one%meas%pressure%qc .LT. request_qc_max ) ) THEN
                                 value = next_one%meas%pressure%data
                                 qc    = next_one%meas%pressure%qc
                              END IF 

                           CASE ( 'PMSLPSFC' ) which_sfc_var
                              !This case means that we want to obtain the sea level pressure 
                              !based on the observed surface pressure
                              !Although ob header stores surface pressure in obs%ground%psfc%data
                              !Obsgrid seems to use the pressure in the surface level of the ob
                              IF ( ( not_missing ( next_one%meas%pressure%data ) )   .AND. &
                                   ( next_one%meas%pressure%qc .LT. request_qc_max ) .AND. &
                                   ( .NOT. eps_equal ( obs%info%elevation , missing_r , 1. ) ) ) THEN
                                 psfc_value = next_one%meas%pressure%data
                                 psfc_qc    = next_one%meas%pressure%qc
                                 elevation_value = obs%info%elevation
                                 IF ( ( not_missing ( next_one%meas%temperature%data ) ) .AND. &
                                                    ( next_one%meas%temperature%qc .LT. request_qc_max ) ) THEN
                                  !If the ob has a surface temperature ob that is not bad
                                  !then use that to define temperature at level of ob and
                                  !standard atmosphere to get to sea level
                                  t_value = next_one%meas%temperature%data
                                  t_value_at_sea_level = t_value - dTdz_std_atm_lt_11000m * elevation_value
                                  t_layer_average = 0.5 * (t_value_at_sea_level + t_value )  
                                 ELSE
                                  !If the ob does NOT have a good surface temperature ob
                                  !Then use first guess field to determine temperature
                                  IF( PRESENT(fg_3d_t) .AND. PRESENT(fg_3d_h) ) THEN
                                   CALL get_fg_t_at_point(fg_3d_t,fg_3d_h,obs%location%latitude, &
                                    obs%location%longitude, 0.5*elevation_value, t_layer_average) 
                                   !If average temperature of layer is missing and the height of the observation is
                                   !below 11000m (the top of the level where the lapse rate is applicable) then use 
                                   !standard atmosphere to find temperature
                                   IF ( ( ABS ( t_layer_average - missing_r ) .LT. 1. ) .AND. &
                                        ( elevation_value .LE. 11000.0 ) ) THEN
                                    t_value_at_sea_level = t_at_sea_level_std_atm
                                    t_value = t_value_at_sea_level + dTdz_std_atm_lt_11000m * elevation_value 
                                    t_layer_average = 0.5 * (t_value_at_sea_level + t_value )  
                                   END IF

                                  ELSE
                                   !If somehow we do not have the first guess field then error out since 
                                   !this should not happen. 
                                   !standard atmosphere
                                   PRINT *,'ERROR: query_ob: First guess temperature and/or height are not available'
                                   STOP
                                   !In theory we could use standard atmosphere to get temperature
                                   t_value_at_sea_level = t_at_sea_level_std_atm
                                   t_value = t_value_at_sea_level + dTdz_std_atm_lt_11000m * elevation_value 
                                   t_layer_average = 0.5 * (t_value_at_sea_level + t_value )  
                                  END IF
                                 END IF
                                 IF ( not_missing ( t_layer_average ) ) THEN
                                  !Value is the sea level pressure calculated from the 
                                  !surface pressure, elevation, and temperatures
                                  value      = psfc_value / (exp( ( -elevation_value * g ) / (gasr * t_layer_average ) ) )
                                  qc         = psfc_qc
                                 END IF
                              END IF 
                           !BPR END
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

                  !BPR BEGIN
                  !Original code assumes that a given ob will be close enough to
                  !at most one of the pressure levels in the background field
                  !If one expands the tolerance (request_p_diff) enough then one
                  !must determine if current pressure level of analysis is the
                  !one closest to this ob
                  curlev=0
                  uselev=-1
                  minpdiff=9999999
                  search_for_levela : DO WHILE ( still_more )
                   curlev=curlev+1
                   curpdiff=ABS ( next_one%meas%pressure%data - request_level * 100 )
                   found_levela : IF ( curpdiff .LT. request_p_diff ) THEN
                     IF(curpdiff.LT.minpdiff) THEN
                      minpdiff=curpdiff
                      uselev=curlev
                     ENDIF
                   END IF found_levela
                   next_one => next_one%next
                   still_more = ASSOCIATED ( next_one )
                  END DO search_for_levela

                  next_one => obs%surface
                  still_more = ASSOCIATED ( next_one )
                  maxlev=curlev
                  !BPR END 
 
                  !  Does the level exist?  Loop through the linked list of the data
                  !  to find a correct pressure.
   
                  !BPR BEGIN
                  curlev=0
                  !BPR END
                  search_for_level : DO WHILE ( still_more ) 
                     !BPR BEGIN
                     curlev=curlev+1
                     curpdiff=ABS ( next_one%meas%pressure%data - request_level * 100 )
                     !found_level : IF ( ABS ( next_one%meas%pressure%data - request_level * 100 ) .LT. request_p_diff ) THEN
                     !found_level : IF ( curlev .EQ. uselev ) THEN
                     !If single-level ob (e.g., ACARS) and the pressure is
                     !within request_p_diff Pa of analysis level we are working on
                     !OR
                     !If multiple-level ob (e.g., sonde) and the pressure is
                     !very close to the analysis level (within 1 Pa) 
                     !AND
                     !It is not a surface ob (excluded by ensuring that the height
                     !of the ob is not identical to the station elevation,
                     !unless the height of the ob is missing)
                     !
                     !Single-level obs need larger tolerances in order to be
                     !included since they cannot be interpolated to analysis pressures
                     !Note that some single-level obs (FM-97 AIREP and FM-88 SATOB) will
                     !be extended to the closest analysis level if possible by
                     !extend_to_closest if use_p_tolerance_one_lev = .FALSE.  
                     !If this has been done they will have
                     !the original level and the new level and so they will be
                     !treated with the multiple-level obs

                     !Multiple-level obs cannot use the broader tolerance, because they 
                     !have been pre-interpolated to the analysis levels.  If we
                     !use the broader tolerance we may include both the
                     !non-interpolated levels close to the analysis level and
                     !the new interpolated value at the analysis level

                     !Surface obs are excluded here since they should be
                     !qc'd/analyzed at the surface level (which is dealt with in
                     !the not_slp_which_level switch under the case where
                     !request_level=1001).  Here we are under the DEFAULT case.
                     !In the original code, surface obs were basically excluded
                     !from this section of code because surface obs would only be
                     !included if their pressure was within 0Pa of analysis
                     !level.  With the pressure tolerance expanded for
                     !single-level obs, surface obs would be much more likely to
                     !be included in another level if surface obs were not
                     !specifically excluded.
                     found_level : IF (((( maxlev .EQ. 1 ).AND.( curlev .EQ. uselev )).OR.&
                                        (( maxlev .NE. 1 ).AND.( curpdiff .LT. 1 ))).AND.&
                                       ((.NOT. eps_equal ( next_one%meas%height%data , obs%info%elevation , 1. ) ) .OR. &
                                        (      eps_equal ( next_one%meas%height%data , missing_r          , 1. ) ))) THEN
                        !If ob qualifies because it is a single-level non-surface 
                        !ob close enough to desired level
                        IF (( maxlev .EQ. 1 ).AND.( curlev .EQ. uselev )) THEN
                         !Note that this ob should be marked with the qc flag
                         !usually used for single-level obs that are assigned
                         !the closest analysis pressure level 
                         !This will allow ob to pass following remove_unverified test on
                         !output in obs_sort_module.F90:
                         !IF ( is_sounding .AND. next%meas%pressure%qc .gt. 4 )  keep_data = .TRUE.
                         mark_as_false_extend = .true.
                        ENDIF
                        IF ( mark_as_false_extend ) THEN
                         !Check if pressure qc flag already contains extend_influence flag
                         isin = contains_2n ( next_one%meas%pressure%qc, extend_influence )
                         !If pressure qc flag does NOT already contain the
                         !extend_influence flag, then add it
                         IF ( .NOT. isin ) THEN
                          next_one%meas%pressure%qc = next_one%meas%pressure%qc + extend_influence
                         ENDIF
                        ENDIF

                     !BPR END
   
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

                           !BPR BEGIN
                           !Added so that we can perform QC on dewpoint to improve QC of RH
                           CASE ( 'DEWPOINT' ) which_var
                              IF ( ( ( not_missing ( next_one%meas%rh%data ) ) .AND. &
                                                   ( next_one%meas%rh%qc .LT. request_qc_max ) ) .AND. &
                                   ( ( not_missing ( next_one%meas%temperature%data ) ) .AND. &
                                                   ( next_one%meas%temperature%qc .LT. request_qc_max ) ) ) THEN
                                 rh_value = next_one%meas%rh%data
                                 rh_qc    = next_one%meas%rh%qc
                                 t_value = next_one%meas%temperature%data
                                 t_qc    = next_one%meas%temperature%qc
                                 !Calculate dewpoint using formula in
                                 !obs_sort_module.F90
                                 !Since log(0) is undefined, if RH is near zero then set it to
                                 !some small value
                                 IF( rh_value .LT. very_small_rh ) THEN
                                  value = 1. / ( 1./t_value - Rv_over_L * LOG ( very_small_rh / 100. ) )
                                 ELSE
                                  value = 1. / ( 1./t_value - Rv_over_L * LOG ( rh_value / 100. ) )
                                 END IF
                                 !Since the purpose is to improve the QC of RH,
                                 !pull RH QC values
                                 qc = rh_qc
                              END IF
                           !BPR END
   
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
!     CALL error_handler ( error_number , error_message ,  &
!     fatal , listing )
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

   !BPR BEGIN
   TYPE ( measurement ) , POINTER :: next_next_one
   LOGICAL                 :: multi_level_ob
   INTEGER                 :: request_p_diff_use
   !BPR END

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

                        !BPR BEGIN
                        CASE ( 'PMSLPSFC' ) which_sfc_var
                           !Do not change the value because this case means that the program
                           !requested a sea level pressure based on the surface pressure
                           !so that quality control of surface pressure could be done since
                           !directly QC'ing surface pressure is difficult due to it's elevation
                           !dependence.  
                           !So although we want to update the surface pressure QC flag, we do not
                           !want to put a derived sea level pressure (which is what is in value)
                           !into surface pressure
                           !Surface pressure is in the header but is also the
                           !pressure of the surface level so we need to adjust both
                           obs%ground%psfc%qc   = qc
                           next_one%meas%pressure%qc   = qc
                        !BPR END

                        !BPR BEGIN
                        !Dewpoint is the moisture variable read in from the
                        !littler file, but RH is used internally and output to
                        !the OBS_DOMAIN file. However, the qc_obs* files use
                        !dewpoint so dewpoint is also important
                        !When we do QC on dewpoint we are trying to check for
                        !bad RH.  So if we are here we are trying to update
                        !the RH QC flags based on the dewpoint QC
                        CASE ( 'DEWPOINT' ) which_sfc_var
                            next_one%meas%dew_point%data = value
                            next_one%meas%dew_point%qc   =qc
                           next_one%meas%rh%qc   =qc
                        !BPR END
            
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

               !BPR BEGIN
               !Determine the pressure tolerance to use when storing observations
               !If this is a multi-level observation then one should effectively
               ! use no tolerance since the ob has already been interpolated to
               ! the pressure of the first-guess field so it should exactly match
               !If this is a single-level (above-surface) observation then we
               ! need to use the same tolerance used in pulling the ob in the
               ! first place.  However, we do not need to worry about determining
               ! whether the current pressure level is the one closest to the ob
               ! since we would not have retrieved this ob in the first place if
               ! the current request_level was not the closest pressure level
               !If the first level exists...
               IF ( still_more ) THEN
                !Make a pointer to the next level
                next_next_one => next_one%next
                !If this pointer is valid then this ob has at least two levels
                multi_level_ob = ASSOCIATED ( next_next_one)
                IF(multi_level_ob) THEN
                 request_p_diff_use = 1
                ELSE
                 request_p_diff_use = request_p_diff
                END IF
               END IF
               !BPR END

               !  Loop through the linked list of the data to find a correct pressure.

               search_for_level : DO WHILE ( still_more ) 
                  !BPR BEGIN
                  !found_level : IF ( ABS ( next_one%meas%pressure%data - request_level * 100 ) .LT. request_p_diff ) THEN
                  found_level : IF ( ABS ( next_one%meas%pressure%data - request_level * 100 ) .LT. request_p_diff_use ) THEN
                  !BPR END

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

                        !BPR BEGIN
                        !Dewpoint is the moisture variable read in from the
                        !littler file, but RH is used internally and output to
                        !the OBS_DOMAIN file. However, the qc_obs* files use
                        !dewpoint so dewpoint is also important
                        !When we do QC on dewpoint we are trying to check for
                        !bad RH.  So if we are here we are trying to update
                        !the RH QC flags based on the dewpoint QC
                        CASE ( 'DEWPOINT' ) which_var
                           next_one%meas%dew_point%data = value
                           next_one%meas%dew_point%qc = qc
                           next_one%meas%rh%qc = qc
                        !BPR END

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
