! ------------------------------------------------------------------------------

SUBROUTINE proc_oa ( t , u , v , rh , slp_x , &
pressure , &
iew_alloc , jns_alloc , kbu_alloc , &
date , time , fdda_loop , &
total_numobs , num_obs_found , total_dups , &
map_projection , obs , dxd , lat_center , &
print_oa , print_found_obs , print_obs_files , use_first_guess , & 
smooth_type              , smooth_sfc_wind          , & 
smooth_sfc_temp          , smooth_sfc_rh            , & 
smooth_sfc_slp           , smooth_upper_wind        , & 
smooth_upper_temp        , smooth_upper_rh          , &
oa_type                  , mqd_minimum_num_obs      , & 
mqd_maximum_num_obs      , oa_max_switch            , &
radius_influence         , oa_min_switch            )

!  This routine is a driver routine for objective analysis.

   USE obj_analysis
   USE observation

   IMPLICIT NONE

   INCLUDE 'error.inc'
   INCLUDE 'first_guess_size.inc'
   INCLUDE 'first_guess.inc'

   INTERFACE
      INCLUDE 'error.int'
      INCLUDE 'proc_ob_access.int'
   END INTERFACE

   REAL    , INTENT ( IN ), DIMENSION ( kbu_alloc )   :: pressure
   LOGICAL                                            :: print_oa           , & 
                                                         print_found_obs    , &
                                                         print_obs_files    , &
                                                         use_first_guess    , &
                                                         oa_min_switch      , &
                                                         oa_max_switch
   INTEGER                                            :: date               , &
                                                         time               , &
                                                         fdda_loop          , &
                                                         total_numobs       , & 
                                                         total_dups         , &
                                                         num_obs_found      , & 
                                                         map_projection     , & 
                                                         smooth_type        , & 
                                                         smooth_sfc_wind    , & 
                                                         smooth_sfc_temp    , & 
                                                         smooth_sfc_rh      , & 
                                                         smooth_sfc_slp     , & 
                                                         smooth_upper_wind  , & 
                                                         smooth_upper_temp  , & 
                                                         smooth_upper_rh    , &
                                                         mqd_minimum_num_obs, &
                                                         mqd_maximum_num_obs
   INTEGER , DIMENSION(10)                            :: radius_influence
   CHARACTER ( LEN = 132 )                            :: oa_type

   TYPE (report) , DIMENSION ( total_numobs )         :: obs
   REAL , INTENT(IN)                                  :: dxd                , &
                                                         lat_center
 
   REAL                  , DIMENSION ( num_obs_found ) :: obs_value          , &
                                                          xob                , &            
                                                          yob
   INTEGER               , DIMENSION ( num_obs_found ) :: array_index        , &
                                                          qc_flag          
   CHARACTER ( LEN =   8 ),DIMENSION ( num_obs_found ) :: station_id
   INTEGER               , PARAMETER                  :: num_var = 5
   CHARACTER ( LEN = 8 ) , DIMENSION ( num_var)       :: name
   INTEGER                                            :: passes             , &
                                                         num_scan
   INTEGER , PARAMETER                                :: max_scan = 10
   INTEGER               , DIMENSION ( num_var)       :: crsdot             , &
                                                         passes_sfc         , &
                                                         passes_upper
   REAL                  , DIMENSION ( num_var)       :: scale   
   REAL , DIMENSION ( iew_alloc , jns_alloc )         :: dum2d
   REAL , DIMENSION ( total_numobs)                   :: diff
   REAL                                               :: dxob               , &
                                                         dyob               , & 
                                                         aob  
   INTEGER                                            :: num_obs_pass       , & 
                                                         ivar               , & 
                                                         kp                 , & 
                                                         num                , & 
                                                         iob                , & 
                                                         job               
   LOGICAL                                            :: connected

   REAL , DIMENSION ( iew_alloc , jns_alloc )         :: u_banana           , &
                                                         v_banana
   CHARACTER (LEN=24)                                 :: date_time_char

   !  Names of all of the variables to objectively analyze.  Also, define each of the
   !  horizontal staggerings for the variables.

   name         = (/ 'U       '        , 'V       '        , 'T       '        , 'RH      '        , 'PSEALVLC'         /)
   crsdot       = (/      1            ,     1             ,      1            ,     1             ,      1             /)
   scale        = (/      1.           ,     1.            ,      1.           ,     1.            ,    100.            /)
   passes_sfc   = (/ smooth_sfc_wind   , smooth_sfc_wind   , smooth_sfc_temp   , smooth_sfc_rh     ,  smooth_sfc_slp    /)
   passes_upper = (/ smooth_upper_wind , smooth_upper_wind , smooth_upper_temp , smooth_upper_rh   ,      0             /)

   !  This is where we are outputting the information about the observations
   !  that went into this analysis, for this level, this variable, this time.
   !  First, check to see if this UNIT has been OPENed yet.  The first thing
   !  that we WRITE is the FORMAT to READ the data.

   IF ( print_obs_files ) THEN
      CALL make_date ( date , time , date_time_char )
      INQUIRE ( UNIT = 4 , OPENED = connected )
      IF ( .NOT. connected ) THEN
         IF ( fdda_loop.EQ.1) THEN
            OPEN ( UNIT   = 4                                      , &
                   FILE   = 'obs_used_for_oa_out_'//date_time_char , &
                   FORM   = 'FORMATTED'                            , &
                   ACCESS = 'SEQUENTIAL'                           , &         
                   STATUS = 'REPLACE'                                )
         ELSE
            OPEN ( UNIT   = 4                                      , &
                   FILE   = 'obs_used_for_oa_out_sfc_fdda_'//date_time_char , &
                   FORM   = 'FORMATTED'                            , &
                   ACCESS = 'SEQUENTIAL'                           , &         
                   STATUS = 'REPLACE'                                )
         END IF
   
         WRITE ( UNIT = 4 , &
         FMT = '( "( 3x,a8,3x,i6,3x,i5,3x,a8,3x,2(g13.6,3x),2(f7.2,3x),i7 )" )' )
      ENDIF
   
      INQUIRE ( UNIT =74 , OPENED = connected )
      IF ( .NOT. connected ) THEN
         IF ( fdda_loop.EQ.1) THEN
            OPEN ( UNIT   =74                                      , &
                   FILE   = 'plotobs_out_'//date_time_char , &
                   FORM   = 'FORMATTED'                            , &
                   ACCESS = 'SEQUENTIAL'                           , &
                   STATUS = 'REPLACE'                                )
         ELSE
            OPEN ( UNIT   =74                                      , &
                   FILE   = 'plotobs_out_sfc_fdda_'//date_time_char , &
                   FORM   = 'FORMATTED'                            , &
                   ACCESS = 'SEQUENTIAL'                           , &
                   STATUS = 'REPLACE'                                )
         END IF
   
         WRITE ( UNIT =74 , &
         FMT = '( "( 3x,a8,3x,i6,3x,i5,3x,a8,3x,g13.6,3x,16x,2(f7.2,3x),i7 )" )' )
      ENDIF
   END IF

   !  Loop through all analysis levels (remember that level kp=1 is the
   !  surface value of the 3-D field).

   vertical_level_loop : DO kp = 1 , kbu_alloc

      !  If we are doing surface FDDA, and this is one of the surface FDDA time
      !  periods, we can pop out of this loop once we are not doing the surface level.

      IF ( ( fdda_loop .GT. 1 ) .AND. ( kp .GT. 1 ) ) THEN
         EXIT vertical_level_loop
      END IF
      
      !  We need the u and v components of the horizontal wind for the Cressman
      !  scheme when we are doing a banana-type analysis.  The
      !  data going into this routine are still on the (y,x) orientation, but
      !  the 2-D fields on output (u_banana and v_banana) are oriented (x,y).
   
      CALL get_background_for_oa ( t , u , v , rh , slp_x , &
      u_banana , kp , name(1) , &
      iew_alloc , jns_alloc , kbu_alloc )
   
      CALL get_background_for_oa ( t , u , v , rh , slp_x , &
      v_banana , kp , name(2) , &
      iew_alloc , jns_alloc , kbu_alloc )

      !  Loop over each variable.
   
      variable_loop : DO ivar = 1 , num_var

         !  For sea level pressure, the analysis is only at one level, so cycle the loop
         !  (which should end this variable loop, but that is not important).


         IF ( ( name(ivar)(1:8) .EQ. 'PSEALVLC' ) .AND. ( kp .GT. 1 ) ) THEN
            CYCLE variable_loop
         END IF

         !  For the smoothing, the number of passes is a function of variable and level.
         !  This information is passed to the objective analysis schemes so that the
         !  perturbation fields may be smoothed before being added back to the base
         !  first-guess field.

         IF ( kp .EQ. 1 ) THEN
            passes = passes_sfc(ivar)
         ELSE
            passes = passes_upper(ivar)
         END IF

         !  Get obs for printing and plotting.

         CALL proc_ob_access ( 'use', name(ivar) , print_found_obs , &
         pressure(kp) , date , time , 1 , &
         200000 , &
         num_obs_found , num_obs_pass , obs , &
         iew_alloc , jns_alloc , kbu_alloc , &
         total_dups , map_projection , &
         get_value=obs_value , get_x_location=xob , get_y_location=yob , &
         get_id=station_id , get_array_index = array_index , get_qc_info = qc_flag )

         !  Write out the qc'ed data to a file.

         IF ( print_obs_files ) THEN
            WRITE ( UNIT =74 , FMT = '( " Number of Observations " / i8.8 )' ) num_obs_pass
            WRITE ( UNIT =74 , FMT = '( "   Variable   Press    Obs     Station   &
            &       Obs                        X          Y        QC    " )' )
            WRITE ( UNIT =74 , FMT = '( "    Name      Level    Number    ID      &
            &      Value                     Location  Location  Value   " )' )
            station_loop_74 : DO num = 1, num_obs_pass
               WRITE ( UNIT =74 , FMT = '( 3x,a8,3x,i6,3x,i5,3x,a8,3x,g13.6,3x,16x,2(f7.2,3x),i7 )' ) &
                 name(ivar) , NINT ( pressure(kp) ) , num , station_id(num) , &
                 obs_value(num) , xob(num) , yob(num) , qc_flag(num)
            END DO station_loop_74
         END IF


         !  Obtain observations for kp level and for variable ivar for this
         !  time period.
 
         CALL proc_ob_access ( 'use', name(ivar) , print_found_obs , &
         pressure(kp) , date , time , 1 , &
         MIN ( fails_error_max , fails_buddy_check ) , &
         num_obs_found , num_obs_pass , obs , &
         iew_alloc , jns_alloc , kbu_alloc , &
         total_dups , map_projection , &
         get_value=obs_value , get_x_location=xob , get_y_location=yob , &
         get_id=station_id , get_array_index = array_index , get_qc_info = qc_flag ) 

         !  To do the print out, objective analysis, clean up, and storage, we need to have
         !  at least a minimum for the MQD technique.  We also need to have NOT too many
         !  obs for MQD (we may switch to Cressman in that case).

         mqd_num_obs : IF ( ( ( num_obs_pass .GE. mqd_minimum_num_obs ) .AND. &
                              ( num_obs_pass .LE. mqd_maximum_num_obs ) )     &
                                             .OR.                             &
                              ( oa_min_switch .AND.                           &
                              ( num_obs_pass .LT. mqd_minimum_num_obs ) )     &
                                             .OR.                             &
                              ( oa_max_switch .AND.                           &
                              ( num_obs_pass .GT. mqd_maximum_num_obs ) )     &
                                             .OR.                             &
                              ( oa_type      .NE. 'MQD'               ) ) THEN

            !  Obtain first guess for kp level and for variable name(ivar).  The
            !  data going into this routine are still on the (y,x) orientation, but
            !  the 2-D field on output (dum2d) is oriented (x,y).
   
            CALL get_background_for_oa ( t , u , v , rh , slp_x , &
            dum2d , kp , name(ivar) , &
            iew_alloc , jns_alloc , kbu_alloc )
   
            !  Calculate the difference between the observations and first guess fields.
            !  Estimate the values of the first guess at each of the observation location.  
            !  This is performed with the four surrounding analysis points, using a linear
            !  "x" and a linear "y" weighting based on distance.  We are either using the
            !  first-guess to create perturbations, or we can do an analysis from the
            !  observations only.
   
            IF ( use_first_guess ) THEN
               station_loop_1fg : DO num = 1, num_obs_pass
      
                  iob  = INT ( xob (num) )
                  job  = INT ( yob (num) )
                  dxob = xob ( num ) - REAL ( iob )
                  dyob = yob ( num ) - REAL ( job )
                  aob  = ( 1.-dxob ) * ( ( 1.-dyob ) * dum2d ( iob   , job   )   + &
                                              dyob   * dum2d ( iob   , job+1 ) ) + &
                              dxob   * ( ( 1.-dyob ) * dum2d ( iob+1 , job   )   + &
                                              dyob   * dum2d ( iob+1 , job+1 ) )
                  diff ( num )  = obs_value(num) / scale(ivar) - aob
         IF ( ( name(ivar)(1:8) .EQ. 'PSEALVLC' ) ) THEN
         if ( iob >= 11 .and. iob <= 29 ) then
         if ( job >= 11 .and. job <= 29 ) then
         ENDIF
         ENDIF
         ENDIF
      
               END DO station_loop_1fg
            ELSE
               station_loop_1nfg : DO num = 1, num_obs_pass
      
                  iob  = INT ( xob (num) )
                  job  = INT ( yob (num) )
                  dxob = xob ( num ) - REAL ( iob )
                  dyob = yob ( num ) - REAL ( job )
                  aob  = ( 1.-dxob ) * ( ( 1.-dyob ) * dum2d ( iob   , job   )   + &
                                              dyob   * dum2d ( iob   , job+1 ) ) + &
                              dxob   * ( ( 1.-dyob ) * dum2d ( iob+1 , job   )   + &
                                              dyob   * dum2d ( iob+1 , job+1 ) )
                  diff ( num )  = obs_value(num) / scale(ivar)
      
               END DO station_loop_1nfg
            END IF

            !  If we requested a print out of every observation used for OA, this is where
            !  she goes.  We know all of the observations to be used, value, name, and
            !  difference from the background analysis.
   
            IF ( print_oa ) THEN
               WRITE ( UNIT = 6 , FMT = '( "   Variable   Press    Obs     Station   &
               &       Obs           Obs-1st      X          Y        QC    " )' )
               WRITE ( UNIT = 6 , FMT = '( "    Name      Level    Number    ID      &
               &      Value           Guess     Location  Location  Value   " )' )
               station_loop_3 : DO num = 1, num_obs_pass
                  WRITE ( UNIT = 6 , FMT = '( 3x,a8,3x,i6,3x,i5,3x,a8,3x,2(g13.6,3x),2(f7.2,3x),i7 )' ) &
                  name(ivar) , NINT ( pressure(kp) ) , num , station_id(num) , &
                  obs_value(num) , diff(num) , xob(num) , yob(num) , qc_flag(num)
               END DO station_loop_3

            END IF

            !  Output data that is used in the objective analysis to a file.

            IF ( print_obs_files ) THEN
               WRITE ( UNIT = 4 , FMT = '( " Number of Observations " / i8.8 )' ) num_obs_pass
               WRITE ( UNIT = 4 , FMT = '( "   Variable   Press    Obs     Station   &
               &       Obs           Obs-1st      X          Y        QC    " )' )
               WRITE ( UNIT = 4 , FMT = '( "    Name      Level    Number    ID      &
               &      Value           Guess     Location  Location  Value   " )' )
               station_loop_4 : DO num = 1, num_obs_pass
                  WRITE ( UNIT = 4 , FMT = '( 3x,a8,3x,i6,3x,i5,3x,a8,3x,2(g13.6,3x),2(f7.2,3x),i7 )' ) &
                  name(ivar) , NINT ( pressure(kp) ) , num , station_id(num) , &
                  obs_value(num) , diff(num) , xob(num) , yob(num) , qc_flag(num)
               END DO station_loop_4
            END IF
   
            !  This is the multiquadric solver.
 
            IF      ( ( oa_type      .EQ. 'MQD'               ) .AND.  &
                      ( num_obs_pass .GE. mqd_minimum_num_obs ) .AND.  &
                      ( num_obs_pass .LE. mqd_maximum_num_obs ) ) THEN
   
               CALL mqd ( diff , xob , yob , num_obs_pass , &
               dum2d , iew_alloc , jns_alloc , &
               crsdot(ivar) , name(ivar) , passes , smooth_type , use_first_guess )
   
            !  This is the Cressman method.  When multiple scans are utilized,
            !  the scheme is called repeatedly, and the "final" analysis is 
            !  used as input for the new difference field that is objectively
            !  analyzed.
      
            ELSE IF ( ( oa_type .EQ. 'Cressman' ) .OR. &
                    ( ( oa_min_switch ) .AND. ( num_obs_pass .LT. mqd_minimum_num_obs ) ) .OR. &
                    ( ( oa_max_switch ) .AND. ( num_obs_pass .GT. mqd_maximum_num_obs ) ) ) THEN

               multi_scan : DO num_scan = 1 , max_scan
   
                  IF ( radius_influence(num_scan) .LT. 0 ) THEN
                     EXIT multi_scan
                  END IF

                  CALL cressman ( obs_value , diff , xob , yob , num_obs_pass , &
                  dum2d , iew_alloc , jns_alloc , &
                  crsdot(ivar) , name(ivar) , radius_influence(num_scan) , &
                  dxd , u_banana , v_banana , pressure(kp) , passes , smooth_type , lat_center , &
                  use_first_guess )

                  more_than_1_scan : IF ( max_scan .GT. 1 ) THEN

                     IF ( use_first_guess ) THEN
                        station_loop_5fg : DO num = 1, num_obs_pass
         
                           iob  = INT ( xob (num) )
                           job  = INT ( yob (num) )
                           dxob = xob ( num ) - REAL ( iob )
                           dyob = yob ( num ) - REAL ( job )
                           aob  = ( 1.-dxob ) * ( ( 1.-dyob ) * dum2d ( iob   , job   )   + &
                                                       dyob   * dum2d ( iob   , job+1 ) ) + &
                                       dxob   * ( ( 1.-dyob ) * dum2d ( iob+1 , job   )   + &
                                                       dyob   * dum2d ( iob+1 , job+1 ) )
                           diff ( num )  = obs_value(num) / scale(ivar) - aob
         
                        END DO station_loop_5fg
                     ELSE
                        station_loop_5nfg : DO num = 1, num_obs_pass
         
                           iob  = INT ( xob (num) )
                           job  = INT ( yob (num) )
                           dxob = xob ( num ) - REAL ( iob )
                           dyob = yob ( num ) - REAL ( job )
                           aob  = ( 1.-dxob ) * ( ( 1.-dyob ) * dum2d ( iob   , job   )   + &
                                                       dyob   * dum2d ( iob   , job+1 ) ) + &
                                       dxob   * ( ( 1.-dyob ) * dum2d ( iob+1 , job   )   + &
                                                       dyob   * dum2d ( iob+1 , job+1 ) )
                           diff ( num )  = obs_value(num) / scale(ivar)
         
                        END DO station_loop_5nfg
                     END IF

                  END IF more_than_1_scan
               END DO multi_scan

            !  This is where we drop through if we are not doing the MQD or the 
            !  Cressman objective analysis methods.  Just to let folks know, we
            !  should mention something to that effect.
   
            ELSE

               error_number = 00381003
               error_message(1:31) = 'proc_oa                        '
               WRITE ( error_message(32:) , FMT = '(A,A,A,I6,A)' ) &
               'Hmmm, not doing either the MQD of Cressman OA for field ',&
               name(ivar)(1:8),' at pressure level ',NINT(pressure(kp)),'.'
               fatal = .false.
               listing = .false.
               call error_handler ( error_number , error_message , fatal , listing )

            END IF

            !  After the OA, we go over the generated RH field to be sure that it is 
            !  constrained between 0% and 100%.
      
            IF ( name(ivar)(1:8) .EQ. 'RH      ' ) THEN
               CALL clean_rh ( dum2d , iew_alloc , jns_alloc , 0. , 100. )
            END IF
   
            !  Store the final analysis for kp level and for variable name(ivar).
   
            CALL put_background_from_oa ( t , u , v , rh , slp_x , &
            dum2d , kp , name(ivar) , &
            iew_alloc , jns_alloc , kbu_alloc )

         ELSE IF ( num_obs_pass .LT. mqd_minimum_num_obs ) THEN mqd_num_obs

            error_number = 00381001
            error_message(1:31) = 'proc_oa                        '
            error_message(32:)  = ' Less than the required number of observations, no MQD done.'
            fatal = .false.
            listing = .false.
            call error_handler ( error_number , error_message , fatal , listing )

         ELSE IF ( num_obs_pass .GT. mqd_maximum_num_obs ) THEN mqd_num_obs

            error_number = 00381002
            error_message(1:31) = 'proc_oa                        '
            error_message(32:)  = ' More than the permitted number of observations, no MQD done.'
            fatal = .false.
            listing = .false.
            call error_handler ( error_number , error_message , fatal , listing )

         END IF mqd_num_obs

      END DO variable_loop

   END DO vertical_level_loop

   !  We need to close our sample print-out file.

   IF ( print_obs_files ) THEN
      CLOSE ( 4 ) 
      CLOSE ( 74 ) 
   END IF

END SUBROUTINE proc_oa
