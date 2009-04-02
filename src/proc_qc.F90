!------------------------------------------------------------------------------
SUBROUTINE proc_qc ( iew_alloc , jns_alloc , kbu_alloc , number_of_obs , &
                     total_dups , map_projection , &
                     qc_test_error_max        , qc_test_buddy          , &
                     qc_test_vert_consistency , qc_test_convective_adj , &
                     max_error_t      , max_error_uv    , &
                     max_error_z      , max_error_rh    , &
                     max_error_p      , print_error_max , &
                     max_buddy_t      , max_buddy_uv    , &
                     max_buddy_z      , max_buddy_rh    , &
                     max_buddy_p      , print_buddy     , &
                     print_found_obs  ,                   &
                     print_vert       , print_dry       , & 
                     pressure , date , time , dx , buddy_weight , &
                     obs , index , max_number_of_obs , & 
                     t , u , v , h , rh , slp_x , sst , tobbox , odis )

! Driver routine for QC
!   
!    This routine is a driver routine for quality control ( or QC )
!       of observations. It includes:
!       1. vertical sounding checks
!       2. error checks by comparing observations and the first 
!          guess fields. This includes:
!          a) error max check, and 
!          b) buddy check.
!    Note that this routine assumes that surface is lowest level
!       in the 3D fields.

   USE qc0
   USE qc1
   USE qc2
   USE qc3

   USE observation

   IMPLICIT NONE

   !  Data from the calling routine.  All of this is input.
  
   !  iew_alloc, jns_alloc: 1st and 2nd dimensions of any horizontal array
   !  kbu_alloc:    number of vertical levels
   !  errmxp:       maximum error allowed for SLP
   !  errmxt:       maximum error allowed for temperature
   !  errmxw:       maximum error allowed for winds
   !  errmxrh:      maximum error allowed for RH
   !  pressure:     array of pressure 
   !  time:         hhmm time (UTC)
   !  dx:           grid distance in km
   !  buddy_weight: weight factor for buddy check

   INCLUDE 'first_guess_size.inc'
   INCLUDE 'first_guess.inc'

   INTEGER , INTENT ( IN )                    :: number_of_obs ,     & 
                                                 total_dups ,        &
                                                 map_projection
   LOGICAL , INTENT ( IN )                    :: qc_test_error_max         , &
                                                 qc_test_buddy             , &
                                                 qc_test_vert_consistency  , &
                                                 qc_test_convective_adj
   REAL    , INTENT ( IN )                    :: max_error_p,             &
                                                 max_error_t,             &
                                                 max_error_uv,            &
                                                 max_error_z,             &
                                                 max_error_rh,            &
                                                 max_buddy_p,             &
                                                 max_buddy_t,             &
                                                 max_buddy_uv,            &
                                                 max_buddy_z,             &
                                                 max_buddy_rh
   LOGICAL , INTENT ( IN )                    :: print_error_max,         &
                                                 print_buddy    ,         &
                                                 print_found_obs ,        &
                                                 print_vert      ,        &
                                                 print_dry
   REAL    , INTENT ( IN ) , DIMENSION (kbu_alloc) :: pressure
   INTEGER , INTENT ( IN )                    :: date , time
   REAL    , INTENT ( IN )                    :: dx , &
                                                 buddy_weight
   INTEGER       , INTENT ( IN )                  :: max_number_of_obs
   TYPE (report) , DIMENSION (max_number_of_obs)  :: obs
   INTEGER       , DIMENSION (max_number_of_obs)  :: index
   REAL          , DIMENSION(jns_alloc,iew_alloc) :: tobbox, odis

   !  Data from the call to the routine to provide all of the information
   !  for the the observations that we will need.

   !  num_obs_found: number of observations that fit the requested criteria
   !  obs:          station observation array
   !  xob,yob:      x, y locations of station obs
   !  lonob:        longitude of the observation, for local time computation
   !  station_id:   station identifier
  
   REAL          , DIMENSION (max_number_of_obs)           :: obs_value
   INTEGER       , DIMENSION (max_number_of_obs)           :: index_criteria
   INTEGER                                                 :: num_obs_found
   REAL          , DIMENSION (max_number_of_obs)           :: xob, yob, lonob
   CHARACTER ( LEN =   8 ) , DIMENSION (max_number_of_obs) :: station_id
   INTEGER       , DIMENSION (max_number_of_obs)           :: qc_flag
   LOGICAL       , DIMENSION (max_number_of_obs)           :: ship_report

   !  Data from the routine that supplies the background field for
   !  QC checking

   !  gridded:      background/first field as input, final analysis as output

   REAL          , DIMENSION (iew_alloc,jns_alloc)        :: gridded

   !  Internally computed value.

   !  error_difference:  initial maximum difference between observations and
   !                   the interpolated value of the analysis at the ob
   !                   location.
   !  ivar:            1=U, 2=V, 3=T, 4=RH, 5=SLP
   !  kp:              vertical index
   !  local_time:      integer time (hhmm), crudely corrected to the local
   !                   time based upon the longitude
   !  name:            variable name, TT, RH, PMSL, UU, or VV

   REAL                                       :: error_difference , &
                                                 buddy_difference
   INTEGER                                    :: ivar , kp , i
   INTEGER        , DIMENSION (max_number_of_obs)        :: local_time
   CHARACTER ( LEN =   8 )                    :: name

   INTERFACE
      INCLUDE 'proc_ob_access.int'
   END INTERFACE

   !  The first quality control (QC) that can be performed is with
   !  data that is vertically stacked.  These reports have the temperature,
   !  speed and direction compared against reasonable benchmarks.
   !  Bad values have flags set, though switching the sign of the temperature
   !  (in degrees C) is allowed.  This test is not performed if the entire
   !  report was discarded, and no error checks are performed on bogus
   !  data types.  Since the vertical consistency check and the dry 
   !  convective adjustment are for vertically stacked data, we can ignore 
   !  these test for data with only one level.

   IF ( qc_test_vert_consistency ) THEN
      loop_all_1 : DO i = 1, number_of_obs
         valid_ob_1 : IF ( ( .NOT. obs(i)%info%discard     )   .AND. &
                           ( .NOT. obs(i)%info%bogus       )   .AND. & 
                           ( obs(i)%info%is_sound          )   .AND. & 
                           ( ASSOCIATED ( obs(i)%surface ) ) ) THEN
            CALL vert_cons_check ( obs ( i ) , i , print_vert )
         END IF valid_ob_1
      END DO loop_all_1
   END IF

   !  The second QC check that is available to be performed on the 
   !  the observations is the convective adjustment.  This is used
   !  only for the upper level temperature reports.

   IF ( qc_test_convective_adj ) THEN
      loop_all_2 : DO i = 1, number_of_obs
         valid_ob_2 : IF ( ( .NOT. obs(i)%info%discard     )   .AND. &
                           ( .NOT. obs(i)%info%bogus       )   .AND. & 
                           ( obs(i)%info%is_sound          )   .AND. & 
                           ( ASSOCIATED ( obs(i)%surface ) ) ) THEN
!           CALL dry_convective_adjustment ( obs ( i ) , i , print_dry )
         END IF valid_ob_2
      END DO loop_all_2
   END IF

   !  This loop is is only processed if the QC tests inside the variable
   !  and level loop are included in the requested QC tests to perform.

   any_checks_to_do : IF ( ( qc_test_error_max ) .OR. ( qc_test_buddy ) ) THEN

      !  Loop through all analysis levels and variable types.  This is the
      !  outer loop over all of the data for the error_max and buddy_check
      !  error flagging.  This loop is is only processed if the QC tests

      vertical_level_loop : DO kp = 1, kbu_alloc
   
         variable_loop : DO ivar = 1, 5
   
            ! do error checks only for ivar = PMSL and for kp = 1
   
            IF ( ( kp .GT. 1 ) .AND. ( ivar .EQ. 5 ) ) EXIT variable_loop
   
            which_variable : SELECT CASE ( ivar )
               CASE ( 1 ) 
                  name = 'UU      '
                  error_difference = max_error_uv
                  buddy_difference = max_buddy_uv
                  CALL yx2xy (     u(1,1,kp) , jns_alloc , iew_alloc , gridded )
               CASE ( 2 ) 
                  name = 'VV      '              
                  error_difference = max_error_uv
                  buddy_difference = max_buddy_uv
                  CALL yx2xy (     v(1,1,kp) , jns_alloc , iew_alloc , gridded )
               CASE ( 3 ) 
                  name = 'TT      '
                  error_difference = max_error_t
                  buddy_difference = max_buddy_t
                  CALL yx2xy (     t(1,1,kp) , jns_alloc , iew_alloc , gridded )
               CASE ( 4 ) 
                  name = 'RH      '
                  IF ( pressure(kp) .EQ. 300 ) THEN
                     error_difference = max_error_rh * 2
                  ELSE IF ( pressure(kp) .LT. 300 ) THEN
                     error_difference = max_error_rh * 3
                  ELSE
                     error_difference = max_error_rh
                  END IF
                  buddy_difference = max_buddy_rh
                  CALL yx2xy (    rh(1,1,kp) , jns_alloc , iew_alloc , gridded )
               CASE ( 5 ) 
                  name = 'PMSL    '
                  error_difference = max_error_p
                  buddy_difference = max_buddy_p
                  CALL yx2xy ( slp_x(1,1   ) , jns_alloc , iew_alloc , gridded )
                  gridded = gridded * 100
            END SELECT which_variable
   
            !  Obtain observations for kp level and for variable ivar for this
            !  time period.  
   
            CALL proc_ob_access ( 'get', name , print_found_obs , &
            pressure(kp) , date , time , 1 , 2**20 , number_of_obs , num_obs_found , obs , &
            iew_alloc , jns_alloc , kbu_alloc , &
            total_dups , map_projection , &
            get_value=obs_value , get_x_location=xob , get_y_location=yob , &
            get_longitude=lonob , get_array_index=index_criteria , &
            get_over_water=ship_report , get_id=station_id , get_qc_info=qc_flag )
   
            !  The local time will be used to modify the acceptable differences
            !  between the observations and the first guess analysis.
   
            CALL local ( time / 100 , lonob , local_time , num_obs_found ) 

            !  Accumulate the number of observations at each grid point for the
            !  surface FDDA.

            IF ( kp .EQ. 1 ) THEN
               CALL ob_density   ( xob , yob , dx , num_obs_found , tobbox , iew_alloc , jns_alloc )
            END IF
   
            !  Perform the maximum error difference QC test with the available
            !  observations and the background field analysis.
       
            IF ( qc_test_error_max ) THEN
               CALL error_max ( station_id , obs_value , ship_report , xob , yob , & 
               index_criteria , qc_flag , num_obs_found , &
               gridded , iew_alloc , jns_alloc , &
               name , error_difference , pressure(kp) , local_time , print_error_max )
            END IF
   
            !  Perform the buddy-check QC test.  The compares each observation
            !  with it's nearest neighbors, using the difference between the 
            !  observation value and the interpolated analysis as the metric.
   
            IF ( qc_test_buddy ) THEN
               CALL buddy_check ( station_id, obs_value , num_obs_found , xob, yob, qc_flag,   &
               gridded, iew_alloc , jns_alloc , name,                      &
               pressure(kp), local_time, dx , buddy_weight , buddy_difference , print_buddy )
            END IF
   
            !  After the observations have been QC'ed, we send back the 
            !  qc_flags information to be stored.

            CALL proc_ob_access ( 'put', name , print_found_obs , &
            pressure(kp) , date , time , 1 , 2**20 , number_of_obs , num_obs_found , obs , &
            iew_alloc , jns_alloc , kbu_alloc , &
            total_dups , map_projection , &
            put_value=obs_value , put_array_index=index_criteria , put_qc_info=qc_flag )
   
         END DO variable_loop
   
      END DO vertical_level_loop

   END IF any_checks_to_do

   !  We have the sum of five surface variables in the tobbox (u, v, t, rh and slp).  So,
   !  do a no-brainer divide by 5 to get an average of the observation density
   !  for the surface FDDA.

   tobbox = tobbox / 5.
   WHERE ( tobbox .LT. 1 ) tobbox = 0

   ! Now get the distance to the obs so we can calculate the confidence in the analysis
   CALL obs_distance ( tobbox , odis , iew_alloc , jns_alloc , dx )

   !  Now that each variable has gone throught the quality control checks individually,
   !  we should make sure that there is consistency between the variables that are
   !  related: u and v, T and relative humidity.

   CALL qc_consistency ( obs , number_of_obs )

END SUBROUTINE proc_qc
