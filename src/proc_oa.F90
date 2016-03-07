! ------------------------------------------------------------------------------

SUBROUTINE proc_oa ( t , u , v , rh , slp_x , &
!BPR BEGIN
!Add 3D pressure of first guess to fields passed in 
!pressure is 2D pressure with surface pressure just set to 1001
!pressure , &
pressure , pres , &
!BPR END
iew_alloc , jns_alloc , kbu_alloc , &
date , time , fdda_loop , mqd_count , mqd_abs_min , &
total_numobs , num_obs_found , total_dups , &
map_projection , obs , dxd , lat_center , &
!BPR BEGIN
use_p_tolerance_one_lev , &
!BPR END
print_oa , print_found_obs , print_obs_files , use_first_guess , & 
smooth_type              , smooth_sfc_wind          , & 
smooth_sfc_temp          , smooth_sfc_rh            , & 
smooth_sfc_slp           , smooth_upper_wind        , & 
smooth_upper_temp        , smooth_upper_rh          , &
oa_type                  , oa_3D_type               , &
mqd_minimum_num_obs      , mqd_maximum_num_obs      , &
oa_max_switch            , radius_influence         , &
oa_min_switch            , oa_3D_option             , &
!BPR BEGIN
!grid_id )
grid_id, terrain, h, scale_cressman_rh_decreases, radius_influence_sfc_mult, &
oa_psfc, max_p_tolerance_one_lev_oa )
!BPR END

!  This routine is a driver routine for objective analysis.

   USE obj_analysis
   USE observation
!BPR BEGIN
   USE final_analysis, only : mixing_ratio
!BPR END

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
   !BPR BEGIN
   !                                                     oa_max_switch
                                                         oa_max_switch      , &
                                                         scale_cressman_rh_decreases , &
                                                         oa_psfc
   REAL , INTENT ( IN )                               :: radius_influence_sfc_mult
   INTEGER , INTENT ( IN )                            :: max_p_tolerance_one_lev_oa
   LOGICAL , INTENT ( IN )                            :: use_p_tolerance_one_lev
   !BPR END
   INTEGER                                            :: oa_3D_option
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
                                                         mqd_maximum_num_obs, &
                                                         grid_id
   INTEGER , DIMENSION(10)                            :: radius_influence
   CHARACTER ( LEN = 132 )                            :: oa_type, oa_3D_type, oa_type_tmp

   TYPE (report) , DIMENSION ( total_numobs )         :: obs
   REAL , INTENT(IN)                                  :: dxd                , &
                                                         lat_center
 
   REAL                  , DIMENSION ( num_obs_found ) :: obs_value          , &
                                                          xob                , &            
                                                          yob
   INTEGER               , DIMENSION ( num_obs_found ) :: array_index        , &
                                                          qc_flag          
   CHARACTER ( LEN =   8 ),DIMENSION ( num_obs_found ) :: station_id
   !BPR BEGIN
   !INTEGER               , PARAMETER                  :: num_var = 5
   INTEGER               , PARAMETER                  :: num_var = 6
   !BPR END
   CHARACTER ( LEN = 8 ) , DIMENSION ( num_var)       :: name
   INTEGER                                            :: passes             , &
                                                         num_scan
   INTEGER , PARAMETER                                :: max_scan = 10
   INTEGER               , DIMENSION ( num_var)       :: crsdot             , &
                                                         passes_sfc         , &
                                                         passes_upper
   REAL                  , DIMENSION ( num_var)       :: scale   
   REAL , DIMENSION ( iew_alloc , jns_alloc )         :: dum2d
   !BPR BEGIN
   REAL , DIMENSION ( iew_alloc-1 , jns_alloc-1 )     :: dum2d_fg, t_first_guess, t_analysis, pressure_2d
   
   REAL                                               :: pressure_t
   REAL, DIMENSION(jns_alloc,iew_alloc,kbu_alloc)     :: qv3d_first_guess_yx
   REAL, DIMENSION(jns_alloc,iew_alloc)               :: psfc_first_guess_model_terrain_yx
   REAL, DIMENSION(iew_alloc,jns_alloc)               :: psfc_first_guess_model_terrain_xy
   REAL                  , DIMENSION ( num_obs_found ) :: fg_value  
   !BPR END
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
   CHARACTER ( LEN = 132 )                            :: tmp_filename

   LOGICAL                                            :: skip_to_cressman, first_time
   INTEGER                                            :: num_mqd_uu , num_mqd_vv , &
                                                         num_mqd_tt , num_mqd_rh , &
                                                         min_mqd , max_mqd       , &
                                                         mqd_count , mqd_abs_min , &
                                                         test_count
   !BPR BEGIN
   REAL                                       :: radius_influence_scale_factor
   REAL                                       :: radius_influence_sfc_min_cells, radius_influence_sfc_max_cells
   REAL                                       :: radius_influence_scaled, radius_influence_scaled_scan1
   LOGICAL, DIMENSION ( kbu_alloc,max_scan )  :: roi_not_printed_yet

   REAL, DIMENSION(jns_alloc,iew_alloc) :: slp_x_pa

   INTEGER                                    :: curx, cury
   INTEGER                                    :: num_last_3d_variable !Index to last of variables under consideration 
                                                                      !that is 3D

   !Maximum difference allowed between pressure of ob and pressure of background
   !data it will be incorporated with for the objective analysis (Pa)
   !Default is 1
   INTEGER                                    :: request_p_diff

   !If the user wants to allow a tolerance between an obs pressure and the
   !pressure level it can be used for an objective anlysis on, then use the user-specified tolerance.
   !If not, then use a tolerance of 1 Pa, which is effectively no tolerance.

   IF( use_p_tolerance_one_lev  ) THEN
    request_p_diff = max_p_tolerance_one_lev_oa
   ELSE
    request_p_diff = 1
   END IF

   roi_not_printed_yet(:,:) = .TRUE. 
   !BPR END

   !BPR BEGIN
   !It seems that we should default to NOT skipping to Cressman unless we find
   !evidence that we need to
   !In most cases what we initialize this to should not matter since:
   !1) If the user chose MQD then the code will set skip_to_cressman as appropriate
   !2) If the user chose Cressman, then skipping to Cressman will be consistent
   !   with what the user chose
   !However, this switch seems to be intended to be used to force a switch to
   !Cressman when needed and thus it is most straightforward to only set it to
   !TRUE when the situtation indicates that it is needed.  Also, if there are
   !any other analysis methodologies that do not explicitly set skip_to_cressman, 
   !those methodologies will be skipped over due to skip_to_cressman being
   !defaulted to true.  This situation may occur for the choice of "none" added in this
   !version of Obsgrid indicating that no objective analysis should be done.
   !skip_to_cressman = .TRUE.
   skip_to_cressman = .FALSE.
   !BPR END

    num_mqd_uu = 0
    num_mqd_vv= 0
    num_mqd_tt = 0
    num_mqd_rh = 0
    min_mqd = 100
    max_mqd = 0
   

   !  Names of all of the variables to objectively analyze.  Also, define each of the
   !  horizontal staggerings for the variables.

   !  BPR BEGIN
   !  Add surface pressure 
   name         = (/ 'UU      '        , 'VV      '        , 'TT      '        , 'RH      '        ,&
                     'PMSL    '        , 'PSFC    '        /)
   crsdot       = (/      1            ,     1             ,      1            ,     1             ,&
                          1            ,     1             /)
   scale        = (/      1.           ,     1.            ,      1.           ,     1.            ,&
                        100.           ,     1.            /)
   passes_sfc   = (/ smooth_sfc_wind   , smooth_sfc_wind   , smooth_sfc_temp   , smooth_sfc_rh     ,&
                     smooth_sfc_slp    , smooth_sfc_slp    /)
   passes_upper = (/ smooth_upper_wind , smooth_upper_wind , smooth_upper_temp , smooth_upper_rh   ,&
                          0            ,     0             /)
   !name         = (/ 'UU      '        , 'VV      '        , 'TT      '        , 'RH      '        , 'PMSL    '         /)
   !crsdot       = (/      1            ,     1             ,      1            ,     1             ,      1             /)
   !scale        = (/      1.           ,     1.            ,      1.           ,     1.            ,    100.            /)
   !passes_sfc   = (/ smooth_sfc_wind   , smooth_sfc_wind   , smooth_sfc_temp   , smooth_sfc_rh     ,  smooth_sfc_slp    /)
   !passes_upper = (/ smooth_upper_wind , smooth_upper_wind , smooth_upper_temp , smooth_upper_rh   ,      0             /)

   !Index to last of the 3D variables
   num_last_3d_variable = 4
   !BPR END

   !  This is where we are outputting the information about the observations
   !  that went into this analysis, for this level, this variable, this time.
   !  First, check to see if this UNIT has been OPENed yet.  The first thing
   !  that we WRITE is the FORMAT to READ the data.

   IF ( print_obs_files ) THEN
      CALL make_date ( date , time , date_time_char )
      !INQUIRE ( UNIT = 4 , OPENED = connected )
      !IF ( .NOT. connected ) THEN
      !   IF ( fdda_loop.EQ.1) THEN
      !      OPEN ( UNIT   = 4                                      , &
      !             FILE   = 'obs_used_for_oa_out_'//date_time_char , &
      !             FORM   = 'FORMATTED'                            , &
      !             ACCESS = 'SEQUENTIAL'                           , &         
      !             STATUS = 'REPLACE'                                )
      !   ELSE
      !      OPEN ( UNIT   = 4                                      , &
      !             FILE   = 'obs_used_for_oa_out_sfc_fdda_'//date_time_char , &
      !             FORM   = 'FORMATTED'                            , &
      !             ACCESS = 'SEQUENTIAL'                           , &         
      !             STATUS = 'REPLACE'                                )
      !   END IF
      !
      !   WRITE ( UNIT = 4 , &
      !   FMT = '( "( 3x,a8,3x,i6,3x,i5,3x,a8,3x,2(g13.6,3x),2(f7.2,3x),i7 )" )' )
      !ENDIF
   
      INQUIRE ( UNIT =74 , OPENED = connected )
      IF ( .NOT. connected ) THEN
         WRITE (tmp_filename,'("plotobs_out.d",i2.2,".")') grid_id
         OPEN ( UNIT   =74                                      , &
                FILE   = trim(tmp_filename)//date_time_char , &
                FORM   = 'FORMATTED'                            , &
                ACCESS = 'SEQUENTIAL'                           , &
                STATUS = 'REPLACE'                                )
  
         WRITE ( UNIT =74 , &
         FMT = '( "( 3x,a8,3x,i6,3x,i5,3x,a8,3x,2(g13.6,3x),2(f7.2,3x),i7 )" )' )
        !FMT = '( "( 3x,a8,3x,i6,3x,i5,3x,a8,3x,g13.6,3x,16x,2(f7.2,3x),i7 )" )' )
      ENDIF
   END IF


   IF ( oa_3D_type .EQ. 'MQD' .AND. fdda_loop == 1 ) THEN
     skip_to_cressman = .FALSE.
     !BPR BEGIN
     !DO ivar = 1 , num_var-1
     !Loop over the 3D variables
     DO ivar = 1 , num_last_3d_variable
     !BPR END
       first_time = .TRUE.
       test_count = 0
       DO kp = 2 , kbu_alloc
         CALL proc_ob_access ( 'use', name(ivar) , print_found_obs , &
         !BPR BEGIN
         !pressure(kp) , date , time , 1 , &
         pressure(kp) , date , time , request_p_diff , &
         !BPR END
         MIN ( fails_error_max , fails_buddy_check ) , &
         num_obs_found , num_obs_pass , obs , &
         iew_alloc , jns_alloc , kbu_alloc , &
         total_dups , map_projection , &
         get_value=obs_value , get_x_location=xob , get_y_location=yob , &
         get_id=station_id , get_array_index = array_index , get_qc_info = qc_flag ) 
         IF ( num_obs_pass .LT. mqd_minimum_num_obs ) THEN
           IF ( first_time .AND. ivar == 1 ) THEN
             WRITE(*,*)
             WRITE(*,'("   WARNING: OA scheme for upper levels are set to MQD, but due to a lack of observations")' )
             WRITE(*,'("            the following variables/levels will revert to Cressman scheme: ")' )
           END IF
           IF ( first_time ) THEN
             WRITE(*,'("   ",A,":   ")',ADVANCE="NO") trim(name(ivar))
             first_time = .FALSE.
           END IF
           skip_to_cressman = .TRUE.
           WRITE(*,'(f5.0,",")',ADVANCE="NO") pressure(kp)
           test_count = test_count + 1
           IF ( test_count == 15 ) THEN
             test_count = 0
             WRITE(*,*)
             WRITE(*,'("         ")',ADVANCE="NO")
           END IF
         END IF
       END DO 
       WRITE(*,*)
     END DO

     IF ( skip_to_cressman .AND. oa_3D_option == 0 .AND. oa_min_switch ) THEN
       write (*,*)
       write (*,'("###########################################################################")')
       write (*,'("  Too few observations have been found to do MQD scheme for all levels      ")')
       write (*,'("  Your options are:                                                        ")')
       write (*,'("    Set oa_3D_option = 1 (if any upper-air level have too few observations,")')
       write (*,'("        revert to Cressman for all upper-air levels for that time)         ")')
       write (*,'("    Set oa_3D_option = 2 (revert to Cressman only for the levels which     ")')
       write (*,'("        have too few observations)                                         ")')
       write (*,'("    Set oa_3D_type = Cressman (revert to Cressman for all levels and times ")')
       write (*,'("        regardless of the number of available observations)                ")')
       write (*,'("    Decrease the min number of required observations (mqd_minimum_num_obs),")')
       write (*,'("        NOTE this may degrade your analysis if set too low                 ")')
       write (*,'("###########################################################################")')
       write (*,*)
       STOP
     END IF

   END IF


   ! BPR BEGIN
   ! Calculate the 3D QV field of the first guess field
   ! To calculate the surface QV, this calculates a surface pressure based on the
   ! first-guess file which includes the terrain field for the target model
   ! configuration.  This means the surface pressure (psfc_first_guess_model_terrain_*) 
   ! is not merely that of the coarse grid model, but a potentially much more
   ! finescale product depending on the difference in horizontal resolution
   ! between the coarse grid model and the target model.

   ! Convert from hPa (mb) to Pa because that is what mixing_ratio expects
   slp_x_pa = slp_x * 100.0
   ! Initialize 3D QV and surface pressure to zero
   qv3d_first_guess_yx=0.0
   psfc_first_guess_model_terrain_yx=0.0
   CALL mixing_ratio ( rh , t , h , terrain , slp_x_pa , pressure , iew_alloc , jns_alloc , &
    kbu_alloc , qv3d_first_guess_yx , psfc_first_guess_model_terrain_yx)
   !Change from (y,x) to (x,y)
   CALL yx2xy( psfc_first_guess_model_terrain_yx, jns_alloc, iew_alloc, psfc_first_guess_model_terrain_xy)
   ! BPR END

   !  Loop through all analysis levels (remember that level kp=1 is the
   !  surface value of the 3-D field).


   vertical_level_loop : DO kp = 1 , kbu_alloc
   
      oa_type_tmp = oa_type
      IF ( kp > 1 ) oa_type_tmp = oa_3D_type
      IF ( kp > 1 .AND. skip_to_cressman .AND. oa_3D_option == 1 ) oa_type_tmp = 'Cressman'

      !  If we are doing surface FDDA, and this is one of the surface FDDA time
      !  periods, we can pop out of this loop once we are not doing the surface level.

      IF ( ( fdda_loop .GT. 1 ) .AND. ( kp .GT. 1 ) ) THEN
         EXIT vertical_level_loop
      END IF
      
      !  We need the u and v components of the horizontal wind for the Cressman
      !  scheme when we are doing a banana-type analysis.  The
      !  data going into this routine are still on the (y,x) orientation, but
      !  the 2-D fields on output (u_banana and v_banana) are oriented (x,y).
   
      !BPR BEGIN
      !CALL get_background_for_oa ( t , u , v , rh , slp_x , &
      CALL get_background_for_oa ( t , u , v , rh , slp_x , pres , &
      !BPR END
      u_banana , kp , name(1) , &
      iew_alloc , jns_alloc , kbu_alloc )
   
      !BPR BEGIN
      !CALL get_background_for_oa ( t , u , v , rh , slp_x , &
      CALL get_background_for_oa ( t , u , v , rh , slp_x , pres , &
      !BPR END
      v_banana , kp , name(2) , &
      iew_alloc , jns_alloc , kbu_alloc )

      !  Loop over each variable.
   
      variable_loop : DO ivar = 1 , num_var

         !BPR BEGIN
         !If user chose not to do objective analysis of surface pressure then go
         !to the next variable
         IF ( (.NOT. oa_psfc) .AND.( name(ivar)(1:8) .EQ. 'PSFC    ' ) )  THEN
          CYCLE variable_loop 
         ENDIF
         !BPR END

         !  For sea level pressure, the analysis is only at one level, so cycle the loop
         !  BPR BEGIN
         !! IGNORE THIS COMMENT AS IT IS NO LONGER TRUE -> (which should end this variable loop, but that is not important).
         !  Same for surface pressure
         !IF ( ( name(ivar)(1:8) .EQ. 'PMSL    ' ) .AND. ( kp .GT. 1 ) ) THEN
         IF ( ( ( name(ivar)(1:8) .EQ. 'PSFC    ' ) .OR. &
                ( name(ivar)(1:8) .EQ. 'PMSL    ' ) ) .AND. ( kp .GT. 1 ) ) THEN
         !  BPR END
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
         !BPR BEGIN
         !pressure(kp) , date , time , 1 , &
         pressure(kp) , date , time , request_p_diff , &
         !BPR END
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
               !WRITE ( UNIT =74 , FMT = '( 3x,a8,3x,i6,3x,i5,3x,a8,3x,g13.6,3x,16x,2(f7.2,3x),i7 )' ) &
               WRITE ( UNIT =74 , FMT = '( 3x,a8,3x,i6,3x,i5,3x,a8,3x,2(g13.6,3x),2(f7.2,3x),i7 )' ) &
                 name(ivar) , NINT ( pressure(kp) ) , num , station_id(num) , &
                  obs_value(num) , diff(num) , xob(num) , yob(num) , qc_flag(num)
                 !obs_value(num) , xob(num) , yob(num) , qc_flag(num)
            END DO station_loop_74
         END IF


         !  Obtain observations for kp level and for variable ivar for this
         !  time period.

         CALL proc_ob_access ( 'use', name(ivar) , print_found_obs , &
         !BPR BEGIN 
         !pressure(kp) , date , time , 1 , &
         pressure(kp) , date , time , request_p_diff , &
         !BPR END
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
                              ( oa_type_tmp  .NE. 'MQD'               ) ) THEN

            !  Obtain first guess for kp level and for variable name(ivar).  The
            !  data going into this routine are still on the (y,x) orientation, but
            !  the 2-D field on output (dum2d) is oriented (x,y).
   
            !BPR BEGIN
            !CALL get_background_for_oa ( t , u , v , rh , slp_x , &
            CALL get_background_for_oa ( t , u , v , rh , slp_x , pres , &
            !BPR END
            dum2d , kp , name(ivar) , &
            iew_alloc , jns_alloc , kbu_alloc )

            !BPR BEGIN
            !Save a copy of the first guess field before it is changed by the
            !analysis
            !Note that although dum2d is dimensioned iew_alloc,jns_alloc, for
            !temperature and qv the values at (:,jns_alloc) and (iew_alloc,:)
            !are meaningless since those points are just there to allow u/v
            !points to fit
            dum2d_fg = dum2d(1:iew_alloc-1,1:jns_alloc-1)
            !BPR END
   
            !  Calculate the difference between the observations and first guess fields.
            !  Estimate the values of the first guess at each of the observation location.  
            !  This is performed with the four surrounding analysis points, using a linear
            !  "x" and a linear "y" weighting based on distance.  We are either using the
            !  first-guess to create perturbations, or we can do an analysis from the
            !  observations only.

!BPR BEGIN
            fg_value(:) = -9999999.0
!BPR END
   
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
!BPR BEGIN
                  fg_value(num) = aob
!BPR END
         IF ( ( name(ivar)(1:8) .EQ. 'PMSL    ' ) ) THEN
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

            !IF ( print_obs_files ) THEN
            !   WRITE ( UNIT = 4 , FMT = '( " Number of Observations " / i8.8 )' ) num_obs_pass
            !   WRITE ( UNIT = 4 , FMT = '( "   Variable   Press    Obs     Station   &
            !   &       Obs           Obs-1st      X          Y        QC    " )' )
            !   WRITE ( UNIT = 4 , FMT = '( "    Name      Level    Number    ID      &
            !   &      Value           Guess     Location  Location  Value   " )' )
            !   station_loop_4 : DO num = 1, num_obs_pass
            !      WRITE ( UNIT = 4 , FMT = '( 3x,a8,3x,i6,3x,i5,3x,a8,3x,2(g13.6,3x),2(f7.2,3x),i7 )' ) &
            !      name(ivar) , NINT ( pressure(kp) ) , num , station_id(num) , &
            !      obs_value(num) , diff(num) , xob(num) , yob(num) , qc_flag(num)
            !   END DO station_loop_4
            !END IF

            !BPR BEGIN
            !Ensure that oa_type stays as None if that is what the user chose
            IF ( oa_type .EQ. 'None' ) THEN
             oa_type_tmp = 'None'
            ENDIF
            !BPR END
   
            IF      ( ( oa_type_tmp  .EQ. 'MQD'               ) .AND.  &
                      ( pressure(kp) .EQ. 1001.               ) .AND.  &
                      ( num_obs_pass .GE. mqd_minimum_num_obs ) .AND.  &
                      ( num_obs_pass .LE. mqd_maximum_num_obs ) ) THEN
                if ( ivar == 1 ) print*," "
                if ( ivar == 1 ) print*," Doing MQD for surface data"
   
            !  This is the multiquadric solver.
               CALL mqd ( diff , xob , yob , num_obs_pass , &
               dum2d , iew_alloc , jns_alloc , &
               crsdot(ivar) , name(ivar) , passes , smooth_type , use_first_guess )
 
            ELSE IF ( ( oa_type_tmp  .EQ. 'MQD'               ) .AND.  &
                      ( pressure(kp) .LE. 1000.               ) .AND.  &
                      ( num_obs_pass .GE. mqd_minimum_num_obs ) .AND.  &
                      ( num_obs_pass .LE. mqd_maximum_num_obs ) ) THEN
                min_mqd = min(num_obs_pass, min_mqd)
                max_mqd = max(num_obs_pass, max_mqd)
                if ( ivar == 1 ) num_mqd_uu = num_mqd_uu + 1
                if ( ivar == 2 ) num_mqd_vv = num_mqd_vv + 1
                if ( ivar == 3 ) num_mqd_tt = num_mqd_tt + 1
                if ( ivar == 4 ) num_mqd_rh = num_mqd_rh + 1
   
            !  This is the multiquadric solver.
               CALL mqd ( diff , xob , yob , num_obs_pass , &
               dum2d , iew_alloc , jns_alloc , &
               crsdot(ivar) , name(ivar) , passes , smooth_type , use_first_guess )
 
             
            !BPR BEGIN 
            !Code to ensure oa_type stays as None if that is what the user chose
            !ELSE IF ( ( oa_type_tmp .EQ. 'Cressman' ) .OR. &
            !        ( ( oa_min_switch ) .AND. ( num_obs_pass .LT. mqd_minimum_num_obs ) ) .OR. &
            !        ( ( oa_max_switch ) .AND. ( num_obs_pass .GT. mqd_maximum_num_obs ) ) ) THEN
            ELSE IF ( ( ( oa_type_tmp .EQ. 'Cressman' ) .OR. &
                      ( ( oa_min_switch ) .AND. ( num_obs_pass .LT. mqd_minimum_num_obs ) ) .OR. &
                      ( ( oa_max_switch ) .AND. ( num_obs_pass .GT. mqd_maximum_num_obs ) ) ) .AND. &
                      ( oa_type .NE. 'None' ) ) THEN
            !BPR END

                  if ( ivar == 1 .and. pressure(kp) .eq. 1001 ) print*," "
                  if ( ivar == 1 .and. pressure(kp) .eq. 1001 ) print*," Doing Cressman for surface data"

                  if ( pressure(kp) .LE. 1000 ) min_mqd = min(num_obs_pass, min_mqd)
                  if ( pressure(kp) .LE. 1000 ) max_mqd = max(num_obs_pass, max_mqd)

               multi_scan : DO num_scan = 1 , max_scan
   
                  IF ( radius_influence(num_scan) .LT. 0 ) THEN
                     EXIT multi_scan
                  END IF


                  !  This is the Cressman method.  When multiple scans are utilized,
                  !  the scheme is called repeatedly, and the "final" analysis is 
                  !  used as input for the new difference field that is objectively
                  !  analyzed.
                  
                  !BPR BEGIN
                  !IF (( name(ivar)(1:8) .EQ. 'RH      ' ).AND.(kp.eq.9)) THEN
                  ! DO cury=1,jns_alloc
                  !  WRITE(110+num_scan,*) dum2d(:,cury)
                  ! ENDDO
                  !ENDIF
                  !BPR END

                  !BPR BEGIN
                  !If we are dealing with the surface then use a smaller radius of influence
                  IF(abs(pressure(kp)-1001.0).lt.0.0001) THEN
                   radius_influence_scale_factor = radius_influence_sfc_mult
                  ELSE
                   radius_influence_scale_factor = 1.00
                  ENDIF
                  radius_influence_scaled = radius_influence(num_scan)*radius_influence_scale_factor
                  radius_influence_scaled_scan1 = radius_influence(1)*radius_influence_scale_factor

                  !IF the user chose to use reduced radii of influence for
                  !surface obs compared to other obs
                  IF(radius_influence_sfc_mult.lt.0.99999) THEN
                   !Ensure that for the surface the radius of influence in terms of model grid
                   !cells stays between radius_influence_sfc_min_cells and radius_influence_sfc_max_cells
                   radius_influence_sfc_max_cells = 100.0;
                   radius_influence_sfc_min_cells = 4.5;
                   IF(abs(pressure(kp)-1001.0).lt.0.0001) THEN
 
                    !Ensure that the radius of influence on the first scan
                    !does not exceed a certain number of grid points. If it does, then adjust the radius
                    !on each scan by the amount necessary so that the first scan has a radius equal to 
                    !the maximum allowed number of grid points.  
                    !This is to avoid everything being washed out on finer
                    !resolution domains even when there is dense data
                    IF(radius_influence_scaled_scan1.gt.radius_influence_sfc_max_cells) THEN
                     radius_influence_scaled = radius_influence_scaled * &
                      ( radius_influence_sfc_max_cells / radius_influence_scaled_scan1 ) 
                    ENDIF
                    !Ensure that the radius of influence on each scan does not go
                    !below a certain number of grid points.  If it does, then
                    !skip the remaining scans.  This is to minimize "spots" in the
                    !data.
                    IF(radius_influence_scaled.lt.radius_influence_sfc_min_cells) THEN
                     IF(num_scan.eq.1) THEN
                      PRINT *,'ERROR: The effective radius of influence for the first Cressman scan ',&
                              'incorporating surface obs into the analysis was too few grid points so ',&
                              'ZERO Cressman scans would have been used and thus the surface data would ',&
                              'not be incorporated into the analysis.  The radius of influence for the first',&
                              ' Cressman scan was set to ',radius_influence_scaled,' model grid points, but ',&
                              'the minimum value is ',radius_influence_sfc_min_cells,'.'
                      STOP 'ERROR: Effective radius of influence for the first Cressman scan of surface obs too small to continue.'
                     ENDIF
                     EXIT multi_scan   
                    ENDIF
                   ENDIF
                  ENDIF
                  IF( roi_not_printed_yet(kp,num_scan) ) THEN
                   IF(abs(pressure(kp)-1001.0).lt.0.0001) THEN
!BPR BEGIN
!                   WRITE(*,'(A,A,I3,A,F6.1,A,F7.2,A,F6.1,A)'), 'ROI for  the surface,',&
                    WRITE(*,'(A,A,I3,A,F6.1,A,F7.2,A,F6.1,A)')  'ROI for  the surface,',&
!BPR END
                           ' Cressman scan number ',num_scan,' is ',radius_influence_scaled,&
                           ' grid points or ',radius_influence_scaled*dxd,'km (',dxd,' km spacing).'
                   ELSE
!BPR BEGIN
!                   WRITE(*,'(A,F8.2,A,I3,A,F6.1,A,F7.2,A,F6.1,A)'), 'ROI for ',pressure(kp),&
                    WRITE(*,'(A,F8.2,A,I3,A,F6.1,A,F7.2,A,F6.1,A)')  'ROI for ',pressure(kp),&
!BPR END
                           ' hPa, Cressman scan number ',num_scan,' is ',radius_influence_scaled,&
                           ' grid points or ',radius_influence_scaled*dxd,'km (',dxd,' km spacing).'
                   ENDIF
                   roi_not_printed_yet(kp,num_scan) = .FALSE.
                  ENDIF
                  !BPR END

                  CALL cressman ( obs_value , diff , xob , yob , num_obs_pass , &
                  dum2d , iew_alloc , jns_alloc , &
!BPR BEGIN
                  crsdot(ivar) , name(ivar) , nint(radius_influence_scaled) , &
!                 crsdot(ivar) , name(ivar) , radius_influence(num_scan) , &
!BPR END
                  dxd , u_banana , v_banana , pressure(kp) , passes , smooth_type , lat_center , &
!BPR BEGIN
                  use_first_guess, fg_value, scale_cressman_rh_decreases )
!                 use_first_guess )
!BPR END

                  !BPR BEGIN
                  !IF (( name(ivar)(1:8) .EQ. 'RH      ' ).AND.(kp.eq.9)) THEN
                  ! DO cury=1,jns_alloc
                  !  WRITE(130+num_scan,*) dum2d(:,cury)
                  ! ENDDO
                  !ENDIF
                  !BPR END

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

            !BPR BEGIN
            !If we are processing temperature, then save temperature for use in
            !calculating mixing ratio once we start processing RH
            IF ( name(ivar)(1:8) .EQ. 'TT      ' )  THEN
             t_first_guess = dum2d_fg
             t_analysis = dum2d(1:iew_alloc-1,1:jns_alloc-1)
             pressure_t = pressure(kp)
            ENDIF
            ! BPR END

   
            !  Store the final analysis for kp level and for variable name(ivar).
   
            !BPR BEGIN
            !CALL put_background_from_oa ( t , u , v , rh , slp_x , &
            CALL put_background_from_oa ( t , u , v , rh , slp_x , pres , &
            !BPR END
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


   IF ( ( oa_type .EQ. 'MQD' ) ) THEN
      mqd_abs_min = min(min_mqd, mqd_abs_min)
      IF ( (num_mqd_uu > 0) .OR. (num_mqd_vv > 0) .OR. (num_mqd_tt > 0) .OR. (num_mqd_rh > 0) ) THEN
        IF ( (num_mqd_uu < kbu_alloc-1) .OR. (num_mqd_vv < kbu_alloc-1) .OR. &
             (num_mqd_tt < kbu_alloc-1) .OR. (num_mqd_rh < kbu_alloc-1) ) THEN
        ! Find the value that will switch off MQD in upper levels
        mqd_count = max(max_mqd, mqd_count)
        WRITE(*,'("  MQD was performed on ", i3, " of ", i3, " upper levels for variable UU" )' ) num_mqd_uu, kbu_alloc-1
        WRITE(*,'("  MQD was performed on ", i3, " of ", i3, " upper levels for variable VV" )' ) num_mqd_vv, kbu_alloc-1
        WRITE(*,'("  MQD was performed on ", i3, " of ", i3, " upper levels for variable TT" )' ) num_mqd_tt, kbu_alloc-1
        WRITE(*,'("  MQD was performed on ", i3, " of ", i3, " upper levels for variable RH" )' ) num_mqd_rh, kbu_alloc-1
        ELSE
          print*," Doing MQD for all upper level data"
        END IF
      ELSE
        IF (kp .eq. kbu_alloc+1 ) print*," Doing Cressman for all upper level data"
      END IF
   !BPR BEGIN
   ELSEIF ( ( oa_type .EQ. 'None' ) ) THEN
    !min_mqd is not actually specific to MQD but is used to store the minimum
    !number of obs when using any method
    !min_mqd is initialized to 100, but if we are doing no analysis it is never
    !updated so it will erroneously print out that there was a minimum of 100
    !obs on each level processed.
    min_mqd = 0
   !BPR END
   ELSE
      IF (kp .eq. kbu_alloc+1 ) print*," Doing Cressman for all upper level data"
   END IF
   print*,"  "
   IF (kp .eq. kbu_alloc+1 ) WRITE(*, '( A, i3, A )' ) &
     "  MIN number of observation on any upper level for this time period is : ", &
     min_mqd, " observations per variable"
   IF (kp .eq. kbu_alloc+1 ) WRITE(*, '( A, i3, A)' ) &
     "  MAX number of observation on any upper level for this time period is : ", &
     max_mqd, " observations per variable"

  
   !  We need to close our sample print-out file.

   IF ( print_obs_files ) THEN
      !CLOSE ( 4 ) 
      CLOSE ( 74 ) 
   END IF

END SUBROUTINE proc_oa
