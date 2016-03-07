!------------------------------------------------------------------------------

!  This module inputs and stores the NAMELIST data to the program.  The
!  simple consistency checks available using only the NAMELIST data
!  are also performed.  

MODULE namelist

   TYPE nml_record_1
      CHARACTER ( LEN = 19 )  :: start_date                   , & ! initial time to generate an
                                                                  ! analysis, process obs
                                 end_date                         ! final time to generate an
                                                                  ! analysis, process obs,
                                                                  ! both times are YYYY-MM-DD_HH:mm:ss
      INTEGER                 :: interval                         ! time (s) between processing
                                                                  ! points
   END TYPE nml_record_1

   INTEGER , PARAMETER :: max_times = 1000

   TYPE nml_record_2
      CHARACTER ( LEN = 132)  :: fg_filename                      ! first-guess filename 
      CHARACTER ( LEN = 132)  :: obs_filename                     ! observation filename 
      LOGICAL :: trim_domain, remove_unverified_data
      INTEGER :: trim_value, grid_id, remove_data_above_qc_flag
   END TYPE nml_record_2

   TYPE nml_record_3
      INTEGER                 :: max_number_of_obs                ! number larger than total number 
                                                                  ! observation reports expected
  
      LOGICAL                 :: fatal_if_exceed_max_obs          ! stop program if the maximum number of 
                                                                  ! observations is exceeded
   END TYPE nml_record_3

   TYPE nml_record_4
      LOGICAL                 :: qc_test_error_max            , & ! perform error-max check
                                 qc_test_buddy                , & ! perform buddy check
                                 qc_test_vert_consistency     , & ! T, u, v vertical consistency
!BPR BEGIN
!                                qc_test_convective_adj           ! lapse rate tests
                                 qc_test_convective_adj       , & ! lapse rate tests
                                 qc_psfc                          ! quality control surface pressure
!BPR END

      REAL                    :: max_error_t                  , &
                                 max_error_uv                 , &
                                 max_error_z                  , &
                                 max_error_rh                 , &
!BPR BEGIN
                                 max_error_dewpoint           , &
!BPR END
                                 max_error_p                  , &
                                 max_buddy_t                  , &
                                 max_buddy_uv                 , &
                                 max_buddy_z                  , &
                                 max_buddy_rh                 , &
!BPR BEGIN
                                 max_buddy_dewpoint           , &
!BPR END
                                 max_buddy_p                  , & 
                                 buddy_weight                 , &  
                                 max_p_extend_t               , &
                                 max_p_extend_w    
!BPR BEGIN
     LOGICAL                  :: use_p_tolerance_one_lev          ! T/F Allow single level obs vertically placed within
                                                                  ! a pressure tolerance of the pressure level for which 
                                                                  ! obs are being sought to be used
     INTEGER                  :: max_p_tolerance_one_lev_qc       ! Pressure tolerance between ob and the first guess 
                                                                  ! pressure level that can be used to QC it (Pa)
                                                                  ! use_p_tolerance_one_lev must be set to true to use this
!BPR END
   END TYPE nml_record_4

   TYPE nml_record_5
      LOGICAL                 :: print_found_obs              , & ! T/F print each observation found
                                 print_header                 , & ! T/F print header information
                                 print_analysis               , & ! T/F print intitial and final analysis means
                                 print_qc_vert                , & ! T/F print QC messages in vertical consistency
                                 print_qc_dry                 , & ! T/F print QC messages from convective adjustment
                                 print_obs_files              , & ! T/F generate additional obs files
                                 print_error_max              , & ! T/F print message if fails error max
                                 print_buddy                  , & ! T/F print message if fails buddy check
                                 print_oa                         ! T/F print obs used in OA

   END TYPE nml_record_5

   TYPE nml_record_7
      LOGICAL                 :: use_first_guess              , & ! T/F add objective analysis perturbations
                                                                  ! to the first-guess fields
                                 f4d                          , & ! T/F processing for sfc FDDA
                                 lagtem                           ! T/F use lagged time for off-time first guess

      INTEGER                 :: intf4d                           ! time (s) between sfc FDDA time periods
   END TYPE nml_record_7

   TYPE nml_record_8
      INTEGER                 :: smooth_type                  , & ! 1 = 5 point stencil; 2 = smoother desmoother
                                 smooth_sfc_wind              , & ! smoothing passes for surface wind
                                 smooth_sfc_temp              , & ! smoothing passes for surface temperature
                                 smooth_sfc_rh                , & ! smoothing passes for surface relative humidity
                                 smooth_sfc_slp               , & ! smoothing passes for sea level pressure
                                 smooth_upper_wind            , & ! smoothing passes for upper level winds
                                 smooth_upper_temp            , & ! smoothing passes for upper  temperature
                                 smooth_upper_rh                  ! smoothing passes for upper relative humidity
 
   END TYPE nml_record_8

   TYPE nml_record_9
      CHARACTER ( LEN = 132 ) :: oa_type                          ! MQD or Cressman
      CHARACTER ( LEN = 132 ) :: oa_3D_type                       ! MQD or Cressman
      INTEGER                 :: mqd_minimum_num_obs          , & ! minimum number of obs to use MQD, revert to Cressman
                                 mqd_maximum_num_obs              ! maximum number of obs to use MQD, revert to Cressman
      INTEGER , DIMENSION(10) :: radius_influence                 ! radius of influence in grid point units for Cressman
      LOGICAL                 :: oa_min_switch                , & ! if MQD has less than number of obs, T/F use Cressman
                                 oa_max_switch                    ! if MQD has more than number of obs, T/F use Cressman
      INTEGER                 :: oa_3D_option  
!BPR BEGIN
      LOGICAL                 :: oa_psfc                          ! Should surface pressure be objectively analyzed
      LOGICAL                 :: scale_cressman_rh_decreases      ! Should Cressman-caused decreases to RH at locations 
                                                                  ! with lower RH than at the ob location be scaled based 
                                                                  ! on the relationship between RH at the ob location and
                                                                  ! at the application location
      REAL                    :: radius_influence_sfc_mult        ! To get the radius of influence used for surface obs 
                                                                  ! what should the radius of influence used by non-surface
                                                                  ! obs be multiplied by  
      INTEGER                 :: max_p_tolerance_one_lev_oa       ! Maximum pressure difference allowed between a single 
                                                                  ! level ob and the pressure level it is analyzed on
                                                                  ! User must enable use_p_tolerance_one_lev for this to
                                                                  ! have an effect and must be less than max_p_tolerance_one_lev_qc
                                                                  
!BPR END
   END TYPE nml_record_9

   TYPE all_nml
      TYPE ( nml_record_1)   :: record_1
      TYPE ( nml_record_2)   :: record_2  
      TYPE ( nml_record_3)   :: record_3
      TYPE ( nml_record_4)   :: record_4
      TYPE ( nml_record_5)   :: record_5
      TYPE ( nml_record_7)   :: record_7
      TYPE ( nml_record_8)   :: record_8
      TYPE ( nml_record_9)   :: record_9
   END TYPE all_nml

CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE check_namelist ( nml ) 
 
!  This routine performs simple error and consistency checks on the NAMELIST
!  input.  Reports are made through the error handling routines.

   IMPLICIT NONE

   INCLUDE 'error.inc'

   TYPE ( all_nml ) :: nml

   INTEGER :: loop_test
   LOGICAL :: exist

   INTERFACE
      INCLUDE 'error.int'
   END INTERFACE

   !  001: Check that bounding dates make chronological sense.

   test_101 : IF ( nml%record_1%start_date .GT. nml%record_1%end_date ) THEN
      error_message = ' '
      error_number = 00013101
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' The starting date is greater than the &
      &ending date in NAMELIST record 1'
      fatal = .true.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_101

   !  002: Check that bounding dates are with the 1900s or 2000s.

   test_102 : IF (( nml%record_1%start_date(1:2) .LT. '19' ) .OR. &
                  ( nml%record_1%start_date(1:2) .GT. '20' ) .OR. &
                  ( nml%record_1%end_date(1:2)   .LT. '19' ) .OR. &
                  ( nml%record_1%end_date(1:2)   .GT. '20' )) THEN
      error_message = ' '
      error_number = 00013102
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' The start and/or ending dates have a &
      &wrong century (' // nml%record_1%start_date(1:2) // ' and ' // &
      nml%record_1%end_date(1:2) // '.'
      fatal = .true.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_102

   !  003: Interval to process has to be .GE. 0 if the
   !  bounding times are not the same.

   test_103 : IF (( nml%record_1%start_date .NE.  &
                    nml%record_1%end_date ) .AND.  &
                  ( nml%record_1%interval .LE. 0 )) THEN
      error_message = ' '
      error_number = 00013103
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' The interval is .LE. 0 for processing'
      fatal = .true.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_103
 

   !  202: Does the observation data file supplied from the NAMELIST file
   !  exist?

   !test_202 : DO loop_test = 1 , max_times
   !   IF ( nml%record_2%obs_filename(loop_test)(1:4) .NE. 'null'     ) THEN
   !      INQUIRE ( EXIST = exist , FILE = nml%record_2%obs_filename(loop_test) ) 
   !      IF ( .NOT. exist ) THEN
   !         error_message = ' '
   !         error_number = 00013202
   !         error_message(1:31) = 'check_namelist                 '
   !         error_message(32:)  = ' The observation data filename supplied &
   !         &in the namelist (' // TRIM ( nml%record_2%obs_filename(loop_test) ) // &
   !         ') does not exist.'
   !         fatal = .true.
   !         listing = .false.
   !         call error_handler ( error_number , error_message , fatal , listing )
   !      END IF
   !   ELSE IF ( nml%record_2%obs_filename(loop_test)(1:4) .EQ. 'null'     ) THEN
   !      EXIT test_202
   !   END IF
   !END DO test_202
 
   !  203: Does the sfc FDDA observation data file supplied from the NAMELIST file
   !  exist?

   !IF ( nml%record_7%f4d ) THEN
   !   test_203 : DO loop_test = 1 , max_times
   !      IF ( nml%record_2%sfc_obs_filename(loop_test)(1:4) .NE. 'null'     ) THEN
   !         INQUIRE ( EXIST = exist , FILE = nml%record_2%sfc_obs_filename(loop_test) ) 
   !         IF ( .NOT. exist ) THEN
   !            error_message = ' '
   !            error_number = 00013203
   !            error_message(1:31) = 'check_namelist                 '
   !            error_message(32:)  = ' The sfc FDDA observation data filename supplied &
   !            &in the namelist (' // TRIM ( nml%record_2%sfc_obs_filename(loop_test) ) // &
   !            ') does not exist.'
   !            fatal = .true.
   !            listing = .false.
   !            call error_handler ( error_number , error_message , fatal , listing )
   !         END IF
   !      ELSE IF ( nml%record_2%sfc_obs_filename(loop_test)(1:4) .EQ. 'null'     ) THEN
   !         EXIT test_203
   !      END IF
   !   END DO test_203
   !END IF

   !  401:  Perform error max QC test?

   test_401 : IF ( .NOT. nml%record_4%qc_test_error_max ) THEN
      error_message = ' '
      error_number = 00013401
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' No error max check is being performed.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_401

   !  402:  Perform buddy check QC test?

   test_402 : IF ( .NOT. nml%record_4%qc_test_buddy ) THEN
      nml%record_4%buddy_weight = 0
      error_message = ' '
      error_number = 00013402
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' No buddy check is being performed.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_402

   !  403:  Perform buddy check QC test?

   test_403 : IF ( nml%record_4%buddy_weight .LT. 0.01 ) THEN
      nml%record_4%qc_test_buddy = .FALSE. 
      error_message = ' '
      error_number = 00013403
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' No buddy check is being performed.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_403

   !  404:  Perform vertical consistency check QC test?

   test_404 : IF ( .NOT. nml%record_4%qc_test_vert_consistency ) THEN
      error_message = ' '
      error_number = 00013404
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' No vertical consistency check is being performed.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_404

   !  405:  Perform convective adjustment check QC test?

   test_405 : IF ( .NOT. nml%record_4%qc_test_convective_adj ) THEN
      error_message = ' '
      error_number = 00013405
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' No convective adjustment check is being performed.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_405

   !BPR BEGIN
   !  406:  Perform quality control on surface pressure?
   !  Tell user whether or not they are quality controlling surface pressure
   !  If they are not, but they are objectively analyzing that field, then throw
   !  a fatal error given the danger of objectively analyzing non-QC'd obs

   error_message = ' '
   error_number = 00013406
   error_message(1:31) = 'check_namelist                 '
   listing = .false.

   test_406 : IF ( .NOT. nml%record_4%qc_psfc ) THEN
      IF ( nml%record_9%oa_psfc ) THEN
       fatal = .true.
       error_message(32:)  = ' Surface pressure is NOT being QCed but is being objectively analyzed.'
      ELSE
       fatal = .false.
       error_message(32:)  = ' Quality control checks will NOT be performed on surface pressure.'
      END IF
   ELSE
      fatal = .false.
      error_message(32:)  = ' Quality control checks WILL be performed on surface pressure.'
   END IF test_406
   call error_handler ( error_number , error_message , fatal , listing )

   !  407:  Allow pressure tolerance when doing QC of single level above surface obs?

   error_message = ' '
   error_number = 00013407
   error_message(1:31) = 'check_namelist                 '
   listing = .false.
   fatal = .false.

   test_407 : IF ( nml%record_4%use_p_tolerance_one_lev ) THEN
      WRITE(error_message(32:),'(A,I9,A)') ' Single-level above-surface obs will be QCed &
       &against the nearest pressure level in the analysis as long as the pressure difference &
       &is no more than',& 
       nml%record_4%max_p_tolerance_one_lev_qc, &
       ' Pa.  No obs will be extended to the nearest pressure level.'
   ELSE
      error_message(32:) = ' Single-level above-surface obs will generally only be QCed if they fall&
       & on an analysis pressure level. However, will attempt to adjust obs marked as FM-88 SATOB or FM-97 AIREP&
       & to the nearest pressure level.'
   END IF test_407
   call error_handler ( error_number , error_message , fatal , listing )

   !  408:  Make user aware that a pressure tolerance of 1 hPa for using a single-level above-surface
   !        ob on the objective analysis at a pressure level, is like having no
   !        tolerance at all.

   test_408 : IF ( ( nml%record_4%use_p_tolerance_one_lev ) .AND. &
                   ( nml%record_4%max_p_tolerance_one_lev_qc .EQ. 1 ) ) THEN
      error_message = ' '
      error_number = 00013408
      error_message(1:31) = 'check_namelist                 '
      listing = .false.
      fatal = .false.

      WRITE(error_message(32:),'(A)') ' Pressure tolerance for single-level above-surface obs being QCed &
       &against an analysis is 1 hPa.  This is equivalent to having no pressure tolerance!'

      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_408


   !BPR END

   !  501:  Print out message for each observation found?

   test_501 : IF ( .NOT. nml%record_5%print_found_obs ) THEN
      error_message = ' '
      error_number = 00013501
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' No message is printed when each observation is found.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_501

   !  502:  Print out header information for initial and final analysis?

   test_502 : IF ( .NOT. nml%record_5%print_header ) THEN
      error_message = ' '
      error_number = 00013502
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' No header information is printed for initial or final analysis.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_502

   !  503:  Print out mean/max/min for initial and final analysis?

   test_503 : IF ( .NOT. nml%record_5%print_analysis ) THEN
      error_message = ' '
      error_number = 00013503
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' No 2D mean/max/min is printed for initial or final analysis.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_503

   !  504:  Print out QC information from vertical consistency check?

   test_504 : IF ( ( .NOT. nml%record_5%print_qc_vert ) .AND. &
                   ( nml%record_4%qc_test_vert_consistency ) ) THEN
      error_message = ' '
      error_number = 00013504
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' No QC information is printed from the vertical consistency check.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_504

   !  505:  Print out QC information from the dry convective adjustment?

   test_505 : IF ( ( .NOT. nml%record_5%print_qc_dry ) .AND. & 
                   ( nml%record_4%qc_test_convective_adj ) ) THEN
      error_message = ' '
      error_number = 00013505
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' No header information is printed for initial or final analysis.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_505

   !  506:  Produce additional files for the observations?

   test_506 : IF ( .NOT. nml%record_5%print_obs_files ) THEN
      error_message = ' '
      error_number = 00013506
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' No additional files for the observations will be produced.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_506

   !  507:  Produce print out when observation fails error max check?

   test_507 : IF ( ( .NOT. nml%record_5%print_error_max ) .AND. &
                   ( nml%record_4%qc_test_error_max ) ) THEN
      error_message = ' '
      error_number = 00013507
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' No print out will be generated when observation fails error max check.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_507

   !  508:  Produce print out when observation fails buddy check?

   test_508 : IF ( ( .NOT. nml%record_5%print_buddy ) .AND. & 
                   ( nml%record_4%qc_test_buddy ) ) THEN
      error_message = ' '
      error_number = 00013508
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' No print out will be generated when observation fails buddy check.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_508

   !  509:  Produce print out of every observation used in OA?

   test_509 : IF ( .NOT. nml%record_5%print_oa ) THEN
      error_message = ' '
      error_number = 00013509
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' No print out will be generated for observations used in OA.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   ELSE
      nml%record_5%print_found_obs = .TRUE.
   END IF test_509

   !  701:  The user may select to have the first-guess fields not used in the 
   !  objective analysis, using only the observations.

   test_701 : IF ( .NOT. nml%record_7%use_first_guess  ) THEN
      error_message = ' '
      error_number = 00013701
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' The final analysis will be from the observations, no first-guess.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_701

   !  702:  The SFC FDDA interval should be less than or equal to the interval specified
   !  for the 3d analysis.

   test_702 : IF ( (  nml%record_7%f4d ) .AND. ( nml%record_7%intf4d .GT. nml%record_1%interval ) ) THEN
      error_message = ' '
      error_number = 00013702
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' The FDDA sfc interval is larger than the 3d interval.'
      fatal = .true.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_702

   !  703:  We want the FDDA frequency in seconds, not hours.

   test_703 : IF ( ( nml%record_7%f4d ) .AND. ( nml%record_7%intf4d .LT. 3600 ) ) THEN
      error_message = ' '
      error_number = 00013703
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' The FDDA sfc interval seems to be in hours, not seconds.'
      fatal = .true.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_703

   !  801:  Which smoother to use on final analysis.

   test_801 : IF ( ( nml%record_8%smooth_type .NE. 1 ) .AND. & 
                   ( nml%record_8%smooth_type .NE. 2 ) ) THEN
      error_message = ' '
      error_number = 00013801
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' The available smoothing options are either (1) 5 point stencil or (2) smoother-desmoother.'
      fatal = .true.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_801

   !  802:  How many smoothing passes on the surface winds?

   test_802 : IF ( nml%record_8%smooth_sfc_wind .LE. 0 ) THEN
      error_message = ' '
      error_number = 00013802
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' There will be no smoothing for the surface winds on the final analysis.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_802

   !  803:  How many smoothing passes on the surface temp?

   test_803 : IF ( nml%record_8%smooth_sfc_temp .LE. 0 ) THEN
      error_message = ' '
      error_number = 00013803
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' There will be no smoothing for the surface temperature on the final analysis.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_803

   !  804:  How many smoothing passes on the surface rh?

   test_804 : IF ( nml%record_8%smooth_sfc_rh   .LE. 0 ) THEN
      error_message = ' '
      error_number = 00013804
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' There will be no smoothing for the surface relative humidity on the final analysis.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_804

   !  805:  How many smoothing passes on the sea level pressure?

   test_805 : IF ( nml%record_8%smooth_sfc_slp  .LE. 0 ) THEN
      error_message = ' '
      error_number = 00013805
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' There will be no smoothing for the sea level pressure on the final analysis.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_805

   !  806:  How many smoothing passes on the upper level winds?

   test_806 : IF ( nml%record_8%smooth_upper_wind .LE. 0 ) THEN
      error_message = ' '
      error_number = 00013806
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' There will be no smoothing for the upper level winds on the final analysis.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_806

   !  807:  How many smoothing passes on the upper level temp?

   test_807 : IF ( nml%record_8%smooth_upper_temp .LE. 0 ) THEN
      error_message = ' '
      error_number = 00013807
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' There will be no smoothing for the upper level temperature on the final analysis.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_807

   !  808:  How many smoothing passes on the upper level rh?

   test_808 : IF ( nml%record_8%smooth_upper_rh   .LE. 0 ) THEN
      error_message = ' '
      error_number = 00013808
      error_message(1:31) = 'check_namelist                 '
      error_message(32:)  = ' There will be no smoothing for the upper level relative humidity on the final analysis.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_808

   !  901:  Which objective analysis technique is to be used?

   !BPR BEGIN
   !Allow "None" option for cases in which we only want to QC data
   !test_901 : IF ( ( nml%record_9%oa_type           .NE. 'MQD'      ) .AND. &
   !                ( nml%record_9%oa_type           .NE. 'Cressman' ) ) THEN
   test_901 : IF ( ( nml%record_9%oa_type           .NE. 'MQD'      ) .AND. &
                   ( nml%record_9%oa_type           .NE. 'Cressman' ) .AND. &
                   ( nml%record_9%oa_type           .NE. 'None' ) ) THEN
   !BPR END
      nml%record_9%oa_type = 'Cressman' 
   END IF test_901

   !  902:  What is the minimum number of observations required to still do MQD?

   test_902 : IF ( ( nml%record_9%oa_type             .EQ. 'MQD'    ) .AND. &
                   ( nml%record_9%mqd_minimum_num_obs .GT. 0        ) ) THEN
      error_message = ' '
      error_number = 00013902
      write (error_message(110:113) , fmt = '(i4)' ) nml%record_9%mqd_minimum_num_obs
      error_message(1:31) = 'check_namelist                 '
      error_message(32:109) = ' The minimum number of observations required for MQD analysis &
      &has been set to '
      error_message(114:) = '.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_902

   !  903:  What is the minimum number of observations required to still do MQD?

   test_903 : IF ( ( nml%record_9%oa_type             .EQ. 'MQD'    ) .AND. &
                   ( nml%record_9%mqd_maximum_num_obs .GT. 0        ) ) THEN
      error_message = ' '
      error_number = 00013903
      write (error_message(111:116) , fmt = '(i6)' ) nml%record_9%mqd_maximum_num_obs
      error_message(1:31) = 'check_namelist                 '
      error_message(32:109) = ' The maximum number of observations permitted for MQD analysis &
      &has been set to '
      error_message(117:) = '.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_903

   !        Are we doing Cressman for the upper-air?

   IF ( ( nml%record_9%oa_3D_type        .NE. 'Cressman' ) ) THEN
          nml%record_9%oa_3D_type = nml%record_9%oa_type            
   END IF
   IF ( ( nml%record_9%oa_type           .EQ. 'MQD'      ) .AND. &
                   ( nml%record_9%oa_3D_type        .EQ. 'Cressman' ) ) THEN
      error_message = ' '
      error_message(32:109) = ' Switching to Cressman for the upper-air analysis'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF 
        IF ( ( nml%record_9%oa_3D_type .EQ. 'MQD' ) .AND. (nml%record_9%oa_3D_option == 0) ) THEN
      error_message = ' '
      error_message(32:112) = ' The code will stop if we do not have enough data to perform MQD on upper levels.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   ELSE IF ( ( nml%record_9%oa_3D_type .EQ. 'MQD' ) .AND. (nml%record_9%oa_3D_option == 1) ) THEN
      error_message = ' '
      error_message(32:110) = ' Revert to Cressman for each time which lack enough data to perform MQD scheme.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   ELSE IF ( ( nml%record_9%oa_3D_type .EQ. 'MQD' ) .AND. (nml%record_9%oa_3D_option == 2) ) THEN
      error_message = ' '
      error_message(32:111) = ' Revert to Cressman for each level which lack enough data to perform MQD scheme.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF


   !  903:  Requested radius of influence for the Cressman scheme.

   test_904 : DO loop_test = 1 , 10
      IF ( ( nml%record_9%radius_influence(loop_test) .GT. 100        ) .AND. &
           ( nml%record_9%oa_type                     .EQ. 'Cressman' ) ) THEN
         error_message = ' '
         error_number = 00013904
         error_message(1:31) = 'check_namelist                 '
         error_message(32:109) = ' The radius_influence variable is in grid units, not km.'
         fatal = .true.
         listing = .false.
         call error_handler ( error_number , error_message , fatal , listing )
      !ELSE IF ( ( nml%record_9%radius_influence(loop_test) .LE. 0          ) .AND. &
      !          ( nml%record_9%oa_type                     .EQ. 'Cressman' ) ) THEN
      !   nml%record_9%radius_influence(loop_test) = -1
      END IF
   END DO test_904

   !  905:  What is the minimum number of observations required to still do MQD?

   test_905 : IF ( ( nml%record_9%oa_type             .EQ. 'MQD'    ) .AND. &
                   ( nml%record_9%mqd_minimum_num_obs .GT. 0        ) .AND. &
                   ( nml%record_9%oa_min_switch                     ) ) THEN
      error_message = ' '
      error_number = 00013905
      error_message(1:31) = 'check_namelist                 '
      error_message(32:109) = ' Will switch from MQD to Cressman scheme when too few observations.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_905

   !  906:  What is the minimum number of observations required to still do MQD?

   test_906 : IF ( ( nml%record_9%oa_type             .EQ. 'MQD'    ) .AND. &
                   ( nml%record_9%mqd_maximum_num_obs .GT. 0        ) .AND. &
                   ( nml%record_9%oa_max_switch                     ) ) THEN
      error_message = ' '
      error_number = 00013906
      error_message(1:31) = 'check_namelist                 '
      error_message(32:109) = ' Will switch from MQD to Cressman scheme when too many observations.'
      fatal = .false.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_906


   !BPR BEGIN
   !  907:  What is the minimum number of observations required to still do MQD?

   error_message = ' '
   error_number = 00013907
   error_message(1:31) = 'check_namelist                 '
   fatal = .false.
   listing = .false.

   test_907 : IF ( nml%record_9%oa_psfc ) THEN
      error_message(32:) = ' Will objectively analyze surface pressure [psfc will NOT be affected by SLP analysis].'
   ELSE
      error_message(32:) = ' Will NOT objectively analyze surface pressure [psfc will be affected by SLP analysis].'
   END IF test_907
   call error_handler ( error_number , error_message , fatal , listing )

   !  908:  Allow pressure tolerance when doing OA of single level above
   !  surface obs?

   error_message = ' '
   error_number = 00013908
   error_message(1:31) = 'check_namelist                 '
   listing = .false.
   fatal = .false.

   test_908 : IF ( nml%record_4%use_p_tolerance_one_lev ) THEN
      WRITE(error_message(32:),'(A,I9,A)') ' Single-level above-surface obs will be included in &
       &the objective analysis on the nearest pressure level as long as the pressure difference &
       &is no more than',& 
       nml%record_9%max_p_tolerance_one_lev_oa, &
       ' Pa.  No obs will be extended to the nearest pressure level.'
   ELSE
      error_message(32:) = ' Single-level above-surface obs will generally only be included in &
       &the objective analysis if they fall on an analysis pressure level. However, will attempt &
       &to adjust obs marked as FM-88 SATOB or FM-97 AIREP to the nearest pressure level.'
   END IF test_908
   call error_handler ( error_number , error_message , fatal , listing )

   !  909:  Make sure pressure tolerance allowed between single-level above-surface obs and the 
   !  pressure level they are analyzed on is less than or equal to the pressure tolerance used to  
   !  do the QC

   test_909 : IF ( ( nml%record_4%use_p_tolerance_one_lev ) .AND. &
                   ( nml%record_9%max_p_tolerance_one_lev_oa .GT. nml%record_4%max_p_tolerance_one_lev_qc ) ) THEN
      error_message = ' '
      error_number = 00013909
      error_message(1:31) = 'check_namelist                 '
      listing = .false.
      fatal = .true.

      WRITE(error_message(32:),'(A,I9,A,I9,A)') ' Pressure tolerance for single-level above-surface obs being included in &
       &the objective analysis is greater than the pressure tolerance for QCing the data.  For QC tolerance is: ',&
       nml%record_4%max_p_tolerance_one_lev_qc, &
       ' Pa, for OA tolerance is: ',&
       nml%record_9%max_p_tolerance_one_lev_oa,&
       ' Pa'

      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_909

   !  910:  Make user aware that a pressure tolerance of 1 hPa for using a single-level above-surface
   !        ob on the objective analysis at a pressure level, is like having no
   !        tolerance at all. 

   test_910 : IF ( ( nml%record_4%use_p_tolerance_one_lev ) .AND. &
                   ( nml%record_9%max_p_tolerance_one_lev_oa .EQ. 1 ) ) THEN
      error_message = ' '
      error_number = 00013910
      error_message(1:31) = 'check_namelist                 '
      listing = .false.
      fatal = .false.

      WRITE(error_message(32:),'(A)') ' Pressure tolerance for single-level above-surface obs being included in &
       &the objective analysis is 1 hPa.  This is equivalent to having no pressure tolerance!'

      call error_handler ( error_number , error_message , fatal , listing )
   END IF test_910
   !BPR END




END SUBROUTINE check_namelist

!------------------------------------------------------------------------------

SUBROUTINE exist_namelist ( unit , filename ) 

!  This routine performs existence checking on the namelist file.

   IMPLICIT NONE

   INCLUDE 'error.inc'

   TYPE ( all_nml ) :: nml

   INTEGER , INTENT ( IN ) :: unit
   CHARACTER *(*) , INTENT ( IN ) :: filename

   LOGICAL :: exist , &
              any_error

   INTERFACE
      INCLUDE 'error.int'
   END INTERFACE

   any_error = .FALSE.

   !  For the file that has been OPENed with the given filename
   !  and unit number, does it exist.

   INQUIRE ( UNIT = unit , EXIST = exist )

   !  001:  If the file does not exist, trap that error.

   test_001 : IF ( .NOT. exist ) THEN
      error_number = 00011001
      error_message(32:)  = ' The filename '// TRIM ( filename) //' does not exist'
      any_error = .TRUE.
   END IF test_001

   any_probs : IF ( any_error ) THEN
      error_message(1:31) = 'exist_namelist                 '
      fatal = .true.
      listing = .false.
      call error_handler ( error_number , error_message , fatal , listing )
   END IF any_probs

END SUBROUTINE exist_namelist

!------------------------------------------------------------------------------

SUBROUTINE store_namelist ( nml )
   
!  This routine stores the NAMELIST variables in the correct locations
!  in the NAMELIST data structure.

   IMPLICIT NONE

   TYPE ( all_nml ) :: nml

   INCLUDE 'namelist.inc'
   INCLUDE 'namelist.common'

   !  Record 1 NAMELIST values:

   WRITE ( nml%record_1%start_date , FMT='(I4.4,"-",I2.2,"-",I2.2,"_",I2.2,":",I2.2,":",I2.2)' ) &
   start_year , start_month , start_day , start_hour , start_minute , start_second
   WRITE ( nml%record_1%end_date , FMT='(I4.4,"-",I2.2,"-",I2.2,"_",I2.2,":",I2.2,":",I2.2)' ) &
   end_year , end_month , end_day , end_hour , end_minute , end_second
   nml%record_1%interval                 = interval    

   !  Record 2 NAMELIST values:

   nml%record_2%fg_filename              = fg_filename    
   nml%record_2%obs_filename             = obs_filename    
   nml%record_2%remove_data_above_qc_flag = remove_data_above_qc_flag    
   nml%record_2%remove_unverified_data   = remove_unverified_data    
   nml%record_2%trim_domain              = trim_domain         
   nml%record_2%trim_value               = trim_value          
   nml%record_2%grid_id                  = grid_id             

   !  Record 3 NAMELIST values:

   nml%record_3%max_number_of_obs        = max_number_of_obs    
   nml%record_3%fatal_if_exceed_max_obs  = fatal_if_exceed_max_obs    

   !  Record 4 NAMELIST values:

   nml%record_4%qc_test_error_max        = qc_test_error_max            
   nml%record_4%qc_test_buddy            = qc_test_buddy                
   nml%record_4%qc_test_vert_consistency = qc_test_vert_consistency     
   nml%record_4%qc_test_convective_adj   = qc_test_convective_adj        
   !BPR BEGIN
   nml%record_4%qc_psfc                  = qc_psfc
   !BPR END

   nml%record_4%max_error_t              = max_error_t                  
   nml%record_4%max_error_uv             = max_error_uv                 
   nml%record_4%max_error_z              = max_error_z                  
   nml%record_4%max_error_rh             = max_error_rh                 
   !BPR BEGIN
   nml%record_4%max_error_dewpoint       = max_error_dewpoint           
   !BPR END
   nml%record_4%max_error_p              = max_error_p                  

   nml%record_4%max_buddy_t              = max_buddy_t                  
   nml%record_4%max_buddy_uv             = max_buddy_uv                 
   nml%record_4%max_buddy_z              = max_buddy_z                  
   nml%record_4%max_buddy_rh             = max_buddy_rh                
   !BPR BEGIN
   nml%record_4%max_buddy_dewpoint       = max_buddy_dewpoint          
   !BPR END
   nml%record_4%max_buddy_p              = max_buddy_p                    

   nml%record_4%buddy_weight             = buddy_weight    

   nml%record_4%max_p_extend_t           = max_p_extend_t    
   nml%record_4%max_p_extend_w           = max_p_extend_w    
   !BPR BEGIN
   nml%record_4%use_p_tolerance_one_lev  = use_p_tolerance_one_lev
   nml%record_4%max_p_tolerance_one_lev_qc  = max_p_tolerance_one_lev_qc
   !BPR END
   
   !  Record 5 NAMELIST values:

   nml%record_5%print_found_obs          = print_found_obs    
   nml%record_5%print_header             = print_header      
   nml%record_5%print_analysis           = print_analysis     
   nml%record_5%print_qc_vert            = print_qc_vert     
   nml%record_5%print_qc_dry             = print_qc_dry     
   nml%record_5%print_obs_files          = print_obs_files    
   nml%record_5%print_error_max          = print_error_max    
   nml%record_5%print_buddy              = print_buddy    
   nml%record_5%print_oa                 = print_oa    

   !  Record 7 NAMELIST values:

   nml%record_7%use_first_guess          = use_first_guess    
   nml%record_7%f4d                      = f4d
   nml%record_7%intf4d                   = intf4d
   nml%record_7%lagtem                   = lagtem

   !  Record 8 NAMELIST values:

   nml%record_8%smooth_type              = smooth_type            
   nml%record_8%smooth_sfc_wind          = smooth_sfc_wind        
   nml%record_8%smooth_sfc_temp          = smooth_sfc_temp        
   nml%record_8%smooth_sfc_rh            = smooth_sfc_rh          
   nml%record_8%smooth_sfc_slp           = smooth_sfc_slp         
   nml%record_8%smooth_upper_wind        = smooth_upper_wind      
   nml%record_8%smooth_upper_temp        = smooth_upper_temp      
   nml%record_8%smooth_upper_rh          = smooth_upper_rh     

   !  Record 9 NAMELIST values:

   nml%record_9%oa_type                  = oa_type    
   nml%record_9%oa_3D_type               = oa_3D_type    
   nml%record_9%mqd_minimum_num_obs      = mqd_minimum_num_obs    
   nml%record_9%mqd_maximum_num_obs      = mqd_maximum_num_obs    
   nml%record_9%radius_influence         = radius_influence    
   nml%record_9%oa_min_switch            = oa_min_switch    
   nml%record_9%oa_max_switch            = oa_max_switch    
   nml%record_9%oa_3D_option             = oa_3D_option    
   !BPR BEGIN
   nml%record_9%oa_psfc                  = oa_psfc
   nml%record_9%scale_cressman_rh_decreases = scale_cressman_rh_decreases    
   nml%record_9%radius_influence_sfc_mult = radius_influence_sfc_mult   
   nml%record_9%max_p_tolerance_one_lev_oa  = max_p_tolerance_one_lev_oa
   !BPR END

END SUBROUTINE store_namelist

!------------------------------------------------------------------------------

!Find pressure ranges where observations cannot be QC'ed against the first-guess
!field or cannot be included in the objective analysis
SUBROUTINE one_lev_coverage_check ( kbu_alloc , pressure ,           &
           max_p_tolerance_one_lev_qc, max_p_tolerance_one_lev_oa )

    IMPLICIT NONE
 
    !Number of vertical levels in first guess field
    INTEGER, INTENT(IN)                           :: kbu_alloc
    !Pressure of each vertical level (Pa).  The first level, the surface just
    !has 100100 since in reality the pressure of the surface varies across the
    !domain
    REAL , INTENT(IN)  , DIMENSION(kbu_alloc)     :: pressure
    !Maximum pressure difference (Pa) between ob and the pressure level it can
    !be QCed against
    INTEGER , INTENT(IN)                          :: max_p_tolerance_one_lev_qc
    !Maximum pressure difference (Pa) between ob and the pressure level on which
    !it can be used in an objective analysis
    INTEGER , INTENT(IN)                          :: max_p_tolerance_one_lev_oa

    INTEGER      :: curk
    INTEGER , DIMENSION(kbu_alloc)                :: top_of_p_range_qc, bot_of_p_range_qc
    INTEGER , DIMENSION(kbu_alloc)                :: top_of_p_range_oa, bot_of_p_range_oa

    !Find the top and bottom of pressure ranges around each pressure level in
    !the first-guess field
    top_of_p_range_qc(1) = 0
    bot_of_p_range_qc(1) = 0
    DO curk = 2,kbu_alloc
     top_of_p_range_qc(curk) = pressure(curk) - max_p_tolerance_one_lev_qc 
     bot_of_p_range_qc(curk) = pressure(curk) + max_p_tolerance_one_lev_qc 
     top_of_p_range_oa(curk) = pressure(curk) - max_p_tolerance_one_lev_oa 
     bot_of_p_range_oa(curk) = pressure(curk) + max_p_tolerance_one_lev_oa 
    ENDDO

    IF( max_p_tolerance_one_lev_qc .GT. 1 ) THEN
     PRINT *,' '
     PRINT *,'*******************************************************************************'
     PRINT *,'Quality control of single-level above-surface observations can use first-guess '
     PRINT *,' fields within ',max_p_tolerance_one_lev_qc,' Pa of observation'
     PRINT *,'Single-level above-surface obs will not be quality controlled against '
     PRINT *,' the first guess field if their pressure falls in the following ranges: '
     WRITE(*,'(A3,I6,A3)') '<= ',top_of_p_range_qc(kbu_alloc),' Pa'
     DO curk = kbu_alloc,3,-1
      !Use ".GE." here because pressure difference must be LESS than the
      !tolerance, not equal to the tolerance
      IF( top_of_p_range_qc(curk-1) .GE. bot_of_p_range_qc(curk) ) THEN
       WRITE(*,'(I6,A3,I6,A3)') top_of_p_range_qc(curk-1),' - ',bot_of_p_range_qc(curk),' Pa'
      END IF
     ENDDO
     WRITE(*,'(A3,I6,A3)') '>= ',bot_of_p_range_qc(2),' Pa' 
     PRINT *,'*******************************************************************************'
     PRINT *,' '
    END IF 

    IF( max_p_tolerance_one_lev_oa .GT. 1 ) THEN
     PRINT *,'*******************************************************************************'
     PRINT *,' '
     PRINT *,'Objective analysis can use single-level above-surface observations within '
     PRINT *,max_p_tolerance_one_lev_oa,' Pa of observation'
     PRINT *,'Single-level above-surface obs will not be included in the objective analysis'
     PRINT *,' at any level if their pressure falls in the following ranges: '
     WRITE(*,'(A3,I6,A3)') '<= ',top_of_p_range_oa(kbu_alloc),' Pa'
     DO curk = kbu_alloc,3,-1
      !Use ".GE." here because pressure difference must be LESS than the
      !tolerance, not equal to the tolerance
      IF( top_of_p_range_oa(curk-1) .GE. bot_of_p_range_oa(curk) ) THEN
       WRITE(*,'(I6,A3,I6,A3)') top_of_p_range_oa(curk-1),' - ',bot_of_p_range_oa(curk),' Pa'
      END IF
     ENDDO
     WRITE(*,'(A3,I6,A3)') '>= ',bot_of_p_range_oa(2),' Pa' 
     PRINT *,'*******************************************************************************'
     PRINT *,' '
    END IF 


END SUBROUTINE one_lev_coverage_check 
!------------------------------------------------------------------------------

END MODULE namelist
