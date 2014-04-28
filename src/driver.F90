!------------------------------------------------------------------------------
SUBROUTINE driver ( filename , filename_out , &
bhi , bhr , nml , iew_alloc , jns_alloc , kbu_alloc , & 
current_date_8 , current_time_6 , date_char , icount , total_count , &
mqd_count , mqd_abs_min )

!  This subroutine is called from the main analysis program, it is called
!  only once per time period.  This subroutine calls all of the routines and 
!  functions that are required to ingest the data, qc the data, perform the 
!  objective analysis, and output the final analysis and observations.

   USE date_pack
   USE input_data
   USE map_utils
   USE map_utils_helper
   USE namelist
   USE tinterp
   USE observation

   IMPLICIT NONE
 
   !  The record header for this domain, for this time, in the original
   !  form.

   INCLUDE 'big_header.inc'

   CHARACTER ( LEN = 132 )  , INTENT ( IN )       :: filename
   CHARACTER ( LEN = 132 )  , INTENT ( INOUT )    :: filename_out
   CHARACTER ( LEN = 132 )                        :: obs_nudge_file
   INTEGER , INTENT ( IN )                        :: current_date_8 , & 
                                                     current_time_6 , & 
                                                     icount , total_count
   CHARACTER (LEN=19) , INTENT(IN)                :: date_char
   CHARACTER (LEN=24)                             :: dt_char
   LOGICAL                                        :: exist

   !  Observation information.

   INTEGER                                        :: number_of_obs, i_num, i_found_time
   INTEGER                                        :: total_dups
   INTEGER       , ALLOCATABLE , DIMENSION ( : )  :: index
   TYPE (report) , ALLOCATABLE , DIMENSION ( : )  :: obs
   TYPE (report)                                  :: obs_sort_tmp
   CHARACTER ( LEN = 132 )                        :: dummy_filename, tmp_filename

   !  The NAMELIST variables, with a small amount of error checks
   !  from the main program.

   TYPE ( all_nml )                 :: nml      ! all namelist information from 
                                                ! all namelist records

   INTEGER                                         :: iew_map , jns_map
   INTEGER                                         :: final_fdda_count
   INTEGER                                         :: i_rad,max_rad, max_rad_used


   INCLUDE 'error.inc'
   INCLUDE 'proc_get_info_header.inc'
   INCLUDE 'first_guess_size.inc'
   INCLUDE 'first_guess.inc'
   INCLUDE 'map.inc'
   REAL                        , DIMENSION (jns_alloc,iew_alloc) :: tobbox, odis
   REAL                        , DIMENSION (jns_alloc,iew_alloc) :: tobbox_ana, odis_ana

   INTERFACE 
      INCLUDE 'error.int'
      INCLUDE 'proc_get_info_header.int'
      INCLUDE 'proc_first_guess.int'
      INCLUDE 'proc_qc.int'
      INCLUDE 'proc_obs_sort.int'
      INCLUDE 'proc_oa.int'
      INCLUDE 'proc_final_analysis.int'
      INCLUDE 'sample.int'
   END INTERFACE

   INTEGER                          :: loop_count
   INTEGER                          :: unit=9
   CHARACTER ( LEN = 132 )          :: filename_output

   !  FDDA counters and such galore!

   INTEGER :: fdda_loop , fdda_loop_max , icount_fdda , icount_1 , icount_2 
   INTEGER :: fdda_date_8 , fdda_time_6
   CHARACTER (LEN=19) :: date_char_fdda

   INTEGER :: obs_file_count
   INTEGER :: mqd_count , mqd_abs_min

   !  If we are doing sfc FDDA, the first time in we do the traditional analysis AND the surface
   !  FDDA fields for that time (fdda_loop_max=2).  All of the subsequent times, we process the 
   !  traditional analysis + the FDDA times upto and including the analyis period (i.e., 12 Z 3d
   !  analysis AND THEN the 03, 06, 09, and 12 Z surface analysis FDDA fields).

   IF      ( ( nml%record_7%f4d ) .AND. ( icount .EQ. 1 ) ) THEN
      fdda_loop_max = 1 + 1
   ELSE IF ( ( nml%record_7%f4d ) .AND. ( icount .GT. 1 ) ) THEN
      fdda_loop_max = 1 + nml%record_1%interval / nml%record_7%intf4d
   ELSE IF ( .NOT. nml%record_7%f4d ) THEN
      fdda_loop_max = 1
   END IF

   final_fdda_count = (total_count-2)*fdda_loop_max + fdda_loop_max - (total_count-2)

   fdda_ish_loop : DO fdda_loop = 1 , fdda_loop_max

      !  If this is an FDDA run, we need to process the intermediate FDDA dates.  Always after
      !  the first loop thorugh fdda_ish_loop, we pick up the FDDA processing.  During the first
      !  icount loop, the FDDA time and the traditional time are the same.   For all of the 
      !  subsequent times, we need to increment forward so that the obs files are placed in an
      !  easy to understand fashion (actually, we are decrementing from the current time, but
      !  let's not get too picky here).
      
      IF      ( ( nml%record_7%f4d ) .AND. ( icount .EQ. 1 ) ) THEN
         IF ( fdda_loop .GT. 1 ) THEN
            date_char_fdda = date_char
         END IF
      ELSE IF ( ( nml%record_7%f4d ) .AND. ( icount .GT. 1 ) ) THEN
         IF ( fdda_loop .GT. 1 ) THEN
            CALL  geth_newdate ( date_char_fdda , date_char , &
            -1 * nml%record_7%intf4d * ( fdda_loop_max - fdda_loop ) )
         END IF
      END IF
   
      !  Initialize the observations per box to be zero.
      !  Initialize the observation density to be zero.

      tobbox = 0
      odis = 0

      !  Print out so that we know the time period that we are processing.

      IF ( (  nml%record_7%f4d ) .AND. ( fdda_loop .GT. 1 ) ) THEN
         IF ( fdda_loop .NE. fdda_loop_max ) &
         WRITE ( UNIT = * , FMT = '( ///,"----------------------------------------------------------------------------",/,&
         &"Time Loop Processing for SFC FDDA, date = ",A,//)' ) date_char_fdda
      ELSE
         WRITE ( UNIT = * , FMT = '( ///,"----------------------------------------------------------------------------",/,&
         &"Time Loop Processing, date = ",A,//)' ) date_char
      END IF

      !  We only do this for the traditional analysis time periods, not for the FDDA processing.  So, we 
      !  either do this when there is no FDDA, or we do it during the first of the fdda_loop iterations.
      !  This is all of the stuff that is related to input of the first guess data.
      
      IF (   (  .NOT. nml%record_7%f4d ) .OR. &
           ( (        nml%record_7%f4d ) .AND. ( fdda_loop .EQ. 1 ) ) ) THEN 

         !  Show what is available as addressable from the record header, this provides a printout file.

         IF ( nml%record_5%print_header ) THEN
            CALL proc_get_info_header ( nml%record_5%print_header , help = help ) 
         END IF

         !  Initialize all of the constants for this data set, at this time slice,
         !  for this domain, with the current vertical resolution.  This routine
         !  is defined with optional arguments, allowing for limited access to
         !  specific values.

         CALL proc_get_info_header ( nml%record_5%print_header , &
         iewd , jnsd , kbu_alloc , grid_id , map_projection , expanded , iewe , jnse , &
         dxc , lat_center , lon_center , cone_factor , true_lat1 , true_lat2 , pole , &
         dxd , ptop )

         IF ( icount .EQ. 1 ) THEN
         IF ( nml%record_9%radius_influence(1) .LE. 0 ) THEN
            max_rad_used = 0
            max_rad = ABS (nml%record_9%radius_influence(1) )
            IF ( nml%record_9%radius_influence(1) == 0 ) max_rad = 4
            ! ASSUME AVG. SPACING OF UPPER AIR STATIONS IS 325 KM
            nml%record_9%radius_influence(1) = MAX (nint( 325 * ( 1.6 / (dxd/1000.) ) + .45 ) , 5 )
            print*," #########################################################################"
            print*," Setting radius of influence for Cressman scheme"
            write(*,'("    The radius of influence of each scan is set to:")')
            write(*,'("     ",i3,",",$)') nml%record_9%radius_influence(1)
            DO i_rad = 2, max_rad
              nml%record_9%radius_influence(i_rad) = nint ( nml%record_9%radius_influence(i_rad-1) * 0.7 + .45 )
              IF ( nml%record_9%radius_influence(i_rad) == nml%record_9%radius_influence(i_rad-1) ) &
                   nml%record_9%radius_influence(i_rad) =  nml%record_9%radius_influence(i_rad-1) - 1
              IF (nml%record_9%radius_influence(i_rad) .LT. 2 ) THEN
                  nml%record_9%radius_influence(i_rad) = -1
              ELSE
                write(*,'(i3,",",$)') nml%record_9%radius_influence(i_rad)
                max_rad_used = max_rad_used + 1
              END IF
            END DO
            write(*,'("  ")') 
            IF (max_rad_used > 5) print*,"  WARNING: Not recommended to use more than 5 scans" 
            print*," #########################################################################"
            print*,"   "
         ELSE
            print*," #########################################################################"
            print*," Settings for the radius of influence for Cressman scheme"
            write(*,'("     ",i3,",",$)') nml%record_9%radius_influence(1)
            DO i_rad = 2, 10
              IF (nml%record_9%radius_influence(i_rad) .NE. -1 ) THEN
                write(*,'(i3,",",$)') nml%record_9%radius_influence(i_rad)
              END IF
            END DO
            write(*,'("  ")') 
            print*," #########################################################################"
            print*,"   "
         ENDIF
         ENDIF

         !  All of the following routines assume distances to be in km, and the pressures
         !  are assumed to be in hPa.

         dxc = dxc * 0.001 
         dxd = dxd * 0.001
         ptop = NINT ( ptop * 0.01 )

         !  These values (*_map) are set for the mapping utilities.  
         !  They need to know how big this domain is, which we supply from the above values. 
         
         iew_map   = iewe
         jns_map   = jnse

         !  Process the first guess data set.

         CALL proc_first_guess ( filename , &
         bhi , bhr , num3d , num2d , num1d , &
         t , u , v , uA , vA , uC , vC , h , rh , pres , terrain , &
         latitude_x , longitude_x , latitude_d , longitude_d , &
         slp_x , slp_C , sst , snow , &
         iew_alloc , jns_alloc , kbu_alloc , pressure , & 
         nml%record_5%print_analysis , & 
         current_date_8 , current_time_6 , date_char , icount )
#ifdef PRESSURE
pressure(:kbu_alloc) = pressure(:kbu_alloc) * 100.
#endif

         !  Initialize the mapping information.  Note we do a quick units switch-a-roo for the grid
         !  distance: from km to m to km.  You hardly notice if you blink.

         IF ( icount .EQ. 1 ) THEN
            dxd = dxd * 1000.
            IF      ( map_projection .EQ. 1 ) THEN
               CALL map_set ( PROJ_LC   , latitude_x(1,1) , longitude_x(1,1) , dxd , lon_center , true_lat1 , true_lat2 , projx )
               CALL map_set ( PROJ_LC   , latitude_d(1,1) , longitude_d(1,1) , dxd , lon_center , true_lat1 , true_lat2 , projd )
            ELSE IF ( map_projection .EQ. 2 ) THEN
               CALL map_set ( PROJ_PS   , latitude_x(1,1) , longitude_x(1,1) , dxd , lon_center , true_lat1 , true_lat2 , projx )
               CALL map_set ( PROJ_PS   , latitude_d(1,1) , longitude_d(1,1) , dxd , lon_center , true_lat1 , true_lat2 , projd )
            ELSE IF ( map_projection .EQ. 3 ) THEN
               CALL map_set ( PROJ_MERC , latitude_x(1,1) , longitude_x(1,1) , dxd , lon_center , true_lat1 , true_lat2 , projx )
               CALL map_set ( PROJ_MERC , latitude_d(1,1) , longitude_d(1,1) , dxd , lon_center , true_lat1 , true_lat2 , projd )
            ELSE
               PRINT '(A)','Whoa there pardner, what projection yew tryin to pull on us?'
               PRINT '(A,I8,A)','We accept 1 through 3, and you tried a ',map_projection,'.'
               STOP 'map_proj_incorrect'
            END IF
            dxd = dxd * 0.001
         END IF

         !  All of the following routines assume this pressure to be in hPa.

         pressure(:kbu_alloc) = NINT ( pressure(:kbu_alloc) * 0.01 )
         slp_x(:jns_alloc-1,:iew_alloc-1) = slp_x(:jns_alloc-1,:iew_alloc-1) * 0.01

      END IF

      !  To keep things ordered here, since we just input all of the first guess data for the
      !  non FDDA sorts of processing, let's mimic that and generate the interpolated fields
      !  for the FDDA processing.  This only happens when f4d is TRUE, and when we are not in the
      !  first loop of the fdda_loop iterations (which is for the traditional processing).

      IF      ( ( nml%record_7%f4d ) .AND. ( fdda_loop .GT. 1 ) .AND. ( icount .EQ. 1 ) ) THEN
         icount_fdda = 1
         icount_1 = 1
         icount_2 = 1
         IF ( nml%record_7%lagtem ) THEN
            slp_x(:jns_alloc-1,:iew_alloc-1) = slp_x(:jns_alloc-1,:iew_alloc-1) * 0.01
         ELSE
            CALL temporal_interp ( t , u , v , uA , vA , uC , vC , h , rh , pres , slp_x , slp_C , snow , pressure , &
            iew_alloc , jns_alloc , kbu_alloc ,  num3d , num2d , &
            icount_fdda , icount_1 , icount_2 )
            slp_x(:jns_alloc-1,:iew_alloc-1) = slp_x(:jns_alloc-1,:iew_alloc-1) * 0.01
         END IF
      ELSE IF ( ( nml%record_7%f4d ) .AND. ( fdda_loop .GT. 1 ) .AND. ( icount .GT. 1 ) ) THEN
         icount_fdda = ( icount - 1 ) * ( nml%record_1%interval / nml%record_7%intf4d ) + 1 - ( fdda_loop_max - fdda_loop ) 
         icount_1 = ( icount - 2 ) * ( nml%record_1%interval / nml%record_7%intf4d ) + 1  
         icount_2 = ( icount - 1 ) * ( nml%record_1%interval / nml%record_7%intf4d ) + 1  
         IF ( nml%record_7%lagtem ) THEN
            CALL lagtem_assign ( t , u , v , uA , vA , uC , vC , h , rh , pres , slp_x , slp_C , snow , &
            iew_alloc , jns_alloc , kbu_alloc ,  num3d , num2d , &
            icount_fdda , icount_1 , icount_2 ) 
            slp_x(:jns_alloc-1,:iew_alloc-1) = slp_x(:jns_alloc-1,:iew_alloc-1) * 0.01
         ELSE
            CALL temporal_interp ( t , u , v , uA , vA , uC , vC , h , rh , pres , slp_x , slp_C , snow , pressure , &
            iew_alloc , jns_alloc , kbu_alloc ,  num3d , num2d , &
            icount_fdda , icount_1 , icount_2 ) 
            slp_x(:jns_alloc-1,:iew_alloc-1) = slp_x(:jns_alloc-1,:iew_alloc-1) * 0.01
         END IF
      END IF

      !  This is an easy way to make sure that a "clean" name gets
      !  passed into the couple of routines that will deal with the
      !  observation file name.  
      
      dummy_filename = '                                                  ' // &
                       '                                                  ' // &
                       '                                '

      !  If fdda time, we need an updated date for the time we are processing.
      IF ( ( nml%record_7%f4d ) .AND. ( fdda_loop .GT. 1 ) ) THEN 
         CALL get_date ( fdda_date_8 , fdda_time_6 , date_char_fdda ) 
      END IF

      IF ( nml%record_7%f4d  .AND. (fdda_loop == fdda_loop_max) ) GOTO 101

      IF (   (  .NOT. nml%record_7%f4d ) .OR. &
           ( (        nml%record_7%f4d ) .AND. ( fdda_loop .EQ. 1 ) ) ) THEN 
         dummy_filename = &
         TRIM ( nml%record_2%obs_filename )//":"//date_char(1:13)
      ELSE
         IF ( icount .EQ. 1 ) THEN
            icount_fdda = ( icount - 1 ) * ( nml%record_1%interval / nml%record_7%intf4d ) + 1
         ELSE
            icount_fdda = ( icount - 1 ) * ( nml%record_1%interval / nml%record_7%intf4d ) + 1 - ( fdda_loop_max - fdda_loop ) 
         END IF
         dummy_filename = &
         TRIM ( nml%record_2%obs_filename )//":"//date_char_fdda(1:13)
      END IF
      INQUIRE ( EXIST = exist , FILE = dummy_filename )
      IF ( .NOT. exist ) THEN
        print*,"WARNING: The following file does not exist - process as if file is empty"
        print*,"        ", trim(dummy_filename)
      ELSE
        print*,"Using ", trim(dummy_filename), " as obs input file"
      END IF

      !  For real time applications, we are only processing the initial times with
      !  observations.  All of the other time periods will simply be "read in
      !  and written out".

      !  Here is the way that we find out if this time period is to be
      !  objectively analyzed, or just passed through the program with
      !  the header modified - the observation file must not be called
      !  'null'.  We have already checked to see if the observation file
      !  exists.

      IF ( dummy_filename(1:4) .NE. 'null' ) THEN
      
         !  We are ready to process the observations.  We need to allocate more
         !  space than is required.  Since only some constant space is 
         !  allocated (not the space for the entire sounding), this is not
         !  too inefficient.
      
         ALLOCATE ( obs   ( nml%record_3%max_number_of_obs ) )
         ALLOCATE ( index ( nml%record_3%max_number_of_obs ) )
       
         !  The observation data can now be ingested, sorted, duplicates
         !  removed.  On return, we have all of the sorted observations 
         !  filled (obs), an index of the sequence of the sorted observations 
         !  (index), and the number of observations with valid data 
         !  (number_of_obs).  This is separated into traditional processing and
         !  the FDDA tasks.  Only the input date is different, all other functionality
         !  and input is the same.
      
         IF ( ( .NOT. nml%record_7%f4d ) .OR. & 
              ( (     nml%record_7%f4d ) .AND. ( fdda_loop .EQ. 1 ) ) ) THEN 
            CALL proc_obs_sort ( dummy_filename  , unit+10 , &
            obs , number_of_obs , nml%record_3%max_number_of_obs , nml%record_3%fatal_if_exceed_max_obs , total_dups , &
            index , nml%record_5%print_found_obs , nml%record_5%print_obs_files , &
            kbu_alloc , pressure , & 
            nml%record_4%max_p_extend_t           , nml%record_4%max_p_extend_w , &
            h , iew_alloc , jns_alloc , map_projection , current_date_8 , current_time_6 , fdda_loop )  
         ELSE
            CALL get_date ( fdda_date_8 , fdda_time_6 , date_char_fdda ) 
            CALL proc_obs_sort ( dummy_filename  , unit+10 , &
            obs , number_of_obs , nml%record_3%max_number_of_obs , nml%record_3%fatal_if_exceed_max_obs , total_dups , &
            index , nml%record_5%print_found_obs , nml%record_5%print_obs_files , &
            kbu_alloc , pressure , & 
            nml%record_4%max_p_extend_t           , nml%record_4%max_p_extend_w , &
            h , iew_alloc , jns_alloc , map_projection , fdda_date_8 , fdda_time_6 , fdda_loop )  
         END IF
         
         !  Run the quality control (QC) procedures on the observations.
      
         IF ( ( .NOT. nml%record_7%f4d ) .OR. & 
              ( (     nml%record_7%f4d ) .AND. ( fdda_loop .EQ. 1 ) ) ) THEN 
            CALL proc_qc ( iew_alloc , jns_alloc , kbu_alloc , number_of_obs , &
            total_dups , map_projection , &
            nml%record_4%qc_test_error_max        , nml%record_4%qc_test_buddy          , &
            nml%record_4%qc_test_vert_consistency , nml%record_4%qc_test_convective_adj , &
            nml%record_4%max_error_t              , nml%record_4%max_error_uv           , &
            nml%record_4%max_error_z              , nml%record_4%max_error_rh           , &
            nml%record_4%max_error_p              , nml%record_5%print_error_max        , &
            nml%record_4%max_buddy_t              , nml%record_4%max_buddy_uv           , &
            nml%record_4%max_buddy_z              , nml%record_4%max_buddy_rh           , &
            nml%record_4%max_buddy_p              , nml%record_5%print_buddy            , &
            nml%record_5%print_found_obs                                                , &
            nml%record_5%print_qc_vert            , nml%record_5%print_qc_dry           , &
            pressure , current_date_8 , current_time_6 , dxd , 1. , &
            obs , index , nml%record_3%max_number_of_obs , &
            t , u , v , h , rh , slp_x , sst , tobbox , odis )
         ELSE
            CALL proc_qc ( iew_alloc , jns_alloc , kbu_alloc , number_of_obs , &
            total_dups , map_projection , &
            nml%record_4%qc_test_error_max        , nml%record_4%qc_test_buddy          , &
            nml%record_4%qc_test_vert_consistency , nml%record_4%qc_test_convective_adj , &
            nml%record_4%max_error_t              , nml%record_4%max_error_uv           , &
            nml%record_4%max_error_z              , nml%record_4%max_error_rh           , &
            nml%record_4%max_error_p              , nml%record_5%print_error_max        , &
            nml%record_4%max_buddy_t              , nml%record_4%max_buddy_uv           , &
            nml%record_4%max_buddy_z              , nml%record_4%max_buddy_rh           , &
            nml%record_4%max_buddy_p              , nml%record_5%print_buddy            , &
            nml%record_5%print_found_obs                                                , &
            nml%record_5%print_qc_vert            , nml%record_5%print_qc_dry           , &
            pressure , fdda_date_8 , fdda_time_6 , dxd , 1. , &
            obs , index , nml%record_3%max_number_of_obs , &
            t , u , v , h , rh , slp_x , sst , tobbox , odis )
         END IF
      
         !  After the QC process, the observations are available for output
         !  in a similar fashion to the non-QC'ed data.  The differences
         !  between this file, "qc_out", and "useful_out" are the QC flags
         !  and any changes from QC routines (vert_consistency_check and/or
         !  dry_conv_adjustment).
      
         IF ( ( .NOT. nml%record_7%f4d ) .OR. & 
              ( (     nml%record_7%f4d ) .AND. ( fdda_loop .EQ. 1 ) ) ) THEN 
            CALL make_date ( current_date_8 , current_time_6 , dt_char )
         ELSE 
            CALL make_date ( fdda_date_8 , fdda_time_6 , dt_char )
         END IF
         IF ( nml%record_5%print_obs_files ) THEN
            !  April 2009 - name changed from qc_out to qc_obs_raw
            !  We also add domain info here
            WRITE (tmp_filename,'("qc_obs_raw.d",i2.2,".")') nml%record_2%grid_id
            CALL output_obs ( obs , 2 , trim(tmp_filename)//dt_char , number_of_obs ,   &
                              1 , .TRUE., .TRUE., 200000, .FALSE., kbu_alloc, pressure, &
                              filename )
         END IF
      
         !  With the observations properly QC'ed and stored, and the first guess
         !  background field already ingested, we can proceed with the objective
         !  analysis.
      
         IF ( ( .NOT. nml%record_7%f4d ) .OR. & 
              ( (     nml%record_7%f4d ) .AND. ( fdda_loop .EQ. 1 ) ) ) THEN 
            CALL proc_oa ( t , u , v , rh , slp_x , pressure , &
            iew_alloc , jns_alloc , kbu_alloc , &
            current_date_8 , current_time_6 , fdda_loop , mqd_count , mqd_abs_min , &
            nml%record_3%max_number_of_obs , number_of_obs , total_dups , &
            map_projection , obs , dxd , lat_center , &
            nml%record_5%print_oa                 , nml%record_5%print_found_obs          , &
            nml%record_5%print_obs_files                                                  , &
            nml%record_7%use_first_guess                                                  , &
            nml%record_8%smooth_type              , nml%record_8%smooth_sfc_wind          , & 
            nml%record_8%smooth_sfc_temp          , nml%record_8%smooth_sfc_rh            , & 
            nml%record_8%smooth_sfc_slp           , nml%record_8%smooth_upper_wind        , & 
            nml%record_8%smooth_upper_temp        , nml%record_8%smooth_upper_rh          , &
            nml%record_9%oa_type                  , nml%record_9%oa_3D_type               , &
            nml%record_9%mqd_minimum_num_obs      , nml%record_9%mqd_maximum_num_obs      , &
            nml%record_9%oa_max_switch            , nml%record_9%radius_influence         , &
            nml%record_9%oa_min_switch            , nml%record_9%oa_3D_option             , &
            nml%record_2%grid_id )

            !  Store the final analysis back into the all_3d and all_2d arrays if we are doing
            !  SFC FDDA.  Why?  So that when we do the LAGTEM or temporal interpolation, we are
            !  using the final analysis as the first guess field for all of the time periods.

            IF ( nml%record_7%f4d ) THEN
               CALL store_fa (  t , u , v , rh , slp_x , iew_alloc , jns_alloc , kbu_alloc , num3d , num2d , icount )
            END IF
         ELSE
            CALL proc_oa ( t , u , v , rh , slp_x , pressure , &
            iew_alloc , jns_alloc , kbu_alloc , &
            fdda_date_8 , fdda_time_6 , fdda_loop , mqd_count , mqd_abs_min , &
            nml%record_3%max_number_of_obs , number_of_obs , total_dups , &
            map_projection , obs , dxd , lat_center , &
            nml%record_5%print_oa                 , nml%record_5%print_found_obs          , &
            nml%record_5%print_obs_files                                                  , &
            nml%record_7%use_first_guess                                                  , &
            nml%record_8%smooth_type              , nml%record_8%smooth_sfc_wind          , & 
            nml%record_8%smooth_sfc_temp          , nml%record_8%smooth_sfc_rh            , & 
            nml%record_8%smooth_sfc_slp           , nml%record_8%smooth_upper_wind        , & 
            nml%record_8%smooth_upper_temp        , nml%record_8%smooth_upper_rh          , &
            nml%record_9%oa_type                  , nml%record_9%oa_3D_type               , &
            nml%record_9%mqd_minimum_num_obs      , nml%record_9%mqd_maximum_num_obs      , &
            nml%record_9%oa_max_switch            , nml%record_9%radius_influence         , &
            nml%record_9%oa_min_switch            , nml%record_9%oa_3D_option             , &
            nml%record_2%grid_id )
         END IF

      END IF

      !  After the final analysis has been generated, it needs to be
      !  output.  OR, if this is a time period that is not the initial
      !  time (so that there was no objective analysis) we write it out
      !  as well.

      IF ( ( .NOT. nml%record_7%f4d ) .OR. & 
           ( (     nml%record_7%f4d ) .AND. ( fdda_loop .EQ. 1 ) ) ) THEN 
         CALL make_date ( current_date_8 , current_time_6 , dt_char )
      ELSE 
         CALL make_date ( fdda_date_8 , fdda_time_6 , dt_char )
      END IF
      IF ( nml%record_5%print_obs_files ) THEN
         WRITE (tmp_filename,'("qc_obs_used.d",i2.2,".")') nml%record_2%grid_id
         CALL output_obs ( obs , 2 , trim(tmp_filename)//dt_char , number_of_obs ,  1 , &
                           .TRUE., .TRUE., nml%record_2%remove_data_above_qc_flag,      &
                           nml%record_2%remove_unverified_data, kbu_alloc, pressure,    &
                           filename )  
      END IF

 101  CONTINUE   ! come here if end of fdda loop, so we don't repeat the analysis time

      IF ( ( .NOT. nml%record_7%f4d ) .OR. & 
           ( (     nml%record_7%f4d ) .AND. ( fdda_loop .EQ. 1 ) ) ) THEN 
         CALL proc_final_analysis ( filename , filename_out , &
         bhi , bhr , t , u , v , uA , vA , uC , vC , h , rh , pres , terrain , &
         latitude_x , longitude_x , latitude_d , longitude_d , &
         slp_x , slp_C , sst , snow , tobbox , odis , &
         pressure , ptop , &
         iew_alloc , jns_alloc , kbu_alloc , iewd , jnsd , & 
         nml%record_5%print_header , nml%record_5%print_analysis , &
         current_date_8 , current_time_6 , fdda_loop , icount_fdda , &
         icount , total_count , nml%record_1%interval , &
         nml%record_4%max_error_t , nml%record_4%max_error_uv           , &
         nml%record_4%max_error_z , nml%record_4%max_error_p/100. , &
         nml%record_4%buddy_weight , date_char , &
         nml%record_2%fg_filename , nml%record_9%oa_3D_option , &
         nml%record_7%intf4d , nml%record_7%lagtem , &
         nml%record_9%oa_type , nml%record_9%oa_3D_type , nml%record_9%radius_influence )
         tobbox_ana = tobbox
         odis_ana = odis
      ELSE 
         IF ( fdda_loop == fdda_loop_max ) THEN
           tobbox = tobbox_ana
           odis = odis_ana
         ENDIF
         CALL proc_final_analysis ( filename , filename_out , &
         bhi , bhr , t , u , v , uA , vA , uC , vC , h , rh , pres , terrain , &
         latitude_x , longitude_x , latitude_d , longitude_d , &
         slp_x , slp_C , sst , snow , tobbox , odis , &
         pressure , ptop , &
         iew_alloc , jns_alloc , kbu_alloc , iewd , jnsd , & 
         nml%record_5%print_header , nml%record_5%print_analysis , &
         fdda_date_8 , fdda_time_6 , fdda_loop , icount_fdda , &
         icount , final_fdda_count , nml%record_7%intf4d , &
         nml%record_4%max_error_t , nml%record_4%max_error_uv           , &
         nml%record_4%max_error_z , nml%record_4%max_error_p/100. , &
         nml%record_4%buddy_weight , date_char , &
         nml%record_2%fg_filename , nml%record_9%oa_3D_option , &
         nml%record_7%intf4d , nml%record_7%lagtem , &
         nml%record_9%oa_type , nml%record_9%oa_3D_type , nml%record_9%radius_influence )
         !nml%record_4%buddy_weight , nml%record_1%start_date )
      END IF

      !!! Let's make sure the output to OBS_DOMAINnxx is sorted in time.
      IF ( dummy_filename(1:4) .NE. 'null' .AND. fdda_loop .NE. fdda_loop_max ) THEN
        i_found_time = 0
        loop_sort_time : DO
        DO i_num = 1, number_of_obs-1
     
          IF ( obs(i_num)%valid_time%date_char > obs(i_num+1)%valid_time%date_char ) THEN
            !print*,"FOUND one - ", i_num, " of ", number_of_obs, " (", i_found_time,")"
            !print* ,"DATES:",obs(i_num)%valid_time%date_char ,"  -  ", obs(i_num+1)%valid_time%date_char 
            obs_sort_tmp = obs(i_num)
            obs(i_num) = obs(i_num+1)
            obs(i_num+1) = obs_sort_tmp
            i_found_time = i_found_time + 1
          END IF
        
        END DO
        if (i_found_time == 0 ) exit loop_sort_time
        if (i_found_time  > 0 ) i_found_time = 0
        END DO loop_sort_time
      ENDIF


      IF ( fdda_loop.EQ.1) THEN
        obs_file_count = (icount-1)*2 + 1
        IF ( .NOT. nml%record_7%f4d ) obs_file_count = icount
        WRITE (obs_nudge_file,'("OBS_DOMAIN",i1,i2.2)') nml%record_2%grid_id, obs_file_count
        INQUIRE ( EXIST = exist , FILE = obs_nudge_file )
        CALL output_obs ( obs , 2 , trim(obs_nudge_file), number_of_obs ,  1 ,      &
                          .TRUE., exist, nml%record_2%remove_data_above_qc_flag,    &
                          nml%record_2%remove_unverified_data, kbu_alloc, pressure, &
                          filename )  
      ELSEIF ( fdda_loop.EQ.2 .AND. fdda_loop.NE.fdda_loop_max) THEN
        obs_file_count = (icount-1)*2 
        WRITE (obs_nudge_file,'("OBS_DOMAIN",i1,i2.2)') nml%record_2%grid_id, obs_file_count
        INQUIRE ( EXIST = exist , FILE = obs_nudge_file )
        CALL output_obs ( obs , 2 , trim(obs_nudge_file), number_of_obs ,  1 ,      &
                          .TRUE., exist, nml%record_2%remove_data_above_qc_flag,    &
                          nml%record_2%remove_unverified_data, kbu_alloc, pressure, &
                          filename )  
      ELSEIF ( fdda_loop.GT.2 .AND. fdda_loop.NE.fdda_loop_max) THEN
        obs_file_count = (icount-1)*2 
        WRITE (obs_nudge_file,'("OBS_DOMAIN",i1,i2.2)') nml%record_2%grid_id, obs_file_count
        CALL output_obs ( obs , 2 , trim(obs_nudge_file), number_of_obs ,  1 ,      &
                          .TRUE., .FALSE., nml%record_2%remove_data_above_qc_flag,  &
                          nml%record_2%remove_unverified_data, kbu_alloc, pressure, &
                          filename )  
      ENDIF


      !  If this is the time that we have observations, we can de-allocate the space.

      IF ( dummy_filename(1:4) .NE. 'null' .AND. fdda_loop .NE. fdda_loop_max ) THEN

         !  We are finished with this time period.  We can remove all of the
         !  allocated space for the observation data.
        
         DEALLOCATE ( index ) 
         zap_space : DO loop_count = 1 , number_of_obs
            IF ( ASSOCIATED (obs(loop_count)%surface ) ) THEN
               CALL dealloc_meas ( obs(loop_count)%surface ) 
            END IF
         END DO zap_space
         DEALLOCATE ( obs   )
      
      END IF

   END DO fdda_ish_loop

   !  DEALLOCATE the 3d, 2d, and 1d space that holds the raw input.  Because we
   !  ALLOCATE the space in the main program (only once, not in a loop), all
   !  we want to DEALLOCATE is the large space associated with the arrays.
   !  Due to FDDA temporal interpolation, we can only DEALLOCATE the first
   !  time.

   IF ( nml%record_7%f4d ) THEN
      IF ( icount .NE. 1 ) THEN
         DO loop_count = 1 , loop3
            DEALLOCATE ( all_3d(loop_count,first_time)%array )
         END DO 
         DO loop_count = 1 , loop2
            DEALLOCATE ( all_2d(loop_count,first_time)%array )
         END DO 
      END IF
   ELSE IF ( .NOT. nml%record_7%f4d ) THEN
      DO loop_count = 1 , loop3
         DEALLOCATE ( all_3d(loop_count,first_time)%array )
      END DO 
      DO loop_count = 1 , loop2
         DEALLOCATE ( all_2d(loop_count,first_time)%array )
      END DO 
   END IF
   DO loop_count = 1 , loop1
      DEALLOCATE ( all_1d(loop_count)%array )
   END DO 

END SUBROUTINE driver
