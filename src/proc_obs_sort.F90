!-------------------------------------------------------------------------------

SUBROUTINE proc_obs_sort ( obs_filename , unit , &
obs , number_of_obs , total_number_of_obs , fatal_if_exceed_max_obs , total_dups , index1 , &
print_out_obs_found , print_out_files , & 
!BPR BEGIN
!levels , pressure , &
levels , pressure , slp_x , temperature , use_p_tolerance_one_lev , &
!BPR END
max_p_extend_t , max_p_extend_w , &
height , iew , jns , map_projection , date , time , fdda_loop ) 

!  This routine handles the processing of all of the observations.

   USE observation

   IMPLICIT NONE

   CHARACTER ( LEN = 132 ) , INTENT ( IN )          :: obs_filename
   INTEGER , INTENT ( IN )                          :: unit
   INTEGER , INTENT ( OUT )                         :: number_of_obs
   INTEGER , INTENT ( IN )                          :: total_number_of_obs
   LOGICAL , INTENT ( IN )                          :: fatal_if_exceed_max_obs
   INTEGER , INTENT ( OUT )                         :: total_dups
   TYPE (report) , DIMENSION (total_number_of_obs ) :: obs
   INTEGER       , DIMENSION (total_number_of_obs ) :: index1
   LOGICAL , INTENT ( IN )                          :: print_out_obs_found , &
                                                       print_out_files
   INTEGER , INTENT ( IN )                          :: levels
   REAL    , INTENT ( IN ) , DIMENSION ( levels )   :: pressure
   REAL                                             :: max_p_extend_t , &
                                                       max_p_extend_w
   INTEGER                                          :: iew , jns , map_projection
   !BPR BEGIN
   REAL    , INTENT ( IN ) , DIMENSION ( jns, iew ) :: slp_x 
   REAL    , INTENT ( IN ) , DIMENSION ( jns , iew , levels ) :: temperature
   LOGICAL , INTENT ( IN )                          :: use_p_tolerance_one_lev
   !BPR END
   REAL , DIMENSION ( jns , iew , levels )          :: height
   INTEGER , INTENT(IN)                             :: date , time , fdda_loop
   CHARACTER (LEN=24)                               :: date_time_char

   INCLUDE 'error.inc'
   INTERFACE
      INCLUDE 'error.int'
   END INTERFACE

   INTEGER                                          :: i
   
   !  Ingest the observations.  On return, we get the number of observations
   !  in the data set, as well as the observation data structure all filled
   !  filled up.  All data types are ingested at the same time, so the
   !  upper air, surface, ship, aircraft, satellite, and bogus data are
   !  all read at this point.  The data is vertically sorted with this
   !  call.

   CALL read_observations ( obs_filename , unit , obs , number_of_obs , & 
   total_number_of_obs , fatal_if_exceed_max_obs , print_out_obs_found , & 
!BPR BEGIN
!  height , pressure , iew , jns , levels , map_projection )
   height , pressure , slp_x , temperature , iew , jns , levels , map_projection )
!BPR END

   !  There should be at least a single observation to make the rest of
   !  this routine worth anyone's time.

   IF ( number_of_obs .GT. 0 ) THEN

      !  Sort the observations so that they can be merged easily.  Really
      !  what happens is that the index1 array is uniquely sorted by location
      !  (except for observations that are from the same "place").  This 
      !  puts duplicate location observations next to each other.
   
      CALL sort_obs ( obs , number_of_obs , index1 ) 
   
      !  Merge the observations to (try to) remove all duplicates and
      !  build composite data.
   
      CALL check_duplicate_ob ( obs , index1 , number_of_obs , total_dups , date , time ) 
   
      !  The final stage of this procedure is to vertically interpolate 
      !  any missing levels that the analysis would like to have available
      !  for the QC and OA.  The requested levels are the analysis levels.
      !  The only requirements are that there are at least 2 levels of data,
      ! BPR BEGIN
      !!  and that the observation has not already been discarded.  While not
      !! really an interpolation, the data that comes in at a specific level
      !! (SATOB and AIREP), we should put them at the nearest analysis 
      !! pressure surface.
      !   and that the observation has not already been discarded.  
      !  For single-level above surface observations:
      !  1) If the user has not chosen to use a pressure tolerance for using
      !  these obs (i.e., the user has set use_p_tolerance_one_lev = .FALSE.) 
      !  then FM-97 and FM-88 observations will be "adjusted" to the
      !  neareast analysis pressure surface.  Other single-level above-surface
      !  observations will not be used unless they happen to fall exactly on an
      !  analysis pressure level.  
      !  2) If the user has chosen to use a pressure tolerance for using these
      !  obs (i.e., the user has set use_p_tolerance_one_lev = .TRUE. )
      !  then all such obs will be treated the same.  Namely, none will be
      !  "adjusted" to the nearest analysis pressure surface.  Rather, as long
      !  as they are close enough to the nearest analysis pressure surface, they
      !  will be directly compared against / used for that surface 
      ! BPR END
   
      DO i = 1 , number_of_obs
         IF ( ( ASSOCIATED ( obs(i)%surface ) ) .AND. &
              ( .NOT. obs(i)%info%discard ) ) THEN

            !  Make sure the surface level is first in this sounding.  We have already
            !  done this for each individual sounding on input, but this is after the
            !  soundings have been combined.

            CALL surf_first ( obs(i)%surface , obs(i)%info%elevation )

            IF  ( ASSOCIATED ( obs(i)%surface%next ) ) THEN
               CALL interp_all ( pressure , levels , obs(i)%surface , obs(i)%location%longitude ) 
            END IF


            !BPR BEGIN
            !If the user has chosen to use a pressure tolerance for comparing
            !single-level above surface obs to the analysis pressure levels then
            !we do not "adjust" the FM-97 and FM-88 data to the nearest analysis
            !pressure surface
            IF ( .NOT. use_p_tolerance_one_lev ) THEN    
            !BPR END
             IF ( ( obs(i)%info%platform(1:11) .EQ. 'FM-97 AIREP' ) .OR. &
                  ( obs(i)%info%platform(1:11) .EQ. 'FM-88 SATOB' ) ) THEN
                CALL extend_to_closest ( obs(i)%surface , pressure , levels , &
                max_p_extend_t , max_p_extend_w )
             END IF
            !BPR BEGIN
            END IF
            !BPR END
   
         END IF
      END DO
   
      !  The processed observations can now be printed out in a few formats.  The
      !  file useful_out can be ingested by this program.  The file result out
      !  has the same information as useful_out, but with imbedded text for
      !  interpretation.  The file discard_out is the data that has been flagged
      !  as bad by these observation ingesting routines.
   
      CALL make_date ( date , time , date_time_char )

      ! no one really wants these files - so lets just clobber them - April 2009
      !IF ( fdda_loop .EQ. 1 ) THEN
      !   IF ( print_out_files ) THEN 
      !      CALL output_obs ( obs , 2 , 'useful_out_'//date_time_char , number_of_obs ,  1 , .TRUE., 100  )  
      !      CALL output_obs ( obs , 2 , 'result_out_'//date_time_char  , number_of_obs ,  0 , .FALSE., 100 )  
      !      CALL output_obs ( obs , 2 , 'discard_out_'//date_time_char , number_of_obs , -1 , .FALSE., 100 )  
      !   END IF 
      !ELSE
      !   IF ( print_out_files ) THEN 
      !      CALL output_obs ( obs , 2 , 'useful_out_sfc_fdda_'//date_time_char , number_of_obs ,  1 , .TRUE., 100  )  
      !      CALL output_obs ( obs , 2 , 'result_out_sfc_fdda_'//date_time_char  , number_of_obs ,  0 , .FALSE., 100 )  
      !      CALL output_obs ( obs , 2 , 'discard_out_sfc_fdda_'//date_time_char , number_of_obs , -1 , .FALSE., 100 )  
      !   END IF 
      !END IF

   !  We did not find any observations, this could be important.

   ELSE
      error_number        = 0330001
      error_message(1:31) = 'proc_obs_sort                  '
      error_message(32:)  = ' There were no valid observations found within the domain &
                              &for this time period.'
      fatal               = .FALSE.
      listing             = .FALSE.
      CALL error_handler ( error_number , error_message , fatal , listing )
   END IF

END SUBROUTINE proc_obs_sort
