!------------------------------------------------------------------------------

SUBROUTINE proc_ob_access ( request_type , request_variable , print_found_obs , &
request_level , date , time , request_p_diff , request_qc_max , total_numobs , &
request_numobs , obs , &
iew_alloc , jns_alloc , kbu_alloc , &
total_dups , map_projection , &
get_value , get_x_location , get_y_location , get_longitude , &
get_array_index , get_over_water , get_id , get_qc_info , &
!BPR BEGIN
get_elevation , get_fg_3d_t , get_fg_3d_h , &
!BPR END
put_value , put_array_index , put_qc_info )

!  This routine accepts user level requests for modifying the observational
!  data.  These values are stored in a defined type structure, requiring
!  an interface so that the defined type MODULE is not required in all of
!  the quality control (QC) and objective analysis (OA) routines.

!  With the exception of request_* variables, all of the routine arguments are
!  optional.  This allows the program to function as both retrieve and
!  store.  It also allows later development to provide additional 
!  capabilities without compromising the present usage.

   USE map_utils
   USE map_utils_helper
   USE observation
   USE ob_access

   IMPLICIT NONE

   !  Required arguments.  These specify the data requested, and if this is
   !  a retrieve or store.

   CHARACTER ( LEN =   3 )                            :: request_type
   CHARACTER ( LEN =   8 )                            :: request_variable
   LOGICAL                                            :: print_found_obs
   REAL                                               :: request_level  
   INTEGER                                            :: date , time ,  & 
                                                         request_p_diff   , &
                                                         request_qc_max
   INTEGER                                            :: request_numobs   , &
                                                         total_numobs
   TYPE (report) , DIMENSION ( total_numobs )         :: obs
   INTEGER , INTENT ( IN )                            :: iew_alloc , &
                                                         jns_alloc , &
                                                         kbu_alloc , &
                                                         total_dups , &
                                                         map_projection

   !  Optional arguments for retrieving data.

   REAL                    , OPTIONAL , DIMENSION (:) :: get_value , &
                                                         get_x_location , &
                                                         get_y_location , &
                                                         get_longitude
   INTEGER                 , OPTIONAL , DIMENSION (:) :: get_array_index , &
                                                         get_qc_info
   LOGICAL                 , OPTIONAL , DIMENSION (:) :: get_over_water
   CHARACTER ( LEN =   8 ) , OPTIONAL , DIMENSION (:) :: get_id
  !BPR BEGIN
   REAL                    , OPTIONAL , DIMENSION (:) :: get_elevation 
   REAL , INTENT(IN)       , OPTIONAL , DIMENSION(:,:,:) :: get_fg_3d_t, &
                                                            get_fg_3d_h
  !BPR END
   
   !  Optional arguments for storing data.

   REAL                    , OPTIONAL , DIMENSION (:) :: put_value
   INTEGER                 , OPTIONAL , DIMENSION (:) :: put_array_index , &
                                                         put_qc_info

   !  Local variables, loop counters, and such.

   INTEGER                                            :: loop_index , &
                                                         counter    , &
                                                         qc
   REAL                                               :: value
   CHARACTER ( LEN =   4 )                            :: level_message
   
   INCLUDE 'error.inc'

   INTERFACE
      INCLUDE 'error.int'
   END INTERFACE

   !BPR BEGIN
   INTERFACE
      INCLUDE 'query_ob.int'
   END INTERFACE
   !BPR END

   !  There are a couple of simple errors to detect.  Make sure that the surface
   !  is requested for sea level pressure.

   oops_01 : IF ( ( request_variable .EQ. 'PMSL    ' ) .AND. &
                  ( NINT ( request_level ) .NE. 1001 ) ) THEN
      error_number = 99999001
      error_message(1:31) = 'proc_ob_access                 '
      error_message(32:)  = ' For sea level pressure, request pressure &
      &level 1001 mb.'
      fatal = .true.
      listing = .false.
      CALL error_handler ( error_number , error_message , &
      fatal , listing )
   END IF oops_01

   ! BPR BEGIN
   !  Make sure that the surface is requested for mean sea level pressure
   !  derived from surface pressure.

   oops_01a : IF ( ( request_variable .EQ. 'PMSLPSFC' ) .AND. &
                  ( NINT ( request_level ) .NE. 1001 ) ) THEN
      error_number = 99999001
      error_message(1:31) = 'proc_ob_access                 '
      error_message(32:)  = ' For sea level pressure from surface pressure, request pressure &
      &level 1001 mb.'
      fatal = .true.
      listing = .false.
      CALL error_handler ( error_number , error_message , &
      fatal , listing )
   END IF oops_01a

   !  Make sure that the surface is requested for surface pressure

   oops_01b : IF ( ( request_variable .EQ. 'PSFC    ' ) .AND. &
                  ( NINT ( request_level ) .NE. 1001 ) ) THEN
      error_number = 99999001
      error_message(1:31) = 'proc_ob_access                 '
      error_message(32:)  = ' For surface pressure, request pressure &
      &level 1001 mb.'
      fatal = .true.
      listing = .false.
      CALL error_handler ( error_number , error_message , &
      fatal , listing )
   END IF oops_01b
   ! BPR END

   !  For QC, we are either retrieving data from the data structure (get), or we are 
   !  storing it there (put).  For OA, we are requesting the data (use), but a bit 
   !  differently from (get).  This routine is separated into those three pieces 
   !  in the following CASE statement.

   get_use_put : SELECT CASE ( request_type )

      CASE ( 'get' , 'GET' ) get_use_put
         
         !  Initialize these two counters.  counter is the array index of the 
         !  returned data, request_numobs is the total number of returned
         !  data values.  These are incremented below.

         request_numobs = 0
         counter = 1 

         !  Loop over each of the reports.

         all_reports_get : DO loop_index = 1 , total_numobs

            !  All of the data we handle are either in the surface linked list
            !  OR it is the sea level pressure.  If this field is not sea level
            !  pressure, we do not need to consider observations without any
            !  surface ASSOCIATED.

            no_upper_air_get : IF ( ( request_variable .NE. 'PMSL    ' )    .AND. & 
                 ( .NOT. ( ASSOCIATED ( obs(loop_index)%surface )  )  )  ) THEN
               CYCLE all_reports_get
            END IF no_upper_air_get

            !  This observation may have the discard flag set.  If so, then
            !  we bypass it.

            discard_data_get : IF ( obs(loop_index)%info%discard ) THEN
               CYCLE all_reports_get
            END IF discard_data_get

            !  Since this is QC only, we can skip any observation with the
            !  bogus flag set.

            bogus_data_get : IF ( obs(loop_index)%info%bogus ) THEN
               CYCLE all_reports_get
            END IF bogus_data_get            

            !  Initialize a flag to missing prior to every call for "success"
            !  detect afterwards.
     
            qc = missing

            !  See if the next report in the array satifies the requested
            !  criteria.  

            CALL query_ob ( obs(loop_index) , date , time , &
            request_variable , request_level , request_qc_max , request_p_diff , &
!BPR BEGIN
!           value , qc )
            value , qc , fg_3d_h = get_fg_3d_h, fg_3d_t = get_fg_3d_t )
!BPR END

            !  If qc is returned .NE. to missing, then this routine found a 
            !  valid observation.

            IF ( qc .NE. missing ) THEN
               IF ( PRESENT ( get_value ) ) get_value(counter) = value
               IF ( PRESENT ( get_x_location ) .AND. PRESENT ( get_y_location ) ) THEN
                  CALL latlon_to_ij(projd , &
                                    obs(loop_index)%location%latitude  , &
                                    obs(loop_index)%location%longitude , &
                                    get_x_location(counter)            , &
                                    get_y_location(counter))
               END IF
               IF ( PRESENT ( get_longitude ) ) &
               get_longitude(counter) = obs(loop_index)%location%longitude
               IF ( PRESENT ( get_array_index ) ) get_array_index(counter) = loop_index
               IF ( PRESENT ( get_qc_info ) ) get_qc_info(counter) = qc
               IF ( PRESENT ( get_id ) ) &
               get_id(counter) = obs(loop_index)%location%id(1:8)

               !BPR BEGIN
               IF ( PRESENT ( get_elevation ) ) get_elevation(counter) = obs(loop_index)%info%elevation 
               !BPR END

               !  Is this observation located in our domain?  We can 
               !  check easily if this has x,y available.

               have_x_and_y_get: IF ( PRESENT ( get_x_location ) .AND. PRESENT ( get_y_location ) ) THEN
                  IF ( ( ( get_x_location(counter) .GE.         1   ) .AND. &
                         ( get_x_location(counter) .LE. iew_alloc   ) .AND. &
                         ( get_y_location(counter) .GE.         1   ) .AND. &
                         ( get_y_location(counter) .LE. jns_alloc   ) .AND. &
                       ( ( request_variable .EQ. 'U       ' ) .OR. ( request_variable .EQ. 'V       ' ) ) ) &
                                                     .OR. &  
                       ( ( get_x_location(counter) .GE.         2   ) .AND. &
                         ( get_x_location(counter) .LE. iew_alloc-1 ) .AND. &
                         ( get_y_location(counter) .GE.         2   ) .AND. &
                         ( get_y_location(counter) .LE. jns_alloc-1 ) ) ) THEN
                     request_numobs = counter
                     counter = counter + 1
                  ELSE
                     get_qc_info(counter) = outside_of_domain
                     obs(loop_index)%info%discard = .TRUE.
                  END IF
               END IF have_x_and_y_get
            END IF

         END DO all_reports_get

         !  There is an important piece of information that we now have, so share it
         !  with the user: total number of observations that will influence the 
         !  the analysis

         IF ( print_found_obs ) THEN
            WRITE ( UNIT = * , FMT = '( // "Of the ",i5," observations found in the &
            &observation input file, after QC, only ",i5," of them are available inside &
            &the domain." /&
            &"This field is ",a8," for pressure level ",i4," hPa (1001 hPa is defined as the surface)."//)' ) &
            total_numobs - total_dups , request_numobs , request_variable , NINT ( request_level )
         END IF

         oops_02 : IF ( ( request_numobs .EQ. 0 ) .AND. ( print_found_obs ) ) THEN
            error_number = 99999002
            error_message(1:31) = 'proc_ob_access                 '
            WRITE ( level_message  , FMT = '(i4)' ) NINT ( request_level )
            error_message(32:)  = ' We found no QC observations for ' // request_variable &
            // ' at level = ' // level_message // ' mb.'
            fatal = .false.
            listing = .false.
            CALL error_handler ( error_number , error_message , &
            fatal , listing )
         END IF oops_02

      CASE ( 'use' , 'USE' ) get_use_put
         
         !  Initialize these two counters.  counter is the array index of the 
         !  returned data, request_numobs is the total number of returned
         !  data values.  These are incremented below.

         request_numobs = 0
         counter = 1 

         !  Loop over each of the reports.

         all_reports_use : DO loop_index = 1 , total_numobs

            !  All of the data we handle is either in the surface linked list
            !  OR it is the sea level pressure.  If this field is not sea level
            !  pressure, we do not need to consider observations without any
            !  surface ASSOCIATED.

            no_upper_air_use : IF ( ( request_variable .NE. 'PMSL    ' )    .AND. & 
                 ( .NOT. ( ASSOCIATED ( obs(loop_index)%surface )  )  )  ) THEN
               CYCLE all_reports_use
            END IF no_upper_air_use

            !  This observation may have the discard flag set.  If so, then
            !  we bypass it.

            discard_data_use : IF ( obs(loop_index)%info%discard ) THEN
               CYCLE all_reports_use
            END IF discard_data_use            

            !  Initialize a flag to missing prior to every call for "success"
            !  detect afterwards.
     
            qc = missing

            !  See if the next report in the array satifies the requested
            !  criteria.  

            CALL query_ob ( obs(loop_index) , date , time , &
            request_variable , request_level , request_qc_max , request_p_diff , &
            value , qc )

            !  If qc is returned .NE. to missing, then this routine found a 
            !  valid observation

            IF ( qc .NE. missing ) THEN
               IF ( PRESENT ( get_value ) ) get_value(counter) = value
               IF ( PRESENT ( get_x_location ) .AND. PRESENT ( get_y_location ) ) THEN
                  CALL latlon_to_ij(projd , &
                                    obs(loop_index)%location%latitude  , &
                                    obs(loop_index)%location%longitude , &
                                    get_x_location(counter)            , &
                                    get_y_location(counter))
               END IF
               IF ( PRESENT ( get_longitude ) ) &
               get_longitude(counter) = obs(loop_index)%location%longitude
               IF ( PRESENT ( get_array_index ) ) get_array_index(counter) = loop_index
               IF ( PRESENT ( get_qc_info ) ) get_qc_info(counter) = qc
               IF ( PRESENT ( get_id ) ) &
               get_id(counter) = obs(loop_index)%location%id(1:8)

               !  Is this observation located in our domain?  We can 
               !  check easily if this has x,y available.

               have_x_and_y_use: IF ( PRESENT ( get_x_location ) .AND. PRESENT ( get_y_location ) ) THEN

                  IF ( ( ( get_x_location(counter) .GE.         1   ) .AND. &
                         ( get_x_location(counter) .LE. iew_alloc   ) .AND. &
                         ( get_y_location(counter) .GE.         1   ) .AND. &
                         ( get_y_location(counter) .LE. jns_alloc   ) .AND. &
                       ( ( request_variable .EQ. 'U       ' ) .OR. ( request_variable .EQ. 'V       ' ) ) ) &
                                                     .OR. &  
                       ( ( get_x_location(counter) .GE.         2   ) .AND. &
                         ( get_x_location(counter) .LE. iew_alloc-1 ) .AND. &
                         ( get_y_location(counter) .GE.         2   ) .AND. &
                         ( get_y_location(counter) .LE. jns_alloc-1 ) ) ) THEN
                     request_numobs = counter
                     counter = counter + 1
                  ELSE
                     get_qc_info(counter) = outside_of_domain
                     obs(loop_index)%info%discard = .TRUE.
                  END IF
               END IF have_x_and_y_use
            END IF

         END DO all_reports_use

         !  There is an important piece of information that we now have, so share it
         !  with the user: total number of observations that will influence the 
         !  the analysis

         IF ( print_found_obs ) THEN
            WRITE ( UNIT = * , FMT = '( // "Of the ",i5," observations found in the &
            &observation input file, for OA, only ",i5," of them are available inside &
            &the domain." /&
            &"This field is ",a8," for pressure level ",i4," hPa (1001 hPa is defined as the surface)."//)' ) &
            total_numobs , request_numobs , request_variable , NINT ( request_level )
         END IF

         oops_03 : IF ( ( request_numobs .EQ. 0 ) .AND. ( print_found_obs ) ) THEN
            error_number = 99999003
            error_message(1:31) = 'proc_ob_access                 '
            WRITE ( level_message  , FMT = '(i4)' ) NINT ( request_level )
            error_message(32:)  = ' We found no OA observations for ' // request_variable &
            // ' at level = ' // level_message // ' mb.'
            fatal = .false.
            listing = .false.
            CALL error_handler ( error_number , error_message , &
            fatal , listing )
         END IF oops_03

      CASE ( 'put' , 'PUT' ) get_use_put

         !  Loop over all of the reports that we will replace.  Only
         !  the "chosen" indices of the original observations
         !  are required.  

         specific_ones : DO loop_index = 1 , request_numobs
    
            !  We know that there are request_numobs valid observations,
            !  and their location is given as the index into the
            !  "put_array_index" array.

            counter = put_array_index(loop_index)

            !  With the array of observations, with known variable
            !  name and level, we store the value and QC information 
            !  back into the original observation database.

            CALL store_ob ( obs(counter) , &
            request_variable , request_level , request_p_diff , & 
            put_value(loop_index) , put_qc_info(loop_index) ) 

         END DO specific_ones

      CASE DEFAULT get_use_put

         error_number = 99999003
         error_message(1:31) = 'proc_ob_access                 '
         error_message(32:)  = ' The selection of ' // request_type // &
         ' is neither *put* nor *get*, the only available options.'
         fatal = .true.
         listing = .false.
         CALL error_handler ( error_number , error_message , &
         fatal , listing )

   END SELECT get_use_put

END SUBROUTINE proc_ob_access
