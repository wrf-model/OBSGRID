!---------------------------------------------------------------------------

!  MODULE: observation

!  DESCRIPTION:  Contains definitions of data structures used to hold data
!  of weather observations.  Contains routines for reading data into data
!  structures, for reading and sorting sounding data in order of increasing
!  height, for sorting observation by location, for deleting or merging
!  duplicate location/time observations.

!  The first section is the PARAMETERS.  These are some handy constants that
!  this module will use heavily.  The next section is the data structure of
!  the observations, which is used to read, store and retrieve data.  The
!  third section describes a couple of interfaces for the logical operators
!  both .EQ. and .LT.) which point to a couple of subsequent function in the
!  CONTAINS fourth section.  The CONTAINS section has an alphabetical list
!  of included routines.

!-------------------------------------------------------------------------

MODULE observation

!--------------------------------------------------------------------------
!                             PARAMETERS
!---------------------------------------------------------------------------

   INCLUDE 'missing.inc'

   !  Define error return codes used by 'read_measurements' routine.

   INTEGER , PARAMETER                            ::  ok       = 0 , &
                                                      eof_err  = 1 , &
                                                      no_data  = 2 , &
                                                      read_err = 3

!  FORMAT STRINGS for input/output of data.
!  These format strings correspond to the data structures in this file
!  and are used for input and output of values in the 'report' structure
!  (first format string) and the 'measurement' structure (second format).
!  Note that report struct contains the first of a linked list of 
!  measurements; this first meas is read using the 'meas_format'.

   CHARACTER ( LEN = 120 ) , PARAMETER :: rpt_format =  &
                ' ( 2f20.5 , 2a40 , ' &             ! format for location_type
             // ' 2a40 , 1f20.5 , 5i10 , 3L10 , ' & ! format for source_info
             // ' 2i10 , a20 , ' &                  ! fmt for valid_time
             // ' 13( f13.5 , i7 ) ) '              ! fmt for 'terrestrial' 

   CHARACTER ( LEN = 120 ) , PARAMETER :: meas_format = & 
                ' ( 10( f13.5 , i7 ) ) '            ! fmt for measurement rcd

   CHARACTER ( LEN = 120 ) , PARAMETER :: end_format = &
                ' ( 3 ( i7 ) ) '                    ! fmt for end record

!
!-------------------------------------------------------------------------
!                          DATA STRUCTURES
!--------------------------------------------------------------------------

   !  These data structures are built to hold all of the required information 
   !  from a single level report.  This includes, but is not limited to,
   !  land based surface observations, ocean based surface observations,
   !  and aircraft data.  All data is assumed to be given a horizontal
   !  location (lat/lon) and a timestamp.  For the data to be useful, there
   !  needs to be a vertical location presribed, or derivable.  For data with
   !  multiple vertical observations, the 'meas' structure is also required.
   !  This includes, but is not limited to, sounding data, satellite derived
   !  winds and satellite derived thickness.

   !  The information in the following two records is usually available 
   !  in every report.  

!-------------------------------------------------------------------------

   TYPE location_type

      !  The fields in this record uniquely identify the source of the 
      !  data, so that duplicates of data can be merged or discarded.
      !  The horizontal location of this report (assumed constant, even
      !  for balloon ascent) is geven by the lat/lon of the site.
      
      REAL                   :: latitude  , &   ! latitude (+ degrees east)
                                longitude       ! longitude (+ degrees north)

      CHARACTER ( LEN = 40 ) :: id , &          ! 5 digit identifier, 
                                                ! consisting of a 2 digit block 
                                                ! number and a 3 digit 
                                                ! identifier (for soundings)
                                                ! for WMO sttn; non digit
                                                ! for other sources
                                name            ! The name corresponds to
                                                ! the id (is obtained from id
                                                ! in the program that is 
                                                ! source of data
   END TYPE location_type


!---------------------------------------------------------------------------

   TYPE source_info

      CHARACTER ( LEN = 40 ) :: platform , &    ! description of the 
                                                ! measurement device
                                source          ! GTS data, NCAR ADP files, 
                                                ! bogus information, etc
      REAL                   :: elevation       ! station elevation

      !  During the decoding process, how many errors (non conforming
      !  codes) were encountered, and how many warnings (this is a subjective
      !  call based on repeated incorrect -- but probably not garbled --
      !  GTS codes).  If a bulletin is completely garbled, the logical
      !  flag to not consider this report is set.

      INTEGER              :: num_vld_fld , & ! number of valid fields in the
                                              ! entire report; used as the
                                              ! first tie-breaker in deciding
                                              ! which conflicting data items
                                              ! to keep if have duplicate rpts
                              num_error , &   ! number of errors 
                                              ! encountered during the
                                              ! decoding process
                              num_warning , & ! number of warnings 
                                              ! encountered during the 
                                              ! decoding process
                              seq_num , &     ! sequence numbers that tell
                                              ! which of 2 reports is more
                                              ! recent.
                              num_dups        ! number of duplicates found of
                                              ! this observation 
      LOGICAL              :: is_sound        ! is-a-sounding tells whether
                                              ! the observation possibly has
                                              ! multiple levels vs having 
                                              ! only one level for srfc ob.
      LOGICAL              :: bogus           ! T/F if this is a bogus 
                                              ! observation
      LOGICAL              :: discard         ! Tells whether this observation
                                              ! has been found to be a dup
                                              ! AND has been discarded or
                                              ! merged.
   END TYPE source_info

!--------------------------------------------------------------------------

   TYPE field

      !  Defines a data type consisting of a paired data value (real) with a
      !  quality control flag that holds a binary-coded combination of error
      !  codes; the codes  identify possible problems with the data.

      REAL                   :: data
      INTEGER                :: qc              !  Quality control flags
                                                !  that are 0 if data is
                                                !  good, or different 
                                                !  integers depending upon
                                                !  what error(s) occured
   END TYPE field

!-------------------------------------------------------------------------

   TYPE terrestrial

      !  The data that will occur, at most, once during a report is 
      !  listed here.  These are typically terrestrial measured values.  The
      !  surface met fields are stored in a separate TYPE, to allow a 
      !  POINTER to the next level (if one exists).  This needs to be a 
      !  separate TYPE so as to allow a POINTER to it 

      TYPE ( field )         :: slp       , &   ! sea level pressure
                                ref_pres  , &   ! reference pres level for
                                                ! the thickness
                                ground_t  , &   ! ground temperature
                                sst       , &   ! sea surface temperature
                                psfc      , &   ! surface pressure
                                precip    , &   ! precipitation accumulation
                                t_max     , &   ! daily temperature max
                                t_min     , &   ! daily temperature min
                                t_min_night , & ! min overnight temperature
                                p_tend03  , &   ! pressure tendency in 3hr
                                p_tend24  , &   ! pressure tendency in 24hr
                                cloud_cvr , &   ! total cloud cover (oktas)
                                ceiling         ! height of lowest cloud base
   END TYPE terrestrial 

!-------------------------------------------------------------------------

   TYPE time_info

      !  GTS report time: the valid time of the report.  The largest INTEGER values 
      !  require only 8 digits, so that this should function properly with 
      !  32-bit INTEGERS.  

      INTEGER                :: sut      , &    ! number of seconds since 1 Jan
                                                ! 0000 UTC 1970
                                julian          ! Julian day
      CHARACTER ( LEN = 14 )    date_char       ! CCYYMMDDHHmmss date

   END TYPE time_info

!--------------------------------------------------------------------------

   TYPE meas_data

      !  The met data involved with this program is defined in this TYPE.  The
      !  standard state variables (wind, temperature, moisture, with pressure
      !  and/or height to fix the vertical location) are stored here.  For 
      !  single level observations, only one of these records is used per     
      !  observation.   For multi-level reports, a linked list of these 
      !  measurement TYPEs is generated.

      TYPE ( field )         :: pressure    , & ! pressure of observation
                                height      , & ! height (above sea level) 
                                temperature , & ! 
                                dew_point   , & ! 
                                speed       , & ! 
                                direction   , & ! 
                                u           , & ! u and v components of wind
                                v           , & ! are derived from spd and dir
                                rh          , & !
                                thickness       ! 

   END TYPE meas_data

!--------------------------------------------------------------------------

   TYPE measurement

      TYPE ( meas_data )               :: meas  ! contains data and qc code
      TYPE ( measurement ) ,  POINTER  :: next  ! the met info is handled
                                                ! as a linked list of the  
                                                ! measurement type

   END TYPE measurement

!-------------------------------------------------------------------------

   TYPE report                                 
                                               ! this is the entire report
      TYPE ( location_type ) :: location       ! for a single time, from a 
      TYPE ( source_info )   :: info           ! single reporting platform,
      TYPE ( time_info )     :: valid_time     ! a sounding, surface, buoy,
      TYPE ( terrestrial )   :: ground         ! aircraft or ship report
      TYPE ( measurement ) , &
               POINTER       :: surface        

   END TYPE report                            
                                             
!
!-------------------------------------------------------------------------
!                            INTERFACES
!------------------------------------------------------------------------

!INTERFACE OPERATOR ( .LT. )
!  !  definition of .LT. operator for use with 'report' data type
!  MODULE PROCEDURE compare           
!END INTERFACE 

!INTERFACE OPERATOR ( .EQ. )
!   !  definition of .EQ. operator for use with field, location, 
!   !  ground, and time
!   MODULE PROCEDURE field_eq , loc_eq , time_eq , ground_eq
!END INTERFACE


!
CONTAINS   

! -------------------------------------------------------------------------
!                            FUNCTIONS
! ---------------------------------------------------------------------------

LOGICAL FUNCTION compare ( a , b )

!  This defines the comparison operator '.LT.' for use with the 'report'
!  data type.  NOTE that the other operators LE, GE, GT are NOT
!  defined at all for the 'report' data type.

   IMPLICIT NONE

   TYPE ( report ) , INTENT ( IN )     :: a  ! the first data item compared
   TYPE ( report ) , INTENT ( IN )     :: b  ! the second data item compared

   compare = .FALSE.
   IF ( a%location%longitude .LT. b%location%longitude ) THEN
      compare = .TRUE.
   ELSE IF ( a%location%longitude .eq. b%location%longitude ) THEN
      IF ( a%location%latitude .LT. b%location%latitude ) THEN
         compare = .TRUE. 
      ELSE IF ( a%location%latitude .EQ. b%location%latitude ) THEN
         IF ( LLT ( a%location%id , b%location%id ) ) THEN
            compare = .TRUE.
         ELSE IF ( a%location%id .EQ. b%location%id ) THEN
            IF ( LLT ( a%location%name , b%location%name ) ) THEN
               compare = .TRUE.
            END IF
         END IF
      END IF
   END IF

END FUNCTION compare

!
! -------------------------------------------------------------------------

LOGICAL FUNCTION eps_equal ( a , b , eps )

!  Compare two real numbers a and b, and return TRUE if they are within
!  parameter 'eps' of one another.  

   IMPLICIT NONE 

   REAL , INTENT ( IN )                     :: a , b , eps

   IF ( ABS ( a - b ) .LT. eps ) THEN
      eps_equal = .TRUE.
   ELSE
      eps_equal = .FALSE.
   END IF

END FUNCTION eps_equal
      
!
! ------------------------------------------------------------------------

LOGICAL FUNCTION field_eq ( a , b )

! This defines operator .EQ. for 'field' data type

   IMPLICIT NONE 

   TYPE ( field ) , INTENT ( IN )                :: a , b
   
   IF ( a%data .EQ. b%data .AND. a%qc .EQ. b%qc ) THEN
      field_eq = .TRUE.
   ELSE
      field_eq = .FALSE.
   END IF
  
END FUNCTION field_eq

!
! -------------------------------------------------------------------------

LOGICAL FUNCTION ground_eq ( a , b )

! This defines operator .EQ. for 'terrestrial' data type

   IMPLICIT NONE 

   TYPE ( terrestrial ) , INTENT ( IN )    :: a , b

   IF ( eps_equal ( a%slp%data         , b%slp%data         , .01 ) .AND. &
        eps_equal ( a%ref_pres%data    , b%ref_pres%data    , .01 ) .AND. &
        eps_equal ( a%ground_t%data    , b%ground_t%data    , .01 ) .AND. &
        eps_equal ( a%sst%data         , b%sst%data         , .01 ) .AND. &
        eps_equal ( a%psfc%data        , b%psfc%data        , .01 ) .AND. &
        eps_equal ( a%precip%data      , b%precip%data      , .01 ) .AND. &
        eps_equal ( a%t_max%data       , b%t_max%data       , .01 ) .AND. &
        eps_equal ( a%t_min%data       , b%t_min%data       , .01 ) .AND. &
        eps_equal ( a%t_min_night%data , b%t_min_night%data , .01 ) .AND. &
        eps_equal ( a%p_tend03%data    , b%p_tend03%data    , .01 ) .AND. &
        eps_equal ( a%p_tend24%data    , b%p_tend24%data    , .01 ) .AND. &
        eps_equal ( a%cloud_cvr%data   , b%cloud_cvr%data   , .01 ) .AND. &
        eps_equal ( a%ceiling%data     , b%ceiling%data     , .01 ) .AND. &
        a%slp%qc         .EQ. b%slp%qc         .AND. &
        a%ref_pres%qc    .EQ. b%ref_pres%qc    .AND. &
        a%ground_t%qc    .EQ. b%ground_t%qc    .AND. &
        a%sst%qc         .EQ. b%sst%qc         .AND. &
        a%psfc%qc        .EQ. b%psfc%qc        .AND. &
        a%precip%qc      .EQ. b%precip%qc      .AND. &
        a%t_max%qc       .EQ. b%t_max%qc       .AND. &
        a%t_min%qc       .EQ. b%t_min%qc       .AND. &
        a%t_min_night%qc .EQ. b%t_min_night%qc .AND. &
        a%p_tend03%qc    .EQ. b%p_tend03%qc    .AND. &
        a%p_tend24%qc    .EQ. b%p_tend24%qc    .AND. &
        a%cloud_cvr%qc   .EQ. b%cloud_cvr%qc   .AND. &
        a%ceiling%qc     .EQ. b%ceiling%qc           ) THEN
      ground_eq = .TRUE.
   ELSE
      ground_eq = .FALSE.
   END IF
  
END FUNCTION ground_eq

!
! -------------------------------------------------------------------------

LOGICAL FUNCTION loc_eq ( a , b )

! This defines operator .EQ. for 'location' data type

   IMPLICIT NONE 

   TYPE ( location_type ) , INTENT ( IN )           :: a , b
   
   IF ( eps_equal ( a%latitude , b%latitude , .01 ) .AND. &
        eps_equal ( a%longitude , b%longitude , .01 ) .AND. &
        a%id .EQ. b%id .AND. a%name .EQ. b%name ) THEN
      loc_eq = .TRUE.
   ELSE
      loc_eq = .FALSE.
   END IF
  
END FUNCTION loc_eq

!
! -------------------------------------------------------------------------

LOGICAL FUNCTION time_eq_old ( a , b )

! This defines operator .EQ. for 'time_info' data type

   IMPLICIT NONE 

   TYPE ( time_info ) , INTENT ( IN )           :: a , b
   
   IF ( ( a%sut       .EQ. b%sut       ) .AND. &
        ( a%julian    .EQ. b%julian    ) .AND. &
        ( a%date_char .EQ. b%date_char ) ) THEN
      time_eq_old = .TRUE.
   ELSE
      time_eq_old = .FALSE.
   END IF
  
END FUNCTION time_eq_old

!
! -------------------------------------------------------------------------

LOGICAL FUNCTION time_eq ( a , b , date , time )

! This defines operator .EQ. for 'time_info' data type

   USE date_pack

   IMPLICIT NONE 

   TYPE ( time_info ) , INTENT ( INOUT )        :: a , b
   INTEGER            , INTENT ( IN )           :: date , time 

   !  Local variables.

   CHARACTER (LEN=19)               :: target_date , a_date , b_date
   INTEGER                          :: diff_seconds , a_diff_seconds , b_diff_seconds

   !  Compute the character string date and time for the current analysis time.

   WRITE ( target_date , '(I4.4,"-",I2.2,"-",I2.2,"_",I2.2,":",I2.2,":",I2.2)' ) &
   date / 10000 ,  &
   ( date - (date / 10000 ) * 10000 ) / 100 , &
   date - ( date / 100 ) * 100 , &
   time / 10000 , &
   ( time - ( time / 10000 ) * 10000 ) / 100, &
   time - ( time / 100 ) * 100

   !  Get the date/time for observations a and b in a YYYY-MM-DD_HH:mm:ss format.

   a_date( 1: 5) = a%date_char( 1: 4) // '-'
   a_date( 6: 8) = a%date_char( 5: 6) // '-'
   a_date( 9:11) = a%date_char( 7: 8) // '_'
   a_date(12:14) = a%date_char( 9:10) // ':'
   a_date(15:17) = a%date_char(11:12) // ':'
   a_date(18:19) = a%date_char(13:14)

   b_date( 1: 5) = b%date_char( 1: 4) // '-'
   b_date( 6: 8) = b%date_char( 5: 6) // '-'
   b_date( 9:11) = b%date_char( 7: 8) // '_'
   b_date(12:14) = b%date_char( 9:10) // ':'
   b_date(15:17) = b%date_char(11:12) // ':'
   b_date(18:19) = b%date_char(13:14)

   !  Compute the time difference between the two observations in seconds.

   CALL geth_idts ( a_date , b_date , diff_seconds )
   
   !  If the times (a and b) are within half an hour of each other, we say that they
   !  are the same time.  

   IF ( ABS ( diff_seconds ) .LT. 1800 ) THEN

      !  Now that we know a and b are the same time, the important question is now
      !  are they the time that we want?  If they are the same (which means either
      !  a or b is within an hour of the target time), we set both of these times
      !  to the target time.

      CALL geth_idts ( target_date , a_date , a_diff_seconds )
      CALL geth_idts ( target_date , b_date , b_diff_seconds )

      IF ( ( ABS ( a_diff_seconds ) .LT. 3600 ) .OR. &
           ( ABS ( b_diff_seconds ) .LT. 3600 ) ) THEN

         a%date_char( 1: 4) = target_date( 1: 4)
         a%date_char( 5: 6) = target_date( 6: 7) 
         a%date_char( 7: 8) = target_date( 9:10) 
         a%date_char( 9:10) = target_date(12:13) 
         a%date_char(11:12) = target_date(15:16) 
         a%date_char(13:14) = target_date(18:19) 

         b%date_char( 1: 4) = target_date( 1: 4)
         b%date_char( 5: 6) = target_date( 6: 7) 
         b%date_char( 7: 8) = target_date( 9:10) 
         b%date_char( 9:10) = target_date(12:13) 
         b%date_char(11:12) = target_date(15:16) 
         b%date_char(13:14) = target_date(18:19) 
      END IF

      time_eq = .TRUE.

   ELSE

      time_eq = .FALSE.

   END IF
  
END FUNCTION time_eq

!
! -------------------------------------------------------------------------
!                            ROUTINES
! -------------------------------------------------------------------------

SUBROUTINE check_duplicate_ob ( obs , index , num_obs , total_dups , date , time )

!  Checks array of reports (obs), which has a sorted index to the reports,
!  to determine if any reports are for the same time/location.  If so,
!  and the data is duplicated exactly in all fields, one is discarded.  If
!  they are from same time/location and data is not identical, data from
!  two reports is merged:  'missing' is replaced by known values; data at
!  different levels is merged into one linked list.

   IMPLICIT NONE

   TYPE ( report ) , INTENT ( INOUT ) , DIMENSION ( : ) :: obs   ! array of observations
   INTEGER         , INTENT ( IN )    , DIMENSION ( : ) :: index ! gives sort order of obs
   INTEGER         , INTENT ( IN )                      :: num_obs 

   INTEGER                                :: current , &
                                             next    , & 
                                             first   , &
                                             second
   INTEGER         , INTENT ( OUT )       :: total_dups
   INTEGER         , INTENT ( IN  )       :: date    , &
                                             time

   INCLUDE 'error.inc'
   INTERFACE
      INCLUDE 'error.int'
   END INTERFACE

   !  Count the total number of duplicate reports.

   total_dups = 0

   obsloop: DO current = 1 , num_obs - 1

      first = index(current)

      !  If this obs has been merged with another obs or discarded, skip it.

      IF ( obs(first)%info%discard ) THEN
         CYCLE obsloop
      END IF

      !  Get second obs to compare with first; compare first obs to second obs 
      !  until next obs does not match.

      compare: DO next = current + 1 , num_obs

         second = index(next)

         ! Sorted by location, so if locations NE, then no chance of any
         ! more matches with first.

! foo
!        IF ( .NOT. ( obs(first)%location .EQ. obs(second)%location ) ) THEN
         IF ( .NOT. loc_eq ( obs(first)%location , obs(second)%location ) ) THEN
            CYCLE obsloop
         END IF

         !  If this obs has been merged with another obs or discarded, skip it.

         IF ( obs(second)%info%discard ) THEN
            CYCLE compare
         END IF

         !  If time fields are not completely identical, go to next observation.
         !  Sort is by location ONLY, not by time; so next+1 may be identical
         !  even though next has different time.

! foo
!        IF ( .NOT. ( obs(first)%valid_time .EQ. obs(second)%valid_time ) ) THEN
         IF ( .NOT. time_eq ( obs(first)%valid_time , obs(second)%valid_time , date , time ) ) THEN
            error_number = 0332001
            error_message(1:31) = 'check_duplicate_ob             '
            error_message(32:)  = ' Found multiple times for ' &
            // TRIM ( obs(first)%location%id ) // ' ' &
            // TRIM ( obs(first)%location%name ) // ', ' &
            // TRIM ( obs(first)%valid_time%date_char )  // ' and ' &
            // TRIM ( obs(second)%valid_time%date_char ) // '.'
            fatal = .false.
            listing = .false.
!           CALL error_handler ( error_number , error_message , &
!           fatal , listing )
            CYCLE compare
         END IF

         !  Observations are from same location and time, so merge them.

         CALL merge_obs ( obs(first) , obs(second) )

         !  Mark second of pair as discarded; data is put in 'first'.  Note that
         !  a duplicate has been found by incrementing the counter.

         obs(second)%info%discard  = .true.
         obs(first)%info%num_dups  = obs(first)%info%num_dups + 1
         total_dups = total_dups + 1

         !  Free up the space for the observation report that is discarded.
         !  Unfortunately, OR NOT!  

! foo
!        CALL dealloc_meas ( obs(second)%surface ) 
         NULLIFY ( obs(second)%surface ) 

      END DO compare

   END DO obsloop

   WRITE ( UNIT = * , FMT = '( "Of the ",i5," observations reported, ", &
   &i5," of them are duplicates, leaving ",i5," merged locations" )' )  &
   num_obs , total_dups , num_obs - total_dups

END SUBROUTINE check_duplicate_ob

!
! -----------------------------------------------------------------------

SUBROUTINE dealloc_meas ( head )

!  This deallocates all nodes in a linked list of measurements.

   IMPLICIT NONE 

   TYPE ( measurement ) , POINTER           :: head     ! head of linked list

   TYPE ( measurement ) , POINTER           :: previous &
                                             , temp
   INTEGER                                  :: status

   !  Start at the head, kill everything that is pointed to.  After no longer 
   !  associated, deallocate the head.

   IF ( ASSOCIATED ( head ) ) THEN

      previous => head
      list_loop : DO WHILE ( ASSOCIATED ( previous%next ) )
         temp => previous
         previous => previous%next
         DEALLOCATE ( temp , STAT = status) 
         IF ( status .NE. 0 ) THEN
            PRINT *,'Error in DEALLOCATE, continuing by stopping DEALLOCATE on this list.'
            EXIT list_loop
         END IF
      END DO list_loop

   END IF

!  NULLIFY ( head ) 
  
END SUBROUTINE dealloc_meas

!
!------------------------------------------------------------------------------

SUBROUTINE diagnostics ( new , longitude )

!  This routine computes the derived diagnostic fields that are to be
!  objectively analyzed (RH from T and Td; u and v (rotated to map) from
!  speed and direction).

   IMPLICIT NONE

   TYPE ( meas_data )    :: new
   REAL                  :: longitude

   REAL                  :: t , td

   REAL                  :: dir_map , lon_dif
                            
   INTEGER               :: iew_map , jns_map

   REAL , PARAMETER :: L_over_Rv = 5418.12
   REAL , PARAMETER :: piover180 = 3.14159265358 / 180.

   INCLUDE 'proc_get_info_header.inc'
   INCLUDE 'map.inc'
   INCLUDE 'error.inc'
   INTERFACE
      INCLUDE 'error.int'
   END INTERFACE

   !  Given a valid temperature and dew point, the relative humidity is
   !  computed.  The relative humidity QC flag is arbitrarily set to the
   !  temperature flag value.

   IF      ( (       eps_equal ( new%temperature%data , 0.        , 1. ) ) .OR. &
             (       eps_equal ( new%dew_point%data   , 0.        , 1. ) ) .OR. &
             (                   new%temperature%data .LT.          0.   ) .OR. &
             (                   new%dew_point%data   .LT.          0.   ) ) THEN
      new%rh%data = missing_r
      new%rh%qc   = zero_t_td
   ELSE IF ( ( .NOT. eps_equal ( new%temperature%data , missing_r , 1. ) ) .AND. &
             ( .NOT. eps_equal ( new%dew_point%data   , missing_r , 1. ) ) ) THEN
      t  = new%temperature%data
      td = new%dew_point%data
      new%rh%data = 100. * exp ( L_over_Rv * (1./T - 1./Td) )
      new%rh%qc   = new%temperature%qc
   END IF

   !  The wind components are computed in two steps.  First we need the
   !  meteorological u and v from the speed and direction.  Second, those
   !  values of u and v are rotated to the model grid.

   IF      ( (       eps_equal ( new%speed%data       , missing_r , 1. ) ) .OR. &
             (       eps_equal ( new%direction%data   , missing_r , 1. ) ) ) THEN
      new%u%data = missing_r
      new%u%qc   = missing
      new%v%data = missing_r
      new%v%qc   = missing
   ELSE IF ( (       eps_equal ( new%speed%data , 0.        , 0.1 ) ) .AND. &
             (       eps_equal ( new%direction%data , 0.    , 0.1 ) ) ) THEN
      new%u%data = missing_r
      new%u%qc   = zero_spd
      new%v%data = missing_r
      new%v%qc   = zero_spd
   ELSE IF   (                   new%speed%data .LT. 0              ) THEN
      new%u%data = missing_r
      new%u%qc   = negative_spd
      new%v%data = missing_r
      new%v%qc   = negative_spd
   ELSE IF ( (                   new%direction%data .LT. 0          ) .OR. &
             (                   new%direction%data .GT. 360        ) ) THEN
      new%u%data = missing_r
      new%u%qc   = wrong_direction
      new%v%data = missing_r
      new%v%qc   = wrong_direction
   ELSE
      lon_dif = longitude - lon_center
      if ( lon_dif .gt. 180. ) lon_dif = lon_dif - 360.
      if ( lon_dif .lt. -180. ) lon_dif = lon_dif + 360.
      dir_map = new%direction%data - ( lon_dif ) * cone_factor * SIGN ( 1. , lat_center )
      if ( dir_map .GT. 360. ) dir_map = dir_map - 360.
      if ( dir_map .LT.   0. ) dir_map = 360     + dir_map
      new%u%data = -1. * new%speed%data * SIN ( dir_map * piover180 ) 
      new%v%data = -1. * new%speed%data * COS ( dir_map * piover180 ) 
      new%u%qc   = new%direction%qc
      new%v%qc   = new%direction%qc
   END IF

END SUBROUTINE diagnostics

!
! ----------------------------------------------------------------------------

SUBROUTINE extend_to_closest ( surface , pressure , levels , &
max_p_extend_t , max_p_extend_w )

!  measurements.

   IMPLICIT NONE 

   INTEGER , INTENT ( IN )                   :: levels
   REAL , INTENT ( IN ) , DIMENSION (levels) :: pressure
   TYPE ( measurement ) , POINTER            :: surface
   REAL                                      :: max_p_extend_t ,  &
                                                max_p_extend_w

   INTEGER                                   :: k
   REAL                                      :: closest , &
                                                pres_ob , &
                                                diff    , & 
                                                dtdp    , &
                                                dp

   TYPE ( measurement ) , POINTER            :: new_level

   !  The only criteria is that we have a pressure for this data.

   IF ( .NOT. eps_equal ( surface%meas%pressure%data , missing_r , 1. ) ) THEN

      pres_ob = surface%meas%pressure%data

      !  Find the closest pressure from the analysis levels.  The second level
      !  skips past the first level, which is designed as the surface.

      closest = pressure(2) * 100.
      diff = ABS ( pressure(2) * 100. - pres_ob )
      pressure_loop_find : DO k = 2 , levels
         IF ( ABS ( pressure(k) * 100. - pres_ob ) .LT. diff ) THEN
            closest = pressure(k) * 100.
            diff = ABS ( pressure(k) * 100. - pres_ob )
         END IF
      END DO pressure_loop_find

      !  Make space for the new observation level.  Make the data the same 
      !  on the two levels.

      ALLOCATE ( new_level )
      new_level = surface
      NULLIFY ( new_level%next )

      !  Modify the pressure to the correct level.

      new_level%meas%pressure%data  = closest

      !  These are the variables that we want, but only if they are within the 
      !  specified delta-p.  The variables are already set from the previous
      !  new_level = surface statement, so this IF test sets things to missing
      !  if the data are not close enough in vertical space.  The majority of the
      !  upper-air data is flight level (200 - 300 hPa), so a standard atmosphere
      !  lapse rate at that level is appropriate.

      IF        ( diff .GT. max_p_extend_t ) THEN
         new_level%meas%temperature%data = missing_r
      ELSE IF ( ( diff .LE. max_p_extend_t ) .AND. &
                ( .NOT. eps_equal ( surface%meas%temperature%data , missing_r , 1. ) ) ) THEN
         !  These constants are computable from CRC, Appendix F, page 171.
         dtdp = 14.6928 * ( ( pres_ob + closest ) * 0.5 * 0.01 ) ** (-0.8097368)
         dp   = ( closest -pres_ob ) * 0.01
         new_level%meas%temperature%data = surface%meas%temperature%data + dtdp * dp
      END IF

      IF ( diff .GT. max_p_extend_w ) THEN
         new_level%meas%speed%data       = missing_r
         new_level%meas%direction%data   = missing_r
         new_level%meas%u%data           = missing_r
         new_level%meas%v%data           = missing_r
      END IF

      !  For the data we want to use, set the QC flags to the value for this
      !  influence extension.

      new_level%meas%pressure%qc    = new_level%meas%pressure%qc    + extend_influence
      new_level%meas%temperature%qc = new_level%meas%temperature%qc + extend_influence
      new_level%meas%speed%qc       = new_level%meas%speed%qc       + extend_influence
      new_level%meas%direction%qc   = new_level%meas%direction%qc   + extend_influence
      new_level%meas%u%qc           = new_level%meas%u%qc           + extend_influence
      new_level%meas%v%qc           = new_level%meas%v%qc           + extend_influence

      !  For the data we do not want to propagate, set both the value and QC flag
      !  to the correct flag values.

      new_level%meas%dew_point%data = missing_r
      new_level%meas%dew_point%qc   = new_level%meas%dew_point%qc   + extend_influence
      new_level%meas%height%data    = missing_r
      new_level%meas%height%qc      = new_level%meas%height%qc      + extend_influence
      new_level%meas%thickness%data = missing_r
      new_level%meas%thickness%qc   = new_level%meas%thickness%qc   + extend_influence

      !  Does the fake level go before or after the current observation in the
      !  linked list?

      IF      ( closest .GT. pres_ob ) THEN
         new_level%next => surface
         surface => new_level
      ELSE IF ( closest .LT. pres_ob ) THEN
         surface%next => new_level
      END IF

   END IF

END SUBROUTINE extend_to_closest

!
!------------------------------------------------------------------------------

SUBROUTINE height_to_pres ( height , pressure , iew , jns , kbu , &
map_projection , latitude , longitude , value ) 

!  This routine is used to compute the pressure of an observation
!  using the lat/lon location to compute the (x,y) in the height
!  array.  This computed height location is used to give an
!  approximate pressure value.  The aternative is to use a standard
!  atmosphere approximation from the height (height_to_pres_old is
!  that routine).
  
   USE map_utils
   USE map_utils_helper

   IMPLICIT NONE

   INTEGER                                 :: iew , jns , kbu , &
                                              map_projection
   REAL    , DIMENSION ( jns , iew , kbu ) :: height
   REAL    , DIMENSION ( kbu )             :: pressure
   REAL                                    :: latitude , & 
                                              longitude
   TYPE ( meas_data )                      :: value

   REAL                                    :: x_location , &
                                              y_location , & 
                                              h_ob       , & 
                                              h_above    , & 
                                              h_below    , & 
                                              p_above    , & 
                                              p_below

   INTEGER                                 :: i , j , k

   LOGICAL                                 :: found

   !  We can zoom directly out of here if the height with which we are 
   !  going to do the approximation is a flag (meaning it is missing).

   missing_height : IF ( .NOT. eps_equal ( value%height%data , missing_r , 1. ) ) THEN

      !  First, we get the (x,y) location of the observation in the 
      !  domain of interest.
      
      CALL latlon_to_ij ( projd , latitude  , longitude , x_location  , y_location )
   
      !  If this observation is not WELL within the domain confines, we need to
      !  exit out of here without modifying the input data set (value).
   
      inside_domain : IF ( ( x_location .GT. 1       ) .AND. &
                           ( x_location .LT. iew - 2 ) .AND. &
                           ( y_location .GT. 1       ) .AND. &
                           ( y_location .LT. jns - 2 ) ) THEN 
   
         i = NINT ( x_location )
         j = NINT ( y_location )
   
         !  Now, we figure out where the pressure level is located.  We know the height
         !  at each (i,j,k), and we know the pressure at each of those k-levels.  We use
         !  a "close enough" value for the (i,j), and vertically search that column for 
         !  the two surrounding heights.  We save the heights and pressure values above 
         !  and below the target height for later interpolation.
   
         h_ob = value%height%data
   
         found = .FALSE.
         find_trapping : DO k = 2 , kbu-1
   
            IF ( ( h_ob .GE. height(j,i,k) ) .AND. ( h_ob .LT. height(j,i,k+1) ) ) THEN
               found = .TRUE.
               h_below =  height(j,i,k  )
               h_above =  height(j,i,k+1)
               p_below =  pressure(k  ) * 100.
               p_above =  pressure(k+1) * 100.
               EXIT find_trapping            
            END IF
         END DO find_trapping
   
         !  If we found the trapping levels, great.  If not, we drop through this
         !  block and out of the routine.  This would happen if an observation was
         !  reported very high in the atmosphere, beyond the analysis levels.  We 
         !  compute the new pressure from linear interpolation.  The QC flag for the
         !  pressure reflects that this data was interpolated.
   
         found_level : IF ( found ) THEN
            value%pressure%data = EXP ( ( LOG ( p_below ) * ( h_ob - h_above) + LOG ( p_above ) * ( h_below - h_ob ) ) /  &
            ( h_below - h_above ) )
            value%pressure%qc   = p_from_h_first_guess
         END IF found_level
   
      END IF inside_domain

   END IF missing_height

END SUBROUTINE height_to_pres

!
!------------------------------------------------------------------------------

SUBROUTINE height_to_pres_old ( new )

!  This takes a measurement in which the pressure is missing, check
!  if there is height data. If so, map the height to a pressure.

   IMPLICIT NONE

   TYPE ( meas_data )    :: new

   REAL                  :: t

   !  Compute the standard atmosphere pressure from the height (CRC F-191).

   height_is_there : IF ( .NOT. eps_equal ( new%height%data , missing_r , 1. ) ) THEN

      IF ( new%height%data .LT. 11000 ) THEN

         t = 288.15 - 0.0065 * new%height%data
         new%pressure%data = 1013.25 * ( 288.15 / t )**(-5.255877) * 100
         new%pressure%qc   = p_std_atm_and_height

      ELSE IF ( new%height%data .LT. 20000 ) THEN

         new%pressure%data = 226.32 * EXP (-0.000156768832 * ( new%height%data - 11000 ) ) * 100
         new%pressure%qc   = p_std_atm_and_height
  
      ELSE IF ( new%height%data .LT. 32000 ) THEN

         t = 216.65 + 0.001 * ( new%height%data - 19999.997 ) 
         new%pressure%data = 54.7487 * ( 216.65 / t ) ** 34.16319 * 100
         new%pressure%qc   = p_std_atm_and_height

      END IF

   END IF height_is_there

END SUBROUTINE height_to_pres_old

!
!---------------------------------------------------------------------------

SUBROUTINE insert_at ( surface , new , elevation )

!  This takes a new measurement (new) and inserts it in a linked list
!  of measurements (surface points to first in list) in decreasing order of
!  pressure value.  If two levels' pressures are 'eps_equal', the levels
!  are merged instead of being linked.

   IMPLICIT NONE

   TYPE ( measurement ) ,  POINTER         :: surface , new
   REAL , INTENT(IN)                       :: elevation

   TYPE ( measurement ) , POINTER          :: current , previous , oldptr
   REAL                                    :: new_pres , new_height
   CHARACTER ( LEN = 32 ) , PARAMETER      :: name = 'insert_at'

   INCLUDE 'error.inc'
   INTERFACE
      INCLUDE 'error.int'
   END INTERFACE

   !  Initialize the variable to test the pressure and the place where the
   !  to-be-inserted measurement points.

   new_pres = new%meas%pressure%data
   new_height = new%meas%pressure%data
   NULLIFY ( new%next )

   !  The first check is to see if we are at the head of the linked list.  This
   !  drops us through to exit the routine.

   IF ( .NOT. ASSOCIATED ( surface ) ) THEN

      surface => new

   !  We are either between a couple of values, after a last value, or we could need
   !  to be merged with a level.  All those tests are handled in this else block.

   ELSE

      !  Initialize some dummy pointers to traverse to where we need to be.

      previous => surface 
      current => surface

      !  Loop to find correct location to link in 'new'.  The pressure is monotonically
      !  decreasing, so as soon as we find one where the current pressure is less than
      !  the new pressure, the new pressure goes just before it (or we run out of data 
      !  looking!).  Additionally, if both of the heights are equal AND the heights are
      !  the same as the input elevation of the station, then these need to be merged
      !  surface observations.

      still_some_data : DO WHILE ( ASSOCIATED ( current ) )
         IF ( current%meas%pressure%data .LT. new_pres ) EXIT still_some_data
         previous => current
         current => current%next
      END DO still_some_data 

      !  There are several cases:
      !  1) the new value has the same pressure as the previous value, or
      !     both heights are equal to the station elevation: merge them
      !  2) ran out of data finding where to insert level: add it to the end
      !  3) the new value has the same pressure as the current pressure value, or
      !     both heights are equal to the station elevation: merge them
      !  4) new pressure is < the previous value: stick it at end of previous
      !  5) new pressure > than previous: put at head of list
      !  ***** THE ORDER OF THE TESTS IS IMPORTANT *****

      IF ( ( eps_equal ( previous%meas%pressure%data , new_pres   , 1. ) ) .OR. &
         ( ( eps_equal ( previous%meas%height%data   , new_height , 1. ) ) .AND. &
           ( eps_equal ( previous%meas%height%data   , elevation  , 1. ) ) ) ) THEN

         CALL merge_measurements ( previous%meas , new%meas , 1 )
         DEALLOCATE ( new )

      ELSE IF ( .NOT. ASSOCIATED ( current ) ) THEN

         previous%next => new

      ELSE IF ( ( eps_equal ( current%meas%pressure%data , new_pres   , 1. ) ) .OR. &
              ( ( eps_equal ( current%meas%height%data   , new_height , 1. ) ) .AND. &
                ( eps_equal ( current%meas%height%data   , elevation  , 1. ) ) ) ) THEN

         CALL merge_measurements ( current%meas , new%meas , 1 )
         DEALLOCATE ( new )

      ELSE IF ( previous%meas%pressure%data .GT. new_pres ) THEN

         oldptr => previous%next
         previous%next => new
         new%next => oldptr

      ELSE IF ( previous%meas%pressure%data .LT. new_pres ) THEN

         ! If we aren't at head of list, have some internal (fatal) error.

         IF ( .NOT. ASSOCIATED ( previous , surface ) ) THEN
            CALL error_handler ( 33341 , name // 'Logic error in IF' , &
                                 .true. , .false. )
         ELSE
            oldptr => surface
            surface => new
            new%next => oldptr
         END IF 

      ELSE

         !  One of those "should never get here" logic errors, fatal.

         CALL error_handler ( 33342 , name &
         // 'Logic error in IF test for where to put the new observation level.' , &
         .true. , .false. )

      END IF

   END IF

END SUBROUTINE insert_at

!
! ----------------------------------------------------------------------------

SUBROUTINE inside_window ( lat , lon , iew , jns , outside_window )

!  This routine determines if an observation with the input latitude and
!  longitude is within the current domain.

   USE map_utils
   USE map_utils_helper

   IMPLICIT NONE

   REAL , INTENT(IN) :: lat , lon
   INTEGER , INTENT(IN) :: iew , jns
   LOGICAL , INTENT(OUT) :: outside_window

   !  Local data

   REAL :: x_location , y_location

   !  common

   REAL ::  XLATC, XLONC, XN, POLE, DS,   &
            TRUE_LAT1, TRUE_LAT2 , &
            dxd 

   INTEGER ::  KPROJ,    &
            IMAX, JMAX, &
            iewe , jnse

   COMMON /MAP_STUFF/ XLATC, XLONC, XN, POLE, DS, KPROJ,    &
                      IMAX, JMAX, TRUE_LAT1, TRUE_LAT2 , &
                      dxd , iewe , jnse


   IF ( ABS(lat) .GT. 90. ) THEN

      outside_window = .TRUE.

   ELSE IF ( ABS(lon) .GT. 360. ) THEN

      outside_window = .TRUE.

   ELSE IF ( ( projd%code .EQ. PROJ_MERC ) .AND. ( ABS(lat) .GT. 89. ) ) THEN

      outside_window = .TRUE.

   ELSE IF ( ( projd%code .EQ. PROJ_PS   ) .AND. ( xlatc .GT. 0 ) .AND. ( lat .LT. -89. ) ) THEN

      outside_window = .TRUE.

   ELSE IF ( ( projd%code .EQ. PROJ_PS   ) .AND. ( xlatc .LT. 0 ) .AND. ( lat .GT.  89. ) ) THEN

      outside_window = .TRUE.

   ELSE IF ( ( projd%code .EQ. PROJ_LC   ) .AND. ( xlatc .GT. 0 ) .AND. ( lat .LT. -89. ) ) THEN

      outside_window = .TRUE.

   ELSE IF ( ( projd%code .EQ. PROJ_LC   ) .AND. ( xlatc .LT. 0 ) .AND. ( lat .GT.  89. ) ) THEN

      outside_window = .TRUE.

   ELSE

      CALL latlon_to_ij ( projd , lat , lon , x_location  , y_location )
   
      IF ( ( x_location .GE. 1 ) .AND. ( x_location .LE. iew ) .AND. & 
           ( y_location .GE. 1 ) .AND. ( y_location .LE. jns ) ) THEN
         outside_window = .FALSE.
      ELSE
         outside_window = .TRUE.
      END IF

   END IF

END SUBROUTINE inside_window

!
! ----------------------------------------------------------------------------

SUBROUTINE interp_all ( pres_mb , levels , surface , lon )

!  This takes a linked list of measurements and an input array ofpressure levels,
!  finds the pressure levels closest to the new pressure, then interpolates
!  all measurement values for the new pressure level.  It links in a new
!  measurement, computed by interpolation, into the linked list of
!  measurements.

   IMPLICIT NONE 

   INCLUDE 'proc_get_info_header.inc'
   INCLUDE 'map.inc'
   INTEGER               :: iew_map , jns_map

   INTEGER , INTENT ( IN )                   :: levels
   REAL , INTENT ( IN ) , DIMENSION (levels) :: pres_mb
   REAL , INTENT(IN)                         :: lon
   TYPE ( measurement ) , POINTER            :: surface

   INTEGER                                   :: k , kstart

   INTEGER                                   :: linear_or_log_p_w , & 
                                                linear_or_log_p_t , & 
                                                linear_or_log_p_r , & 
                                                linear_or_log_p_z

   LOGICAL                                   :: start , found_level

   TYPE ( field )                            :: val1_u , val2_u , p1_w , p2_w , &
                                                val1_v , val2_v               , &
                                                val1_t , val2_t , p1_t , p2_t , &
                                                val1_r , val2_r , p1_r , p2_r , &
                                                val1_z , val2_z , p1_z , p2_z

   TYPE ( measurement ) , POINTER            :: previous_lower_w , previous_upper_w , &
                                                previous_lower_t , previous_upper_t , &
                                                previous_lower_r , previous_upper_r , &
                                                previous_lower_z , previous_upper_z

   TYPE ( measurement ) , POINTER            :: next_w , &
                                                next_t , &
                                                next_r , &
                                                next_z

   LOGICAL                                   :: found_lower_w , found_upper_w , & 
                                                found_lower_t , found_upper_t , & 
                                                found_lower_r , found_upper_r , & 
                                                found_lower_z , found_upper_z

   TYPE ( measurement ) , POINTER            :: previous , & 
                                                next     , & 
                                                new_level
   REAL , PARAMETER :: L_over_Rv = 5418.12
   REAL , PARAMETER :: Rv_over_L = 1. / L_over_Rv

   REAL :: delta_lon

   !  Assign the linear in p or ln-p interpolation based upon the
   !  variables: 1=linear in pressure interpolation; 0=linear in 
   !  ln(p) interpolation.

   linear_or_log_p_w = 1 
   linear_or_log_p_t = 0 
   linear_or_log_p_r = 1 
   linear_or_log_p_z = 0 

   !  Initialize the variable that says we have pressure levels to interpolate
   !  to that are possible, given the first available pressure from the 
   !  linked list of data.  For example, no interpolation to 850 mb if the
   !  linked list starts at 300 mb.  If all of the pressure levels are larger
   !  than the first pressure surface in the linked list, start remains false.

   start = .FALSE.

   !  If we do not have at least two levels of data, then there is no
   !  need to waste any time.  The surface POINTER is required to be
   !  ASSOCIATED on input.

   IF ( ASSOCIATED ( surface%next ) ) THEN

      !  Skip over level (1) if it is the "surface" level.  This is not an
      !  isobaric surface, but a level with a flag value for pressure.
   
      IF ( eps_equal ( pres_mb(1) , 1001. , 0.01 ) ) THEN
         kstart = 2
      ELSE
         kstart = 1
      END IF

      !  Loop over each of the requested pressure levels to find where we can begin.
   
      p_level_search : DO k = kstart , levels
   
         IF ( pres_mb(k) * 100 .GT. surface%meas%pressure%data ) THEN
            CYCLE p_level_search
         ELSE 
            kstart = k
            start = .TRUE.
            EXIT p_level_search
         END IF

      END DO p_level_search

   END IF

   !  A quicky test to see if this would be a valid vertical interpolation
   !  is that there are at least 2 levels, and one of them is physically 
   !  below at least one of the requested pressure surfaces.

   IF ( ASSOCIATED ( surface%next ) .AND. start ) THEN

      !  Initialize the pointers to the first two locations in the linked 
      !  list.  This allows us to keep our place when looping through for the 
      !  different levels.  Identifiers: _w is wind, _t is temperature, _r is 
      !  relative humidity, and _z is geopotential height.  Since u and v are 
      !  tightly related, their location is the same ( _w ).
      
      previous_lower_w => surface
      previous_upper_w => surface%next
   
      previous_lower_t => surface
      previous_upper_t => surface%next
   
      previous_lower_r => surface
      previous_upper_r => surface%next
   
      previous_lower_z => surface
      previous_upper_z => surface%next

      !  Loop over all of the allowable pressure levels.

      pressure_loop : DO k = kstart , levels

         !  All data for the interpolation is initialized to missing.  The 
         !  naming convention is as above, though _u and _v refer to the 
         !  specific horizontal wind components.

         val1_u%data = missing_r
         val1_u%qc   = missing  
         p1_w%data   = missing_r
         p1_w%qc     = missing  

         val1_v%data = missing_r
         val1_v%qc   = missing  

         val1_t%data = missing_r
         val1_t%qc   = missing  
         p1_t%data   = missing_r
         p1_t%qc     = missing  

         val1_r%data = missing_r
         val1_r%qc   = missing  
         p1_r%data   = missing_r
         p1_r%qc     = missing  

         val1_z%data = missing_r
         val1_z%qc   = missing  
         p1_z%data   = missing_r
         p1_z%qc     = missing  

         next_w => previous_lower_w
         next_t => previous_lower_t
         next_r => previous_lower_r
         next_z => previous_lower_z

         !  Find the lower level (larger pressure) for the horizontal winds ( _w ).

         loop_w1 : DO 
            IF      ( ( eps_equal ( pres_mb(k) * 100 , next_w%meas%pressure%data , 1. ) ) .AND. &
                      ( .NOT. eps_equal ( next_w%meas%u%data , missing_r , 1. ) ) .AND. &
                      ( .NOT. eps_equal ( next_w%meas%v%data , missing_r , 1. ) ) ) THEN
               val1_u%data = missing_r
               val1_u%qc   = missing  
               p1_w%data   = missing_r
               p1_w%qc     = missing  
               val1_v%data = missing_r
               val1_v%qc   = missing  
               EXIT loop_w1
            ELSE IF ( ( eps_equal ( pres_mb(k) * 100 , next_w%meas%pressure%data , 1. ) ) .AND. &
                    ( ( eps_equal ( next_w%meas%u%data , missing_r , 1. ) ) .OR. &
                      ( eps_equal ( next_w%meas%v%data , missing_r , 1. ) ) ) ) THEN
               EXIT loop_w1
            ELSE IF ( ( pres_mb(k) * 100 .LT. next_w%meas%pressure%data ) .AND. &
                      ( next_w%meas%u%qc .NE. vert_interpolated ) .AND. & 
                      ( next_w%meas%v%qc .NE. vert_interpolated ) .AND. & 
                      ( .NOT. eps_equal ( next_w%meas%u%data , missing_r , 1. ) ) .AND. &
                      ( .NOT. eps_equal ( next_w%meas%v%data , missing_r , 1. ) ) ) THEN
               val1_u = next_w%meas%u
               val1_v = next_w%meas%v
               p1_w   = next_w%meas%pressure
               previous_lower_w => next_w
               next_w => next_w%next
            ELSE IF ( ( pres_mb(k) * 100 .LT. next_w%meas%pressure%data ) .AND. &
                    ( ( next_w%meas%u%qc .EQ. vert_interpolated ) .OR. & 
                      ( next_w%meas%v%qc .EQ. vert_interpolated ) .OR. & 
                      ( eps_equal ( next_w%meas%u%data , missing_r , 1. ) ) .OR. &
                      ( eps_equal ( next_w%meas%v%data , missing_r , 1. ) ) ) ) THEN
               next_w => next_w%next
            ELSE IF ( pres_mb(k) * 100 .GT. next_w%meas%pressure%data ) THEN
               EXIT loop_w1
            END IF
            IF ( .NOT. ASSOCIATED ( next_w ) )  EXIT loop_w1
         END DO loop_w1

         IF ( ( eps_equal ( val1_u%data , missing_r , 1. ) ) .OR. &
              ( eps_equal ( p1_w%data   , missing_r , 1. ) ) ) THEN
            found_lower_w = .FALSE.
         ELSE
            found_lower_w = .TRUE.
         END IF

         !  Find the lower level (larger pressure) for the temperature ( _t ).

         loop_t1 : DO 
            IF      ( ( eps_equal ( pres_mb(k) * 100 , next_t%meas%pressure%data , 1. ) ) .AND. &
                      ( .NOT. eps_equal ( next_t%meas%temperature%data , missing_r , 1. ) ) ) THEN
               val1_t%data = missing_r
               val1_t%qc   = missing  
               p1_t%data   = missing_r
               p1_t%qc     = missing  
               EXIT loop_t1
            ELSE IF ( ( eps_equal ( pres_mb(k) * 100 , next_t%meas%pressure%data , 1. ) ) .AND. &
                      ( eps_equal ( next_t%meas%temperature%data , missing_r , 1. ) ) ) THEN
               EXIT loop_t1
            ELSE IF ( ( pres_mb(k) * 100 .LT. next_t%meas%pressure%data ) .AND. &
                      ( next_t%meas%temperature%qc .NE. vert_interpolated ) .AND. & 
                      ( .NOT. eps_equal ( next_t%meas%temperature%data , missing_r , 1. ) ) ) THEN
               val1_t = next_t%meas%temperature
               p1_t   = next_t%meas%pressure
               previous_lower_t => next_t
               next_t => next_t%next
            ELSE IF ( ( pres_mb(k) * 100 .LT. next_t%meas%pressure%data ) .AND. &
                    ( ( next_t%meas%temperature%qc .EQ. vert_interpolated ) .OR. & 
                      ( eps_equal ( next_t%meas%temperature%data , missing_r , 1. ) ) ) ) THEN
               next_t => next_t%next
            ELSE IF ( pres_mb(k) * 100 .GT. next_t%meas%pressure%data ) THEN
               EXIT loop_t1
            END IF
            IF ( .NOT. ASSOCIATED ( next_t ) )  EXIT loop_t1
         END DO loop_t1

         IF ( ( eps_equal ( val1_t%data , missing_r , 1. ) ) .OR. &
              ( eps_equal ( p1_t%data   , missing_r , 1. ) ) ) THEN
            found_lower_t = .FALSE.
         ELSE
            found_lower_t = .TRUE.
         END IF

         !  Find the lower level (larger pressure) for the relative humidity ( _r ).

         loop_r1 : DO 
            IF      ( ( eps_equal ( pres_mb(k) * 100 , next_r%meas%pressure%data , 1. ) ) .AND. &
                      ( .NOT. eps_equal ( next_r%meas%rh%data , missing_r , 1. ) ) .AND. &
                      ( .NOT. eps_equal ( next_r%meas%temperature%data , missing_r , 1. ) ) ) THEN
               val1_r%data = missing_r
               val1_r%qc   = missing  
               p1_r%data   = missing_r
               p1_r%qc     = missing  
               EXIT loop_r1
            ELSE IF ( ( eps_equal ( pres_mb(k) * 100 , next_r%meas%pressure%data , 1. ) ) .AND. &
                    ( ( eps_equal ( next_r%meas%rh%data , missing_r , 1. ) ) .OR. &
                      ( eps_equal ( next_r%meas%temperature%data , missing_r , 1. ) ) ) ) THEN
               EXIT loop_r1
            ELSE IF ( ( pres_mb(k) * 100 .LT. next_r%meas%pressure%data ) .AND. &
                      ( next_r%meas%rh%qc .NE. vert_interpolated ) .AND. & 
                      ( next_r%meas%temperature%qc .NE. vert_interpolated ) .AND. & 
                      ( .NOT. eps_equal ( next_r%meas%rh%data , missing_r , 1. ) ) .AND. &
                      ( .NOT. eps_equal ( next_r%meas%temperature%data , missing_r , 1. ) ) ) THEN
               val1_r = next_r%meas%rh
               p1_r   = next_r%meas%pressure
               previous_lower_r => next_r
               next_r => next_r%next
            ELSE IF ( ( pres_mb(k) * 100 .LT. next_r%meas%pressure%data ) .AND. &
                    ( ( next_r%meas%rh%qc .EQ. vert_interpolated ) .OR. & 
                      ( next_r%meas%temperature%qc .EQ. vert_interpolated ) .OR. & 
                      ( eps_equal ( next_r%meas%rh%data , missing_r , 1. ) ) .OR. &
                      ( eps_equal ( next_r%meas%temperature%data , missing_r , 1. ) ) ) ) THEN
               next_r => next_r%next
            ELSE IF ( pres_mb(k) * 100 .GT. next_r%meas%pressure%data ) THEN
               EXIT loop_r1
            END IF
            IF ( .NOT. ASSOCIATED ( next_r ) )  EXIT loop_r1
         END DO loop_r1

         IF ( ( eps_equal ( val1_r%data , missing_r , 1. ) ) .OR. &
              ( eps_equal ( p1_r%data   , missing_r , 1. ) ) ) THEN
            found_lower_r = .FALSE.
         ELSE
            found_lower_r = .TRUE.
         END IF

         !  Find the lower level (larger pressure) for the geopotential height ( _z ).

         loop_z1 : DO 
            IF      ( ( eps_equal ( pres_mb(k) * 100 , next_z%meas%pressure%data , 1. ) ) .AND. &
                      ( .NOT. eps_equal ( next_z%meas%height%data , missing_r , 1. ) ) ) THEN
               val1_z%data = missing_r
               val1_z%qc   = missing  
               p1_z%data   = missing_r
               p1_z%qc     = missing  
               EXIT loop_z1
            ELSE IF ( ( eps_equal ( pres_mb(k) * 100 , next_z%meas%pressure%data , 1. ) ) .AND. &
                      ( eps_equal ( next_z%meas%height%data , missing_r , 1. ) ) ) THEN
               EXIT loop_z1
            ELSE IF ( ( pres_mb(k) * 100 .LT. next_z%meas%pressure%data ) .AND. &
                      ( next_z%meas%height%qc .NE. vert_interpolated ) .AND. & 
                      ( .NOT. eps_equal ( next_z%meas%height%data , missing_r , 1. ) ) ) THEN
               val1_z = next_z%meas%height
               p1_z   = next_z%meas%pressure
               previous_lower_z => next_z
               next_z => next_z%next
            ELSE IF ( ( pres_mb(k) * 100 .LT. next_z%meas%pressure%data ) .AND. &
                    ( ( next_z%meas%height%qc .EQ. vert_interpolated ) .OR. & 
                      ( eps_equal ( next_z%meas%height%data , missing_r , 1. ) ) ) ) THEN
               next_z => next_z%next
            ELSE IF ( pres_mb(k) * 100 .GT. next_z%meas%pressure%data ) THEN
               EXIT loop_z1
            END IF
            IF ( .NOT. ASSOCIATED ( next_z ) )  EXIT loop_z1
         END DO loop_z1

         IF ( ( eps_equal ( val1_z%data , missing_r , 1. ) ) .OR. &
              ( eps_equal ( p1_z%data   , missing_r , 1. ) ) ) THEN
            found_lower_z = .FALSE.
         ELSE
            found_lower_z = .TRUE.
         END IF

         !  Initialize the upper-level interpolation variables to 
         !  missing.

         val2_u%data = missing_r
         val2_u%qc   = missing  
         p2_w%data   = missing_r
         p2_w%qc     = missing  

         val2_v%data = missing_r
         val2_v%qc   = missing  

         val2_t%data = missing_r
         val2_t%qc   = missing  
         p2_t%data   = missing_r
         p2_t%qc     = missing  

         val2_r%data = missing_r
         val2_r%qc   = missing  
         p2_r%data   = missing_r
         p2_r%qc     = missing  

         val2_z%data = missing_r
         val2_z%qc   = missing  
         p2_z%data   = missing_r
         p2_z%qc     = missing  

         next_w => previous_upper_w
         next_t => previous_upper_t
         next_r => previous_upper_r
         next_z => previous_upper_z

         !  Find the upper level (smaller pressure) for the horizontal winds ( _w ).

         IF ( found_lower_w ) THEN
            loop_w2 : DO 
               IF ( ( pres_mb(k) * 100 .GT. next_w%meas%pressure%data ) .AND. &
                    ( .NOT. eps_equal ( next_w%meas%u%data , missing_r , 1. ) ) .AND. &
                    ( .NOT. eps_equal ( next_w%meas%v%data , missing_r , 1. ) ) .AND. &
                    ( next_w%meas%u%qc .NE. vert_interpolated ) .AND. &
                    ( next_w%meas%v%qc .NE. vert_interpolated ) ) THEN
                  val2_u = next_w%meas%u
                  val2_v = next_w%meas%v
                  p2_w   = next_w%meas%pressure
                  previous_upper_w => next_w
                  EXIT loop_w2
               ELSE
                  next_w => next_w%next
               END IF
               IF ( .NOT. ASSOCIATED ( next_w ) )  EXIT loop_w2
            END DO loop_w2
         END IF

         IF ( ( eps_equal ( val2_u%data , missing_r , 1. ) ) .OR. &
              ( eps_equal ( p2_w%data   , missing_r , 1. ) ) ) THEN
            found_upper_w = .FALSE.
         ELSE
            found_upper_w = .TRUE.
         END IF

         !  Find the upper level (smaller pressure) for the temperature ( _t ).

         IF ( found_lower_t ) THEN
            loop_t2 : DO 
               IF ( ( pres_mb(k) * 100 .GT. next_t%meas%pressure%data ) .AND. &
                    ( .NOT. eps_equal ( next_t%meas%temperature%data , missing_r , 1. ) ) .AND. &
                    ( next_t%meas%temperature%qc .NE. vert_interpolated ) ) THEN
                  val2_t = next_t%meas%temperature
                  p2_t   = next_t%meas%pressure
                  previous_upper_t => next_t
                  EXIT loop_t2
               ELSE
                  next_t => next_t%next
               END IF
               IF ( .NOT. ASSOCIATED ( next_t ) )  EXIT loop_t2
            END DO loop_t2
         END IF

         IF ( ( eps_equal ( val2_t%data , missing_r , 1. ) ) .OR. &
              ( eps_equal ( p2_t%data   , missing_r , 1. ) ) ) THEN
            found_upper_t = .FALSE.
         ELSE
            found_upper_t = .TRUE.
         END IF

         !  Find the upper level (smaller pressure) for the relative humidity ( _r ).

         IF ( found_lower_r ) THEN
            loop_r2 : DO 
               IF ( ( pres_mb(k) * 100 .GT. next_r%meas%pressure%data ) .AND. &
                    ( .NOT. eps_equal ( next_r%meas%rh%data , missing_r , 1. ) ) .AND. &
                    ( .NOT. eps_equal ( next_r%meas%temperature%data , missing_r , 1. ) ) .AND. &
                    ( next_r%meas%rh%qc .NE. vert_interpolated ) .AND. &
                    ( next_r%meas%temperature%qc .NE. vert_interpolated ) ) THEN
                  val2_r = next_r%meas%rh
                  p2_r   = next_r%meas%pressure
                  previous_upper_r => next_r
                  EXIT loop_r2
               ELSE
                  next_r => next_r%next
               END IF
               IF ( .NOT. ASSOCIATED ( next_r ) )  EXIT loop_r2
            END DO loop_r2
         END IF

         IF ( ( eps_equal ( val2_r%data , missing_r , 1. ) ) .OR. &
              ( eps_equal ( p2_r%data   , missing_r , 1. ) ) ) THEN
            found_upper_r = .FALSE.
         ELSE
            found_upper_r = .TRUE.
         END IF

         !  Find the upper level (smaller pressure) for the geopotential height ( _z ).

         IF ( found_lower_z ) THEN
            loop_z2 : DO 
               IF ( ( pres_mb(k) * 100 .GT. next_z%meas%pressure%data ) .AND. &
                    ( .NOT. eps_equal ( next_z%meas%height%data , missing_r , 1. ) ) .AND. &
                    ( next_z%meas%height%qc .NE. vert_interpolated ) ) THEN
                  val2_z = next_z%meas%height
                  p2_z   = next_z%meas%pressure
                  previous_upper_z => next_z
                  EXIT loop_z2
               ELSE
                  next_z => next_z%next
               END IF
               IF ( .NOT. ASSOCIATED ( next_z ) )  EXIT loop_z2
            END DO loop_z2
         END IF

         IF ( ( eps_equal ( val2_z%data , missing_r , 1. ) ) .OR. &
              ( eps_equal ( p2_z%data   , missing_r , 1. ) ) ) THEN
            found_upper_z = .FALSE.
         ELSE
            found_upper_z = .TRUE.
         END IF

         !  If we did not find the upper level (by implication the lower level
         !  as well) of any of the variables, we do not want to add a level.  There
         !  wouldn't be any data to add, just missing values.

         IF ( ( .NOT. found_upper_w ) .AND. & 
              ( .NOT. found_upper_t ) .AND. & 
              ( .NOT. found_upper_r ) .AND. & 
              ( .NOT. found_upper_z ) ) THEN
            CYCLE pressure_loop
         END IF

         !  We would like to add a new level, but it is possible that it already
         !  exists, and just needs the new variable to replace the undefined
         !  fields.  Start checking to see if we already have a level with the
         !  correct pressure.  Because of the vertical interpolation, we know 
         !  that we are going to insert a new level, never put it at the front or 
         !  back of the list.  We either find the correct level, we have yet to 
         !  find it, or we have passed where it would be.

         next => surface
         found_level = .FALSE.

         make_this_level : DO
            IF ( eps_equal ( pres_mb(k) * 100 , next%meas%pressure%data , 1. ) ) THEN
               found_level = .TRUE.
               new_level => next
               EXIT make_this_level
            ELSE IF ( pres_mb(k) * 100 .LT. next%meas%pressure%data ) THEN
               previous => next
               next => next%next
            ELSE IF ( pres_mb(k) * 100 .GT. next%meas%pressure%data ) THEN
               ALLOCATE ( new_level )
               new_level%meas%pressure%data = pres_mb(k) * 100
               new_level%meas%pressure%qc   = vert_interpolated
               new_level%next => next
               previous%next => new_level

               IF ( found_upper_w ) THEN
                  new_level%meas%u%qc             = 0
                  new_level%meas%v%qc             = 0
                  new_level%meas%direction%qc     = 0
                  new_level%meas%speed%qc         = 0
               ELSE
                  new_level%meas%u%data           = missing_r
                  new_level%meas%v%data           = missing_r
                  new_level%meas%direction%data   = missing_r
                  new_level%meas%speed%data       = missing_r
                  new_level%meas%u%qc             = missing
                  new_level%meas%v%qc             = missing
                  new_level%meas%direction%qc     = missing
                  new_level%meas%speed%qc         = missing
               END IF

               IF ( found_upper_t ) THEN
                  new_level%meas%temperature%qc   = 0
               ELSE
                  new_level%meas%temperature%data = missing_r
                  new_level%meas%temperature%qc   = missing
               END IF

               IF ( found_upper_r ) THEN
                  new_level%meas%dew_point%qc     = 0
                  new_level%meas%rh%qc            = 0
               ELSE
                  new_level%meas%dew_point%data   = missing_r
                  new_level%meas%rh%data          = missing_r
                  new_level%meas%dew_point%qc     = missing
                  new_level%meas%rh%qc            = missing
               END IF

               IF ( found_upper_z ) THEN
                  new_level%meas%height%qc        = 0
               ELSE
                  new_level%meas%height%data      = missing_r
                  new_level%meas%height%qc        = missing
               END IF

               new_level%meas%thickness%data      = missing_r
               new_level%meas%thickness%qc        = missing

               EXIT make_this_level
            END IF
         END DO make_this_level

         !  Vertically interpolate the horizontal winds.  We need speed and
         !  direction as well.

         IF ( ( found_upper_w ) .AND. ( ABS ( p1_w%data - p2_w%data ) .LE. 15000 ) ) THEN
            CALL interp_level ( val1_u , val2_u , new_level%meas%u , &
                                 p1_w   , p2_w   , new_level%meas%pressure , & 
                                 linear_or_log_p_w ) 
   
            CALL interp_level ( val1_v , val2_v , new_level%meas%v , &
                                p1_w   , p2_w   , new_level%meas%pressure , & 
                                linear_or_log_p_w ) 
            new_level%meas%speed%data     = SQRT ( new_level%meas%u%data**2 + &
                                                   new_level%meas%v%data**2 ) 
            new_level%meas%speed%qc       = new_level%meas%u%qc
            IF      ( ( eps_equal ( new_level%meas%v%data , 0. , .01 ) ) .AND. &
                 ( new_level%meas%u%data .GT. 0 ) ) THEN
               new_level%meas%direction%data = 270.
            ELSE IF ( ( eps_equal ( new_level%meas%v%data , 0. , .01 ) ) .AND. &
                 ( new_level%meas%u%data .LT. 0 ) ) THEN
               new_level%meas%direction%data = 90.
            ELSE IF ( ( eps_equal ( new_level%meas%u%data , 0. , .01 ) ) .AND. &
                 ( new_level%meas%v%data .GT. 0 ) ) THEN
               new_level%meas%direction%data = 180.
            ELSE IF ( ( eps_equal ( new_level%meas%u%data , 0. , .01 ) ) .AND. &
                 ( new_level%meas%v%data .LT. 0 ) ) THEN
               new_level%meas%direction%data = 360
            ELSE IF ( eps_equal ( new_level%meas%speed%data , 0. , .01 ) ) THEN
               new_level%meas%direction%data = 0.
            ELSE
               new_level%meas%direction%data = 270. - ATAN2((new_level%meas%v%data),(new_level%meas%u%data)) * (180./3.14159265358)
               IF      ( new_level%meas%direction%data  .GT. 360 ) THEN
                  new_level%meas%direction%data = new_level%meas%direction%data - 360.
               ELSE IF ( new_level%meas%direction%data  .LT.   0 ) THEN
                  new_level%meas%direction%data = new_level%meas%direction%data + 360.
               END IF
            END IF
            IF      ( lon-lon_center .GT.  180 ) THEN
               delta_lon = lon-lon_center - 360.
            ELSE IF ( lon-lon_center .LT. -180 ) THEN
               delta_lon = lon-lon_center + 360.
            ELSE
               delta_lon = lon-lon_center
            END IF
            new_level%meas%direction%data = new_level%meas%direction%data  + delta_lon * cone_factor * SIGN ( 1. , lat_center )
            IF      ( new_level%meas%direction%data .LT.   0 ) THEN
               new_level%meas%direction%data = new_level%meas%direction%data + 360.
            ELSE IF ( new_level%meas%direction%data .GT. 360 ) THEN
               new_level%meas%direction%data = new_level%meas%direction%data - 360.
            END IF

            new_level%meas%direction%qc   = new_level%meas%u%qc
         ELSE IF ( .NOT. found_level ) THEN
            new_level%meas%u%data         = missing_r
            new_level%meas%u%qc           = vert_interpolated
            new_level%meas%v%data         = missing_r
            new_level%meas%v%qc           = vert_interpolated
            new_level%meas%speed%data     = missing_r
            new_level%meas%speed%qc       = vert_interpolated
            new_level%meas%direction%data = missing_r
            new_level%meas%direction%qc   = vert_interpolated
         END IF

         !  Vertically interpolate the temperature.

         IF ( ( found_upper_t ) .AND. ( ABS ( p1_t%data - p2_t%data ) .LE. 15000 ) ) THEN
            CALL interp_level ( val1_t , val2_t , new_level%meas%temperature , &
                                 p1_t   , p2_t   , new_level%meas%pressure , & 
                                 linear_or_log_p_t ) 
         ELSE IF ( .NOT. found_level ) THEN
            new_level%meas%temperature%data         = missing_r
            new_level%meas%temperature%qc           = vert_interpolated
         END IF

         !  Vertically interpolate the relative humidity.  Diagnose the
         !  dew point temperature, store that as well.

         IF ( ( found_upper_r ) .AND. ( ABS ( p1_r%data - p2_r%data ) .LE. 15000 ) .AND. &
              (  .NOT. eps_equal ( new_level%meas%temperature%data , missing_r , 1. ) ) ) THEN
            CALL interp_level ( val1_r , val2_r , new_level%meas%rh , &
                                 p1_r   , p2_r   , new_level%meas%pressure , & 
                                 linear_or_log_p_r ) 
            new_level%meas%dew_point%data     = 1. / &
            ( 1./new_level%meas%temperature%data - &
              Rv_over_L * LOG ( new_level%meas%rh%data / 100. ) )  
            new_level%meas%dew_point%qc       = new_level%meas%rh%qc
         ELSE IF ( .NOT. found_level ) THEN
            new_level%meas%rh%data                = missing_r
            new_level%meas%rh%qc                  = vert_interpolated
            new_level%meas%dew_point%data         = missing_r
            new_level%meas%dew_point%qc           = vert_interpolated
         END IF

         !  Vertically interpolate the geopotential height.

         IF ( ( found_upper_z ) .AND. ( ABS ( p1_z%data - p2_z%data ) .LE. 15000 ) ) THEN
            CALL interp_level ( val1_z , val2_z , new_level%meas%height , &
                                 p1_z   , p2_z   , new_level%meas%pressure , & 
                                 linear_or_log_p_z ) 
         ELSE IF ( .NOT. found_level ) THEN
            new_level%meas%height%data                = missing_r
            new_level%meas%height%qc                  = vert_interpolated
         END IF

      END DO pressure_loop
                  
   END IF

END SUBROUTINE interp_all
       
!
!---------------------------------------------------------------------------

SUBROUTINE interp_level ( first_met , second_met , new_met , &
                          first_p   , second_p   , new_p   , &
                          linear_or_log_p ) 

   IMPLICIT NONE

   TYPE ( field ) , INTENT ( IN )    :: first_p   , second_p  , &
                                        first_met , second_met
   TYPE ( field ) , INTENT ( OUT )   :: new_met
   TYPE ( field ) , INTENT ( IN )    :: new_p 
   INTEGER                           :: linear_or_log_p

   REAL                              :: p1 , p2 , pn
   
   INTEGER                           :: qc_small , & 
                                        qc_large

   INCLUDE 'constants.inc'

   !  The QC flag of the computed variable is the QC value of the same
   !  variable at the closest level + the value for vertical interpolation.

   IF ( ABS ( first_p%data - new_p%data ) .LE. ABS ( second_p%data - new_p%data ) ) THEN
      qc_small = mod ( first_met%qc , vert_interpolated ) 
      qc_large = ( first_met%qc / ( vert_interpolated * 2 ) ) * 2
      new_met%qc = qc_large*vert_interpolated + qc_small + vert_interpolated
   ELSE
      qc_small = mod ( second_met%qc , vert_interpolated ) 
      qc_large = ( second_met%qc / ( vert_interpolated * 2 ) ) * 2
      new_met%qc = qc_large*vert_interpolated + qc_small + vert_interpolated
   END IF

   !  Simple linear interpolation between the two known values.  If
   !  linear_or_log_p = 1, then linear in pressure.  If linear_or_log_p = 0,
   !  linear in ln(pressure).

   IF ( linear_or_log_p .EQ. 1 ) THEN
      p1 = first_p%data
      p2 = second_p%data
      pn = new_p%data
   ELSE IF ( linear_or_log_p .EQ. 0 ) THEN
      p1 = LOG ( first_p%data  ) 
      p2 = LOG ( second_p%data ) 
      pn = LOG ( new_p%data    ) 
   ELSE IF ( linear_or_log_p .EQ. 2 ) THEN
      p1 = ( first_p%data  ) ** rcp
      p2 = ( second_p%data ) ** rcp
      pn = ( new_p%data    ) ** rcp
   END IF

   new_met%data = ( first_met%data * ( pn - p2 ) + second_met%data * ( p1 - pn ) ) / ( p1 - p2 ) 

END SUBROUTINE interp_level
       
!
!---------------------------------------------------------------------------

SUBROUTINE keep_best ( field1 , field2 , best )

!  Use quality control (qc) info and keep the best value; if one is missing
!  keep the one that is known; if both present keep the one with better qc
!  flag;  if qc flags the same, keep the one chosen 'best'.

   IMPLICIT none

   TYPE ( field ) , INTENT ( INOUT )      :: field1
   TYPE ( field ) , INTENT ( IN )         :: field2
   INTEGER        , INTENT ( IN )         :: best

   CHARACTER ( LEN = 32 ) , PARAMETER     :: sub_name = 'keep_best'
   CHARACTER ( LEN = 80 )                 :: msg

   INCLUDE 'error.inc'
   INTERFACE
      INCLUDE 'error.int'
   END INTERFACE

! foo
!  IF ( field1 .EQ. field2 ) THEN 
   IF ( field_eq ( field1 , field2 ) ) THEN 

      ! If there is no difference, have nothing to do.

   ELSE IF ( (       eps_equal ( field1%data , missing_r , 1. ) ) .AND. &
             ( .NOT. eps_equal ( field2%data , missing_r , 1. ) ) ) THEN

      !  Copy both data and quality control flag.

      field1 = field2

   ELSE IF ( ( .NOT. eps_equal ( field1%data , missing_r , 1. ) ) .AND. &
             (       eps_equal ( field2%data , missing_r , 1. ) ) ) THEN

      !  Already have data in report1, so do nothing.

   ELSE IF ( (       eps_equal ( field1%data , missing_r , 1. ) ) .AND. &
             (       eps_equal ( field2%data , missing_r , 1. ) ) ) THEN

      !  When both data fields are empty, do nothing.

   ELSE IF ( ( .NOT. eps_equal ( field1%data , missing_r , 1. ) ) .AND. &
             ( .NOT. eps_equal ( field2%data , missing_r , 1. ) ) ) THEN

      !  There are several cases to consider if both fields have data:

      !  First,  use quality control flags to differentiate the differences.

      IF ( field1%qc .LT. field2%qc ) THEN

         !  Since the data is already in field1, do nothing.

      ELSE IF ( field1%qc .GT. field2%qc ) THEN

         field1 = field2

      ELSE IF ( field1%qc .EQ. field2%qc ) THEN

         !  Second, if they have the same quality control values, use data 
         !  that was chosen 'best'

         IF ( best .EQ. 1 ) THEN

            !  Again, since the data is already in field1, do nothing.

         ELSE IF ( best .EQ. 2 ) THEN

            field1 = field2

         ELSE

            ! Should never execute this part; have invalid 'best' integer.
            msg = 'Internal logic error.  Invalid value of ''best'''
            CALL error_handler ( 331103 , sub_name // msg , .TRUE. , .FALSE. )

         END IF

      ELSE  
 
         ! Should never execute this; if so, have fatal error
         msg = 'Internal logic error.  Either the QCs are different or the same.'
         CALL error_handler ( 331101 , sub_name // msg , .TRUE. , .FALSE. )

      END IF

   ELSE  

      ! should never execute this; if so have fatal error
      msg = 'Internal logic error.  Only four combinations of fields missing are possible.'
      CALL error_handler ( 331102 , sub_name // msg , .TRUE. , .FALSE. )

   END IF

END SUBROUTINE keep_best

!
!---------------------------------------------------------------------------

SUBROUTINE link_levels ( list1 , list2 , info1 , info2 , best )

!  Starting at the surface level, link levels into one list if pressure levels
!  are different;  if have two levels at the same pressure level (within 
!  epsilon) then keep the best data.  The resulting (output) linked list
!  starts from list1; list2 contains nothing useful on return.

   IMPLICIT NONE

   TYPE ( measurement ) , POINTER           :: list1 , list2
   INTEGER , INTENT ( IN )                  :: best

   TYPE ( measurement ) , POINTER           :: next1 , &
                                               next2 , &
                                               current , &
                                               delete_it

   TYPE ( source_info )                     :: info1 , info2

   !  Initialize both traversal pointers.

   next1 => list1
   next2 => list2
   NULLIFY ( current )

   !  Merge until the end of either list1 or list2 is reached.

   still_associated : DO WHILE ( ASSOCIATED ( next1 ) .AND. ASSOCIATED ( next2 ) )

      IF (    ( eps_equal ( next1%meas%pressure%data , & 
                            next2%meas%pressure%data , 1. ) ) &
                             .OR.  &
           (  ( eps_equal ( info1%elevation        , & 
                            next1%meas%height%data , .1 ) ) .AND. &
              ( eps_equal ( info2%elevation        , & 
                            next2%meas%height%data , .1 ) ) .AND. &
              ( .NOT. eps_equal ( info1%elevation , missing_r , 1. ) ) .AND. &
              ( eps_equal ( next1%meas%height%data , & 
                            next2%meas%height%data , .1 ) ) ) ) THEN

         !  There are two ways that cause us to merge the data into one level:
         !  1) Both levels are at same pressure level within precision of pressure,
         !  so merge data from both levels into one measurement.
         !  2) If both of the observations are surface reports, then the pressure
         !  may be different, but the height = terrain elevation, and the two
         !  heights are equal.

         CALL merge_measurements ( next1%meas , next2%meas , best )

         !  Advance the pointers.

         IF ( .NOT. ASSOCIATED ( current ) ) THEN
            ! are at the head of the output linked list;
            ! already have list1 => next1 so do nothing
         ELSE
            current%next => next1
         END IF 
         current => next1         !  set so current points to next output node 
         next1 => next1%next      !  get next node in list1
         delete_it => next2       !  record location of next2 to delete it
         next2 => next2%next      !  get next node in list2

         !  The bypassed observation can be deleted.

         DEALLOCATE ( delete_it )

         !  Because of the way that the data is merged (allowing the surface data
         !  to be recognized through the elevation = height), there may arise
         !  conditions that allow replicated pressure surfaces.  Those instances
         !  are what we check for in the next two IF blocks.

         duplicates_list1 : DO WHILE ( ASSOCIATED ( next1 ) ) 

            IF      (    ( eps_equal ( current%meas%pressure%data , next1%meas%pressure%data , 1. ) ) &
                                        .OR.  &
                      (  ( eps_equal ( current%meas%height%data   , next1%meas%height%data , .1 ) ) .AND. &
                         ( .NOT. eps_equal ( current%meas%height%data   , missing_r , 1. ) ) .AND. &
                         ( eps_equal ( info1%elevation            , next1%meas%height%data , .1 ) ) ) ) THEN
      
               CALL merge_measurements ( current%meas , next1%meas , best )
      
               !  Advance the next1 pointer, kill the old location.
      
               delete_it => next1       !  record location of next1 to delete it
               next1 => next1%next      !  get next node in list1
      
               !  The bypassed observation can be deleted.
      
               DEALLOCATE ( delete_it )
           
               !  We need to continue checking for more duplicates.
   
               CYCLE duplicates_list1
   
            ELSE
   
               !  There are no more duplicates.
   
               EXIT duplicates_list1 
   
            END IF
         
         END DO duplicates_list1

         duplicates_list2 : DO WHILE ( ASSOCIATED ( next2 ) ) 

            IF      (    ( eps_equal ( current%meas%pressure%data , next2%meas%pressure%data , 1. ) ) &
                                        .OR.  &
                      (  ( eps_equal ( current%meas%height%data   , next2%meas%height%data , .1 ) ) .AND. &
                         ( .NOT. eps_equal ( current%meas%height%data   , missing_r , 1. ) ) .AND. &
                         ( eps_equal ( info2%elevation            , next2%meas%height%data , .1 ) ) ) ) THEN
      
               CALL merge_measurements ( current%meas , next2%meas , best )
      
               !  Advance the next2 pointer, kill the old location.
      
               delete_it => next2       !  record location of next2 to delete it
               next2 => next2%next      !  get next node in list2
      
               !  The bypassed observation can be deleted.
      
               DEALLOCATE ( delete_it )

            ELSE
   
               !  There are no more duplicates.
   
               EXIT duplicates_list2 
   
            END IF
         
         END DO duplicates_list2

      ELSE IF ( next1%meas%pressure%data .LT. next2%meas%pressure%data ) THEN

         ! Link node from list2 in current list.                  

         IF ( .NOT. ASSOCIATED ( current ) ) THEN
            ! are at the head of the output list
            list1 => next2
         ELSE
            current%next => next2
         END IF
         current => next2
         next2 => next2%next

      ELSE

         ! Link node from list1 into the current list.

         IF ( .NOT. ASSOCIATED ( current ) ) THEN
            ! are at the head of the output list; list1 already points to next1
            !  have list1 => next1 so do nothing
         ELSE
            current%next => next1
         END IF
         current => next1
         next1 => next1%next

      END IF

   END DO still_associated

   !  The end of either list1 or list2 was reached.  The list that is still
   !  associated (not finished), still has data in the list tail.  Have the
   !  current list include that tail.  If both lists are exhausted, nullify 
   !  the last pointer.
   
   IF      ( ASSOCIATED ( next2 ) ) THEN
      current%next => next2
   ELSE IF ( ASSOCIATED ( next1 ) ) THEN
      current%next => next1
   ELSE
      NULLIFY ( current%next )
   END IF

END SUBROUTINE link_levels

!
! ------------------------------------------------------------------------

SUBROUTINE make_date ( date , time , date_time_char ) 

!  This routine takes an 8-digit date (YYYYMMDD) and a
!  6-digit time (HHmmss) and converts it to a 24-digit
!  string (YYYY-MM-DD_HH:mm:ss.ffff).

   INTEGER , INTENT(IN) :: date , &
                           time

   CHARACTER (LEN=24) , INTENT(OUT) :: date_time_char


   !  Local data.

   INTEGER :: year , month , day , hour , minute , second , fraction

   year = date / 10000
   month = ( date - year*10000 ) / 100
   day   =  date - year*10000 - month*100 

   hour = time / 10000
   minute = ( time - hour*10000 ) / 100
   second = time - hour*10000 - minute*100

   fraction = 0

   WRITE ( date_time_char , &
           FMT = '(I4.4,"-",I2.2,"-",I2.2,"_",I2.2,":",I2.2,":",I2.2,".",I4.4) ' ) &
                   year , month ,    day ,    hour  , minute , second , fraction

END SUBROUTINE make_date

!
! ------------------------------------------------------------------------

SUBROUTINE get_date ( date , time , date_time_char ) 

!  This routine takes a 19-digit date string  (YYYY-MM-DD_HH:mm:ss) and converts
!  it to an 8-digit date (YYYYMMDD) and a 6-digit time (HHmmss).

   CHARACTER (LEN=19) , INTENT(IN) :: date_time_char

   INTEGER , INTENT(OUT) :: date , &
                            time

   !  Local data.

   INTEGER :: year , month , day , hour , minute , second

   READ ( date_time_char , &
           FMT = '(I4.4,1x,I2.2,1x,I2.2,1x,I2.2,1x,I2.2,1x,I2.2) ' ) &
                   year , month , day , hour , minute , second

   date = year * 10000 + month  * 100 + day
   time = hour * 10000 + minute * 100 + second

END SUBROUTINE get_date

!
! ------------------------------------------------------------------------

SUBROUTINE merge_measurements ( first ,second , best )

!  This takes two measurements that have been found to be at the same
!  pressure level and takes the best data from each.  Criterion for 
!  determining which to keep is which has better quality control value.

   IMPLICIT NONE 

   TYPE ( meas_data ) , INTENT ( INOUT )       :: first
   TYPE ( meas_data ) , INTENT ( IN )          :: second
   INTEGER , INTENT ( IN )                     :: best

   CALL keep_best ( first%pressure     , second%pressure     , best ) 
   CALL keep_best ( first%height       , second%height       , best ) 
   CALL keep_best ( first%temperature  , second%temperature  , best ) 
   CALL keep_best ( first%dew_point    , second%dew_point    , best ) 
   CALL keep_best ( first%speed        , second%speed        , best ) 
   CALL keep_best ( first%direction    , second%direction    , best ) 
   CALL keep_best ( first%u            , second%u            , best ) 
   CALL keep_best ( first%v            , second%v            , best ) 
   CALL keep_best ( first%rh           , second%rh           , best ) 
   CALL keep_best ( first%thickness    , second%thickness    , best ) 
   
END SUBROUTINE merge_measurements

!
! --------------------------------------------------------------------------

SUBROUTINE merge_obs ( first , second )

!  Reports 'first' and 'second' have been found to have same location and
!  time, therefore they must be merged and one of them discarded.
!  The result of merge is put in 'first' report; second is discarded.
!  If either has data and other has 'missing', keep data; if both have data, take
!  data from one with greatest num_vld_fld or fewest 'num_error'.

   IMPLICIT NONE 

   TYPE ( report ) , INTENT ( INOUT )            :: first , &
                                                    second  
   INTEGER                                       :: best

   INCLUDE 'error.inc'
   INTERFACE
      INCLUDE 'error.int'
   END INTERFACE

   IF      ( first%info%num_vld_fld .GT. second%info%num_vld_fld ) THEN
      best = 1
   ELSE IF ( first%info%num_vld_fld .LT. second%info%num_vld_fld ) THEN
      best = 2
   ELSE IF ( first%info%num_error   .LT. second%info%num_error   ) THEN
      best = 1
   ELSE IF ( first%info%num_error   .GT. second%info%num_error   ) THEN
      best = 2
   ELSE IF ( first%info%num_warning .LT. second%info%num_warning ) THEN
      best = 1
   ELSE IF ( first%info%num_warning .GT. second%info%num_warning ) THEN
      best = 2
   ELSE IF ( first%info%seq_num     .GT. second%info%seq_num     ) THEN
      best = 1
   ELSE IF ( first%info%seq_num     .LT. second%info%seq_num     ) THEN
      best = 2
   ELSE
      best = 1
      error_number =3321001
      error_message(1:31) = 'merge_obs                      '
      error_message(32:)  = ' Arbitrarily assuming "first" obs &
      &is better than "second" for ' // & 
      TRIM ( first%location%name ) // '  ' // &
      TRIM ( first%location%id ) // '.'
      fatal = .false.
      listing = .false.
!     CALL error_handler ( error_number , error_message , &
!     fatal , listing )
   END IF

   !  Will put all useful information in first report; discard second
   !  report.  Update report being kept with the num_vld_fld, num_error,
   !  num_warnings, etc, from best.

   IF ( best .EQ. 2 ) THEN
      first%info = second%info
   END IF

   ! Now look at all terrestrial fields, keeping the best values.

! foo
!  IF ( .NOT. ( first%ground .EQ. second%ground ) ) THEN
   IF ( .NOT. ground_eq ( first%ground , second%ground ) ) THEN
      CALL keep_best ( first%ground%slp         , second%ground%slp         , best )
      CALL keep_best ( first%ground%ref_pres    , second%ground%ref_pres    , best )
      CALL keep_best ( first%ground%ground_t    , second%ground%ground_t    , best )
      CALL keep_best ( first%ground%sst         , second%ground%sst         , best )
      CALL keep_best ( first%ground%psfc        , second%ground%psfc        , best )
      CALL keep_best ( first%ground%precip      , second%ground%precip      , best )
      CALL keep_best ( first%ground%t_max       , second%ground%t_max       , best )
      CALL keep_best ( first%ground%t_min       , second%ground%t_min       , best )
      CALL keep_best ( first%ground%t_min_night , second%ground%t_min_night , best )
      CALL keep_best ( first%ground%p_tend03    , second%ground%p_tend03    , best )
      CALL keep_best ( first%ground%p_tend24    , second%ground%p_tend24    , best )
      CALL keep_best ( first%ground%cloud_cvr   , second%ground%cloud_cvr   , best )
      CALL keep_best ( first%ground%ceiling     , second%ground%ceiling     , best )
   END IF

   !  Merge data at different levels, starting at ground (or lowest level).  On return
   !  all of linked list from second is deallocated, so not much additional memory
   !  is used by keeping the absorbed observation.  The info types are provided so 
   !  that it can be determined if these are both surface observations.

   CALL link_levels ( first%surface , second%surface , &
   first%info , second%info , best )

END SUBROUTINE merge_obs

!
! -----------------------------------------------------------------------------

RECURSIVE SUBROUTINE merge_sort ( obs , index , low , high )

!  This is the recursive part of the merge sort routine.  This routine recurses
!  all the way down to a list of 2 items, which it sorts with no more
!  recursion.  It would be faster to stop at list of 10 or so and use n^2
!  sort like selection-sort; but this is easier and plenty fast.  This does 
!  not actually sort the observations, it creates an index to a sorted list
!  list of ovbservations.

   IMPLICIT NONE

   TYPE ( report )        , INTENT ( IN    ) , DIMENSION ( : ) :: obs 
   INTEGER                , INTENT ( INOUT ) , DIMENSION ( : ) :: index
   INTEGER                , INTENT ( IN    )                   :: low        , &
                                                                  high

   INTEGER , ALLOCATABLE                     , DIMENSION ( : ) :: tmp_old
   INTEGER                                                        mid        , &
                                                                  current    , & 
                                                                  first_ndx  , & 
                                                                  second_ndx , &
                                                                  temp

   !  The list is either small (2 items), or too big.  If it is too large, it
   !  is recursively broken in half.

   break_it_down : IF ( high - low .GE. 2 ) THEN

      !  Half-way point for the list that is too large to handle.

      mid = ( low + high ) / 2

      !  Have merge_sort work on both halves of the list.

      CALL merge_sort ( obs , index , low , mid )
      CALL merge_sort ( obs , index , mid + 1 , high )

      !  Have now got two sorted lists.  They are now to be merged into one list.
      !  Allocate  additional temporary space for this list merging.

      ALLOCATE ( tmp_old ( low : mid ) )

      !  The lower half of the list is stored in the new space.

      tmp_old = index ( low : mid )

      !  Initialize the list sorting counters.

      first_ndx = low
      second_ndx = mid + 1
      current = low

      !  The lists are sorted until one of the indices is out of bounds.

      sort_two_groups : DO WHILE ( first_ndx .LE. mid .AND. second_ndx .LE. high )

         !  Either the first index or the second index is the lowest.  The
         !  current counter keeps the correct choice, and the used index is
         !  incremented.

! foo
!        IF ( obs(tmp_old(first_ndx)) .LT. obs(index(second_ndx)) ) THEN
         IF ( compare ( obs(tmp_old(first_ndx)) , obs(index(second_ndx)) ) ) THEN
            index(current) = tmp_old(first_ndx)
            first_ndx = first_ndx + 1
         ELSE 
            index(current) = index(second_ndx)
            second_ndx = second_ndx + 1
         ENDIF

         !  Increment the cureent counter, which is the sorted list of observations.

         current = current + 1

      END DO sort_two_groups

      !  After the above sort, must still copy the tail end of tmp_old (if any is 
      !  left over from above).  There is NO need to move the tail end if the tail 
      !  is in index -- it is already there.

      tail_of_first_group : DO WHILE ( first_ndx .LE. mid )
         index(current) = tmp_old(first_ndx)
         current = current + 1
         first_ndx = first_ndx + 1
      END DO tail_of_first_group

      !  Free up the temporary memory before continuing.

      DEALLOCATE ( tmp_old )

   ELSE break_it_down

      !  Now we have a list that is EASY to sort, just 1 or 2 elements.  Either the
      !  indices need to be swapped (this if test), or they are already ordered.

! foo
!     small_enough : IF ( obs(index(high)) .LT. obs(index(low)) ) THEN
      small_enough : IF ( compare ( obs(index(high)) , obs(index(low)) ) ) THEN
         temp = index(low)
         index(low) = index(high)
         index(high) = temp       
      END IF small_enough

   END IF break_it_down

END SUBROUTINE merge_sort

!
! -------------------------------------------------------------------------

SUBROUTINE output_obs ( obs , unit , file_name , num_obs , out_opt, forinput, &
                        new_file, qc_flag_keep, remove_unverified, num_klvls, pressure, &
                        met_file)

!  Take the array of observations and write them including measurements
!  at all levels.  The two options (out_opt and forinput) are described
!  below.

   !  If ( out_opt is 0 ) , write everything
   !                > 0   , write only non-discard data
   !                < 0   , write only discarded data  
   
   !  If ( forinput is true ) output can be pipe back for input.


   IMPLICIT NONE

   INCLUDE 'netcdf.inc'

   TYPE ( report ) , INTENT (INOUT), DIMENSION ( : ) :: obs
   INTEGER , INTENT ( IN )                           :: num_obs
   INTEGER , INTENT ( IN )                           :: out_opt   
   INTEGER , INTENT ( IN )                           :: unit
   CHARACTER ( LEN = * ) , INTENT ( IN )             :: file_name
   LOGICAL , INTENT ( IN )                           :: forinput, new_file

   INTEGER                                           :: i , iout, num_klvls
   INTEGER                                           :: qc_flag_keep
   LOGICAL                                           :: remove_unverified        
   TYPE ( measurement ) , POINTER                    :: next
   TYPE ( meas_data   )                              :: end_meas

   REAL, DIMENSION(num_klvls)                        :: pressure, pres_hPA
   LOGICAL                                           :: keep_data
   LOGICAL                                           :: OBS_data=.FALSE.
   LOGICAL                                           :: is_sounding      
   LOGICAL                                           :: no_qc_done      
   INTEGER                                           :: true_num_obs
   INTEGER                                           :: track_surface_data

   INTEGER                                           :: ncid, dimid, varid, vardims
   INTEGER, DIMENSION(2)                             :: ishape, start, count
   INTEGER                                           :: ilev
   CHARACTER (LEN=132)                               :: nfile
   REAL, ALLOCATABLE, DIMENSION(:)                   :: press_nc, z_nc, temp_nc, td_nc
   REAL, ALLOCATABLE, DIMENSION(:)                   :: spd_nc, dir_nc, u_nc, v_nc, rh_nc 

   INTEGER, ALLOCATABLE, DIMENSION(:)                :: press_qc_nc, z_qc_nc, temp_qc_nc, td_qc_nc
   INTEGER, ALLOCATABLE, DIMENSION(:)                :: spd_qc_nc, dir_qc_nc, u_qc_nc, v_qc_nc, rh_qc_nc             
   INTEGER                                           :: int_sounding, int_bogus, int_discard, iout_nc

   CHARACTER ( LEN = 132 )                           :: met_file
   INTEGER                                           :: met_ncid
   INTEGER                                           :: idummy
   REAL                                              :: rdummy
   

   end_meas%pressure%data    = end_data_r
   end_meas%height%data      = end_data_r
   end_meas%temperature%data = end_data_r
   end_meas%dew_point%data   = end_data_r
   end_meas%speed%data       = end_data_r
   end_meas%direction%data   = end_data_r
   end_meas%u%data           = end_data_r
   end_meas%v%data           = end_data_r
   end_meas%rh%data          = end_data_r
   end_meas%thickness%data   = end_data_r
   end_meas%pressure%qc      = end_data  
   end_meas%height%qc        = end_data  
   end_meas%temperature%qc   = end_data  
   end_meas%dew_point%qc     = end_data  
   end_meas%speed%qc         = end_data  
   end_meas%direction%qc     = end_data  
   end_meas%u%qc             = end_data  
   end_meas%v%qc             = end_data  
   end_meas%rh%qc            = end_data  
   end_meas%thickness%qc     = end_data  

   pres_hPA = pressure * 100.


   IF ( new_file ) THEN
     OPEN ( UNIT = unit , FILE = file_name ,  ACTION = 'write' , FORM = 'formatted' , &
            STATUS = 'replace' )
   ELSE
     OPEN ( UNIT = unit , FILE = file_name ,  ACTION = 'write' , FORM = 'formatted' , &
            STATUS = 'unknown', POSITION = 'append' )
   ENDIF
   OBS_data = .FALSE.
   IF ( file_name(1:10) == "OBS_DOMAIN" ) OBS_data = .TRUE.

   IF ( .not. OBS_data ) THEN
     nfile = trim(adjustl(file_name))//".nc"
     CALL check(nf_create(nfile, NF_CLOBBER, ncid))
     dimid = 1
     CALL check (nf_def_dim(ncid, "report", NF_UNLIMITED, dimid))
     dimid = 2
     CALL check(nf_def_dim(ncid, "lev", 300, dimid))
     dimid = 3
     CALL check(nf_def_dim(ncid, "date_len", 14, dimid))
     dimid = 4
     CALL check(nf_def_dim(ncid, "char_len", 40, dimid))

     vardims = 2
     ishape(1)=3
     ishape(2)=1
     varid = 1
     CALL check(nf_def_var(ncid, "date" , NF_CHAR, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",33,"Observation Date - YYYYMMDDHHmmss"))
     ishape(1)=4
     varid = varid + 1
     CALL check(nf_def_var(ncid, "stationID" , NF_CHAR, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",13,"ID of Station"))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "stationNAME" , NF_CHAR, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",15,"Name of Station"))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "platform" , NF_CHAR, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",33,"Description of Measurement Device"))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "source" , NF_CHAR, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",21,"Source of Observation"))

     vardims = 1
     ishape(1)=1
     varid = varid + 1
     CALL check(nf_def_var(ncid, "lon" , NF_FLOAT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",36,"Station Longitude - east is positive"))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "lat" , NF_FLOAT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",36,"Station Latitude - north is positive"))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "elevation" , NF_FLOAT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",17,"Station Elevation"))
       CALL check(nf_put_att_text(ncid,varid,"units",1,"m"))
       CALL check(nf_put_att_real(ncid,varid,"_FillValue",NF_REAL,1,-888888.00))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "slp" , NF_FLOAT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",18,"Sea Level Pressure"))
       CALL check(nf_put_att_text(ncid,varid,"units",2,"Pa"))
       CALL check(nf_put_att_real(ncid,varid,"_FillValue",NF_REAL,1,-888888.00))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "slp_QC" , NF_INT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",27,"Sea Level Pressure QC Flag"))
       CALL check(nf_put_att_int(ncid,varid,"_FillValue",NF_INT,1,-888888))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "ref_pressure" , NF_FLOAT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",24,"Reference Pressure Level"))
       CALL check(nf_put_att_text(ncid,varid,"units",2,"Pa"))
       CALL check(nf_put_att_real(ncid,varid,"_FillValue",NF_REAL,1,-888888.00))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "ref_pressure_QC" , NF_INT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",32,"Reference Pressure Level QC Flag"))
       CALL check(nf_put_att_int(ncid,varid,"_FillValue",NF_INT,1,-888888))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "levels" , NF_INT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",36,"Number of Valid Levels in the Report"))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "sounding" , NF_INT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",55,"1/0 Multiple Level Sounding or Single Level Observation"))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "bogus" , NF_INT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",33,"1/0 Bogus Report or Normal Report"))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "discard" , NF_INT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",31,"1/0 Keep or Discard Observation"))

     vardims = 2
     ishape(1)=2
     ishape(2)=1
     varid = varid + 1
     CALL check(nf_def_var(ncid, "pressure" , NF_FLOAT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",8,"Pressure"))
       CALL check(nf_put_att_text(ncid,varid,"units",2,"Pa"))
       CALL check(nf_put_att_real(ncid,varid,"_FillValue",NF_REAL,1,-888888.00))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "pressure_QC" , NF_INT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",16,"Pressure QC Flag"))
       CALL check(nf_put_att_int(ncid,varid,"_FillValue",NF_INT,1,-888888))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "height" , NF_FLOAT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",12,"Height (MSL)"))
       CALL check(nf_put_att_text(ncid,varid,"units",1,"m"))
       CALL check(nf_put_att_real(ncid,varid,"_FillValue",NF_REAL,1,-888888.00))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "height_QC" , NF_INT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",14,"Height QC Flag"))
       CALL check(nf_put_att_int(ncid,varid,"_FillValue",NF_INT,1,-888888))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "temperature" , NF_FLOAT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",11,"Temperature"))
       CALL check(nf_put_att_text(ncid,varid,"units",1,"K"))
       CALL check(nf_put_att_real(ncid,varid,"_FillValue",NF_REAL,1,-888888.00))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "temperature_QC" , NF_INT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",19,"Temperature QC Flag"))
       CALL check(nf_put_att_int(ncid,varid,"_FillValue",NF_INT,1,-888888))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "dew_point" , NF_FLOAT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",20,"Dewpoint Temperature"))
       CALL check(nf_put_att_text(ncid,varid,"units",1,"K"))
       CALL check(nf_put_att_real(ncid,varid,"_FillValue",NF_REAL,1,-888888.00))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "dew_point_QC" , NF_INT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",28,"Dewpoint Temperature QC Flag"))
       CALL check(nf_put_att_int(ncid,varid,"_FillValue",NF_INT,1,-888888))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "speed" , NF_FLOAT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",10,"Wind Speed"))
       CALL check(nf_put_att_text(ncid,varid,"units",3,"m/s"))
       CALL check(nf_put_att_real(ncid,varid,"_FillValue",NF_REAL,1,-888888.00))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "speed_QC" , NF_INT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",18,"Wind Speed QC Flag"))
       CALL check(nf_put_att_int(ncid,varid,"_FillValue",NF_INT,1,-888888))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "direction" , NF_FLOAT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",14,"Wind Direction"))
       CALL check(nf_put_att_text(ncid,varid,"units",7,"degrees"))
       CALL check(nf_put_att_real(ncid,varid,"_FillValue",NF_REAL,1,-888888.00))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "direction_QC" , NF_INT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",22,"Wind Direction QC Flag"))
       CALL check(nf_put_att_int(ncid,varid,"_FillValue",NF_INT,1,-888888))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "u" , NF_FLOAT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",19,"U Component of Wind"))
       CALL check(nf_put_att_text(ncid,varid,"units",3,"m/s"))
       CALL check(nf_put_att_real(ncid,varid,"_FillValue",NF_REAL,1,-888888.00))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "u_QC" , NF_INT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",27,"U Component of Wind QC Flag"))
       CALL check(nf_put_att_int(ncid,varid,"_FillValue",NF_INT,1,-888888))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "v" , NF_FLOAT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",19,"V Component of Wind"))
       CALL check(nf_put_att_text(ncid,varid,"units",3,"m/s"))
       CALL check(nf_put_att_real(ncid,varid,"_FillValue",NF_REAL,1,-888888.00))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "v_QC" , NF_INT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",27,"V Component of Wind QC Flag"))
       CALL check(nf_put_att_int(ncid,varid,"_FillValue",NF_INT,1,-888888))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "rh" , NF_FLOAT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",17,"Relative Humidity"))
       CALL check(nf_put_att_text(ncid,varid,"units",1,"%"))
       CALL check(nf_put_att_real(ncid,varid,"_FillValue",NF_REAL,1,-888888.00))
     varid = varid + 1
     CALL check(nf_def_var(ncid, "rh_QC" , NF_INT, vardims, ishape, varid))
       CALL check(nf_put_att_text(ncid,varid,"description",25,"Relative Humidity QC Flag"))
       CALL check(nf_put_att_int(ncid,varid,"_FillValue",NF_INT,1,-888888))

     IF ( file_name(1:10) == "qc_obs_raw" ) &
       CALL check(nf_put_att_text(ncid,NF_GLOBAL,"TITLE",24,"Observational Data - raw"))
     IF ( file_name(1:10) == "qc_obs_use" ) &
       CALL check(nf_put_att_text(ncid,NF_GLOBAL,"TITLE",25,"Observational Data - used"))

     CALL check(nf_open(met_file, NF_NOCLOBBER, met_ncid))
     CALL check(nf_get_att_int(met_ncid, NF_GLOBAL, "WEST-EAST_GRID_DIMENSION", idummy))
       CALL check(nf_put_att_int(ncid, NF_GLOBAL, "WEST-EAST_GRID_DIMENSION", NF_INT,1,idummy))
     CALL check(nf_get_att_int(met_ncid, NF_GLOBAL, "SOUTH-NORTH_GRID_DIMENSION", idummy))
       CALL check(nf_put_att_int(ncid, NF_GLOBAL, "SOUTH-NORTH_GRID_DIMENSION", NF_INT,1,idummy))
     CALL check(nf_get_att_real(met_ncid, NF_GLOBAL, "DX", rdummy))
       CALL check(nf_put_att_real(ncid, NF_GLOBAL, "DX", NF_REAL,1,rdummy))
     CALL check(nf_get_att_real(met_ncid, NF_GLOBAL, "DY", rdummy))
       CALL check(nf_put_att_real(ncid, NF_GLOBAL, "DY", NF_REAL,1,rdummy))
     CALL check(nf_get_att_real(met_ncid, NF_GLOBAL, "CEN_LAT", rdummy))
       CALL check(nf_put_att_real(ncid, NF_GLOBAL, "CEN_LAT", NF_REAL,1,rdummy))
     CALL check(nf_get_att_real(met_ncid, NF_GLOBAL, "CEN_LON", rdummy))
       CALL check(nf_put_att_real(ncid, NF_GLOBAL, "CEN_LON", NF_REAL,1,rdummy))
     CALL check(nf_get_att_real(met_ncid, NF_GLOBAL, "TRUELAT1", rdummy))
       CALL check(nf_put_att_real(ncid, NF_GLOBAL, "TRUELAT1", NF_REAL,1,rdummy))
     CALL check(nf_get_att_real(met_ncid, NF_GLOBAL, "TRUELAT2", rdummy))
       CALL check(nf_put_att_real(ncid, NF_GLOBAL, "TRUELAT2", NF_REAL,1,rdummy))
     CALL check(nf_get_att_real(met_ncid, NF_GLOBAL, "MOAD_CEN_LAT", rdummy))
       CALL check(nf_put_att_real(ncid, NF_GLOBAL, "MOAD_CEN_LAT", NF_REAL,1,rdummy))
     CALL check(nf_get_att_real(met_ncid, NF_GLOBAL, "STAND_LON", rdummy))
       CALL check(nf_put_att_real(ncid, NF_GLOBAL, "STAND_LON", NF_REAL,1,rdummy))
     CALL check(nf_get_att_real(met_ncid, NF_GLOBAL, "POLE_LAT", rdummy))
       CALL check(nf_put_att_real(ncid, NF_GLOBAL, "POLE_LAT", NF_REAL,1,rdummy))
     CALL check(nf_get_att_real(met_ncid, NF_GLOBAL, "POLE_LON", rdummy))
       CALL check(nf_put_att_real(ncid, NF_GLOBAL, "POLE_LON", NF_REAL,1,rdummy))
     CALL check(nf_get_att_int(met_ncid, NF_GLOBAL, "MAP_PROJ", idummy))
       CALL check(nf_put_att_int(ncid, NF_GLOBAL, "MAP_PROJ", NF_INT,1,idummy))
     CALL check(nf_get_att_int(met_ncid, NF_GLOBAL, "grid_id", idummy))
       CALL check(nf_put_att_int(ncid, NF_GLOBAL, "grid_id", NF_INT,1,idummy))
     CALL check(nf_get_att_int(met_ncid, NF_GLOBAL, "parent_id", idummy))
       CALL check(nf_put_att_int(ncid, NF_GLOBAL, "parent_id", NF_INT,1,idummy))
     CALL check(nf_get_att_int(met_ncid, NF_GLOBAL, "i_parent_start", idummy))
       CALL check(nf_put_att_int(ncid, NF_GLOBAL, "i_parent_start", NF_INT,1,idummy))
     CALL check(nf_get_att_int(met_ncid, NF_GLOBAL, "j_parent_start", idummy))
       CALL check(nf_put_att_int(ncid, NF_GLOBAL, "j_parent_start", NF_INT,1,idummy))
     CALL check(nf_get_att_int(met_ncid, NF_GLOBAL, "i_parent_end", idummy))
       CALL check(nf_put_att_int(ncid, NF_GLOBAL, "i_parent_end", NF_INT,1,idummy))
     CALL check(nf_get_att_int(met_ncid, NF_GLOBAL, "j_parent_end", idummy))
       CALL check(nf_put_att_int(ncid, NF_GLOBAL, "j_parent_end", NF_INT,1,idummy))
     CALL check(nf_get_att_int(met_ncid, NF_GLOBAL, "parent_grid_ratio", idummy))
       CALL check(nf_put_att_int(ncid, NF_GLOBAL, "parent_grid_ratio", NF_INT,1,idummy))

     CALL check(nf_close(met_ncid))

     CALL check(nf_enddef(ncid))
   ENDIF

   iout = 0
   iout_nc = 0

   DO i = 1 , num_obs
      true_num_obs = 0
      no_qc_done = .FALSE.

      IF (   out_opt .EQ. 0                                   .OR. &
           ( out_opt .GT. 0 .AND. .NOT. obs(i)%info%discard ) .OR. &
           ( out_opt .LT. 0 .AND.       obs(i)%info%discard ) ) THEN

         iout = iout + 1
         IF ( .NOT. forinput ) write(unit,*) '**************** Next Observation *******************'

         next => obs(i)%surface
         is_sounding = obs(i)%info%is_sound
         IF ( is_sounding ) THEN
           DO WHILE ( ASSOCIATED ( next ) )
            IF ( next%meas%pressure%qc .eq. no_qc_possible ) no_qc_done = .TRUE.
            next => next%next
           END DO
         ELSE
           true_num_obs = 1
         ENDIF

         next => obs(i)%surface
         IF ((obs(i)%info%platform(1:11) == 'FM-88 SATOB') .OR. &
                  (obs(i)%info%platform(1:11) == 'FM-97 AIREP')) THEN
           is_sounding = .TRUE.
           obs(i)%info%is_sound = .TRUE.
         ENDIF
         !!!!IF ( obs(i)%info%num_vld_fld == 1 .AND. (obs(i)%info%elevation .ne. missing) .AND. &
         IF ( obs(i)%info%num_vld_fld == 1  .AND. &
            ( next%meas%height%data .eq. obs(i)%info%elevation ) ) THEN
           is_sounding = .FALSE.
           obs(i)%info%is_sound = .FALSE.
         ENDIF
         IF ( is_sounding ) THEN
           DO WHILE ( ASSOCIATED ( next ) )
              IF ( no_qc_done ) THEN
                next%meas%pressure%qc     = no_qc_possible
                next%meas%height%qc       = no_qc_possible
                next%meas%temperature%qc  = no_qc_possible
                next%meas%dew_point%qc    = no_qc_possible
                next%meas%u%qc            = no_qc_possible
                next%meas%v%qc            = no_qc_possible
                next%meas%rh%qc           = no_qc_possible
                next%meas%speed%qc        = no_qc_possible
                next%meas%direction%qc    = no_qc_possible
                next%meas%thickness%qc    = no_qc_possible
              ENDIF
              IF ( remove_unverified ) THEN
                keep_data = any ( ( pres_hPA .eq. next%meas%pressure%data ) )
                IF ( (obs(i)%info%elevation .ne. missing) .AND. & 
                      (next%meas%height%data .eq. obs(i)%info%elevation) ) keep_data = .TRUE.
                !!IF ( is_sounding .AND. next%meas%pressure%qc .gt. 4 )  keep_data = .TRUE.
                IF ( keep_data ) true_num_obs = true_num_obs + 1
              ELSE
                true_num_obs = true_num_obs + 1
              ENDIF
              next => next%next
           END DO
         ELSE
           true_num_obs = 1
         ENDIF

        !! 2012-12-17 cB - If surface and we don't want data above a set qc value - discard
        !!                 Also discard soundings with no data
        if ( .NOT. is_sounding .and. (obs(i)%surface%meas%pressure%qc >= qc_flag_keep) ) & 
                      obs(i)%info%discard = .TRUE.
        if (       is_sounding .and. true_num_obs.eq.0) obs(i)%info%discard = .TRUE.

         if ( obs(i)%ground%slp%qc   >= qc_flag_keep ) then
           obs(i)%ground%slp%data = missing_r
           obs(i)%ground%slp%qc = missing
         endif
         if ( obs(i)%ground%slp%data   == missing_r ) then
           obs(i)%ground%slp%qc = missing
         endif

         IF ( OBS_data ) THEN
           IF ( .NOT. obs(i)%info%discard ) THEN
           IF ( is_sounding ) THEN
             WRITE ( UNIT = unit , FMT='(1x,A14)' ) obs(i)%valid_time%date_char(1:14)
             WRITE ( UNIT = unit , FMT='(2x,2(F9.4,1x))' ) obs(i)%location%latitude, obs(i)%location%longitude
             IF ( obs(i)%location%id(1:5) == "US un" ) THEN
               WRITE ( UNIT = unit , FMT='(2x,2(A40,3x))' )    &
                  "-----                                   " , &
                  "Unknown Station                         "
             ELSE
               WRITE ( UNIT = unit , FMT='(2x,2(A40,3x))' ) obs(i)%location%id, obs(i)%location%name
             ENDIF
             WRITE ( UNIT = unit , FMT='(2x,2(A16,2x),F8.0,2x,2(L4,2x),I5)' )       &
               obs(i)%info%platform, obs(i)%info%source, obs(i)%info%elevation, &
               is_sounding, obs(i)%info%bogus, true_num_obs
           ENDIF
           ENDIF
         ELSE
           IF ( .NOT. obs(i)%info%discard ) THEN
             WRITE ( UNIT = unit , FMT = rpt_format ) &
                obs(i)%location , obs(i)%info , obs(i)%valid_time , obs(i)%ground
           ENDIF
         ENDIF
    
         track_surface_data = 0

         next => obs(i)%surface

       IF ( .not. OBS_data ) THEN
         ilev = 0
         if(allocated(press_nc)) deallocate(press_nc)
         allocate(press_nc(true_num_obs))
         if(allocated(press_qc_nc)) deallocate(press_qc_nc)
         allocate(press_qc_nc(true_num_obs))
         if(allocated(z_nc)) deallocate(z_nc)
         allocate(z_nc(true_num_obs))
         if(allocated(z_qc_nc)) deallocate(z_qc_nc)
         allocate(z_qc_nc(true_num_obs))
         if(allocated(temp_nc)) deallocate(temp_nc)
         allocate(temp_nc(true_num_obs))
         if(allocated(temp_qc_nc)) deallocate(temp_qc_nc)
         allocate(temp_qc_nc(true_num_obs))
         if(allocated(td_nc)) deallocate(td_nc)
         allocate(td_nc(true_num_obs))
         if(allocated(td_qc_nc)) deallocate(td_qc_nc)
         allocate(td_qc_nc(true_num_obs))
         if(allocated(spd_nc)) deallocate(spd_nc)
         allocate(spd_nc(true_num_obs))
         if(allocated(spd_qc_nc)) deallocate(spd_qc_nc)
         allocate(spd_qc_nc(true_num_obs))
         if(allocated(dir_nc)) deallocate(dir_nc)
         allocate(dir_nc(true_num_obs))
         if(allocated(dir_qc_nc)) deallocate(dir_qc_nc)
         allocate(dir_qc_nc(true_num_obs))
         if(allocated(u_nc)) deallocate(u_nc)
         allocate(u_nc(true_num_obs))
         if(allocated(u_qc_nc)) deallocate(u_qc_nc)
         allocate(u_qc_nc(true_num_obs))
         if(allocated(v_nc)) deallocate(v_nc)
         allocate(v_nc(true_num_obs))
         if(allocated(v_qc_nc)) deallocate(v_qc_nc)
         allocate(v_qc_nc(true_num_obs))
         if(allocated(rh_nc)) deallocate(rh_nc)
         allocate(rh_nc(true_num_obs))
         if(allocated(rh_qc_nc)) deallocate(rh_qc_nc)
         allocate(rh_qc_nc(true_num_obs))

         press_nc    = missing_r
         press_qc_nc = missing
         z_nc        = missing_r
         z_qc_nc     = missing
         temp_nc     = missing_r
         temp_qc_nc  = missing
         td_nc       = missing_r
         td_qc_nc    = missing
         spd_nc      = missing_r
         spd_qc_nc   = missing
         dir_nc      = missing_r
         dir_qc_nc   = missing
         u_nc        = missing_r
         u_qc_nc     = missing
         v_nc        = missing_r
         v_qc_nc     = missing
         rh_nc       = missing_r
         rh_qc_nc    = missing
       ENDIF  

         DO WHILE ( ASSOCIATED ( next ) )
            if ( obs(i)%info%discard ) exit 
            keep_data = any ( ( pres_hPA .eq. next%meas%pressure%data ) )
            IF ( (obs(i)%info%elevation .ne. missing) .AND. (next%meas%height%data .eq. obs(i)%info%elevation) ) keep_data = .TRUE.
            !!IF ( is_sounding .AND. next%meas%pressure%qc .gt. 4 )  keep_data = .TRUE.

            !!! Make sure no surface ob has more than one entry
            track_surface_data = track_surface_data + 1
            IF ( .not. is_sounding .AND. track_surface_data .gt. 1 ) exit

            IF ( .NOT. keep_data .AND. next%meas%pressure%qc == 0 ) THEN
              next%meas%pressure%qc     = missing
              IF (next%meas%height%data       .NE. missing_r .AND. next%meas%height%qc       ==  0 ) &
                  next%meas%height%qc       = missing
              IF (next%meas%temperature%data  .NE. missing_r .AND. next%meas%temperature%qc  ==  0 ) &
                  next%meas%temperature%qc  = missing
              IF (next%meas%dew_point%data    .NE. missing_r .AND. next%meas%dew_point%qc    ==  0 ) &
                  next%meas%dew_point%qc    = missing
              IF (next%meas%u%data            .NE. missing_r .AND. next%meas%u%qc            ==  0 ) &
                  next%meas%u%qc            = missing
              IF (next%meas%v%data            .NE. missing_r .AND. next%meas%v%qc            ==  0 ) &
                  next%meas%v%qc            = missing
              IF (next%meas%rh%data           .NE. missing_r .AND. next%meas%rh%qc           ==  0 ) &
                  next%meas%rh%qc           = missing
              IF (obs(i)%ground%slp%data      .NE. missing_r .AND. obs(i)%ground%slp%qc      ==  0 ) &
                  obs(i)%ground%slp%qc      = missing
              IF (obs(i)%ground%ref_pres%data .NE. missing_r .AND. obs(i)%ground%ref_pres%qc ==  0 ) &
                  obs(i)%ground%ref_pres%qc = missing
              IF (next%meas%speed%data        .NE. missing_r .AND. next%meas%speed%qc        ==  0 ) &
                  next%meas%speed%qc        = missing
              IF (next%meas%direction%data    .NE. missing_r .AND. next%meas%direction%qc    ==  0 ) &
                  next%meas%direction%qc    = missing
              IF (obs(i)%ground%precip%data   .NE. missing_r .AND. obs(i)%ground%precip%qc   ==  0 ) &
                  obs(i)%ground%precip%qc   = missing
            ENDIF

            if ( next%meas%pressure%qc    >= qc_flag_keep ) then
              next%meas%pressure%data = missing_r
              next%meas%pressure%qc = missing
            endif
            if ( next%meas%height%qc      >= qc_flag_keep ) then
              next%meas%height%data = missing_r
              next%meas%height%qc = missing
            endif
            if ( next%meas%temperature%qc >= qc_flag_keep ) then
              next%meas%temperature%data = missing_r
              next%meas%temperature%qc = missing
            endif
            if ( next%meas%dew_point%qc   >= qc_flag_keep ) then
              next%meas%dew_point%data = missing_r
              next%meas%dew_point%qc = missing
            endif
            if ( next%meas%speed%qc       >= qc_flag_keep ) then
              next%meas%speed%data = missing_r
              next%meas%speed%qc = missing
            endif
            if ( next%meas%direction%qc   >= qc_flag_keep ) then
              next%meas%direction%data = missing_r
              next%meas%direction%qc = missing
            endif
            if ( next%meas%u%qc           >= qc_flag_keep ) then
              next%meas%u%data = missing_r
              next%meas%u%qc = missing
            endif
            if ( next%meas%v%qc           >= qc_flag_keep ) then
              next%meas%v%data = missing_r
              next%meas%v%qc = missing
            endif
            if ( next%meas%rh%qc          >= qc_flag_keep ) then
              next%meas%rh%data = missing_r
              next%meas%rh%qc = missing
            endif
            if ( next%meas%thickness%qc   >= qc_flag_keep ) then
              next%meas%thickness%data = missing_r
              next%meas%thickness%qc = missing
            endif

            IF (next%meas%pressure%data     == missing_r) next%meas%pressure%qc     = missing
            IF (next%meas%height%data       == missing_r) next%meas%height%qc       = missing
            IF (next%meas%temperature%data  == missing_r) next%meas%temperature%qc  = missing
            IF (next%meas%dew_point%data    == missing_r) next%meas%dew_point%qc    = missing
            IF (next%meas%u%data            == missing_r) next%meas%u%qc            = missing
            IF (next%meas%v%data            == missing_r) next%meas%v%qc            = missing
            IF (next%meas%rh%data           == missing_r) next%meas%rh%qc           = missing
            IF (obs(i)%ground%slp%data      == missing_r) obs(i)%ground%slp%qc      = missing
            IF (obs(i)%ground%ref_pres%data == missing_r) obs(i)%ground%ref_pres%qc = missing
            IF (next%meas%speed%data        == missing_r) next%meas%speed%qc        = missing
            IF (next%meas%direction%data    == missing_r) next%meas%direction%qc    = missing
            IF (obs(i)%ground%precip%data   == missing_r) obs(i)%ground%precip%qc   = missing
            IF (next%meas%thickness%data    == missing_r) next%meas%thickness%qc    = missing

            IF ( is_sounding ) THEN
              IF ( (keep_data .AND. remove_unverified) .OR. (.NOT. remove_unverified) ) THEN
                IF ( OBS_data ) THEN
                  WRITE ( UNIT = unit , FMT='(1x,6(F11.3,1x,F11.3,1x))' )        &
                    next%meas%pressure%data,    real(next%meas%pressure%qc),     &
                    next%meas%height%data,      real(next%meas%height%qc),       &
                    next%meas%temperature%data, real(next%meas%temperature%qc),  &
                    next%meas%u%data,           real(next%meas%u%qc),            &
                    next%meas%v%data,           real(next%meas%v%qc),            &
                    next%meas%rh%data,          real(next%meas%rh%qc)
                ELSE
                  WRITE ( UNIT = unit , FMT = meas_format )  next%meas
                ENDIF
                IF ( file_name(1:10) == "qc_obs_use" ) THEN
                  ilev = ilev + 1
                  press_nc(ilev)    = next%meas%pressure%data
                  press_qc_nc(ilev) = next%meas%pressure%qc
                  z_nc(ilev)        = next%meas%height%data
                  z_qc_nc(ilev)     = next%meas%height%qc
                  temp_nc(ilev)     = next%meas%temperature%data
                  temp_qc_nc(ilev)  = next%meas%temperature%qc
                  td_nc(ilev)       = next%meas%dew_point%data
                  td_qc_nc(ilev)    = next%meas%dew_point%qc
                  spd_nc(ilev)      = next%meas%speed%data
                  spd_qc_nc(ilev)   = next%meas%speed%qc
                  dir_nc(ilev)      = next%meas%direction%data
                  dir_qc_nc(ilev)   = next%meas%direction%qc
                  u_nc(ilev)        = next%meas%u%data
                  u_qc_nc(ilev)     = next%meas%u%qc
                  v_nc(ilev)        = next%meas%v%data
                  v_qc_nc(ilev)     = next%meas%v%qc
                  rh_nc(ilev)       = next%meas%rh%data
                  rh_qc_nc(ilev)    = next%meas%rh%qc
                ENDIF
              ENDIF
            ELSE
              IF ( (keep_data .AND. remove_unverified) .OR. (.NOT. remove_unverified) ) THEN
                IF ( OBS_data ) THEN
                  WRITE ( UNIT = unit , FMT='(1x,A14)' ) obs(i)%valid_time%date_char(1:14)
                  WRITE ( UNIT = unit , FMT='(2x,2(F9.4,1x))' ) obs(i)%location%latitude, obs(i)%location%longitude
                  IF ( obs(i)%location%id(1:5) == "US un" ) THEN
                    WRITE ( UNIT = unit , FMT='(2x,2(A40,3x))' )    &
                       "-----                                   " , &
                       "Unknown Station                         "
                  ELSE
                    WRITE ( UNIT = unit , FMT='(2x,2(A40,3x))' ) obs(i)%location%id, obs(i)%location%name
                  ENDIF
                  WRITE ( UNIT = unit , FMT='(2x,2(A16,2x),F8.0,2x,2(L4,2x),I5)' )       &
                    obs(i)%info%platform, obs(i)%info%source, obs(i)%info%elevation, &
                    is_sounding, obs(i)%info%bogus, true_num_obs
                  WRITE ( UNIT = unit , FMT='(1x,9(F11.3,1x,F11.3,1x))' )        &
                    obs(i)%ground%slp%data,      real(obs(i)%ground%slp%qc),     &
                    obs(i)%ground%ref_pres%data, real(obs(i)%ground%ref_pres%qc),&
                    next%meas%height%data,       real(next%meas%height%qc),      &
                    next%meas%temperature%data,  real(next%meas%temperature%qc), &
                    next%meas%u%data,            real(next%meas%u%qc),           &
                    next%meas%v%data,            real(next%meas%v%qc),           &
                    next%meas%rh%data,           real(next%meas%rh%qc),          &
                    next%meas%pressure%data,     real(next%meas%pressure%qc),    &
                    obs(i)%ground%precip%data,   real(obs(i)%ground%precip%qc)
                ELSE
                  WRITE ( UNIT = unit , FMT = meas_format )  next%meas
                ENDIF
                IF ( file_name(1:10) == "qc_obs_use" ) THEN
                  ilev = ilev + 1
                  press_nc(ilev)    = next%meas%pressure%data
                  press_qc_nc(ilev) = next%meas%pressure%qc
                  z_nc(ilev)        = next%meas%height%data
                  z_qc_nc(ilev)     = next%meas%height%qc
                  temp_nc(ilev)     = next%meas%temperature%data
                  temp_qc_nc(ilev)  = next%meas%temperature%qc
                  td_nc(ilev)       = next%meas%dew_point%data
                  td_qc_nc(ilev)    = next%meas%dew_point%qc
                  spd_nc(ilev)      = next%meas%speed%data
                  spd_qc_nc(ilev)   = next%meas%speed%qc
                  dir_nc(ilev)      = next%meas%direction%data
                  dir_qc_nc(ilev)   = next%meas%direction%qc
                  u_nc(ilev)        = next%meas%u%data
                  u_qc_nc(ilev)     = next%meas%u%qc
                  v_nc(ilev)        = next%meas%v%data
                  v_qc_nc(ilev)     = next%meas%v%qc
                  rh_nc(ilev)       = next%meas%rh%data
                  rh_qc_nc(ilev)    = next%meas%rh%qc
                ENDIF
              ENDIF
            ENDIF

            
            IF ( file_name(1:10) == "qc_obs_raw" ) THEN
              ilev = ilev + 1
              press_nc(ilev)    = next%meas%pressure%data
              press_qc_nc(ilev) = next%meas%pressure%qc
              z_nc(ilev)        = next%meas%height%data
              z_qc_nc(ilev)     = next%meas%height%qc
              temp_nc(ilev)     = next%meas%temperature%data
              temp_qc_nc(ilev)  = next%meas%temperature%qc
              td_nc(ilev)       = next%meas%dew_point%data
              td_qc_nc(ilev)    = next%meas%dew_point%qc
              spd_nc(ilev)      = next%meas%speed%data
              spd_qc_nc(ilev)   = next%meas%speed%qc
              dir_nc(ilev)      = next%meas%direction%data
              dir_qc_nc(ilev)   = next%meas%direction%qc
              u_nc(ilev)        = next%meas%u%data
              u_qc_nc(ilev)     = next%meas%u%qc
              v_nc(ilev)        = next%meas%v%data
              v_qc_nc(ilev)     = next%meas%v%qc
              rh_nc(ilev)       = next%meas%rh%data
              rh_qc_nc(ilev)    = next%meas%rh%qc
            ENDIF
            next => next%next
         END DO

         IF ( .not. OBS_data .and. .not. obs(i)%info%discard ) THEN

           iout_nc = iout_nc + 1
           int_sounding = 0
           IF (is_sounding) int_sounding = 1
           int_bogus = 0
           IF (obs(i)%info%bogus) int_bogus = 1
           int_discard = 0
           IF (obs(i)%info%discard) int_discard = 1

           start = 1
           count = 1
           start(2) = iout_nc
           count(1) = 14
           CALL check (nf_put_vara_text(ncid, 1,start,count,obs(i)%valid_time%date_char(1:14)))
           count(1) = 40
           CALL check (nf_put_vara_text(ncid, 2,start,count,obs(i)%location%id))
           CALL check (nf_put_vara_text(ncid, 3,start,count,obs(i)%location%name))
           CALL check (nf_put_vara_text(ncid, 4,start,count,obs(i)%info%platform))
           CALL check (nf_put_vara_text(ncid, 5,start,count,obs(i)%info%source))

           start = 1
           count = 1
           start(1) = iout_nc
           CALL check (nf_put_vara_real(ncid, 6,start,count,obs(i)%location%longitude))
           CALL check (nf_put_vara_real(ncid, 7,start,count,obs(i)%location%latitude))
           CALL check (nf_put_vara_real(ncid, 8,start,count,obs(i)%info%elevation))
           CALL check (nf_put_vara_real(ncid, 9,start,count,obs(i)%ground%slp%data))
           CALL check (nf_put_vara_int(ncid,10,start,count,obs(i)%ground%slp%qc))
           CALL check (nf_put_vara_real(ncid,11,start,count,obs(i)%ground%ref_pres%data))
           CALL check (nf_put_vara_int(ncid,12,start,count,obs(i)%ground%ref_pres%qc))
           CALL check (nf_put_vara_int(ncid,13,start,count,true_num_obs))
           CALL check (nf_put_vara_int(ncid, 14,start,count,int_sounding))
           CALL check (nf_put_vara_int(ncid, 15,start,count,int_bogus))
           CALL check (nf_put_vara_int(ncid, 16,start,count,int_discard))

           start = 1
           count = 1
           start(2) = iout_nc
           count(1) = true_num_obs
           CALL check (nf_put_vara_real(ncid,17,start,count,press_nc))
           CALL check (nf_put_vara_int(ncid,18,start,count,press_qc_nc))
           CALL check (nf_put_vara_real(ncid,19,start,count,z_nc))
           CALL check (nf_put_vara_int(ncid,20,start,count,z_qc_nc))
           CALL check (nf_put_vara_real(ncid,21,start,count,temp_nc))
           CALL check (nf_put_vara_int(ncid,22,start,count,temp_qc_nc))
           CALL check (nf_put_vara_real(ncid,23,start,count,td_nc))
           CALL check (nf_put_vara_int(ncid,24,start,count,td_qc_nc))
           CALL check (nf_put_vara_real(ncid,25,start,count,spd_nc))
           CALL check (nf_put_vara_int(ncid,26,start,count,spd_qc_nc))
           CALL check (nf_put_vara_real(ncid,27,start,count,dir_nc))
           CALL check (nf_put_vara_int(ncid,28,start,count,dir_qc_nc))
           CALL check (nf_put_vara_real(ncid,29,start,count,u_nc))
           CALL check (nf_put_vara_int(ncid,30,start,count,u_qc_nc))
           CALL check (nf_put_vara_real(ncid,31,start,count,v_nc))
           CALL check (nf_put_vara_int(ncid,32,start,count,v_qc_nc))
           CALL check (nf_put_vara_real(ncid,33,start,count,rh_nc))
           CALL check (nf_put_vara_int(ncid,34,start,count,rh_qc_nc))

           WRITE ( UNIT = unit , FMT = meas_format ) end_meas
           WRITE ( UNIT = unit , FMT = end_format ) obs(i)%info%num_vld_fld, &
              obs(i)%info%num_error, obs(i)%info%num_warning
         ENDIF
         IF ( .NOT. forinput ) &
            write(unit,*) 'End of measurements for observation ' , i

         END IF 

   END DO

   IF ( .NOT. forinput ) THEN
      write(unit,*) '======================================================='
      write(unit,*) 'Total Number of Measurements output ' , iout
   ENDIF

   !  This routine may be called again, with the same unit number, so CLOSE
   !  up the file so everything is handled cleanly.

   CLOSE ( unit )
   IF ( .not. OBS_data ) THEN
     CALL check(nf_close(ncid))
   ENDIF

END SUBROUTINE output_obs

!
!---------------------------------------------------------------------------

SUBROUTINE read_measurements ( file_num , surface , location , bad_data , error , &
height , pressure , iew , jns , levels , map_projection , elevation )

!  This routine reads in 'measurements' at as many levels as there are in
!  the report, then stops and returns when an end-of-measurements flag is
!  found.  If any reads produce error, return error code which causes entire
!  observation to be discarded (ob is not discarded on eof error).

   IMPLICIT NONE 

   INTEGER , INTENT ( IN )                      :: file_num   ! file to read  
   TYPE ( measurement ) , POINTER               :: surface    ! ptr to 1st msmt
   TYPE ( location_type ) , INTENT ( IN )       :: location   ! 5 digit ID, name, etc
   LOGICAL , INTENT ( IN )                      :: bad_data   ! read, but do not store
   INTEGER , INTENT ( OUT )                     :: error      ! err and type 

   INTEGER , INTENT ( IN )                          :: levels
   REAL    , INTENT ( IN ) , DIMENSION ( levels )   :: pressure
   INTEGER                                          :: iew , jns , map_projection
   REAL , DIMENSION ( jns , iew , levels )          :: height

   CHARACTER ( LEN = 32 ) , PARAMETER    :: sub_name = 'read_measurements'
   INTEGER                                      :: meas_count
   INTEGER                                      :: io_error
   TYPE ( measurement ) , POINTER               :: current

   CHARACTER ( LEN = 40 )                       :: location_id , &
                                                   location_name
   REAL , INTENT(IN)                            :: elevation

   INCLUDE 'error.inc'
   INTERFACE
      INCLUDE 'error.int'
   END INTERFACE

   !  Initialize dummy pointers and counters and observation names, and such.

   ALLOCATE ( current )
   NULLIFY ( current%next )
   NULLIFY ( surface )
   error = ok
   meas_count = 0
   location_id   = TRIM ( location%id )
   location_name = TRIM ( location%name )

   !  This loop continues until either an error occurs, or until the end of
   !  the measurement tag is found (the graceful exit).

   read_meas: DO 

      !  Currently, this read puts in 12 pairs of data, a real observation
      !  value and the accompanying QC flag. 

      READ ( file_num , IOSTAT = io_error , FMT = meas_format )  &
             current%meas
    
      !  An error < 0 means the end of the file (usually), and an error > 0
      !  is just a broken read.  Describe the read error so that the calling
      !  routine knows what happened, then exit this loop (which is exiting
      !  this routine, basically).

      IF ( io_error .GT. 0 ) THEN
         error = read_err
!        CLOSE ( file_num ) 
         EXIT read_meas
      ELSE IF ( io_error .LT. 0 ) THEN
         error = eof_err
         CLOSE ( file_num ) 
         EXIT read_meas
      END IF

      !  If we know a priori that this data is bad, no tests are necessary on
      !  the various flags values.

      bad_loop_1 : IF ( .NOT. bad_data ) THEN
   
         !  A successful read, yahoo!  As the data may not have the flags set up the
         !  way we want, go through directly after this read and make sure that
         !  any special values are all set to missing.
   
         IF ( ( current%meas%pressure%data    .GT. ( undefined1_r - 1. ) )  .OR. &
              ( current%meas%pressure%data    .LT. ( undefined2_r + 1. ) ) ) THEN
            current%meas%pressure%data    = missing_r
         END IF
         IF ( ( current%meas%height%data      .GT. ( undefined1_r - 1. ) )  .OR. &
              ( current%meas%height%data      .LT. ( undefined2_r + 1. ) ) ) THEN
            current%meas%height%data      = missing_r
         END IF
         IF ( ( current%meas%temperature%data .GT. ( undefined1_r - 1. ) )  .OR. &
              ( current%meas%temperature%data .LT. ( undefined2_r + 1. ) ) ) THEN
            current%meas%temperature%data = missing_r
         END IF
         IF (   current%meas%temperature%data .GT. (    99999.0   - 1. ) )   THEN
            current%meas%temperature%data = missing_r
         END IF
         IF ( ( current%meas%dew_point%data   .GT. ( undefined1_r - 1. ) )  .OR. &
              ( current%meas%dew_point%data   .LT. ( undefined2_r + 1. ) ) ) THEN
            current%meas%dew_point%data   = missing_r
         END IF
         IF ( ( current%meas%speed%data       .GT. ( undefined1_r - 1. ) )  .OR. &
              ( current%meas%speed%data       .LT. ( undefined2_r + 1. ) ) ) THEN
            current%meas%speed%data       = missing_r
         END IF
         IF ( ( current%meas%direction%data   .GT. ( undefined1_r - 1. ) )  .OR. &
              ( current%meas%direction%data   .LT. ( undefined2_r + 1. ) ) ) THEN
            current%meas%direction%data   = missing_r
         END IF
         IF ( ( current%meas%u%data           .GT. ( undefined1_r - 1. ) )  .OR. &
              ( current%meas%u%data           .LT. ( undefined2_r + 1. ) ) ) THEN
            current%meas%u%data           = missing_r
         END IF
         IF ( ( current%meas%v%data           .GT. ( undefined1_r - 1. ) )  .OR. &
              ( current%meas%v%data           .LT. ( undefined2_r + 1. ) ) ) THEN
            current%meas%v%data           = missing_r
         END IF
         IF ( ( current%meas%rh%data          .GT. ( undefined1_r - 1. ) )  .OR. &
              ( current%meas%rh%data          .LT. ( undefined2_r + 1. ) ) ) THEN
            current%meas%rh%data          = missing_r
         END IF
         IF ( ( current%meas%thickness%data   .GT. ( undefined1_r - 1. ) )  .OR. &
              ( current%meas%thickness%data   .LT. ( undefined2_r + 1. ) ) ) THEN
            current%meas%thickness%data   = missing_r
         END IF

      END IF bad_loop_1

      !  The data we just read in could have been the flag for the end of the measurement.
      !  This is the graceful way to exit this routine.  If this is the end of the
      !  measurement section for this observation, set all of the data to the same end of
      !  measurement value, just in case there were some stray unset values in the
      !  generating program.

      IF ( eps_equal ( current%meas%pressure%data , end_data_r , 1. ) .OR. &
           eps_equal ( current%meas%height%data   , end_data_r , 1. ) ) THEN
         current%meas%pressure%data    = end_data_r
         current%meas%height%data      = end_data_r
         current%meas%temperature%data = end_data_r
         current%meas%dew_point%data   = end_data_r
         current%meas%speed%data       = end_data_r
         current%meas%direction%data   = end_data_r
         current%meas%u%data           = end_data_r
         current%meas%v%data           = end_data_r
         current%meas%rh%data          = end_data_r
         current%meas%thickness%data   = end_data_r
         current%meas%pressure%qc      = end_data  
         current%meas%height%qc        = end_data  
         current%meas%temperature%qc   = end_data  
         current%meas%dew_point%qc     = end_data  
         current%meas%speed%qc         = end_data  
         current%meas%direction%qc     = end_data  
         current%meas%u%qc             = end_data  
         current%meas%v%qc             = end_data  
         current%meas%rh%qc            = end_data  
         current%meas%thickness%qc     = end_data  
         error = ok
         EXIT read_meas
      END IF

      !  If this is bad data, we needed to make sure that the ending measurement
      !  is the famous end_data flag so that we hit a correct exit from this
      !  loop and routine.  We can just cycle the read loop again.

      IF ( bad_data ) THEN
         CYCLE read_meas
      END IF

      !  Later analysis will want to use the pressure as the vertical coordinate.
      !  If the pressure is missing, this causes a problem.  We can compute it 
      !  through a couple of simple ways.  

      !  1)  With the available geopotential height field, we can compute a 
      !  fairly accurate pressure if the height is available from the observation.

      IF ( ( eps_equal ( current%meas%pressure%data , missing_r , 1. ) ) .OR. &
           ( eps_equal ( current%meas%pressure%data ,        0. , 1. ) ) ) THEN
         CALL height_to_pres ( height , pressure , iew , jns , levels , map_projection , &
         location%latitude , location%longitude , current%meas )
      END IF 

      !  2)  If the pressure is missing, first try to compute it through a standard
      !  atmosphere calculation, using the available height.  This should only happen 
      !  if the previous computation failed - the target being when the missing
      !  pressure is the surface pressure, so that we cannot do a vertical interpolation.

      IF ( ( eps_equal ( current%meas%pressure%data , missing_r , 1. ) ) .OR. &
           ( eps_equal ( current%meas%pressure%data ,        0. , 1. ) ) ) THEN
         CALL height_to_pres_old ( current%meas )
      END IF

      !  3) If the pressure is still missing (which implies that the height was
      !  missing as well), we can still go after the pressure from the standard
      !  atmosphere, but using the temperature (only low levels, though).

!     IF ( ( eps_equal ( current%meas%pressure%data , missing_r , 1. ) ) .OR. &
!          ( eps_equal ( current%meas%pressure%data ,        0. , 1. ) ) ) THEN
!        CALL temp_to_pres ( current%meas )
!     END IF

      !  If the pressure is still missing after trying our darnedest to get 
      !  something from the height and temperature, we admit failure, and discard
      !  this measurement (cycle to the next measurement level).

      IF ( ( eps_equal ( current%meas%pressure%data , missing_r , 1. ) ) .OR. &
           ( eps_equal ( current%meas%pressure%data ,        0. , 1. ) ) ) THEN
         current%meas%pressure%data = missing_r
         CYCLE read_meas
      END IF 

      !  Compute diagnostics of derived variables (u, v, RH).  At this point, where
      !  we first compute u,v, they are rotated to the map grid.  Since the
      !  analysis fields of horizontal wind components are similarly rotated, this
      !  is a consistency that is maintained throughout the code.  The speed and
      !  direction are the original met-data; u and v are particular to the
      !  domain.

      CALL diagnostics ( current%meas , location%longitude )

      !  Increment count of measurements correctly read in (these are the levels).

      meas_count = meas_count + 1

      !  Since it seems that everything went ok, insert this measurement ordered 
      !  by pressure.

      CALL insert_at ( surface , current , elevation )

      !  Allocate space for another measurement, so we can go try and read another level
      !  in this routine.  Initialize it to pointing to nothing.

      ALLOCATE ( current )
      NULLIFY ( current%next )

   END DO read_meas  

   !  The last allocated measurement is not used (no matter how loop is exited)
   !  so deallocate space.

   DEALLOCATE ( current )

   !  If unable to read in at least one measurement, return error so that
   !  entire observation is discarded.  If this was bad data, we forced it
   !  to skip over the observations without storing any data.  That will be
   !  handled in the calling routine.

   IF ( ( meas_count .LT. 1  ) .AND. &
        ( error      .EQ. ok ) .AND. &
        ( .NOT. bad_data     ) ) THEN
      error = no_data
   END IF 

   !  This is some diagnostic print-out to state the problems encountered.  Anything
   !  that is not expected issues a message.  Any read errors mean to throw away
   !  the observation, though if there was no data, there was nothing to toss
   !  out anyways.  If the error condition is not recognized, the program will
   !  stop from this routine.

   SELECT CASE ( error )

      CASE ( eof_err )
         CALL error_handler ( 33221 , sub_name // ' Found EOF, expected ' &
                            // 'measurement.  Continuing.  '  &
                            // TRIM(location_id) // ' ' // TRIM(location_name) , &
                            .false. , .false. )

      CASE ( read_err )
         CALL error_handler ( 33222 , sub_name // ' Error in measurement read.  ' &
                           // ' Discarding entire observation and continuing. ' &
                            // TRIM(location_id) // ' ' // TRIM(location_name) , &
                            .false. , .false. )
         CALL dealloc_meas ( surface )

      CASE ( no_data , ok )

      CASE DEFAULT

         CALL error_handler ( 33225 , sub_name // ' Internal error:  bad ' &
                           // 'error number. ' , .true. , .false. )

   END SELECT

END SUBROUTINE read_measurements

!
!-----------------------------------------------------------------------------

SUBROUTINE read_observations ( file_name , file_num , obs , n_obs , &
total_number_of_obs , fatal_if_exceed_max_obs , print_out_found_obs , &
height , pressure , iew , jns , levels , map_projection )

!  This routine opens file 'file_name' and reads all observations, and 
!  measurements at all levels, from the file into the 'obs' observation array.
!  For each observation, calls read_measurements to read all measurements for
!  that one observation.

   USE date_pack

   IMPLICIT NONE

   CHARACTER ( LEN=132 ) , INTENT ( IN ):: file_name  ! file name string
   INTEGER , INTENT ( IN )              :: file_num   ! file unit number
   INTEGER , INTENT ( OUT )             :: n_obs
   INTEGER , INTENT ( IN )              :: total_number_of_obs
   LOGICAL , INTENT ( IN )              :: fatal_if_exceed_max_obs
   LOGICAL , INTENT ( IN )              :: print_out_found_obs
   TYPE ( report ), INTENT ( OUT ) , &
                    DIMENSION ( : )     :: obs        ! array of 'report'

   INTEGER , INTENT ( IN )                          :: levels
   REAL    , INTENT ( IN ) , DIMENSION ( levels )   :: pressure
   INTEGER                                          :: iew , jns , map_projection
   REAL , DIMENSION ( jns , iew , levels )          :: height

   CHARACTER ( LEN = 32 ) , PARAMETER   :: proc_name = 'read_observations '
   INTEGER                              :: io_error   ! holds ret-val from IO
   INTEGER                              :: obs_num
   INTEGER                              :: error_ret  ! type of error         

   integer  :: iid
   INTEGER                              :: num_empty , num_outside , bad_count
   LOGICAL                              :: outside_window
   TYPE(meas_data)                      :: dummy_middle

   INCLUDE 'error.inc'
   INTERFACE
      INCLUDE 'error.int'
   END INTERFACE

   !  Initialize a couple of counters for how many observations are either
   !  empty or outside of the domain.

   num_empty = 0
   num_outside = 0

   !  Open file for reading; handle any errors in open by quitting since
   !  this is probably a simple-to-fix user mistake.

   OPEN ( UNIT = file_num , FILE = file_name , FORM = 'FORMATTED'  , &
          ACTION = 'READ' , IOSTAT = io_error ) 

   !  If there were any troubles on the previous open, let's just fail straight
   !  away to reduce the anxiety.

   IF ( io_error .NE. 0 ) THEN
      CALL error_handler ( 331001, proc_name(1:32) //  & 
             'Unable to open observations file. ' , .true. , .false. )
   ENDIF

   !  While we are not at the end of the observation file, keep reading data.

   obs_num = 1
   read_obs : DO                                 

      !  This is an array that we are filling.  Are we beyond that limit yet?

      IF      ( ( obs_num .GT. total_number_of_obs ) .AND. ( fatal_if_exceed_max_obs ) ) THEN
         error_number = 0332111
         error_message(1:31) = 'read_observations              '
         error_message(32:)  = ' Too many obs for the NAMELIST value of max_number_of_obs_nml.'
         fatal = .TRUE.
         listing = .FALSE.
         CALL error_handler ( error_number , error_message , &
         fatal , listing )
      ELSE IF ( ( obs_num .GT. total_number_of_obs ) .AND. ( .NOT. fatal_if_exceed_max_obs ) ) THEN
         error_number = 0332111
         error_message(1:31) = 'read_observations              '
         error_message(32:)  = ' Too many obs for the NAMELIST value of max_number_of_obs_nml.'
         fatal = .FALSE.
         listing = .FALSE.
         CALL error_handler ( error_number , error_message , &
         fatal , listing )
         CLOSE ( file_num ) 
         EXIT read_obs
      END IF

      !  The first read is the "once only" type of information.

      READ ( file_num , IOSTAT = io_error , FMT = rpt_format ) &
             obs(obs_num)%location , &
             obs(obs_num)%info , &
             obs(obs_num)%valid_time , &
             obs(obs_num)%ground 

!          The station name for buoys, metars, ships are in a weird place in the
!          NCAR files, so reset the id.
             iid = index(obs(obs_num)%location%name,'ICAO ')
             if (iid .ne. 0) then
               obs(obs_num)%location%id = obs(obs_num)%location%name(iid+5:iid+9)
             else
               iid = index(obs(obs_num)%location%name,' >>>')
               if (iid .ne. 0) then
                 obs(obs_num)%location%id = obs(obs_num)%location%name(iid+5:iid+9)
               endif
             endif

      !  If there are troubles with the "once only" type of data, we keep trying
      !  until we either come to the end of this report (and cycle) or we come
      !  to the end of all of the data (exit).

      IF ( io_error .GT. 0 ) THEN

         PRINT '(A,A,A,A)','Troubles with first line ', TRIM ( obs(obs_num)%location%id ) , &
         ' ', TRIM ( obs(obs_num)%location%name ) 

         !  Keep track of how many loops we are taking so this is not infinite.

         bad_count = 0

         DO WHILE ( io_error .GE. 0 )
            bad_count = bad_count + 1
            IF ( bad_count .LT. 1000 ) THEN
               dummy_middle%pressure%data = 0
               READ ( file_num , IOSTAT = io_error , FMT = meas_format )  dummy_middle
               IF ( eps_equal ( dummy_middle%pressure%data , end_data_r , 1. ) ) THEN
                  READ ( file_num , IOSTAT = io_error , FMT = end_format ) &
                         obs(obs_num)%info%num_vld_fld , &
                         obs(obs_num)%info%num_error , &  
                         obs(obs_num)%info%num_warning    
                  PRINT '(A)','Starting to READ a new report.'
                  CYCLE read_obs
               END IF
            ELSE
               PRINT '(A)','Too many attempts to read the data correctly.  Exiting read loop.'
               EXIT read_obs
            END IF
         END DO 

         !  While trying to find a good read, we came to the end of the file.  It
         !  happens to the best of us.

         IF ( io_error .LT. 0 ) THEN
            !!!PRINT *, 'Have reached end of observations file. '
            CLOSE ( file_num ) 
            EXIT read_obs
         END IF 

      ELSE IF ( io_error .LT. 0 ) THEN

         !  No errors.  This is the intended way to find the end of data mark.

         !!!PRINT *, 'Have reached end of observations file. '
         CLOSE ( file_num ) 
         EXIT read_obs

      ELSE IF ( io_error .EQ. 0 ) THEN 

         CALL inside_window ( obs(obs_num)%location%latitude  , &
                              obs(obs_num)%location%longitude , &
                              iew , jns , outside_window )

         !  The previous read ("once only" data) was valid.  If any of the data is suspicious,
         !  the easiest place to clean it up is as soon as we read it in, so we do not have
         !  to track through the array or accidently hit undefined values that are not exactly
         !  "our" undefined values.

         !  Sometimes the date and time are ingested in a silly way.  Set the valid time to a 
         !  time guaranteed not to be within the time window.

         IF      ( INDEX ( obs(obs_num)%valid_time%date_char , '*' ) .GT. 0 ) THEN
            obs(obs_num)%valid_time%date_char = '19000101000000'
         ELSE IF ( ( obs(obs_num)%valid_time%date_char( 1: 2) .NE. '19' ) .AND. &
                   ( obs(obs_num)%valid_time%date_char( 1: 2) .NE. '20' ) ) THEN
            obs(obs_num)%valid_time%date_char = '19000101000000'
         ELSE IF ( ( obs(obs_num)%valid_time%date_char( 5: 6) .LT. '01' ) .OR. &
                   ( obs(obs_num)%valid_time%date_char( 5: 6) .GT. '12' ) ) THEN
            obs(obs_num)%valid_time%date_char = '19000101000000'
         ELSE IF ( ( obs(obs_num)%valid_time%date_char( 7: 8) .LT. '01' ) .OR. &
                   ( obs(obs_num)%valid_time%date_char( 7: 8) .GT. '31' ) ) THEN
            obs(obs_num)%valid_time%date_char = '19000101000000'
         ELSE IF ( ( obs(obs_num)%valid_time%date_char( 9:10) .LT. '00' ) .OR. &
                   ( obs(obs_num)%valid_time%date_char( 9:10) .GT. '23' ) ) THEN
            obs(obs_num)%valid_time%date_char = '19000101000000'
         ELSE IF ( ( obs(obs_num)%valid_time%date_char(11:12) .LT. '00' ) .OR. &
                   ( obs(obs_num)%valid_time%date_char(11:12) .GT. '59' ) ) THEN
            obs(obs_num)%valid_time%date_char = '19000101000000'
         ELSE IF ( ( obs(obs_num)%valid_time%date_char(13:14) .LT. '00' ) .OR. &
                   ( obs(obs_num)%valid_time%date_char(13:14) .GT. '59' ) ) THEN
            obs(obs_num)%valid_time%date_char = '19000101000000'
         ELSE IF ((( obs(obs_num)%valid_time%date_char( 5: 6) .EQ. '04' ) .OR. &
                   ( obs(obs_num)%valid_time%date_char( 5: 6) .EQ. '06' ) .OR. &
                   ( obs(obs_num)%valid_time%date_char( 5: 6) .EQ. '09' ) .OR. &
                   ( obs(obs_num)%valid_time%date_char( 5: 6) .EQ. '11' ) ) .AND. &
                   ( obs(obs_num)%valid_time%date_char( 7: 8) .GT. '30' ) ) THEN
            obs(obs_num)%valid_time%date_char = '19000101000000'
         ELSE IF ( ( obs(obs_num)%valid_time%date_char( 5: 6) .EQ. '02' ) .AND. &
                   ( nfeb_ch ( obs(obs_num)%valid_time%date_char( 1: 4) ) .LT.  &
                   obs(obs_num)%valid_time%date_char( 7: 8) ) ) THEN
            obs(obs_num)%valid_time%date_char = '19000101000000'
         END IF

         !  These tests are for the ground type data.  Missing data is OK, but sometimes it
         !  comes in as undefined, which meant bad data.  Set it all to missing.

         IF ( ( obs(obs_num)%ground%slp%data         .GT. ( undefined1_r - 1. ) )  .OR. &
              ( obs(obs_num)%ground%slp%data         .LT. ( undefined2_r + 1. ) ) ) THEN
            obs(obs_num)%ground%slp%data         = missing_r
         END IF
         IF ( ( obs(obs_num)%ground%ref_pres%data    .GT. ( undefined1_r - 1. ) )  .OR. &
              ( obs(obs_num)%ground%ref_pres%data    .LT. ( undefined2_r + 1. ) ) ) THEN
            obs(obs_num)%ground%ref_pres%data    = missing_r
         END IF
         IF ( ( obs(obs_num)%ground%ground_t%data    .GT. ( undefined1_r - 1. ) )  .OR. &
              ( obs(obs_num)%ground%ground_t%data    .LT. ( undefined2_r + 1. ) ) ) THEN
            obs(obs_num)%ground%ground_t%data    = missing_r
         END IF
         IF ( ( obs(obs_num)%ground%sst%data         .GT. ( undefined1_r - 1. ) )  .OR. &
              ( obs(obs_num)%ground%sst%data         .LT. ( undefined2_r + 1. ) ) ) THEN
            obs(obs_num)%ground%sst%data         = missing_r
         END IF
         IF ( ( obs(obs_num)%ground%psfc%data        .GT. ( undefined1_r - 1. ) )  .OR. &
              ( obs(obs_num)%ground%psfc%data        .LT. ( undefined2_r + 1. ) ) ) THEN
            obs(obs_num)%ground%psfc%data        = missing_r
         END IF
         IF ( ( obs(obs_num)%ground%precip%data      .GT. ( undefined1_r - 1. ) )  .OR. &
              ( obs(obs_num)%ground%precip%data      .LT. ( undefined2_r + 1. ) ) ) THEN
            obs(obs_num)%ground%precip%data      = missing_r
         END IF
         IF ( ( obs(obs_num)%ground%t_max%data       .GT. ( undefined1_r - 1. ) )  .OR. &
              ( obs(obs_num)%ground%t_max%data       .LT. ( undefined2_r + 1. ) ) ) THEN
            obs(obs_num)%ground%t_max%data       = missing_r
         END IF
         IF ( ( obs(obs_num)%ground%t_min%data       .GT. ( undefined1_r - 1. ) )  .OR. &
              ( obs(obs_num)%ground%t_min%data       .LT. ( undefined2_r + 1. ) ) ) THEN
            obs(obs_num)%ground%t_min%data       = missing_r
         END IF
         IF ( ( obs(obs_num)%ground%t_min_night%data .GT. ( undefined1_r - 1. ) )  .OR. &
              ( obs(obs_num)%ground%t_min_night%data .LT. ( undefined2_r + 1. ) ) ) THEN
            obs(obs_num)%ground%t_min_night%data = missing_r
         END IF
         IF ( ( obs(obs_num)%ground%p_tend03%data    .GT. ( undefined1_r - 1. ) )  .OR. &
              ( obs(obs_num)%ground%p_tend03%data    .LT. ( undefined2_r + 1. ) ) ) THEN
            obs(obs_num)%ground%p_tend03%data    = missing_r
         END IF
         IF ( ( obs(obs_num)%ground%p_tend24%data    .GT. ( undefined1_r - 1. ) )  .OR. &
              ( obs(obs_num)%ground%p_tend24%data    .LT. ( undefined2_r + 1. ) ) ) THEN
            obs(obs_num)%ground%p_tend24%data    = missing_r
         END IF
         IF ( ( obs(obs_num)%ground%cloud_cvr%data   .GT. ( undefined1_r - 1. ) )  .OR. &
              ( obs(obs_num)%ground%cloud_cvr%data   .LT. ( undefined2_r + 1. ) ) ) THEN
            obs(obs_num)%ground%cloud_cvr%data   = missing_r
         END IF
         IF ( ( obs(obs_num)%ground%ceiling%data     .GT. ( undefined1_r - 1. ) )  .OR. &
              ( obs(obs_num)%ground%ceiling%data     .LT. ( undefined2_r + 1. ) ) ) THEN
            obs(obs_num)%ground%ceiling%data     = missing_r
         END IF

         !  Since no I/O errors, read 1 or more measurements.
         !  Note that obs(obs_num)%surface is pointer to first node in linked list,
         !  so it is initially not pointing to anything.  

         NULLIFY ( obs(obs_num)%surface )
         CALL read_measurements( file_num , obs(obs_num)%surface , &
         obs(obs_num)%location , outside_window , error_ret , & 
         height , pressure , iew , jns , levels , map_projection , obs(obs_num)%info%elevation )

         !  An error in the measurements read is handled in a couple of ways.  A 
         !  flat out error in the read requires the process to start again (cycle
         !  to read_obs).  If there was no data, we need to clean up a bit of stuff,
         !  and read the famous last three integers that have some QC information.

         IF ( error_ret .EQ. read_err ) THEN
            IF ( ASSOCIATED ( obs(obs_num)%surface ) ) THEN
               !  dealloc entire linked list if it exists
               CALL dealloc_meas ( obs(obs_num)%surface)
            END IF
            PRINT '(A,A,A,A)','Troubles with measurement lines ', TRIM ( obs(obs_num)%location%id ) , &
            ' ', TRIM ( obs(obs_num)%location%name ) 
            bad_count = 0
            io_error = 0
            DO WHILE ( io_error .GE. 0 )
               bad_count = bad_count + 1
               IF ( bad_count .LT. 1000 ) THEN
                  READ ( file_num , IOSTAT = io_error , FMT = meas_format ) dummy_middle
                  IF ( eps_equal ( dummy_middle%pressure%data , end_data_r , 1. ) ) THEN
                     READ ( file_num , IOSTAT = io_error , FMT = end_format ) &
                            obs(obs_num)%info%num_vld_fld , &
                            obs(obs_num)%info%num_error , &  
                            obs(obs_num)%info%num_warning    
                     PRINT '(A)','Starting to READ a new report.'
                     CYCLE read_obs
                  END IF
               ELSE
                  PRINT '(A)','Too many attempts to read the measurement data correctly.  Exiting read loop.'
                  EXIT read_obs
               END IF
            END DO 
            IF ( io_error .LT. 0 ) THEN
               EXIT read_obs
            END IF
         ELSE IF ( error_ret .EQ. no_data ) THEN
            READ ( file_num , IOSTAT = io_error , FMT = end_format ) &
            obs(obs_num)%info%num_vld_fld , &
            obs(obs_num)%info%num_error , &  
            obs(obs_num)%info%num_warning    
            IF ( ASSOCIATED ( obs(obs_num)%surface ) ) THEN
               !  dealloc entire linked list if it exists
               CALL dealloc_meas ( obs(obs_num)%surface)
            END IF
            IF ( print_out_found_obs ) THEN
               WRITE ( UNIT = * , FMT = '(A,A,A,A,2F9.3,A)' ) 'Found Name and ID = ' , &
               TRIM ( obs(obs_num)%location%id ) , ' ' , &
               TRIM ( obs(obs_num)%location%name ) , &
               obs(obs_num)%location%latitude , obs(obs_num)%location%longitude , &
               '  *** NO SFC OR UPPER AIR DATA ***'
            END IF 
            num_empty = num_empty + 1
            CYCLE read_obs
         END IF

! We need to skip over obs types that obsgrid doesn't understand. These data
! are present in the NCAR data files that are used by OBSPROC and WRF-VAR.

         IF ( ( obs(obs_num)%info%platform(1:12) .EQ. 'FM-116 GPSRF' ) .OR. &
              ( obs(obs_num)%info%platform(1:11) .EQ. 'FM-86 SATEM' ) .OR. &
              ( obs(obs_num)%info%platform(1:12) .EQ. 'FM-111 GPSPW' ) ) THEN
            READ ( file_num , IOSTAT = io_error , FMT = end_format ) &
            obs(obs_num)%info%num_vld_fld , &
            obs(obs_num)%info%num_error , &
            obs(obs_num)%info%num_warning
            IF ( ASSOCIATED ( obs(obs_num)%surface ) ) THEN
               !  dealloc entire linked list if it exists
               CALL dealloc_meas ( obs(obs_num)%surface)
            END IF
            IF ( print_out_found_obs ) THEN
               WRITE ( UNIT = * , FMT = '(A,A,A,A,2F9.3,A)' ) 'Found Name and ID = ' , &
               TRIM ( obs(obs_num)%location%id ) , ' ' , &
               TRIM ( obs(obs_num)%location%name ) , &
               obs(obs_num)%location%latitude , obs(obs_num)%location%longitude , &
               '  *** INVALID OBSERVATION TYPE ***'
            END IF
            CYCLE read_obs
         END IF

         !  We can compare the observation location with the window that
         !  the analysis will require.  If we are significantly outside
         !  of the analysis window, we toss out the observation.

         IF ( outside_window ) THEN
            READ ( file_num , IOSTAT = io_error , FMT = end_format ) &
            obs(obs_num)%info%num_vld_fld , &
            obs(obs_num)%info%num_error , &  
            obs(obs_num)%info%num_warning    
            IF ( ASSOCIATED ( obs(obs_num)%surface ) ) THEN
               !  dealloc entire linked list if it exists
               CALL dealloc_meas ( obs(obs_num)%surface)
            END IF
            IF ( print_out_found_obs ) THEN
               WRITE ( UNIT = * , FMT = '(A,A,A,A,2F9.3,A)' ) 'Found Name and ID = ' , &
               TRIM ( obs(obs_num)%location%id ) , ' ' , &
               TRIM ( obs(obs_num)%location%name ) , &
               obs(obs_num)%location%latitude , obs(obs_num)%location%longitude , &
               '  *** OUTSIDE ANALYSIS DOMAIN ***'
            END IF 
            num_outside = num_outside + 1
            CYCLE read_obs
         END IF

         !  If this is a ship observation, we need to define the elevation as identical
         !  to the geopotential height, and set the height QC flag to ok.  This is the
         !  only way to get SHIP data into the surface analysis.  Since we are at sea level,
         !  we also set the pressure to equal to the sea level pressure.

         IF ( ( obs(obs_num)%info%platform(1:10) .EQ. 'FM-13 SHIP' ) .AND. &
              ( ASSOCIATED ( obs(obs_num)%surface ) ) ) THEN
            obs(obs_num)%info%elevation             = 0.01
            obs(obs_num)%surface%meas%height%data   = 0.01
            obs(obs_num)%surface%meas%height%qc     = 0
            obs(obs_num)%surface%meas%pressure%data = obs(obs_num)%ground%slp%data
            obs(obs_num)%surface%meas%pressure%qc   = 0
         END IF

         !  This may be wasted print-out, but it is comforting to see.

         IF ( print_out_found_obs ) THEN
            WRITE ( UNIT = * , FMT = '(A,A,A,A,2F9.3)' ) 'Found Name and ID = ' , &
            TRIM ( obs(obs_num)%location%id ) , ' ' , &
            TRIM ( obs(obs_num)%location%name ) , &
            obs(obs_num)%location%latitude , obs(obs_num)%location%longitude
         END IF
      END IF

      !  We have now read ground info and soundings, what follows in the
      !  standard format are three integers describing information gleaned
      !  from the program that generated the observational data.

      READ ( file_num , IOSTAT = io_error , FMT = end_format ) &
             obs(obs_num)%info%num_vld_fld , &
             obs(obs_num)%info%num_error , &  
             obs(obs_num)%info%num_warning    

      !  Once again, after a read, was it successful.  If not toss the whole thing
      !  out (this is handled through the dealloc_meas routine if any upper-air
      !  data was encountered).  Discarding all of the ingested data may be a bit much,
      !  which is why the error print-out is provided.  After the error is processed, 
      !  the reading process is re-started.

      IF ( io_error .NE. 0 ) THEN
         CALL error_handler ( 331003, proc_name // & 
              'Error trying to read last 3 integers in observation ' &
              // TRIM ( obs(obs_num)%location%id ) &               
              // TRIM ( obs(obs_num)%location%name ) // '.' &               
              // ' Discarding entire and continuing.' , .false. , .false. )
         IF ( ASSOCIATED ( obs(obs_num)%surface ) ) THEN
            CALL dealloc_meas ( obs(obs_num)%surface)
         END IF
         CYCLE read_obs
      END IF
   
      !  Before we leave this loop, we make sure the surface level is the first
      !  level in the sounding.

      CALL surf_first ( obs(obs_num)%surface , obs(obs_num)%info%elevation )

      !  Have read observation and its measurement(s) without error
      !  so continue to next observation.

      obs_num = obs_num + 1

   END DO read_obs

   !  The end of the observation file has been reached.  Decrement the counter to
   !  get the total number of observations successfully read by the program.  Output
   !  this information to the outside world.  We can also provide the information
   !  on the observations that are NOT included in the analysis.

   obs_num = obs_num - 1
   WRITE ( UNIT = * , FMT = '(//"Number of observations successfully ingested:       ",i8,"."/     &
                               &"Number of empty observations discarded:             ",i8,"."/     &
                               &"Number of observations discarded outside of domain: ",i8,"."/)' ) &
                               obs_num , num_empty , num_outside
   n_obs = obs_num


END SUBROUTINE read_observations
       
!
! -----------------------------------------------------------------------------

SUBROUTINE sort_obs ( obs , array_size , index )

!  This takes an array of observations ('report' type) and sorts them,
!  or more precisely, creates an index giving their order.  This is the
!  driver that calls recursive merge_sort.

   IMPLICIT NONE

   INTEGER , INTENT ( IN )                            :: array_size
   TYPE ( report ) , INTENT ( IN ) , DIMENSION ( : )  :: obs
   INTEGER , INTENT ( OUT )        , DIMENSION ( : )  :: index
   
   INTEGER                                            :: i  

   !  Iinitialize index for sorting.  These indices are sorted, not the
   !  array of observations.

   index(1:array_size) = (/ ( i , i = 1 , array_size ) /)

   !  On entry, index is numerically sequential (1, 2, ... , array_size).  On
   !  return, index is the array that holds the obs sequential order by location
   !  (function compare).  The "ordered" observations may then be tested for
   !  duplicates (which is the purpose for this ordering).

   CALL merge_sort ( obs , index , 1 , array_size )

END SUBROUTINE sort_obs

!
!---------------------------------------------------------------------------

SUBROUTINE surf_first ( surface , elevation )


!  This routine takes the sounding and makes sure that if a surface
!  level exists, that it is the first level.

   IMPLICIT NONE

   TYPE ( measurement ) ,  POINTER         :: surface
   REAL , INTENT(IN)                       :: elevation

   TYPE ( measurement ) , POINTER          :: current

   !  Um, is there any data at all?

   IF ( ASSOCIATED ( surface ) ) THEN

      !  Alrighty, we have data, so loop through the sounding to see if their is 
      !  a surface observation.  We can't very well have the surface be the first
      !  level if we don't have one.  Also, start looking at location #2 (surface%next)
      !  so that we don't "fix" the ones that aren't broken.

      current  => surface%next

      find_sfc : DO WHILE ( ASSOCIATED ( current ) ) 

         IF ( eps_equal ( current%meas%height%data , elevation , 1. ) ) THEN
            surface => current
            EXIT find_sfc
         END IF

         current => current%next

      END DO find_sfc

   END IF

END SUBROUTINE surf_first

!
!---------------------------------------------------------------------------

SUBROUTINE temp_to_pres ( new )

!  This takes a meas_data in which the pressure is missing, check
!  if there is temperature data. If so, map the temperature to a pressure.

   IMPLICIT NONE

   TYPE ( meas_data )    :: new

   REAL                  :: h

   !  Compute the standard lapse rate pressure from the temperature (CRC F-191). 
   !  We are only in here if both pressure and height are missing, so this 
   !  is a crude approximation.

   temp_is_there : IF ( .NOT. eps_equal ( new%temperature%data , missing_r , 1. ) ) THEN

      IF ( new%temperature%data .GT. 216.65 ) THEN

         new%pressure%data = 1013.25 * ( 288.15 / new%temperature%data )**(-5.255877) * 100
         new%pressure%qc   = p_std_atm_and_temp

      END IF

   END IF temp_is_there

END SUBROUTINE temp_to_pres

!
! ----------------------------------------------------------------------------
  subroutine check(status)

   INCLUDE 'netcdf.inc'
   integer, intent ( in) :: status

    if(status /= nf_noerr) then
      print *, trim(nf_strerror(status))
      stop "Stopped"
    end if
  end subroutine check

END MODULE observation
