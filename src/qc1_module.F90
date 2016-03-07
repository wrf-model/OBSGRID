MODULE qc1

   INCLUDE 'missing.inc'

CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE error_max ( station_id , obs , ship_report , xob , yob ,  &
!BPR BEGIN
station_elevation , &
!BPR END
index , qc_flag , numobs , &
gridded,  iew , jns ,  &
name , max_difference , pressure , local_time , print_error_max )

   USE qc0

!  This routine performs error checks by comparing first guess field
!  and observations at the observing locations (xob, yob).
!  The error-check criteria are set in user-defined namelist,
!  and modified in this routine depending on vertical levels,
!  and local time.  This routine accepts observations and a first
!  guess field for a single level, a single time, and a single
!  variable to be checked.
!
!  After checking, the routine defines an error message for those
!  observations that exceed the error maximum allowance, and
!  'attatch' it to the observation records.

!  Note that x and y are in natural coordinates now.  Need to verify that
!  the x and y are used correctly, not in the traditional MM sense.

   IMPLICIT NONE

   INTEGER                                        :: numobs
   CHARACTER ( LEN =   8 ) , DIMENSION ( : )      :: station_id
   REAL                    , DIMENSION ( : )      :: obs, xob, yob
   !BPR BEGIN
   REAL                    , DIMENSION ( : )      :: station_elevation
   !BPR END
   LOGICAL                 , DIMENSION ( : )      :: ship_report
   INTEGER                 , DIMENSION ( : )      :: index
   INTEGER                 , DIMENSION ( : )      :: qc_flag
   REAL                    , DIMENSION( : , : )   :: gridded
   INTEGER                                        :: iew, jns
   CHARACTER ( LEN =   8 )                        :: name
   INTEGER                 , DIMENSION ( : )      :: local_time
   REAL                                           :: max_difference
   REAL                                           :: pressure
   LOGICAL                                        :: print_error_max

   !  factor:      multiplying factor for error tolerance
   !  err:         difference between observations and first guess fields

   INTEGER                                    :: num,                &
                                                 iob, job
   REAL                                       :: dxob, dyob
   REAL  , DIMENSION (numobs)                 :: err,                &
                                                 factor ,            &
                                                 factored_difference , &
                                                 aob
   REAL                                       :: plevel
   !BPR BEGIN
   !CHARACTER ( LEN = 100 )                    :: message
   CHARACTER ( LEN = 125 )                    :: message
   !BPR END
   REAL                                       :: scale

   CHARACTER ( LEN =   8 ) , DIMENSION ( numobs ) :: station_id_small
   REAL                    , DIMENSION ( numobs ) :: obs_small , xob_small , yob_small
   !BPR BEGIN
   REAL                    , DIMENSION ( numobs ) :: station_elevation_small
   !BPR END
   LOGICAL                 , DIMENSION ( numobs ) :: ship_report_small
   INTEGER                 , DIMENSION ( numobs ) :: index_small
   INTEGER                 , DIMENSION ( numobs ) :: qc_flag_small
   INTEGER                 , DIMENSION ( numobs ) :: local_time_small

   INCLUDE 'error.inc'

   INTERFACE
      INCLUDE 'error.int'
   END INTERFACE

   !  The data coming in is dimensioned as the full size in some previous
   !  PARAMETER statement.  To make the arrays conformable, we make new
   !  small arrays.

   station_id_small(:numobs)  = station_id(:numobs)
   obs_small(:numobs)         = obs(:numobs)
   xob_small(:numobs)         = xob(:numobs)
   yob_small(:numobs)         = yob(:numobs)
   !BPR BEGIN
   station_elevation_small(:numobs) = station_elevation(:numobs)
   !BPR END
   ship_report_small(:numobs) = ship_report(:numobs)
   index_small(:numobs)       = index(:numobs)
   qc_flag_small(:numobs)     = qc_flag(:numobs)
   local_time_small(:numobs)  = local_time(:numobs)

   plevel = pressure

   !   The factor for the error tolerance begins at 1.0.  This value
   !   is modified based on level, local time of day,
   !   variable, and platform type.  

   factor = 1.0
   error_factor_by_field : SELECT CASE ( name )
      CASE ( 'UU      ' , 'VV      ' ) error_factor_by_field 
         IF      ( plevel .LE. 499.9 ) THEN
            factor = 1.5
         ELSE IF ( plevel .LE. 1000.1 ) THEN
            factor = 1.25
         END IF

      CASE ( 'TT      ' ) error_factor_by_field 
         level : IF ( plevel .GT. 1000.1 ) THEN
            WHERE ( ( local_time_small .GE. 1100 ) .AND. &
                    ( local_time_small .LE. 2000 ) ) 
               factor = 1.5
            ELSEWHERE
               factor = 1.25
            ENDWHERE
         ELSE IF ( ( plevel .LE. 1000.1 ) .AND. &
                   ( plevel .GT.  700.0 ) ) THEN
            factor = 0.85 
         ELSE 
            factor = 0.65
         END IF level

         !  If the data came from a ship, provide a tighter tolerance.
         !  This is assuming that the quality of the data from a ship
         !  based platform is less than land based.  This contradicts
         !  the notion that the first-guess field over the water is
         !  not as well defined.

!!!!!!!  WHERE ( shipreport_small ) factor = 0.75

!BPR BEGIN
      CASE ('PMSLPSFC') error_factor_by_field
        !  Loop through each of the station locations (numobs).  
        station_loop_0a : DO num = 1, numobs
         !Increase the maximum allowed error in sea level pressure derived from
         !surface pressure as station elevation increases.  Since assumptions must
         !be made in this calculation above regarding the atmosphere below ground,
         !the thicker the layer through which these assumptions are made the
         !larger the potential error caused by the assumptions
         factor(num)=1.0+(ABS(station_elevation(num))/2000.0)
        END DO station_loop_0a

      CASE ( 'DEWPOINT' ) error_factor_by_field
        !  Loop through each of the station locations (numobs).  
        !  Make the tolerance larger for lower dewpoints since the same dewpoint
        !  difference represents less moisture difference at low dewpoints
        station_loop_0b : DO num = 1, numobs
         IF((obs_small(num).lt.255.0).and.(obs_small(num).ge.235.0)) THEN
          factor(num)=1.25
         ELSEIF((obs_small(num).lt.235.0).and.(obs_small(num).ge.215.0)) THEN
          factor(num)=1.50
         ELSEIF(obs_small(num).lt.215.0) THEN
          factor(num)=1.75
         ELSE
          factor(num)=1.00
         END IF
         
        END DO station_loop_0b

!BPR END
      CASE ( 'RH      ' , 'PMSL    ' ) error_factor_by_field

   END SELECT error_factor_by_field

   !  Loop through each of the station locations (numobs).  

   station_loop_1 : DO num = 1, numobs

      !  Calculate values of first guess at the observation location:
      !  note that i is in west-east (or x) direction, and j in south-north 
      !  (or y) direction.  This is a weighted mean in both the x and y 
      !  directions.
      !
      !
      !                    (iob,job+1)   (iob+1,job+1)
      !               -        -------------
      !        1-dyob |        |           |
      !               |        |           |
      !               -        |  *        |
      !          dyob |        | (xob,yob) |
      !               |        |           |
      !               -        -------------
      !                    (iob,job)     (iob+1,job)
      !
      !                        |--|--------|
      !                       dxob  1-dxob

      iob       = xob (num) 
      job       = yob (num)
      dxob      = xob ( num ) - iob
      dyob      = yob ( num ) - job
      aob(num)  = ( 1.-dxob ) * ( ( 1.-dyob ) * gridded ( iob   , job   )   &
                  +                    dyob   * gridded ( iob   , job+1 ) ) &
                  +    dxob *   ( ( 1.-dyob ) * gridded ( iob+1 , job   )   &
                  +                    dyob   * gridded ( iob+1 , job+1 ) )

   END DO station_loop_1

   !  The difference is based linearly upon the variable (RH is 
   !  treated as sqrt(1-rh/100) in this computation).
      
! foo
!  IF      ( name .NE. 'RH      ') THEN
      err = obs_small  - aob
!  ELSE IF ( name .EQ. 'RH      ' ) THEN
!     err = 100. * ( aob * aob - obs_small * obs_small )
!  END IF

   !  Define the error max for a particular variable and pressure level.

   maximum_error : SELECT CASE ( name )

      CASE ( 'UU      ' , 'VV      ' )  

         !  For U and V, factored difference is the larger of 
         !  1) the factored error maximum, or 
         !  2) 15% of the analysis value.

         factored_difference  = max ( factor * max_difference , 0.15 * aob )

!BPR BEGIN
!     CASE ( 'TT      ' , 'RH      ' , 'PMSL    ' ) 
      CASE ( 'TT      ' , 'RH      ' , 'PMSL    ', 'PMSLPSFC', 'DEWPOINT' ) 
!BPR END

         !  The maximum difference is scaled by the factor for the 
         !  temperature, sea level pressure and the relative
         !  humidity.

         factored_difference  = factor * max_difference

   END SELECT maximum_error

   !  If an observation has been flagged as failing the maximum
   !  difference criteria (error [err] larger than the allowable 
   !  difference [factored_difference]),then qc_flag for this variable, 
   !  level, time has to be modified.  

   CALL modify_qc_flag ( qc_flag_small , numobs , err , &
   factored_difference , fails_error_max )

   !  This notifies the user which observations have been flagged as
   !  suspect by the error_max routine.

   IF ( print_error_max ) THEN
      station_loop_2 : DO num = 1 , numobs
         !BPR BEGIN
         !IF ( name(1:8) .EQ. 'PMSL    ' ) THEN
         IF ( ( name(1:8) .EQ. 'PMSL    ' ) .OR. ( name(1:8) .EQ. 'PMSLPSFC' ) )  THEN
         !BPR END 
            scale = 100.
         ELSE
            scale = 1.
         END IF
         IF ( ABS ( err(num) ) .GT. factored_difference(num) ) THEN
            WRITE ( message , FMT = '(&
                  &  " ID=",a8,&    
                  &  ",NAME="  ,a8, &
                  &  ",ERRMX="     ,f5.1,&
                  &  ",LOC=(" ,f6.1,",",f6.1,")",&
                  !BPR BEGIN
                  &  ",PLEVEL=" ,f6.1,&
                  !BPR END
                  &  ",OBS="       ,f6.0,&
                  &  ",1ST GUESS=" ,f6.0,&
                  &  ",DIFF="      ,f5.0)' ) &
                 station_id(num),name,factored_difference(num)/scale, &
                 !BPR BEGIN
                 !Also add pressure level being examined (plevel) -- note that this is not
                 !necessarily the pressure of the ob, but it is the pressure level for which we
                 !were doing QC when we found this ob so it must be at least close to this
                 !pressure
                 !xob(num),yob(num),obs(num)/scale,aob(num)/scale,err(num)/scale
                 xob(num),yob(num),plevel,obs(num)/scale,aob(num)/scale,err(num)/scale
                 !BPR END
            error_number = 00361001
            error_message(1:31) = 'error_max                      '
            error_message(32:)  = message
            fatal = .false.
            listing = .false.
            CALL error_handler ( error_number , error_message ,  &
            fatal , listing )
         END IF
      END DO station_loop_2
   END IF

   !  Put the qc_flag back into the permanent space.

   qc_flag(:numobs) = qc_flag_small(:numobs)     

END SUBROUTINE error_max

END MODULE qc1
