MODULE qc2

   INCLUDE 'missing.inc'

CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE buddy_check ( station_id, obs, numobs, xob, yob, qc_flag,   &
gridded, iew, jns, name,                      &
pressure, local_time, dx , buddy_weight , buddy_difference , print_buddy )

   USE qc0

!   Imported from BUDDY in RAWINS

!   Notes:
!      This routine performs error checks by comparing the difference
!         between the first guess field
!         and observations at the observing locations (xob, yob)
!         with the averaged difference at the nearby stations within
!         a given distance. If the difference value at the observation
!         point differs from the average of the nearby difference values
!         by more than a certain maximum (dependent on variable types,
!         pressure levels, and BUDWGT parameter in namelist), the
!         observation is flagged as suspect, and add the flag message
!         to the observation records.

!   *** Note that x and y are in natural coordinates now.
!         Need to check if the x and y are used correctly...

!   other variables:
!    err:       difference between observations and first guess fields

   IMPLICIT NONE
   INTEGER                                        :: numobs
   CHARACTER ( LEN =   8 ) , DIMENSION ( : )      :: station_id
   REAL                    , DIMENSION ( : )      :: obs, xob, yob
   INTEGER                 , DIMENSION ( : )      :: qc_flag
   REAL                    , DIMENSION( : , : )   :: gridded
   INTEGER                                        :: iew, jns
   CHARACTER ( LEN =   8 )                        :: name
   REAL                                           :: buddy_weight,   &
                                                     buddy_difference , &
                                                     dx          , &
                                                     pressure
   INTEGER                 , DIMENSION ( : )      :: local_time
   LOGICAL                                        :: print_buddy

   !  err:         difference between observations and first guess fields

   INTEGER                                    :: num, num3, numj,    &
                                                 iob, job,           &
                                                 kob
   REAL                                       :: dxob, dyob,         &
                                                 range,              &
                                                 distance,           &
                                                 sum,                &
                                                 average,            &
                                                 x, y,               &
                                                 diff_numj,          &
                                                 diff_check_1,       &
                                                 diff_check_2
   REAL                                          plevel
   REAL  , DIMENSION ( numobs )               :: err,                &
                                                 diff,               &
                                                 difmax,             &
                                                 diff_at_obs,        &
                                                 obs2,               &
                                                 aob
   INTEGER , DIMENSION ( numobs )             :: buddy_num,          &
                                                 buddy_num_final
   CHARACTER ( LEN = 100 )                    :: message
   REAL                                       :: scale

   CHARACTER ( LEN =   8 ) , DIMENSION ( numobs ) :: station_id_small
   REAL                    , DIMENSION ( numobs ) :: obs_small , xob_small , yob_small
   INTEGER                 , DIMENSION ( numobs ) :: qc_flag_small
   INTEGER                 , DIMENSION ( numobs ) :: local_time_small

   INCLUDE 'error.inc'

   INTERFACE
      INCLUDE 'error.int'
   END INTERFACE

   station_id_small(:numobs)  = station_id(:numobs)
   obs_small(:numobs)         = obs(:numobs)
   xob_small(:numobs)         = xob(:numobs)
   yob_small(:numobs)         = yob(:numobs)
   qc_flag_small(:numobs)     = qc_flag(:numobs)
   local_time_small(:numobs)  = local_time(:numobs)

   plevel = pressure

   !  Define distance within which buddy check is to find neighboring stations:

   IF ( plevel .LE. 201.0 ) THEN
      range = 650.
   ELSE IF ( plevel .LE. 1000.1 ) THEN
      range = 300.
   ELSE
      range = 500.
   END IF

   !  Define tolerance value for comparison:

   error_allowance_by_field : SELECT CASE ( name )

      CASE ( 'UU      ' , 'VV      ' )
         IF      ( plevel .LT. 401. ) THEN
            difmax = buddy_difference * 1.7
         ELSE IF ( plevel .LT. 701. ) THEN
            difmax = buddy_difference * 1.4
         ELSE IF ( plevel .LT. 851. ) THEN
            difmax = buddy_difference * 1.1
         ELSE
            difmax = buddy_difference
         END IF

      CASE ( 'PMSL    ' )
         difmax = buddy_difference

      CASE ( 'TT      ' )
         IF      ( plevel .LT. 401. ) THEN
            difmax = buddy_difference * 0.6
         ELSE IF ( plevel .LT. 701. ) THEN
            difmax = buddy_difference * 0.8
         ELSE
            difmax = buddy_difference
         END IF

      CASE ( 'RH      ' )
         difmax = buddy_difference

   END SELECT error_allowance_by_field

   difmax = difmax * buddy_weight

   !  Loop through station locations (numobs) to compute first-guess value
   !  at each observation locations.

   station_loop_1 : DO num = 1, numobs

      !  Calculate values of first guess at the observation location:
      !  note that i is in west-east (or x) direction, and j in south-north 
      !  (or y) direction.

      iob  = xob(num)
      job  = yob(num)
      dxob = xob(num) - iob
      dyob = yob(num) - job
      aob(num)  =  ( 1.-dxob ) * ( ( 1.-dyob ) * gridded(iob,  job)     &
             +                          dyob   * gridded(iob,  job+1) ) &
             +          dxob   * ( ( 1.-dyob ) * gridded(iob+1,job)     &
             +                          dyob   * gridded(iob+1,job+1) )

   END DO station_loop_1

! foo
!  IF ( name .EQ. 'RH      ' ) THEN
!     obs2 = 100. * ( 1. - obs_small * obs_small )
!     aob  = 100. * ( 1. - aob * aob )
!  ELSE
      obs2 = obs_small
!  END IF

   station_loop_2 : DO num = 1, numobs

      !  Compute distance between all pairs of observations, 
      !  find number of observations within a given distance: buddy_num,
      !  and compute the difference between first guess and obs.

      buddy_num(num) = 0
      station_loop_3 : DO num3 = 1, numobs

         IF ( num3 .NE. num ) THEN
            x = xob(num) - xob(num3)
            y = yob(num) - yob(num3)
            distance = SQRT ( x*x + y*y ) * dx
            IF ( distance .LE. range .AND. distance .GE. 0.1 ) THEN
               buddy_num(num) = buddy_num(num) + 1
               diff(buddy_num(num)) = obs2(num3) - aob(num3)
            END IF
         END IF

      END DO station_loop_3

      diff_at_obs(num) = obs(num) - aob(num)
      sum = 0.
      DO numj = 1, buddy_num(num)
         sum = sum + diff(numj)
      END DO

      !  Check to see if there are any buddies.

      IF ( buddy_num(num) .GT. 0 ) average = sum / buddy_num(num)

      !  Check if there is any bad obs among the obs located within the 
      !  the radius surrounding the test ob.

      diff_check_1 = difmax (1) * 1.25
      diff_check_2 = difmax (1)

      check_bad_ob : IF ( buddy_num(num) .GE. 2 ) THEN

         kob = buddy_num(num)
         remove_bad_ob : DO numj = 1, buddy_num(num)
            diff_numj = ABS ( diff ( numj ) - average )
            IF ( diff ( numj ) .GT. diff_check_1  &
                .AND. diff_numj .GT. diff_check_2 ) THEN
               kob = kob - 1
               sum = sum - diff ( numj )
            END IF
         END DO remove_bad_ob
         buddy_num_final(num) = kob

         !  We may have removed too many observations.

         IF ( kob .GE. 2 ) THEN
            average = sum / kob
            err(num) = diff_at_obs(num) - average
         ELSE
            err(num)     = 0.
            qc_flag_small(num) = qc_flag_small(num) + no_buddies
         END IF

      ELSE check_bad_ob

         err(num)     = 0.
         qc_flag_small(num) = qc_flag_small(num) + no_buddies

      END IF check_bad_ob

      IF ( buddy_num_final(num) .EQ. 2 ) difmax(num) = 1.2 * difmax(num)

   END DO station_loop_2

   !  If an observation has been flagged as failing the buddy
   !  check criteria (error [err] larger than the allowable 
   !  difference [difmax]), then qc_flag for this variable, 
   !  level, time has to be modified.  

   CALL modify_qc_flag ( qc_flag_small , numobs , err, &
   difmax , fails_buddy_check )

   !  This notifies the user which observations have been flagged as
   !  suspect by the buddy check routine.

   IF ( print_buddy ) THEN
      station_loop_4 : DO num = 1 , numobs
         IF ( name(1:8) .EQ. 'PMSL    ' ) THEN
            scale = 100.
         ELSE
            scale = 1.
         END IF
         IF ( ABS ( err(num) ) .GT. difmax(num) ) THEN
            WRITE ( message , FMT = '(&
                  &  " ID=",a8,&    
                  &  ",NAME="  ,a8, &
                  &  ",BUDDY="     ,f5.1,&
                  &  ",LOC=(" ,f6.1,",",f6.1,")",&
                  &  ",OBS="       ,f6.0,&
                  &  ",1ST GUESS=" ,f6.0,&
                  &  ",DIFF="      ,f5.0)' ) &
                 station_id(num),name,difmax(num)/scale, &
                 xob(num),yob(num),obs(num)/scale,aob(num)/scale,err(num)
            error_number = 00364001
            error_message(1:31) = 'buddy_check                    '
            error_message(32:)  = message
            fatal = .false.
            listing = .false.
            CALL error_handler ( error_number , error_message ,  &
            fatal , listing )
         END IF
      END DO station_loop_4
   END IF

   !  Put the qc_flag back into the permanent space.

   qc_flag(:numobs) = qc_flag_small(:numobs)     

END SUBROUTINE buddy_check

!-------------------------------------------------------------------------------

END MODULE qc2
