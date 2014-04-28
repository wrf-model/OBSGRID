MODULE qc0

   INCLUDE 'missing.inc'

CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE modify_qc_flag ( qc_flag , numobs , err ,  &
factored_difference , test_number ) 

   IMPLICIT NONE

   INTEGER                                        :: numobs , test_number
   INTEGER                 , DIMENSION ( numobs ) :: qc_flag
   REAL                    , DIMENSION ( numobs ) :: err,                &
                                                     factored_difference 

   INTEGER                 , DIMENSION ( numobs ) :: qc_small , & 
                                                     qc_large

   !  If the difference between observation and analysis at obs location
   !  (err) exceeds errmx, mark those obs as flagged by the error checking
   !  routines.

   WHERE ( ABS ( err ) .GT. factored_difference ) 
      
      qc_small = mod ( qc_flag , test_number ) 
      qc_large = ( qc_flag / ( test_number * 2 ) ) * 2
      qc_flag = qc_large*test_number + qc_small + test_number
      
   END WHERE
      
END SUBROUTINE modify_qc_flag

!------------------------------------------------------------------------------

SUBROUTINE ob_density ( xob , yob , grid_dist , numobs , tobbox , iew , jns )

!  Compute a guess at the observation density for the surface FDDA.

   IMPLICIT NONE

   REAL    , DIMENSION ( : )      :: xob, yob
   REAL                           :: grid_dist
   INTEGER                        :: numobs
   REAL    , DIMENSION ( : , : )  :: tobbox
   INTEGER                        :: iew, jns

   INTEGER                        :: num, iob, job, i, j , &
                                     iobs , iobe , jobs , jobe
   REAL                           :: dist
                            
   !  Loop through each of the station locations (numobs).  

   DO num = 1, numobs

      iob       = xob (num) 
      job       = yob (num)
      jobs = MAX ( job - NINT( 250. / grid_dist ) - 1 ,   1   )
      jobe = MIN ( job + NINT( 250. / grid_dist ) + 1 , jns-1 )
      iobs = MAX ( iob - NINT( 250. / grid_dist ) - 1 ,   1   )
      iobe = MIN ( iob + NINT( 250. / grid_dist ) + 1 , iew-1 )
      DO j = jobs , jobe
         DO i = iobs , iobe
            dist = SQRT ( ( REAL(i) - xob(num) ) **2 + ( REAL(j) - yob(num) ) **2 ) 
            IF ( dist .LT. 250. / grid_dist ) THEN
               tobbox(j,i) = tobbox(j,i) + 1
            END IF
         END DO
      END DO

   END DO

END SUBROUTINE ob_density
                 
!------------------------------------------------------------------------------

  SUBROUTINE obs_distance( tobbox, distance, iew, jns, dx )
!
! This routine computes the distance (array) to the closest obs for each 
! grid point based on the obs density (array) that is used in
! MM5 surface analysis nudging in determining the horizontal weighting
! function (or confidence in the analysis) based on observation density 
! at each grid point.  In WRF surface analysis nudging codes, the confidence 
! is determined by the distance to the nearest obs.  The data density is 
! no longer needed.
!
!                      Aijun Deng     06/11/2008 at Penn State
!
    IMPLICIT NONE

    INTEGER                               :: jns, iew     ! Array x and y dimensions
    REAL, INTENT(IN)                      :: dx           ! Grid distance in meters
    REAL, DIMENSION( : , : ), INTENT(IN)  :: tobbox       ! Incoming density array
    REAL, DIMENSION( : , : ), INTENT(OUT) :: distance     ! Outgoing distance array

    INTEGER :: i, j, ii, jj
    REAL    :: dijmin, di, dj, dij

    distance(:,:) = sqrt( (iew-1)*dx*(iew-1)*dx + (jns-1)*dx*(jns-1)*dx )   ! Initialize 
                    ! the distance array to be the domain diagonal size (something big).
    dijmin = 0.0

    DO i = 1, iew
    DO j = 1, jns
        IF( tobbox(j,i) > 0.5 ) THEN
          distance(j,i) = 0.0
        ELSE
          dijmin = MAXVAL( distance/dx )
          DO ii = 1, iew
          DO jj = 1, jns
            IF( tobbox(jj,ii) > 0.5 ) THEN
              di     = ABS( FLOAT(ii-i) )
              dj     = ABS( FLOAT(jj-j) )
              dij    = SQRT( di*di+dj*dj )
              dijmin = MIN( dijmin,dij )
            ENDIF
          ENDDO
          ENDDO

          distance(j,i) = dijmin * dx

        ENDIF
    ENDDO
    ENDDO

  END SUBROUTINE obs_distance

!------------------------------------------------------------------------------

SUBROUTINE local ( time , lonob , local_time , num_obs )

!  The "local" time is an approximation for the maximum difference 
!  QC check.  This allows warmer times of the day and cooler night times
!  to relax the maximum difference.

   IMPLICIT NONE

   INTEGER                                 :: time , num_obs
   INTEGER                 , DIMENSION (:) :: local_time
   REAL                    , DIMENSION (:) :: lonob

   INTEGER                                 :: num

   INCLUDE 'error.inc'

   INTERFACE
      INCLUDE 'error.int'
   END INTERFACE

   !  The UTC time in hours (time is hhmm), plus a correction for the
   !  longitude (+ is east of prime meridian, - is west).  We do not
   !  care what "day" it is, only the approximate time based on the 
   !  sun's elevation.

   local_time(:num_obs) = time/100 + lonob(:num_obs)/15.

   !  Make the local time between 0 and 24 h.

   WHERE ( local_time(:num_obs) .LT. 0 ) local_time(:num_obs) = 24 + local_time(:num_obs)

   local_time(:num_obs) = mod ( local_time(:num_obs) , 24 ) * 100

   !  Reset the unused values of longitude and the local time.

   lonob(num_obs+1:)      = missing_r
   local_time(num_obs+1:) = missing_r

   !  If there are any values outside the acepted limits, we need to know 
   !  about it.

   DO num = 1 , num_obs

      wrong_time : IF ( local_time(num) .GE. 2400 ) THEN
         WRITE ( UNIT = * , FMT = * ) 'Index = ',num, &
         '  local_time = ',local_time(num), & 
         '  UTC time = ',time, & 
         '  longitude = ',lonob(num)
         error_number = 00363001
         error_message(1:31) = 'local                          '
         error_message(32:)  = ' The value of the computed local_time &
         &is larger than permissible.'
         fatal = .true.
         listing = .false.
         CALL error_handler ( error_number , error_message ,  &
         fatal , listing )
      ELSE IF ( local_time(num) .LT.    0 ) THEN
         WRITE ( UNIT = * , FMT = * ) 'Index = ',num, &
         '  local_time = ',local_time(num), & 
         '  UTC time = ',time, & 
         '  longitude = ',lonob(num)
         error_number = 00363002
         error_message(1:31) = 'local                          '
         error_message(32:)  = ' The value of the computed local_time &
         &is negative.'
         fatal = .true.
         listing = .false.
         CALL error_handler ( error_number , error_message ,  &
         fatal , listing )
      END IF wrong_time

   END DO 

END SUBROUTINE local

!-------------------------------------------------------------------------------

SUBROUTINE qc_consistency ( obs , num_obs )

!  This routine forces related variables to have the same QC flags when one of
!  the variables will be removed from usage in the objective analysis phase.

   USE observation

   IMPLICIT NONE

   TYPE ( report ) , INTENT ( INOUT ) , DIMENSION ( : ) :: obs
   INTEGER , INTENT(IN) :: num_obs

   !  Local data.

   TYPE ( measurement ) , POINTER :: current , next

   INTEGER :: obs_loop

   INTEGER :: qc_u , qc_v
   INTEGER :: qc_t , qc_rh

   !  Loop over all of the observations.

   obs_loop = 1
   loop_all_obs : DO obs_loop = 1 , num_obs

      !  We are now looking at a single observation.  Let's look at all of the levels.
      !  The first level is named the surface in this data structure.

      current => obs(obs_loop)%surface
      loop_all_levels : DO WHILE ( ASSOCIATED ( current ) ) 


         !  If u is bad, that means v has to be bad.  If v is bad, that means
         !  that u has to be bad.  The "badness" is measured as a value that has a 
         !  QC flag that implies that it won't be used in the objective analsyis.

         qc_u = current%meas%u%qc
         qc_v = current%meas%v%qc

         !IF      ( ( qc_u .GE. fails_error_max   ) .AND. ( qc_v .LT. fails_error_max   ) ) THEN 
         !   current%meas%v%qc         = current%meas%v%qc         + fails_error_max
         !   current%meas%speed%qc     = current%meas%speed%qc     + fails_error_max
         !   current%meas%direction%qc = current%meas%direction%qc + fails_error_max
         !ELSE IF ( ( qc_u .GE. fails_buddy_check ) .AND. ( qc_v .LT. fails_buddy_check ) ) THEN 
         !   current%meas%v%qc         = current%meas%v%qc         + fails_buddy_check
         !   current%meas%speed%qc     = current%meas%speed%qc     + fails_buddy_check
         !   current%meas%direction%qc = current%meas%direction%qc + fails_buddy_check
         !ELSE IF ( ( qc_v .GE. fails_error_max   ) .AND. ( qc_u .LT. fails_error_max   ) ) THEN 
         !   current%meas%u%qc         = current%meas%u%qc         + fails_error_max
         !   current%meas%speed%qc     = current%meas%speed%qc     + fails_error_max
         !   current%meas%direction%qc = current%meas%direction%qc + fails_error_max
         !ELSE IF ( ( qc_v .GE. fails_buddy_check ) .AND. ( qc_u .LT. fails_buddy_check ) ) THEN 
         !   current%meas%u%qc         = current%meas%u%qc         + fails_buddy_check
         !   current%meas%speed%qc     = current%meas%speed%qc     + fails_buddy_check
         !   current%meas%direction%qc = current%meas%direction%qc + fails_buddy_check
         !END IF

         IF ( qc_u .GT. qc_v ) current%meas%v%qc = qc_u
         IF ( qc_v .GT. qc_u ) current%meas%u%qc = qc_v
         current%meas%speed%qc     = current%meas%u%qc
         current%meas%direction%qc = current%meas%u%qc


         !  If t is bad, that means td has to be bad.  If td is bad, that DOES NOT
         !  mean that t has to be bad.  The "badness" is measured as a value that has a 
         !  QC flag that implies that it won't be used in the objective analsyis.

         qc_t = current%meas%temperature%qc
         qc_rh = current%meas%rh%qc

         IF      ( ( qc_t .GE. fails_error_max   ) .AND. ( qc_rh .LT. fails_error_max   ) ) THEN 
            current%meas%rh%qc        = current%meas%rh%qc        + fails_error_max
            current%meas%dew_point%qc = current%meas%dew_point%qc + fails_error_max
         ELSE IF ( ( qc_t .GE. fails_buddy_check ) .AND. ( qc_rh .LT. fails_buddy_check ) ) THEN 
            current%meas%rh%qc        = current%meas%rh%qc        + fails_buddy_check
            current%meas%dew_point%qc = current%meas%dew_point%qc + fails_buddy_check
         END IF

         !  If td > t, then td and RH have to be bad.

         IF ( ( current%meas%dew_point%data .GT. current%meas%temperature%data   ) .AND. &
              ( qc_rh                       .LT. fails_buddy_check ) ) THEN 
            current%meas%rh%qc        = current%meas%rh%qc        + fails_buddy_check
            current%meas%dew_point%qc = current%meas%dew_point%qc + fails_buddy_check
         END IF

         !  Point to the next level

         current => current%next

      END DO loop_all_levels

   END DO loop_all_obs

END SUBROUTINE qc_consistency

!-------------------------------------------------------------------------------

END MODULE qc0
