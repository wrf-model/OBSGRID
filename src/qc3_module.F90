MODULE qc3

   INCLUDE 'missing.inc'

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE dry_convective_adjustment ( obs , counter , print_dry ) 

!  This subroutine performs the dry convective adjustment for a sounding 
!  and prints out the points and levels where the adjustment have been
!  invoked. The adjustment removes any super-adiabatic lapse rate in the
!  sounding by making use of the conservation of dry static energy.

   USE observation
   USE qc0, only: add_to_qc_flag

   IMPLICIT NONE

   INTEGER        :: counter
   INTEGER        :: k, n1, kk, numb, num1, i, klop,  &  ! parameters for do loops
                     kst,     & !  kst is the starting level of adjustment
                     ked        !  ked is the ending level of adjustment

   REAL              :: esum   , & !  sum of static energy between two adjacent levels
                        eav    , & !  mean static energy between two adjacent levels
                        psfc       !  surface pressure
   REAL , DIMENSION ( : ) , ALLOCATABLE  :: p , &  !  pressure 
                                            t , &  !  temperature
                                            td, &  !  dew point
                                            dp, &  !  dew point depression
                                            z , &  !  height
                                            zt, &  !  height calculated from temperature
                                            e      !  dry static energy
   REAL , DIMENSION ( : ) , ALLOCATABLE  :: theta       !  potential temperature
   INTEGER                               :: index       !  number of vertical levels
   TYPE ( measurement ) , POINTER        :: current
   TYPE ( report )                       :: obs
   CHARACTER ( LEN = 100 )               :: message
   LOGICAL                               :: print_dry

   INCLUDE 'constants.inc'
   INCLUDE 'error.inc'

   INTERFACE
      INCLUDE 'error.int'
   END INTERFACE

   ! Try to find out the surface level first

   psfc = missing_r
   NULLIFY ( current )
   current => obs%surface
   DO WHILE ( ASSOCIATED ( current ) )
      IF ( eps_equal ( current%meas%height%data, obs%info%elevation, 0.1 ) ) THEN
         psfc = current%meas%pressure%data
      ELSE
         ! Do nothing
      END IF
      current => current%next
   END DO

   !  Find out the number of vertical levels, then ALLOCATE arrays.

   NULLIFY ( current )
   index = 1
   current => obs%surface
   DO WHILE ( ASSOCIATED ( current ) )
      IF ( ( .NOT. eps_equal ( current%meas%pressure%data,    missing_r, 1. ) ) .AND. & 
           ( .NOT. eps_equal ( psfc                      ,    missing_r, 1. ) ) .AND. & 
           ( current%meas%pressure%data .LE. psfc ) .AND. & 
           ( .NOT. eps_equal ( current%meas%temperature%data, missing_r, 1. ) ) ) THEN
         index = index + 1
      ELSE
         ! Do nothing, missing or bad data
      END IF
      current => current%next
   END DO
   index = index - 1

   IF ( ( index .LE. 3 ) .OR. ( eps_equal ( psfc , missing_r, 1. ) ) ) THEN
      !   The vertical levels are less than 3, do nothing
   ELSE ! Make dry adjustment for temperature

      ALLOCATE ( P ( index ) )
      ALLOCATE ( t ( index ) )
      ALLOCATE ( td( index ) )
      ALLOCATE ( dp( index ) )
      ALLOCATE ( z ( index ) )
      ALLOCATE ( zt( index ) )
      ALLOCATE ( e ( index ) )
      ALLOCATE ( theta ( index ) )
      NULLIFY  ( current )

      !  Read data set into array and perform quality check.

      index = 1
      current => obs%surface
      DO WHILE ( ASSOCIATED ( current ) )
         IF ( ( .NOT. eps_equal ( current%meas%pressure%data,    missing_r, 1. ) ) .AND. & 
              ( current%meas%pressure%data .LE. psfc ) .AND. & 
              ( .NOT. eps_equal ( current%meas%temperature%data, missing_r, 1. ) ) ) THEN
            p ( index ) = current%meas%pressure%data
            t ( index ) = current%meas%temperature%data
            td( index ) = current%meas%dew_point%data
            IF ( .NOT. eps_equal ( current%meas%dew_point%data, missing_r, 1. ) ) THEN
               dp ( index ) = t ( index ) - td ( index )
            ELSE
               dp ( index ) = missing_r
            ENDIF
            z ( index ) = current%meas%height%data
            index = index + 1
         ELSE
            ! Do nothing, missing or bad data
         END IF
         current => current%next
      END DO
      index = index - 1
      n1    = index - 1

      IF ( .NOT. eps_equal ( z ( 1 ) , missing_r, 1. ) ) THEN

         !  We have the first level height, go ahead make dry adjustment.

         zt ( 1 ) = z ( 1 )

         !  Compute potential temperature.

         DO k = 1, index
            theta ( k ) = t ( k ) * ( 100000.0 / p ( k ) ) ** rcp
         END DO

         loop_adjust: DO k = 1 , index - 2
      
            IF ( theta ( k ) .LE. theta ( k + 1 ) )  THEN

               !  Stable layer, do nothing.

               CYCLE loop_adjust 

            ELSE

               !  Compute the height z at all levels by taking a mean temperature of the layer.
              
               DO kk = 2 , n1
                  zt ( kk ) = zt ( kk - 1 ) + 0.5 * ( t ( kk - 1 ) + t ( kk ) ) *   &
                             ( gasr / g ) * LOG ( p ( kk -1 ) / p ( kk ) )
               END DO

               !  Compute the dry static energy.

               DO kk = 1 , n1
                  e ( kk ) = g * zt ( kk ) + cp * t ( kk )
               END DO

               esum = e ( k ) + e ( k + 1 )
               eav  = esum * 0.5

               !  Downward checking the layer in which superadiabat exists.

               numb = 2

               loop_down: DO klop = k - 1 , 1 , -1
                  IF ( klop .LT. 1 ) THEN
                     EXIT loop_down
                  ELSE
                     IF ( eav  .GE. e ( klop ) ) THEN
                        EXIT loop_down
                     ELSE
                        esum = esum + e ( klop )
                        numb = numb + 1
                        eav = esum / REAL ( numb )
                     END IF
                  END IF
               END DO loop_down
               num1 = numb
                      
               !  Upward checking the layer in which superadiabat exists.

               loop_up: DO klop = k + 2 , index , 1
                  IF ( klop .GT. n1 ) THEN
                     EXIT loop_up
                  ELSE
                     IF ( eav .LE. e ( klop ) ) THEN
                        EXIT loop_up
                     ELSE
                        esum = esum + e ( klop )
                        numb = numb + 1
                        eav = esum / REAL ( numb )
                     END IF
                  END IF
               END DO  loop_up

               !  Define the starting and ending levels.

               kst = k + 2 - num1
               ked = kst + numb -1

               WRITE ( message , FMT = '(&
               &  "DRY COVECTIVE ADJUSTMENT FROM LEVELS", I5, "  TO  ", I5, " FROM P &
               &  = ",F10.2," TO ",F10.2," Pa." )' ) kst, ked, p(kst), p(ked)
               error_number = 00364001
               error_message(1:31) = 'dry_convective_adjustment      '
               error_message(32:)  = message
               fatal = .false.
               listing = .false.
               IF ( print_dry ) CALL error_handler ( error_number , error_message ,  &
               fatal , listing )

               !  Evaluates the new temperature profile according to the mean static 
               !  energy of the layer.

               DO kk = kst , ked
                  t ( kk ) = ( ( eav - g * zt ( kst ) ) / cp ) *             &
                             ( p ( kk ) / p ( kst ) ) ** rcp
               END DO

            END IF

         END DO loop_adjust

      ELSE
         ! We do not have the first level height, do nothing.
      END IF   

      current => obs%surface
      looptt: DO WHILE ( ASSOCIATED ( current ) )
         DO i = 1, index
            IF ( ( eps_equal ( p ( i ), current%meas%pressure%data, 0.1) ) .AND. &
                 ( current%meas%temperature%data .NE. t ( i ) ) ) THEN
               WRITE ( message , FMT = '(&
               &  " T dry conv adj: ID=",a8,&
               &  ", NAME="  ,a8, " T old=" ,f5.1," (K),  T new=" ,f5.1," (K), P=",f6.1," (hPa)")' )&
               obs%location%id , obs%location%name , current%meas%temperature%data , t(i) , p(i)/100.
               error_number = 00364001
               error_message(1:31) = 'dry_convective_adjustment      '
               error_message(32:)  = message
               fatal = .false.
               listing = .false.
               IF ( print_dry ) CALL error_handler ( error_number , error_message ,  &
               fatal , listing )
               current%meas%temperature%data = t ( i )
               IF ( .NOT. eps_equal ( current%meas%dew_point%data, missing_r, 1. ) ) THEN
                  current%meas%dew_point%data = t ( i ) - dp ( i)
                  current%meas%rh%data = 100. * exp ( 5418.12 * ( 1./t(i) - 1./(t(i)-dp(i))))
               ENDIF
               IF ( current%meas%temperature%qc .NE. missing ) THEN
                  !BPR BEGIN
                  !Previous code would add convective_adjustment flag even it
                  !was already extant in the QC flag
                  !current%meas%temperature%qc   = current%meas%temperature%qc + convective_adjustment
                  CALL add_to_qc_flag( current%meas%temperature%qc, convective_adjustment )
                  !BPR END
               END IF
               IF ( current%meas%dew_point%qc .NE. missing ) THEN
                  !BPR BEGIN
                  !Previous code would add convective_adjustment flag even it
                  !was already extant in the QC flag
                  !current%meas%dew_point%qc   = current%meas%dew_point%qc + convective_adjustment
                  CALL add_to_qc_flag( current%meas%dew_point%qc, convective_adjustment )
                  !BPR END
               END IF
               IF ( current%meas%rh%qc .NE. missing ) THEN
                  !BPR BEGIN
                  !Previous code would add convective_adjustment flag even it
                  !was already extant in the QC flag
                  !current%meas%rh%qc   = current%meas%rh%qc + convective_adjustment
                  CALL add_to_qc_flag( current%meas%rh%qc, convective_adjustment )
                  !BPR END
               END IF
            ELSE
               ! Do nothing
            END IF
         END DO
         current => current%next
      END DO looptt

      DEALLOCATE ( P     )
      DEALLOCATE ( t     )
      DEALLOCATE ( td    )
      DEALLOCATE ( dp    )
      DEALLOCATE ( z     )
      DEALLOCATE ( zt    )
      DEALLOCATE ( e     )
      DEALLOCATE ( theta )
      NULLIFY    ( current )

   END IF

END SUBROUTINE dry_convective_adjustment 

!-------------------------------------------------------------------------------

SUBROUTINE vert_cons_check ( obs , counter , print_vert ) 

!  This subroutine is called to perform vertical 
!  consistency check for a single sounding.     

   USE observation

   IMPLICIT NONE

   INTEGER                                :: i, index, iflag , counter
   REAL   , DIMENSION ( : ), ALLOCATABLE  :: p, t, td, ws, wd, u, v, rh, work1, work2
   INTEGER, DIMENSION ( : ), ALLOCATABLE  :: tqc, tdqc, wsqc, wdqc, uqc, vqc, rhqc, iindex
   TYPE ( measurement ), POINTER          :: current
   TYPE ( report )                        :: obs
   REAL                                   :: d1, d2, che1, che2, th, th1, th2, psfc, depression
   CHARACTER ( LEN = 100 )                :: message
   LOGICAL                                :: print_vert
   REAL                                   :: p1 , p2 , h1 , h2
   LOGICAL                                :: found

   INCLUDE 'error.inc'
   INCLUDE 'constants.inc'

   INTERFACE
      INCLUDE 'error.int'
   END INTERFACE

   !  Check for spikes in temperature curve. 
   !  Both the lapse rate check and potential temperature check will be performed

   !  Is there a duplicate surface value?  We can check the pressure and height values.
   !  If we find two (or more) levels with the height equal to the surface elevation,
   !  then we set all of the information in the duplicate level to missing.

   NULLIFY ( current )
   current => obs%surface
   IF ( ASSOCIATED ( current ) ) THEN
      p1 = current%meas%pressure%data
      h1 = current%meas%height%data
   END IF
   current => current%next
   found = .FALSE.
   DO WHILE ( ASSOCIATED ( current ) )
      IF ( eps_equal ( current%meas%height%data, obs%info%elevation, 0.1 ) ) THEN
         p2 = current%meas%pressure%data
         h2 = current%meas%height%data
         IF ( ( eps_equal ( h1 , h2 , 0.1 ) ) .AND. &
              ( .NOT. eps_equal ( p1 , p2 , 0.1 ) ) ) THEN
            WRITE ( message , FMT = '(" Duplicate surface found at ",a8,a8)' ) &
            TRIM ( obs%location%id ) , TRIM ( obs%location%name )
            error_message(1:31) = 'vert_cons_check                '
            error_message(32:)  = message
            fatal = .false.
            listing = .false.
            CALL error_handler ( error_number , error_message ,  &
            fatal , listing )
            current%meas%pressure%data     = missing_r
            current%meas%height%data       = missing_r
            current%meas%temperature%data  = missing_r
            current%meas%dew_point%data    = missing_r
            current%meas%speed%data        = missing_r
            current%meas%direction%data    = missing_r
            current%meas%rh%data           = missing_r
            current%meas%u%data            = missing_r
            current%meas%v%data            = missing_r
            current%meas%pressure%qc       = missing
            current%meas%height%qc         = missing
            current%meas%temperature%qc    = missing
            current%meas%dew_point%qc      = missing
            current%meas%speed%qc          = missing
            current%meas%direction%qc      = missing
            current%meas%rh%qc             = missing
            current%meas%u%qc              = missing
            current%meas%v%qc              = missing
         END IF
      END IF
      current => current%next
   END DO

   !  Try to find out the surface level first

   psfc = missing_r
   NULLIFY ( current )
   current => obs%surface
   DO WHILE ( ASSOCIATED ( current ) )
      IF ( eps_equal ( current%meas%height%data, obs%info%elevation, 0.1 ) ) THEN
         psfc = current%meas%pressure%data
      ELSE
         ! Do nothing
      END IF
      current => current%next
   END DO

   !  Find out the vertical level number, then ALLOCATE arrays

   NULLIFY ( current )
   index = 1
   current => obs%surface
   DO WHILE ( ASSOCIATED ( current ) )
      IF ( ( .NOT. eps_equal ( current%meas%pressure%data,    missing_r, 1. ) ) .AND. & 
           ( .NOT. eps_equal ( psfc                      ,    missing_r, 1. ) ) .AND. & 
           ( current%meas%pressure%data .LE. psfc ) .AND. & 
           ( .NOT. eps_equal ( current%meas%temperature%data, missing_r, 1. ) ) ) THEN
         index = index + 1
      ELSE
         ! Do nothing, missing or bad data
      END IF
      current => current%next
   END DO
   index = index - 1

   IF ( ( index .LE. 3 ) .OR. ( eps_equal ( psfc , missing_r, 1. ) ) ) THEN
      !   The vertical levels are less than 3, do nothing
   ELSE ! Check superadiabatic and inversion in temperature cure

      ALLOCATE (      P ( index ) )
      ALLOCATE (      t ( index ) )
      ALLOCATE (     td ( index ) )
      ALLOCATE (     rh ( index ) )
      ALLOCATE (    tqc ( index ) )
      ALLOCATE (   tdqc ( index ) )
      ALLOCATE (   rhqc ( index ) )
      ALLOCATE ( iindex ( index ) )
      NULLIFY  ( current )

      !  Read data set into array and perform quality check

      index = 1
      current => obs%surface
      DO WHILE ( ASSOCIATED ( current ) )
         IF ( ( .NOT. eps_equal ( current%meas%pressure%data,    missing_r, 1. ) ) .AND. & 
              ( current%meas%pressure%data .LE. psfc ) .AND. & 
              ( .NOT. eps_equal ( current%meas%temperature%data, missing_r, 1. ) ) ) THEN
            p    ( index ) = current%meas%pressure%data
            t    ( index ) = current%meas%temperature%data
            td   ( index ) = current%meas%dew_point%data
            rh   ( index ) = current%meas%rh%data
            tqc  ( index ) = 0
            tdqc ( index ) = 0
            rhqc ( index ) = 0
            index = index + 1
         ELSE
            ! Do nothing, missing or bad data
         END IF
         current => current%next
      END DO
      index = index - 1

      loop1: DO i = 1, index
         iindex ( i ) = 0
      END DO loop1
         
      ! First level
      che1 = LOG ( t ( 2 ) / t ( 1 ) ) / LOG ( p ( 2 ) / p ( 1 ) ) 
      che2 = LOG ( t ( 3 ) / t ( 2 ) ) / LOG ( p ( 3 ) / p ( 2 ) ) 

      IF ( ( ( che1 .LT. -1.5 ) .OR. ( che1 .GT. 0.9 ) ) .AND.               &
           ( ( che2 .GT. -1.5 ) .AND. ( che2 .LT. 0.9 ) ) ) iindex ( 1 ) = 1
!     IF ( ( ( che1 .LT. -1.5 ) .OR. ( che1 .GT. 0.9 ) ) .AND.               &
!          ( ( che2 .LT. -1.5 ) .OR. ( che2 .GT. 0.9 ) ) )  iindex ( 1 ) = 0

      ! Last level
      che1 = LOG ( t ( index - 1 ) / t ( index - 2 ) ) / LOG ( p ( index - 1 ) / p ( index - 2 ) ) 
      che2 = LOG ( t ( index     ) / t ( index - 1 ) ) / LOG ( p ( index     ) / p ( index - 1 ) ) 

      IF ( ( ( che2 .LT. -0.9 ) .OR. ( che2 .GT. 0.9 ) ) .AND.                   &
           ( ( che1 .GT. -0.9 ) .AND. ( che1 .LT. 0.9 ) ) ) iindex ( index ) = 1
!     IF ( ( ( che2 .LT. -0.9 ) .OR. ( che2 .GT. 0.9 ) ) .AND.                   &
!          ( ( che1 .LT. -0.9 ) .OR. ( che1 .GT. 0.9 ) ) )  iindex ( index ) = 0

      loop2: DO i = 2, index - 1
         che1 = LOG ( t ( i + 1 ) / t ( i     ) ) / LOG ( p ( i + 1 ) / p ( i     ) ) 
         che2 = LOG ( t ( i     ) / t ( i - 1 ) ) / LOG ( p ( i     ) / p ( i - 1 ) ) 
         th1 = t ( i - 1 ) * ( 100000.0 / p ( i - 1 ) ) ** RCP
         th  = t ( i     ) * ( 100000.0 / p ( i     ) ) ** RCP
         th2 = t ( i + 1 ) * ( 100000.0 / p ( i + 1 ) ) ** RCP

         IF ( ( ( ( che1 .LT. -0.9 ) .OR. ( che1 .GT. 0.9 ) )  .AND.                &
                ( ( che2 .LT. -0.9 ) .OR. ( che2 .GT. 0.9 ) ) ) .OR.                &
                ( ( th .LT. ( th1 - 10.0 ) ) .AND. ( th .LT. ( th2 - 10.0) ) ) .OR. &
                ( ( th .GT. ( th1 + 10.0 ) ) .AND. ( th .GT. ( th2 + 10.0) ) ) )    &
            iindex ( i ) = 1

         IF ( ( i .GE. 3 ) .AND. ( ( che1 .GT. -0.6 ) .AND. ( che1 .LT. 0.6 ) ) .AND.   &
                ( ( che2 .GT. -0.6 ) .AND. ( che2 .LT. 0.6 ) ) .AND.                    &
                ( ( che1 * che2 ) .LT. 0.0 ) )                                          &
            iindex ( i ) = 2
             
      END DO loop2
   
      loop3: DO i = 2, index - 1
         IF ( iindex ( i ) .EQ. 0 ) THEN
            ! Do nothing
         ELSE
            ! Change the sign of temperature ( in Celsius )
            che1 = LOG ( t ( i + 1 ) / ( -t ( i ) + 2. * 273.15 ) ) / LOG ( p ( i + 1 ) / p ( i ) ) 
            che2 = LOG ( ( -t ( i ) + 2. * 273.15 ) / t ( i - 1 ) ) / LOG ( p ( i ) / p ( i - 1 ) ) 

            th1 = t ( i - 1 ) * ( 100000.0 / p ( i - 1 ) ) ** RCP
            th  = ( -t ( i ) + 2. * 273.15 ) * ( 100000.0 / p ( i ) ) ** RCP
            th2 = t ( i + 1 ) * ( 100000.0 / p ( i + 1 ) ) ** RCP

            IF ( iindex ( i ) .EQ. 1 ) THEN
               IF ( ( ( ( che1 .LT. -0.6 ) .OR. ( che1 .GT. 0.6 ) )  .OR.               &
                      ( ( che2 .LT. -0.6 ) .OR. ( che2 .GT. 0.6 ) ) ) .OR.              &
                      ( ( th .LT. ( th1 - 10.0 ) ) .AND. ( th .LT. ( th2 - 10.0 ) ) ) .OR. &    
                      ( ( th .GT. ( th1 + 10.0 ) ) .AND. ( th .GT. ( th2 + 10.0 ) ) ) ) THEN    
                   ! Do nothing
               ELSE
                  iindex ( i ) = 0
                  tqc  ( i ) = wrong_t_sign
                  tdqc ( i ) = wrong_t_sign
                  rhqc ( i ) = wrong_t_sign
                  IF ( .NOT. eps_equal ( td ( i ) , missing_r , 1. ) ) THEN
                     depression = t(i) - td(i)
                  END IF
                  t  ( i ) = -t  ( i ) + 2. * 273.15
                  IF ( .NOT. eps_equal ( td ( i ) , missing_r , 1. ) ) THEN
                     td(i) = t(i) - depression
                     rh ( i ) = 100. * exp ( 5418.12 * ( 1./t(i) - 1./td(i)))
                  END IF

                  WRITE ( message , FMT = '(" Temp wrong sign at",a8,a8)' ) &
                  TRIM ( obs%location%id ) , TRIM ( obs%location%name )
                  error_number = 00363001
                  error_message(1:31) = 'vert_cons_check                '
                  error_message(32:)  = message
                  fatal = .false.
                  listing = .false.
                  IF ( print_vert ) CALL error_handler ( error_number , error_message , &
                  fatal , listing )
               END IF
            ELSE
               IF ( ( ( che1 .GT. -0.6 ) .AND. ( che1 .LT. 0.6 ) ) .AND.                 &
                    ( ( che2 .GT. -0.6 ) .AND. ( che2 .LT. 0.6 ) ) .AND.                 &
                    ( ( che1 * che2 ) .GT. 0.0 ) )  THEN
                  iindex ( i ) = 0
                  tqc  ( i ) = wrong_t_sign
                  tdqc ( i ) = wrong_t_sign
                  rhqc ( i ) = wrong_t_sign
                  IF ( .NOT. eps_equal ( td ( i ) , missing_r , 1. ) ) THEN
                     depression = t(i) - td(i)
                  END IF
                  t  ( i ) = -t ( i ) + 2. * 273.15
                  IF ( .NOT. eps_equal ( td ( i ) , missing_r , 1. ) ) THEN
                     td(i) = t(i) - depression
                     rh ( i ) = 100. * exp ( 5418.12 * ( 1./t(i) - 1./td(i)))
                  END IF
                  WRITE ( message , FMT = '(" Temp wrong sign at",a8,a8)' ) &
                  TRIM ( obs%location%id ) , TRIM ( obs%location%name )
                  error_number = 00363001
                  error_message(1:31) = 'vert_cons_check                '
                  error_message(32:)  = message
                  fatal = .false.
                  listing = .false.
                  IF ( print_vert ) CALL error_handler ( error_number , error_message ,  &
                  fatal , listing )
               ELSE
                  iindex ( i ) = 0
               END IF
            END IF
         END IF
      END DO loop3

      loop4: DO i = 1, index
         IF ( iindex ( i ) .EQ. 1 ) THEN
            WRITE ( message , FMT = '(&
            &  " Failed super-adiabatic check: ID=",a8,&    
            &  ", NAME="  ,a8, ", T=" ,f5.1," (K), P=",f6.1," (hPa)")' ) &
            obs%location%id , obs%location%name , t(i) , p(i)/100.
            error_number = 00363001
            error_message(1:31) = 'vert_cons_check                '
            error_message(32:)  = message
            fatal = .false.
            listing = .false.
            IF ( print_vert ) CALL error_handler ( error_number , error_message ,  &
            fatal , listing )
            tqc  ( i ) = t_fail_supa_inver
            tdqc ( i ) = t_fail_supa_inver
            rhqc ( i ) = t_fail_supa_inver
            t  ( i ) = missing_r
            td ( i ) = missing_r
            rh ( i ) = missing_r
         END IF
      END DO loop4
      NULLIFY ( current )

     !  Put the "Quality Flag" and changed data into "obs" data

      current => obs%surface
      loopthermal : DO WHILE ( ASSOCIATED ( current ) )
         DO i = 1, index
            IF ( eps_equal ( p ( i ), current%meas%pressure%data, 0.1 ) ) THEN 
               current%meas%temperature%data = t ( i )
               current%meas%dew_point%data   = td ( i )
               current%meas%rh%data          = rh ( i )
               IF ( current%meas%temperature%qc .NE. missing ) THEN
                  current%meas%temperature%qc = current%meas%temperature%qc + tqc ( i )
               END IF
               IF ( current%meas%dew_point%qc .NE. missing ) THEN
                  current%meas%dew_point%qc = current%meas%dew_point%qc + tdqc ( i )
               END IF
               IF ( current%meas%rh%qc .NE. missing ) THEN
                  current%meas%rh%qc = current%meas%rh%qc + rhqc ( i )
               END IF
            ELSE
               ! Do nothing
            END IF
         END DO
         current => current%next
      END DO loopthermal
      DEALLOCATE ( p )
      DEALLOCATE ( t )
      DEALLOCATE ( td )
      DEALLOCATE ( rh )
      DEALLOCATE ( tqc )
      DEALLOCATE ( tdqc )
      DEALLOCATE ( rhqc )
      DEALLOCATE ( iindex )
      NULLIFY ( current )

   END IF

   !  End checking thermal data

   !  Check for spikes in hodograph and mark them in wind data,
   !  if wind speed is too strong or wind direction is wrong. 
   
   !  Find out the vertical level number, then ALLOCATE arrays

   index = 1
   current => obs%surface
   DO WHILE ( ASSOCIATED ( current ) )
      IF ( ( .NOT. eps_equal ( current%meas%pressure%data,  missing_r, 1. ) ) .AND. & 
           ( .NOT. eps_equal ( psfc                      ,  missing_r, 1. ) ) .AND. & 
           ( current%meas%pressure%data .LE. psfc ) .AND. & 
           ( .NOT. eps_equal ( current%meas%speed%data,     missing_r, 1. ) ) .AND. & 
           ( .NOT. eps_equal ( current%meas%direction%data, missing_r, 1. ) ) ) THEN
         index = index + 1
      ELSE
         ! Do nothing, missing or bad data
      END IF
      current => current%next
   END DO
   index = index - 1

   IF ( ( index .LE. 3 ) .OR. ( eps_equal ( psfc , missing_r, 1. ) ) ) THEN
      !    The vertical levels are less than 3, do nothing
   ELSE !  Check wind data

      ALLOCATE ( p     ( index ) )
      ALLOCATE ( u     ( index ) )
      ALLOCATE ( v     ( index ) )
      ALLOCATE ( ws    ( index ) )
      ALLOCATE ( wd    ( index ) )
      ALLOCATE ( wsqc  ( index ) )
      ALLOCATE ( wdqc  ( index ) )
      ALLOCATE ( uqc   ( index ) )
      ALLOCATE ( vqc   ( index ) )
      ALLOCATE ( work1 ( index ) )
      ALLOCATE ( work2 ( index ) )
      NULLIFY  ( current )

      !  Read data set into array and perform quality check

      index = 1
      current => obs%surface
      DO WHILE ( ASSOCIATED ( current ) )
         IF ( ( .NOT. eps_equal ( current%meas%pressure%data,  missing_r, 1. ) ) .AND. & 
              ( current%meas%pressure%data .LE. psfc ) .AND. & 
              ( .NOT. eps_equal ( current%meas%speed%data,     missing_r, 1. ) ) .AND. & 
              ( .NOT. eps_equal ( current%meas%direction%data, missing_r, 1. ) ) ) THEN
            p    ( index ) = current%meas%pressure%data
            ws   ( index ) = current%meas%speed%data
            wd   ( index ) = current%meas%direction%data
            u    ( index ) = current%meas%u%data
            v    ( index ) = current%meas%v%data
            wsqc ( index ) = 0
            wdqc ( index ) = 0
            uqc  ( index ) = 0
            vqc  ( index ) = 0
            index = index + 1
         ELSE
            ! Do nothing, missing or bad data
         END IF
         current => current%next
      END DO
      index = index - 1
      NULLIFY ( current )

      loop5: DO i = 2, index - 1
         d1 = wd ( i - 1 )
         d2 = wd ( i + 1 )
         IF ( ( d2 - d1 ) .GT.  180. ) d1 = d1 + 360.
         IF ( ( d2 - d1 ) .LT. -180. ) d2 = d2 + 360.
         work1 ( i ) = wd ( i ) - 0.5 * ( d1 + d2 ) + 0.0001
         IF ( work1 ( i ) .LT. -180. ) work1 ( i ) = work1 ( i ) + 360.
         IF ( work1 ( i ) .GT.  180. ) work1 ( i ) = work1 ( i ) - 360.
         work2 ( i ) = ws ( i ) - 0.5 * ( ws ( i - 1 ) + ws ( i + 1 ) ) + 0.0001
      END DO loop5

      loop6: DO i = 3, index - 2
         IF ( ( ( ABS ( work2 ( i ) ) .LT. 50. ) .OR.               & ! The wind speed at this
            ( work2 ( i - 1 ) * work2 ( i + 1 ) .LE. 0.0 ) .OR.     & ! level is reasonable
            ( MAX ( work2 ( i + 1 ) / work2 ( i - 1 ),            & ! compared to up and 
              work2 ( i - 1 ) / work2 ( i + 1 ) ) .GT. 3. ) )       & ! down levels;
             .AND. ( ( ( wd ( i ) .LE. 360. ) .AND. ( ws ( i ) *    & ! The wind direction
              ABS ( work1 ( i ) ) .LT. 1600. ) ) .OR.               & ! at the current level
            ( ( work1 ( i - 1 ) * work1 ( i + 1 ) ) .LE. 0.0 ) .OR. & ! is reasonable
            ( MAX ( work1 ( i + 1 ) / work1 ( i - 1 ),            & ! compared to up
              work1 ( i - 1 ) / work1 ( i + 1 ) ) .GT. 3. ) ) )     & ! and down levels.
         THEN
            ! wind data is good, do nothing
         ELSE
            WRITE ( message , FMT = '(&
            &  " Spikes in wind detected: ID=",a8,&    
            &  ",NAME="  ,a8, ", Speed=" ,f6.1," (m/s), Direction=",f5.1," (degrees)")' )&
            obs%location%id , obs%location%name , ws(i) , wd(i)
            error_number = 00363002
            error_message(1:31) = 'vert_cons_check                '
            error_message(32:)  = message
            fatal = .false.
            listing = .false.
            IF ( print_vert ) CALL error_handler ( error_number , error_message ,  &
            fatal , listing )
            wsqc ( i ) = wrong_wind_data
            wdqc ( i ) = wrong_wind_data
            uqc  ( i ) = wrong_wind_data
            vqc  ( i ) = wrong_wind_data
            ws   ( i ) = missing_r
            wd   ( i ) = missing_r
            u    ( i ) = missing_r
            v    ( i ) = missing_r
         END IF
      END DO loop6

      !  Put the "Quality Flag" into "obs" data

      current => obs%surface
      loopwind : DO WHILE ( ASSOCIATED ( current ) )
         DO i = 1, index
            IF ( eps_equal ( p ( i ), current%meas%pressure%data, 0.1) ) THEN 
               current%meas%speed%data = ws ( i )
               current%meas%direction%data = wd ( i )
               current%meas%u%data = u ( i )
               current%meas%v%data = v ( i )
               IF ( current%meas%speed%qc .NE. missing ) THEN
                  current%meas%speed%qc     = current%meas%speed%qc + wsqc ( i )
               END IF
               IF ( current%meas%direction%qc .NE. missing ) THEN
                  current%meas%direction%qc = current%meas%direction%qc + wdqc ( i )      
               END IF
               IF ( ( current%meas%u%qc .NE. missing ) .AND. &
                    ( current%meas%v%qc .NE. missing ) ) THEN
                  current%meas%u%qc = current%meas%u%qc + uqc ( i )                    
                  current%meas%v%qc = current%meas%v%qc + vqc ( i )                    
               END IF
            ELSE
               ! Do nothing
            END IF
         END DO
         current => current%next
      END DO loopwind
      DEALLOCATE ( p       )
      DEALLOCATE ( u       )
      DEALLOCATE ( v       )
      DEALLOCATE ( ws      )
      DEALLOCATE ( wd      )
      DEALLOCATE ( wsqc    )
      DEALLOCATE ( wdqc    )
      DEALLOCATE ( uqc     )
      DEALLOCATE ( vqc     )
      DEALLOCATE ( work1   )
      DEALLOCATE ( work2   )
      NULLIFY    ( current )

   END IF

   !  End checking wind data

END SUBROUTINE vert_cons_check

!-------------------------------------------------------------------------------

END MODULE qc3
