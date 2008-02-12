!------------------------------------------------------------------------------
SUBROUTINE sample ( name , istagu , istagv , &
field1 , field2 , field3 , &
iewd , jnsd , kbu )

!  This routine calls the correct routines to provide printout
!  information to the calling routine, from the gridded data
!  on input.

   USE printout
 
   IMPLICIT NONE

   !  The required arguments through the SUBROUTINE call.

   CHARACTER *(*) , INTENT ( IN )                 :: name

   !  The OPTIONAL arguments, though there are rules for
   !  being present.  One and only one of the fields must
   !  be present.  The correct number of dimensions for that
   !  field must be present.

   REAL , DIMENSION ( : , : , : ) , &
                         OPTIONAL , INTENT ( IN ) :: field3
   REAL , DIMENSION ( : , : ) , &
                         OPTIONAL , INTENT ( IN ) :: field2
   REAL , DIMENSION ( : ) , &
                         OPTIONAL , INTENT ( IN ) :: field1

   INTEGER , OPTIONAL , INTENT ( IN )             :: iewd , &
                                                     jnsd , &
                                                     kbu
   INTEGER , OPTIONAL , INTENT ( IN )             :: istagu , istagv

   INCLUDE 'error.inc'
   INTERFACE
      INCLUDE 'error.int'
   END INTERFACE

   !  Easy error conditions to test for the gridded data.

   IF ( ( PRESENT (field1) .AND. PRESENT (field2) ) .OR. &
        ( PRESENT (field1) .AND. PRESENT (field3) ) .OR. &
        ( PRESENT (field2) .AND. PRESENT (field3) ) ) THEN
   
      !  001:  You can only request either field1, field2, or field3,
      !  exclusively.

      error_number = 001
      error_message(1:31) = 'sample                         '
      error_message(32:)  = ' Only one field may be used per call &
      &with sample for field ' // TRIM ( name )
      fatal = .false.
      listing = .false.
      CALL error_handler ( error_number , error_message , &
      fatal , listing )
   ELSE IF ( .NOT. ( PRESENT (field1) .OR. &
                     PRESENT (field2) .OR. &
                     PRESENT (field3) ) ) THEN

      !  002: You must request at least one field.

      error_number = 002
      error_message(1:31) = 'sample                         '
      error_message(32:)  = ' No field was given to sample &
      &for field ' // TRIM ( name )
      fatal = .false.
      listing = .false.
      CALL error_handler ( error_number , error_message , &
      fatal , listing )
   ELSE IF ( ( PRESENT (field3) .OR. PRESENT (field2) ) .AND. &
   ( .NOT. ( PRESENT (istagu) ) ) .AND. ( .NOT. ( PRESENT (istagv) ) ) ) THEN

      !  003:  For 2D or 3D fields, the staggering 
      !  configuration is required to be known.

      error_number = 003
      error_message(1:31) = 'sample                         '
      error_message(32:)  = ' Cross/dot point configuration is not &
      &specified for field ' // TRIM ( name )
      fatal = .true.
      listing = .false.
      CALL error_handler ( error_number , error_message , &
      fatal , listing )
   ENDIF

   !  Which of the printout samplers to call depends on the 
   !  number of indices in the array.  Only 1D, 2D, and 3D
   !  arrays are expected.  The three tests in each IF test
   !  are 1) does the array have the right dimension for the
   !  specification, 2) are all of the required dimensions
   !  available, 3) is the data array present.

   IF      ( ( PRESENT (iewd) .AND. PRESENT (jnsd) .AND. PRESENT (kbu) ) .AND. &
   ( PRESENT (field3) ) ) THEN
      CALL sample3 ( field3 , name , iewd , jnsd , kbu , istagu , istagv )
   ELSE IF ( ( PRESENT (iewd) .AND. PRESENT (jnsd) ) .AND. &
   ( PRESENT (field2) ) ) THEN
      CALL sample2 ( field2 , name , iewd , jnsd ,       istagu , istagv )
   ELSE IF ( (  PRESENT (kbu) ) .AND. &
   ( PRESENT (field1) ) ) THEN
      CALL sample1 ( field1 , name ,               kbu          )
   ELSE

      !  010: Error can be several things, let them sort it out.

      error_number = 010
      error_message(1:31) = 'sample                         '
      error_message(32:)  = ' The combination of the field, number of &
      &arguments do not match for field ' // TRIM ( name )
      fatal = .false.
      listing = .false.
      CALL error_handler ( error_number , error_message , &
      fatal , listing )
   ENDIF

END SUBROUTINE sample
