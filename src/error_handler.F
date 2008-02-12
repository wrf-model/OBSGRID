!------------------------------------------------------------------------------
SUBROUTINE error_handler ( error_number , error_message , fatal , listing )

!  This subroutine handles all error exits and warnings from
!  the program.  The error number and message are passed
!  via the calling arguments to this routine, directly
!  from the error producing portion of code.  The calling
!  routine's name is contained in the first 31 characters of
!  the error message.  The actual error message starts at
!  character #32 in the error_message variable.

!  This routine will terminate the program with a fatal error
!  condition if the fatal variable is .TRUE. on input.  For an 
!  uncorrectable error, the program should exit the program from 
!  this routine, otherwise it returns to the calling routine 
!  with the info/warning/error noted in the standard print output.

   USE error

   IMPLICIT NONE

   INTEGER , INTENT ( IN )                    :: error_number 
   CHARACTER *(*) , INTENT ( IN )             :: error_message
   LOGICAL , INTENT ( IN )                    :: fatal 
   LOGICAL , INTENT ( IN )                    :: listing

   LOGICAL                                    :: any_errors    = .false.

   TYPE ( error_info ) , POINTER              :: current      , &
                                                 head
   SAVE current , head , any_errors

   !  They may just want a listing of the errors thus far.  If this
   !  is a fatal error, we are going to get a listing anyways, so 
   !  no need to do it twice.

   listing_only_1 : IF (  listing .AND. ( .NOT. fatal ) ) THEN

      !  We should have been in here at least once for this to work.

      IF ( .NOT. any_errors ) THEN
         WRITE ( UNIT = * , FMT = '(//,"NO INFORMATION/WARNINGS/ERRORS RECORDED IN PROGRAM")' )
      ELSE
         WRITE ( UNIT = * , FMT = '(//,"Error Listing",/)' )
         CALL print_all_error ( head ) 
      ENDIF

   ENDIF listing_only_1

   !  Write error information to the standard output.  If this is a
   !  fatal error, note that before exit.

   IF ( fatal ) THEN
      WRITE ( UNIT = * , FMT = '(/,a60,/,a20)' )  &
      '------------------------------------------------------------', &
      'ERROR EXIT'
   ELSE IF( error_number .GT. 0 ) THEN
!     WRITE ( UNIT = * , FMT = '(/,a60,/,a20)' )  &
!     '------------------------------------------------------------', &
!     'ERROR DETECTION'
   ENDIF
      
   !  Check to see that we should be processing some error, else we can
   !  assume that this was just a sanity check for an error listing.
   !  If the error number <= 0, this is not an error, so we can skip
   !  the rest of the routine.
   
   listing_only_2 : IF ( error_number .GT. 0 ) THEN
!     WRITE ( UNIT = * , FMT = '("Error identifier ",i8,",")' , &
!     ADVANCE = 'YES' ) error_number
!     WRITE ( UNIT = * , FMT = '("Calling Routine: ",A,/,"Error message: ",A)' , &
!     ADVANCE = 'YES' ) TRIM (error_message(1:31)),TRIM(error_message(32:)) 
      WRITE ( UNIT = * , FMT = '(A)' ) TRIM(error_message(32:)) 

      !  Initialize the linked list during the first time inside.
   
      IF ( .NOT. any_errors ) THEN
         CALL initialize_error ( head , current ) 
      ENDIF
   
      !  Store this error information every time we are in this routine.
   
      CALL store_error ( current , error_number , error_message ) 
   
      IF ( fatal ) THEN
         
         !  If this is a fatal error, print out all messages, then exit.
   
!        WRITE ( UNIT = * , &
!        FMT = '(//,"FATAL ERROR detected, INFO/WARNING/ERROR history follows:",/)' )
!        CALL print_all_error ( head ) 
         STOP 'ERROR_EXIT in error_handler'

      ELSE
   
         !  Non-fatal error, so we just return to the calling routine.  We
         !  must keep track that we have scored a hit on an error for the
         !  error history.
   
!        WRITE ( UNIT = * , FMT = '(a,a,a,/,a,///)' ) &
!        'Returning from error_handler to ', &
!        TRIM ( error_message(1:LEN_TRIM(error_message(1:31)))), &
!        '.','------------------------------------------------------------'
         any_errors    = .true.
   
      ENDIF
   ENDIF listing_only_2
 
END SUBROUTINE error_handler
