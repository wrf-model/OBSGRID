!------------------------------------------------------------------------------

!  This module is used only by the error_handler routine.  These utility
!  routines and the error data structure allow the program to record 
!  errors, retrieve the error history, and exit the program when a 
!  fatal error is detected.

MODULE error

   TYPE error_info
      INTEGER                                 :: num
      CHARACTER ( LEN=132 )                   :: msg
      TYPE ( error_info ) , POINTER           :: next
   END TYPE error_info


CONTAINS

   !------------------------------------------------------------------------------
   
   SUBROUTINE initialize_error ( head , current ) 
   
   !  Initialize the linked list of errors reported from the calling
   !  routines.  Send back out the two required pointers to remember
   !  where the beginning is and the current pointer location.
   
      IMPLICIT NONE
      
      TYPE ( error_info ) , POINTER                  :: head , &
                                                        current
      TYPE ( error_info ) , POINTER                  :: error_report
   
      ALLOCATE ( error_report ) 
      error_report%num = -1
      error_report%msg = '                                Head of List'
      NULLIFY ( error_report%next )
      head         => error_report
      current      => error_report
   
   END SUBROUTINE initialize_error
   
   !------------------------------------------------------------------------------
   
   SUBROUTINE print_all_error ( head ) 
   
   !  Print all of the stored error messages.  Head is where the
   !  linked list begins, current is what is used to traverse.
   
      IMPLICIT NONE
   
      TYPE ( error_info ) , POINTER                 :: head 
      TYPE ( error_info ) , POINTER                 :: temp
   
      INTEGER                                       :: loop
      
      loop = 0
      temp => head
      print_loop : DO 
         IF ( .NOT. ASSOCIATED ( temp ) ) EXIT print_loop
         WRITE ( UNIT = * , &
         FMT = '(a,i8,/,a,i8,/,a,a,/,a,a,/,a,/)' ) &
         ' Error Number : '   , loop , &
         ' Error Identifier: ', temp%num , &
         ' Error Issued By: ' , temp%msg(1:31) , &
         ' Error Message: '   , temp%msg(32:LEN_TRIM(temp%msg)) , &
         '------------------------------------------------------------'
         temp => temp%next
         loop = loop + 1
      END DO print_loop
   END SUBROUTINE print_all_error
   
   !------------------------------------------------------------------------------
   
   SUBROUTINE print_single_error ( current )
   
   !  Print a single stored error message.  Current is where the
   !  data is located to provide information.  This is useful only
   !  for debugging.
   
      IMPLICIT NONE
   
      TYPE ( error_info ) , POINTER                 :: current
      
      IF ( .NOT. ASSOCIATED ( current ) ) THEN
         WRITE ( UNIT = * , FMT = * ) 'No error data at this location'
      ELSE
         WRITE ( UNIT = * , &
         FMT = '(28x    ,a20,i8,a20,a31,a20,/,a101,/,a60,/)' ) &
         ADJUSTL ( 'Error Identifier: ' ) , current%num , &
         ADJUSTL ( 'Error Issued By: ' ) , ADJUSTL ( current%msg(1:31) ), &
         ADJUSTL ( 'Error Message: ' ) , &
         ADJUSTL ( current%msg(32:LEN_TRIM(current%msg)) ), &
         '------------------------------------------------------------'
      ENDIF
   END SUBROUTINE print_single_error
   
   !------------------------------------------------------------------------------
   
   SUBROUTINE store_error ( current , error_number , error_message ) 
   
   !  Store this information in the already initialized linked list.
   !  Return the new pointer location.
   
      IMPLICIT NONE
      
      TYPE ( error_info ) , POINTER           :: current 
   
      INTEGER                                 :: error_number
      CHARACTER *(*)                          :: error_message
   
      TYPE ( error_info ) , POINTER           :: error_report
   
      ALLOCATE ( error_report )
      error_report%num = error_number
      error_report%msg = error_message
      NULLIFY ( error_report%next )
      current%next => error_report
      current => error_report
   
   END SUBROUTINE store_error

END MODULE error
