!------------------------------------------------------------------------------
PROGRAM main

!  This program starts with a background analysis, and performs
!  an objective analysis with available observations.  The quality
!  control procedure on the observations makes them suitable for inclusion
!  into data assimilation schemes.  

   USE namelist
   USE date_pack
   USE input_data

   IMPLICIT NONE
 
   INCLUDE 'error.inc'
   INCLUDE 'big_header.inc'

   TYPE ( all_nml )  :: nml ! all NAMELIST information from all NAMELIST records
   CHARACTER ( LEN = 132 ) :: filename, filename_out
   INTEGER                 :: unit

   INTERFACE

      SUBROUTINE proc_namelist ( unit , filename , nml )
         USE namelist
         INTEGER , INTENT ( IN )           :: unit
         CHARACTER *(*) , INTENT ( IN )    :: filename
         TYPE ( all_nml ) , INTENT ( OUT ) :: nml
      END SUBROUTINE proc_namelist

      INCLUDE 'proc_header.int'

      SUBROUTINE driver ( filename , filename_out , &
      bhi , bhr , nml , iew_alloc , jns_alloc, kbu_alloc , &
      current_date_8 , current_time_6 , date_char , icount , total_count , &
      mqd_count , mqd_abs_min )
         USE namelist
         INCLUDE 'big_header.inc'
         CHARACTER ( LEN = 132 ) ,    INTENT ( IN )    :: filename
         CHARACTER ( LEN = 132 ) ,    INTENT ( INOUT ) :: filename_out
         TYPE ( all_nml ) , INTENT (IN )               :: nml
         INTEGER , INTENT ( IN )                       :: iew_alloc , &
                                                          jns_alloc , &
                                                          kbu_alloc
         INTEGER , INTENT ( IN )                       :: current_date_8 , & 
                                                          current_time_6 , &
                                                          icount , total_count
         CHARACTER (LEN=19) , INTENT(IN)               :: date_char
         INTEGER                                       :: mqd_count , mqd_abs_min
      END SUBROUTINE driver

      INCLUDE 'error.int'
      INCLUDE 'proc_get_info_header.int'

   END INTERFACE

   INTEGER                                :: iew_alloc , &
                                             jns_alloc , &
                                             kbu_alloc

   INTEGER                                :: icount , &
                                             total_count    , &
                                             current_date_8 , & 
                                             current_time_6

   CHARACTER (LEN=19)                     :: current_date , next_date
   INTEGER                                :: century_year , &
                                             month , &
                                             day , &
                                             hour , &
                                             minute , &
                                             second , &
                                             interval, ilen
   LOGICAL                                :: is_used
   LOGICAL                                :: does_exist
   INTEGER                                :: mqd_count = 0, mqd_abs_min = 100
#ifdef NCARG
call opngks
#endif

      WRITE ( UNIT = * , FMT = '("                                 ")' ) 
      WRITE ( UNIT = * , FMT = '("################################ ")' ) 
      WRITE ( UNIT = * , FMT = '("          WRF OBSGRID            ")' ) 
      WRITE ( UNIT = * , FMT = '("          Version 3.7.0          ")' )   
      WRITE ( UNIT = * , FMT = '("          July 2014          ")' )  
      !!WRITE ( UNIT = * , FMT = '("     pre-release - 02/-9/10      ")' )  
      WRITE ( UNIT = * , FMT = '("################################ ")' ) 
      WRITE ( UNIT = * , FMT = '("                                 ")' ) 

   !  Read in the NAMELIST information.  This file is connected to the given
   !  compile-time name. All error processing on the NAMELIST data are 
   !  handled by this routine (all simple errors that allow consistency 
   !  checks, etc).  The namelist MODULE is required for the nml data TYPE.  
   !  The NAMELIST file is CLOSED during the routine, and the nml structure 
   !  is filled.

   DO unit=10,100
      INQUIRE(unit=unit, opened=is_used)
      IF (.not. is_used) EXIT
   END DO
   filename = 'namelist.oa'
   CALL proc_namelist ( unit , filename , nml )

   !  Compute the time perids that are to be processed.  This is specified
   !  in the NAMELIST.  The two dates are 19 digit long character strings
   !  of the form YYYY-MM-DD_HH:mm:ss, where:
   !   YYYY = year (1900 - 2099 are valid)
   !     MM = month of the year (01-12)
   !     DD = day of the month (01-31)
   !     HH = UTC hour of the day (00-23)
   !     mm = minute of the hour (00-59)
   !     ss = second of the minute (00-59)
   !  The time interval between the starting and ending times is an integer
   !  specified in seconds.

   interval = nml%record_1%interval
   current_date = nml%record_1%start_date

   CALL geth_newdate ( next_date , current_date , interval )
   
   !  With the starting and ending time, print out the time periods this 
   !  program will process.

   IF ( nml%record_7%f4d) THEN
      WRITE ( UNIT = * , FMT = * ) '3d analysis dates to be processed by this program (excluding SFC FDDA times):'
   ELSE
      WRITE ( UNIT = * , FMT = * ) 'Dates to be processed by this program:'
   END IF

   total_count = 0
   time_loop_1 : DO 

      !  Print out the loop counter increment and the computed date.

      total_count = total_count + 1 
      WRITE ( UNIT = * , FMT = '("      Time period #",i5.5," is for date ",A)' ) &
      total_count , current_date

      !  The next date is the current date plus a time interval.

      CALL geth_newdate ( next_date , current_date , interval )
      current_date = next_date

      !  Exit the loop if we have passed the last requested time period.  Exit with
      !  a fatal error if the NAMELIST request too many time periods (this is
      !  probably a mistake in setting up the NAMELIST).

      IF ( next_date .GT. nml%record_1%end_date ) THEN
         EXIT time_loop_1
      ELSE IF ( total_count .GT. max_times ) THEN
         error_number        = 1
         error_message(1:31) = 'main                           '
         error_message(32:)  = ' Too many time periods for processing have been specified.'
         fatal               = .TRUE.
         listing             = .TRUE.
         CALL error_handler ( error_number , error_message ,  &
         fatal , listing )
      END IF
   END DO time_loop_1

   !  ALLOCATE space for the temporary holders for the 3d, 2d and the 1d
   !  data.  Each of these has a large array space associated with it,
   !  but since we are not ALLOCATing that really big stuff, this is
   !  not a problem.   

   ALLOCATE ( all_3d( 20,2) )
   ALLOCATE ( all_2d(100,2) )
   ALLOCATE ( all_1d( 10) )


   !  Pass all of the data from the NAMELIST and record header to the
   !  driver routine.  Pick off the first guess field and the logical unit
   !  number.  Re-initialize the current date to the first time period
   !  requested.

   current_date = nml%record_1%start_date
   filename = nml%record_2%fg_filename 
   filename = trim(filename)//"."//current_date//".nc"
   INQUIRE ( EXIST = does_exist , FILE = filename )
   IF ( .NOT. does_exist ) THEN
      WRITE ( UNIT = * , FMT = '("   ")' ) 
      WRITE ( UNIT = * , FMT = '("###   Could not find file: ",A )' ) trim(filename)
      STOP '      STOP: Missing input file'
   ENDIF

   !  Now that the NAMELIST has been input, the other source of initial
   !  data is from the record header.  This file provides the information
   !  concerning the specific anlaysis data to be ingested.  This file is
   !  the same as the analysis file to be used later, so this routine CLOSEs
   !  the file after the initial time of the record header is input.  The
   !  NAMELIST structure is passed for error checking purposes.
   !  Initialize the domain size constants, pass them through to the 
   !  driver routine to allow the data arrays to be allocated.

   CALL proc_header ( filename , bhi , bhr , nml )

   !  We will always do the analysis on the incoming domain size - iewe & jnse
   !  We will only need the cut down size for the output

   CALL proc_get_info_header ( nml%record_5%print_header , &
   iewe=iew_alloc , jnse=jns_alloc ,  &
   kbu=kbu_alloc )

   icount = 0

   time_loop_2 : DO

      filename = nml%record_2%fg_filename 
      filename = trim(filename)//"."//current_date//".nc"
      WRITE(filename_out,'("./metoa_em.d",i2.2,".")') nml%record_2%grid_id
      filename_out = trim(filename_out)//current_date//".nc"

      !  Which time period are we processing.

      icount = icount + 1 

      !  Set up the counting for the large arrays for FDDA.  Instead of copying
      !  arrays around, we are just going to change the last index and effectively
      !  just change the variable pointers.

      IF ( nml%record_7%f4d ) THEN
         IF       ( icount      .EQ. 1      ) THEN
            first_time  = 1
            !BPR BEGIN
            !If we are at the first time there is no second time
            !If we set second_time=2 it will try to use this in interpolation
            !even though it does not yet exist
            !second_time = 2
            second_time = 1
            !BPR END
            initial_time = .TRUE.
         ELSE IF  ( icount      .EQ. 2      ) THEN
            first_time  = 1
            second_time = 2
            initial_time = .FALSE.
         ELSE IF ( (icount/2)*2 .NE. icount ) THEN
            first_time  = 2
            second_time = 1
         ELSE IF ( (icount/2)*2 .EQ. icount ) THEN
            first_time  = 1
            second_time = 2
         END IF
      ELSE IF ( .NOT. nml%record_7%f4d ) THEN
         first_time  = 1
         second_time = 1
         initial_time = .TRUE.
      END IF
     

      !  Compute the integer date (YYYYMMDD) and the integer time (HHmmss).

      CALL split_date_char ( current_date , century_year , month , day , hour , minute , second )

      current_date_8 = century_year * 10000 + month  * 100 + day
      current_time_6 = hour         * 10000 + minute * 100 + second

      !  This routine is called once for each time period to be processed.  This
      !  is the main driver routine for the program.


      CALL driver ( filename , filename_out , & 
      bhi , bhr , nml , iew_alloc , jns_alloc , kbu_alloc , &
      current_date_8 , current_time_6 , current_date, icount , total_count , &
      mqd_count , mqd_abs_min )

      !  Increment to the next time and check if we should try to process the
      !  data at that time.  The only way to not process the data is if we have
      !  gone past the requested ending time.
     
      CALL geth_newdate ( next_date , current_date , interval )
      current_date = next_date

      IF ( current_date .GT. nml%record_1%end_date ) THEN
         EXIT time_loop_2
      END IF

   END DO time_loop_2

   IF ( mqd_count > 0 ) THEN
        mqd_count = mqd_count + 1
        WRITE ( UNIT = * , FMT = '( /,"----------------------------------------------------------------------------",/ )' )
        WRITE(*, '( "  ########################################################################" )' )
        WRITE(*, '( "  NOTE: MQD was requested for all upper-levels.                           " )' )
        WRITE(*, '( "  But due to lack of data, CRESSMAN scheme was performed on some levels.  " )' )
        WRITE(*, '( "  To switch all upper levels to Cressman scheme, set either :             " )' )
        WRITE(*, '( "     oa_3D_switch = .TRUE. , or                                           " )' )
        WRITE(*, '( "     mqd_minimum_num_obs = " ,  i3, " and oa_3D_option = 1 or 2           " )' ) mqd_count
        WRITE(*, '( "                                                                          " )' ) 
        WRITE(*, '( "  You could decrease the value of mqd_minimum_num_obs, to perform         " )' )
        WRITE(*, '( "     MQD on all upper-levels, but this may result in a bad analysis.      " )' )
        WRITE(*, '( "     For this run you would need to set this parameter to " , i3            )' ) mqd_abs_min
        WRITE(*, '( "  ########################################################################" )' )
   END IF

   WRITE ( UNIT = * , FMT = '( /,"----------------------------------------------------------------------------",/,&
         &"Successful completion of obsgrid")' )
   WRITE ( UNIT = * , FMT = '("    ")' )

#ifdef NCARG
call clsgks
#endif

END PROGRAM main 

