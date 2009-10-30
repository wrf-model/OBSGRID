!------------------------------------------------------------------------------
SUBROUTINE proc_namelist ( unit , filename , nml )
 
!  This routine reads in the NAMELIST file from the specified unit.
!  Data is stored in the NAMELIST data structure.  The NAMELIST data are 
!  compared for internal consistency, and simple errors.

   USE namelist

   IMPLICIT NONE

   INTEGER , INTENT ( IN )          :: unit
   CHARACTER *(*) , INTENT ( IN )   :: filename
   TYPE ( all_nml )                 :: nml

   INTEGER , DIMENSION ( 100 )      :: nml_read_errors
   INTEGER                          :: loop , i

   INCLUDE 'error.inc'
   INCLUDE 'namelist.inc'
   INCLUDE 'namelist.nml'
   INCLUDE 'namelist.common'

   INTERFACE 
      INCLUDE 'error.int'
   END INTERFACE

   !  Initialize the NAMELIST read error holder to no errors.

   nml_read_errors = 0

   !  Opening the NAMELIST file constitutes an important step.

   OPEN ( FILE = filename  , &
          UNIT = unit , &
          STATUS = 'OLD' , &
          ACCESS = 'SEQUENTIAL' , &
          FORM = 'FORMATTED' , &
          ACTION = 'READ' , &
          IOSTAT = error_number )

   !  Was there an error on the NAMELIST OPEN.

   IF ( error_number .NE. 0 ) then
      error_number = error_number + 00001000
      error_message(1:31) = 'proc_namelist                  '
      error_message(32:)  = ' Error in NAMELIST OPEN'
      fatal = .true.
      listing = .false.
      CALL error_handler ( error_number , error_message ,  &
      fatal , listing )
   ENDIF

   !  With an apparently OPENed NAMELIST file, let us 
   !  check some obvious errors with the file I/O.

   CALL exist_namelist ( unit , filename(1:LEN_TRIM(filename)) )

   ! Default grid_id 
   grid_id = 1

   remove_data_above_qc_flag = 200000
   remove_unverified_data = .FALSE.

   !  Initialize the array of radius of influence scans.  This permits an easy way
   !  to determine the number of requested scans without an additional input value.

   radius_influence = -1

   !  Initialize the array of observation file names to 'null'.

   obs_filename     = 'null                                            &
   &                                                                                    '

   !  Unless explicitly requested, 
   !     we are not going down the FDDA path.
   !     we are not triming down the domain.

   f4d = .FALSE.
   trim_domain = .FALSE.
   oa_3D_option = 0

   !  Input the NAMELIST data, from either the specified
   !  input unit, or the default standard in.

   READ ( UNIT = unit , NML = record1 , IOSTAT = nml_read_errors(1) )
   IF ( nml_read_errors(1) .NE. 0 ) WRITE ( UNIT = * , FMT = * ) 'Error in NAMELIST record 1'
   READ ( UNIT = unit , NML = record2 , IOSTAT = nml_read_errors(2) )
   IF ( nml_read_errors(2) .NE. 0 ) WRITE ( UNIT = * , FMT = * ) 'Error in NAMELIST record 2'
   READ ( UNIT = unit , NML = record3 , IOSTAT = nml_read_errors(3) )
   IF ( nml_read_errors(3) .NE. 0 ) WRITE ( UNIT = * , FMT = * ) 'Error in NAMELIST record 3'
   READ ( UNIT = unit , NML = record4 , IOSTAT = nml_read_errors(4) )
   IF ( nml_read_errors(4) .NE. 0 ) WRITE ( UNIT = * , FMT = * ) 'Error in NAMELIST record 4'
   READ ( UNIT = unit , NML = record5 , IOSTAT = nml_read_errors(5) )
   IF ( nml_read_errors(5) .NE. 0 ) WRITE ( UNIT = * , FMT = * ) 'Error in NAMELIST record 5'
   READ ( UNIT = unit , NML = record7 , IOSTAT = nml_read_errors(7) )
   IF ( nml_read_errors(7) .NE. 0 ) WRITE ( UNIT = * , FMT = * ) 'Error in NAMELIST record 7'
   READ ( UNIT = unit , NML = record8 , IOSTAT = nml_read_errors(8) )
   IF ( nml_read_errors(8) .NE. 0 ) WRITE ( UNIT = * , FMT = * ) 'Error in NAMELIST record 8'
   READ ( UNIT = unit , NML = record9 , IOSTAT = nml_read_errors(9) )
   IF ( nml_read_errors(9) .NE. 0 ) WRITE ( UNIT = * , FMT = * ) 'Error in NAMELIST record 9'
   error_number = SUM ( nml_read_errors )

   !  Was there an error on the NAMELIST read.  If so, then
   !  process the error.

   IF ( error_number .NE. 0 ) then
      error_number = error_number + 00001000
      error_message(1:31) = 'proc_namelist                  '
      error_message(32:)  = ' Error in NAMELIST read, unit&
      & specified by call.'
      fatal = .true.
      listing = .false.
      CALL error_handler ( error_number , error_message ,  &
      fatal , listing )
   END IF

   !  After the NAMELIST file is input, the unit needs to
   !  be properly closed.

   CLOSE ( UNIT = unit )

   !  Default input met_em file names 
   !  Leave the namelist as a placehold, but we are not using it 
   !  This overwrite anything from the namelist
   
   WRITE(fg_filename,'("./met_em.d",i2.2)') grid_id

   !  For completeness we output the NAMELIST record.

   WRITE ( UNIT = * , FMT = * ) 'Record 1 NAMELIST, no checks'
   WRITE ( UNIT = * , NML = record1 )
   WRITE ( UNIT = * , FMT = * ) 'Record 2 NAMELIST, no checks'
   WRITE ( UNIT = * , FMT = '(" &RECORD2")')
   WRITE ( UNIT = * , FMT = '(" FG_FILENAME = ",A)' ) TRIM(fg_filename)
   WRITE ( UNIT = * , FMT = '(" OBS_FILENAMES = ",A)' ) TRIM(obs_filename)
   WRITE ( UNIT = * , FMT = '(" /")' )
   WRITE ( UNIT = * , FMT = * ) 'Record 3 NAMELIST, no checks'
   WRITE ( UNIT = * , NML = record3 )
   WRITE ( UNIT = * , FMT = * ) 'Record 4 NAMELIST, no checks'
   WRITE ( UNIT = * , NML = record4 )
   WRITE ( UNIT = * , FMT = * ) 'Record 5 NAMELIST, no checks'
   WRITE ( UNIT = * , NML = record5 )
   WRITE ( UNIT = * , FMT = * ) 'Record 7 NAMELIST, no checks'
   WRITE ( UNIT = * , NML = record7 )
   WRITE ( UNIT = * , FMT = * ) 'Record 8 NAMELIST, no checks'
   WRITE ( UNIT = * , NML = record8 )
   WRITE ( UNIT = * , FMT = * ) 'Record 9 NAMELIST, no checks'
   WRITE ( UNIT = * , NML = record9 )

   !  Store the NAMELIST in the correct locations

   CALL store_namelist ( nml ) 

   !  Do simple-minded checks of the NAMELIST input with
   !  itself.  Report modifications to the error handler
   !  routines.

   CALL check_namelist ( nml ) 

END SUBROUTINE proc_namelist
