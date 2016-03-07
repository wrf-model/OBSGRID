!------------------------------------------------------------------------------
SUBROUTINE proc_final_analysis ( filename , filename_out , &
bhid , bhrd , t , u , v , uA , vA , uC , vC , h , rh , pres , terrain , &
latitude_x , longitude_x , latitude_d , longitude_d , &
slp_x , slp_C , sst , snow , tobbox , odis , &
pressure , ptop , &
iew_alloc , jns_alloc , kbu_alloc , iewd , jnsd , & 
print_header , print_analysis , & 
current_date_8 , current_time_6 , fdda_loop , icount_fdda , &
icount , total_count , interval , &
max_error_t , max_error_uv           , &
max_error_z , max_error_p , &
buddy_weight , date_char , root_filename , oa_3D_option , &
!BPR BEGIN
!intf4d , lagtem , oa_type , oa_3D_type , rad_influence )
intf4d , lagtem , oa_type , oa_3D_type , rad_influence, oa_psfc )
!BPR END

!  This routine is a driver for the required utilities to output the 
!  final analysis of this program.  The input values are the objectively 
!  analyzed data (the met fields) and the constant values (terrestrial
!  fields).  

   USE final_analysis
   USE observation

   IMPLICIT NONE 
   INCLUDE 'netcdf.inc'

   CHARACTER *(*)              , INTENT ( IN    ) :: filename
   CHARACTER *(*)              , INTENT ( INOUT ) :: filename_out
   CHARACTER *(*)              , INTENT ( IN    ) :: root_filename
   INTEGER , DIMENSION(50)                        :: bhid
   REAL    , DIMENSION(20)                        :: bhrd
   LOGICAL                                        :: print_header , &
                                                     print_analysis
   INTEGER                                        :: current_date_8 , &
                                                     current_time_6 , &
                                                     fdda_loop
   INTEGER                                        :: icount , total_count , &
                                                     interval , icount_fdda
   REAL                                           :: max_error_t  , &
                                                     max_error_uv , &
                                                     max_error_z  , &
                                                     max_error_p  , &
                                                     buddy_weight
   CHARACTER (LEN=19)                             :: date_char
   CHARACTER (LEN=24)                             :: fdda_date_24
   INTEGER , SAVE                                 :: interval_analysis
  
   INTEGER                                        :: rcode, met_ncid, oa_ncid, ilen
   INTEGER , SAVE                                 :: sfc_ncid
   CHARACTER (LEN=50)                             :: sfcfile

   INTEGER                                        :: intf4d
   LOGICAL                                        :: lagtem 
   INTEGER                                        :: oa_3D_option
   CHARACTER *(*)                                 :: oa_type , oa_3D_type
   INTEGER , DIMENSION(10)                        :: rad_influence
   !BPR BEGIN
   LOGICAL                                        :: oa_psfc
   !BPR END


   !  Temporary holding arrays for the header information.  What is in the header
   !  on input is the data from the first guess file.  To output the data, we need
   !  to store the final analysis header in the bhi, bhr, arrays.  So that
   !  we do not clobber the input data, we use this temporary storage.  

   INTEGER , DIMENSION (50) :: bhi_hold
   REAL    , DIMENSION (20) :: bhr_hold

   !  Temporary date integers

   INTEGER :: yy , mm , dd , hh , date8, time6

   !  The size of the data to be processed.

   INCLUDE 'first_guess_size.inc'
   INCLUDE 'first_guess.inc'
   INCLUDE 'proc_get_info_header.inc'
   REAL , DIMENSION(jns_alloc,iew_alloc)    :: tobbox , odis

   INTERFACE
      INCLUDE 'proc_get_info_header.int'
      INCLUDE 'sample.int'
   END INTERFACE
   
   !  The record header is initialized with the call to this
   !  routine.  That data is assigned to these arrays that
   !  are accessible from a COMMON block in the proc_get_info_header
   !  routine.  

   INCLUDE 'header.common'

   !  Save the headers!

   bhi_hold  = bhi
   bhr_hold  = bhr

   !  Since this routine is called repeatedly by time, we
   !  only OPEN the SFC file the first time with a particular UNIT number.

   IF ( ( icount .EQ. 1 ) .AND. ( fdda_loop .EQ. 2 ) ) THEN
      ilen = len_trim(root_filename)
      sfcfile = 'wrfsfdda_'//root_filename(ilen-2:ilen)
      rcode = nf_create( sfcfile, 0, sfc_ncid )
   END IF

   !  Now that the record header is modified, we put the data into 
   !  these storage locations, which are the common variables.  This
   !  allows them to be used to print out information from the
   !  record header from the proc_get_info_header routine.

   bhi  = bhid
   bhr  = bhrd

   IF      ( fdda_loop .EQ. 1 ) THEN

      ! for each time we need to create a new met_em output file
      rcode = nf_open  (filename, NF_CLOBBER, met_ncid)
      rcode = nf_create(filename_out, 0, oa_ncid )

      !  Put pressure values back on to Pa.
      slp_x(:jns_alloc-1,:iew_alloc-1) = slp_x(:jns_alloc-1,:iew_alloc-1) * 100.

      CALL write_analysis ( met_ncid , oa_ncid , &
      t , u , v , uA , vA , uC , vC , h , rh , pres , terrain , &
      latitude_x , longitude_x , latitude_d , longitude_d , &
      slp_x , slp_C , sst , tobbox , odis , &
      iew_alloc , jns_alloc , kbu_alloc , iewd , jnsd , date_char, print_analysis , &
!BPR BEGIN
!     oa_type , oa_3D_type , oa_3D_option , rad_influence )
      oa_type , oa_3D_type , oa_3D_option , rad_influence , oa_psfc )
!BPR END
   ELSE
      CALL make_date ( current_date_8 , current_time_6 , fdda_date_24 )
      slp_x(:jns_alloc-1,:iew_alloc-1) = slp_x(:jns_alloc-1,:iew_alloc-1) * 100.
      IF ( icount == 1 ) rcode = nf_open(filename_out, NF_NOWRITE, oa_ncid)
      CALL write_analysis_fdda ( oa_ncid , sfc_ncid , total_count, icount_fdda , &
      t , u , v , uA , vA , uC , vC , h , rh , pres , &
      terrain , slp_x , slp_C , snow , tobbox , odis , pressure , &
      iew_alloc , jns_alloc , kbu_alloc , iewd , jnsd , fdda_date_24 , ptop , &
      intf4d , lagtem , oa_type , rad_influence )
      IF ( icount == 1 ) rcode = nf_close(oa_ncid)
   END IF

   !  With all of the final analysis contained here, we can print
   !  out some data.

   IF      ( fdda_loop .EQ. 1 ) THEN
      CALL proc_get_info_header ( print_header , &
      iewd , jnsd , kbu , grid_id , map_projection , expanded , iewe , jnse , &
      dxc , lat_center , lon_center , cone_factor , true_lat1 , true_lat2 , pole , &
      dxd , ptop )
   END IF

   !  Now that we have all of this data input, we can do a fast run
   !  through of the values.

   IF ( ( print_analysis ) .AND. ( fdda_loop .EQ. 1 ) ) THEN

      CALL sample ( 'TT      ' , field3=t            , &
      iewd=iew_alloc , jnsd=jns_alloc , kbu=kbu_alloc , istagu=1 , istagv=1 )
      CALL sample ( 'UU      ' , field3=u            , &
      iewd=iew_alloc , jnsd=jns_alloc , kbu=kbu_alloc , istagu=1 , istagv=1 )
      CALL sample ( 'VV      ' , field3=v            , &
      iewd=iew_alloc , jnsd=jns_alloc , kbu=kbu_alloc , istagu=1 , istagv=1 )
      CALL sample ( 'GHT     ' , field3=h            , &
      iewd=iew_alloc , jnsd=jns_alloc , kbu=kbu_alloc , istagu=1 , istagv=1 )
      CALL sample ( 'RH      ' , field3=rh           , &
      iewd=iew_alloc , jnsd=jns_alloc , kbu=kbu_alloc , istagu=1 , istagv=1 )
   
      CALL sample ( 'HGT_M   ' , field2=terrain      , &
      iewd=iew_alloc , jnsd=jns_alloc ,                 istagu=1 , istagv=1 )
      CALL sample ( 'XLAT_M  ' , field2=latitude_x   , &
      iewd=iew_alloc , jnsd=jns_alloc ,                 istagu=1 , istagv=1 )
      CALL sample ( 'XLONG_M ' , field2=longitude_x  , &
      iewd=iew_alloc , jnsd=jns_alloc ,                 istagu=1 , istagv=1 )
      CALL sample ( 'XLAT_D  ' , field2=latitude_d   , &
      iewd=iew_alloc , jnsd=jns_alloc ,                 istagu=1 , istagv=1 )
      CALL sample ( 'XLONG_D ' , field2=longitude_d  , &
      iewd=iew_alloc , jnsd=jns_alloc ,                 istagu=1 , istagv=1 )
      CALL sample ( 'PMSL    ' , field2=slp_x        , &
      iewd=iew_alloc , jnsd=jns_alloc ,                 istagu=1 , istagv=1 )
   END IF

   !  Restore the headers!

   bhi  = bhi_hold  
   bhr  = bhr_hold 

END SUBROUTINE proc_final_analysis
