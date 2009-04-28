!-------------------------------------------------------------------------------
SUBROUTINE proc_first_guess ( filename , &
bhi , bhr , num3d , num2d , num1d , &
t , u , v , uA , vA , uC , vC , h , rh , pres , terrain , &
latitude_x , longitude_x , latitude_d , longitude_d , &
slp_x , slp_C , sst , snow , &
iew_alloc , jns_alloc , kbu_alloc , pressure , &
print_analysis , &
current_date_8 , current_time_6 , date_char , icount )

!  This routine processes the modeling system first guess data.  This data is
!  typically:
!  1) interpolated from a global analysis (NMC , ECMWF)
!  2) interpolated from a global model (MRF, AVN)
!  3) from a previous mesoscale model forecast (ETA)

   USE input_data
   USE first_guess

   IMPLICIT NONE
   INCLUDE 'netcdf.inc'

   CHARACTER ( LEN = 132 ) ,INTENT ( IN )           :: filename
   INTEGER                                          :: num3d , &
                                                       num2d , &
                                                       num1d
   INTEGER                , INTENT ( IN )           :: iew_alloc  , &
                                                       jns_alloc  , &
                                                       kbu_alloc

   REAL                 , DIMENSION(kbu_alloc)      :: pressure

   REAL                 , DIMENSION(jns_alloc,iew_alloc,kbu_alloc) &
                                                    :: t  , &
                                                       u  , &
                                                       v  , &
                                                       uA , &
                                                       vA , &
                                                       uC , &
                                                       vC , &
                                                       h  , &
                                                       rh , &
                                                       pres

   REAL                  , DIMENSION(jns_alloc,iew_alloc)  &
                                                    :: terrain      , &
                                                       latitude_x   , &
                                                       longitude_x  , &
                                                       latitude_d   , &
                                                       longitude_d  , &
                                                       slp_x        , &
                                                       slp_C        , &
                                                       snow         , &
                                                       sst

   LOGICAL                                          :: print_analysis

   INTEGER , INTENT ( IN )                          :: current_date_8 , & 
                                                       current_time_6 , & 
                                                       icount

   INTEGER                                          :: loop , icu, icv

   CHARACTER (LEN=19)                               :: date_char
   INTEGER                                          :: met_ncid
   INTEGER                                          :: rcode

   INCLUDE 'big_header.inc'

   INTERFACE
      INCLUDE 'sample.int'
   END INTERFACE

   !  We need to keep track of where we are sticking the data for the FDDA option.

   IF ( initial_time ) THEN
      tt = first_time
   ELSE
      tt = second_time
   END IF

   !  Each date time is in a new file - so open every time     
      rcode = nf_open(filename, 0, met_ncid)
      IF ( rcode /= 0 ) THEN
         print*,"  ERROR opening file: ", trim(filename)
         STOP
      ENDIF

   !  Read in the analysis data for this time period.

   CALL read_first_guess ( met_ncid , &
   bhi , bhr , num3d , num2d , num1d , &
   t , u , v , uA , vA , uC , vC , h , rh , pres , terrain , &
   latitude_x , longitude_x , latitude_d , longitude_d , &
   slp_x , slp_C , sst , snow , pressure , &
   iew_alloc , jns_alloc , kbu_alloc , & 
   current_date_8 , current_time_6 , date_char , icount , print_analysis ) 

   !  Now that we have all of this data input, we can do a fast run
   !  through of the values.

   IF ( print_analysis ) THEN

      three_d : DO loop = 1 , num3d
         IF      ( all_3d(loop,tt)%small_header%staggering(1:1) .EQ. 'M' ) THEN
            icu = 1
            icv = 1
         ELSE IF ( all_3d(loop,tt)%small_header%staggering(1:1) .EQ. 'X' ) THEN
            icu = 0
            icv = 1
         ELSE IF ( all_3d(loop,tt)%small_header%staggering(1:1) .EQ. 'Y' ) THEN
            icu = 1
            icv = 0
         END IF
         CALL sample ( all_3d(loop,tt)%small_header%name(1:8) , field3=all_3d(loop,tt)%array , &
                       iewd=iew_alloc , jnsd=jns_alloc , kbu=kbu_alloc , istagu=icu , istagv=icv )
      END DO three_d

      two_d : DO loop = 1 , num2d
         IF      ( all_2d(loop,tt)%small_header%staggering(1:1) .EQ. 'M' ) THEN
            icu = 1
            icv = 1
         ELSE IF ( all_2d(loop,tt)%small_header%staggering(1:1) .EQ. 'X' ) THEN
            icu = 0
            icv = 1
         ELSE IF ( all_2d(loop,tt)%small_header%staggering(1:1) .EQ. 'Y' ) THEN
            icu = 1
            icv = 0
         END IF
         CALL sample ( all_2d(loop,tt)%small_header%name(1:8) , field2=all_2d(loop,tt)%array , &
                       iewd=iew_alloc , jnsd=jns_alloc , istagu=icu , istagv=icv )
      END DO two_d
   
   END IF

END SUBROUTINE proc_first_guess
