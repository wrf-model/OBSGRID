PROGRAM window_obs

   IMPLICIT NONE

   CHARACTER (LEN=20) :: date_char
   CHARACTER (LEN=40) :: string1, string2 , string3 , string4
   CHARACTER (LEN=84) :: rpt_format 
   CHARACTER (LEN=22) :: meas_format 
   CHARACTER (LEN=14) :: end_format
  
   REAL :: xlat, xlon, xlatmin, xlatmax, xlonmin, xlonmax
   REAL :: ter , slp , val1, val2, val3, val4, val5, val6, val7, val8, val9, val10, val11, val12
   REAL :: p, z, t, td, spd, dir, u, v, rh, th

   INTEGER :: iunit , num_info , kx, iflag1,iflag2, iseq_num,iflag3, &
              int1, int2, int3, int11, int12, int13, int14, int15, &
              int16, int17, int18, int19, int20, int21, int22              
   INTEGER :: pf, zf, tf, tdf, spdf, dirf, uf, vf, rhf, thf
   INTEGER :: foundend

   LOGICAL :: bogus, keep , lflag1, lflag2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   xlatmin = 20
   xlatmax = 60
   xlonmin = -130
   xlonmax = -60
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   iunit=10

   rpt_format =  ' ( 2F20.5 , 2A40 ,  2A40 , 1F20.5 , 5I10 , 3L10 ,  2I10 , A20 ,  13( F13.5 , I7 ) ) '
   meas_format =  ' ( 10( F13.5 , I7 ) ) '
   end_format = ' ( 3 ( I7 ) ) ' 

   still_more_obs : DO
      READ ( UNIT = iunit , IOSTAT=foundend, FMT = rpt_format )                                  &
             xlat,xlon, string1 , string2 ,                                                & 
             string3 , string4 , ter, num_info, iflag1,iflag2,                             &
             iseq_num,iflag3,lflag1,bogus,lflag2,int1,int2, date_char ,                    &
             slp,int3,val1,int11, val2,int12, val3,int13,val4,int14,val5,int15,val6,int16, &
                      val7,int17, val8,int18, val9,int19,val10,int20,val11,int21, val12,int22
      IF ( foundend .LT. 0 ) THEN
         EXIT still_more_obs
      END IF

      
      IF ( ( xlat .GE. xlatmin ) .AND. ( xlat .LE. xlatmax ) .AND.    &
           ( xlon .GE. xlonmin ) .AND. ( xlon .LE. xlonmax ) ) THEN   
         keep = .TRUE.
      ELSE
         keep = .FALSE.
      END IF

      IF ( keep ) THEN
         WRITE( UNIT = iunit+1, FMT = rpt_format )                                            &
                xlat,xlon, string1 , string2 ,                                                & 
                string3 , string4 , ter, num_info, iflag1,iflag2,                             &
                iseq_num,iflag3,lflag1,bogus,lflag2,int1,int2, date_char ,                    &
                slp,int3,val1,int11, val2,int12, val3,int13,val4,int14,val5,int15,val6,int16, &
                         val7,int17, val8,int18, val9,int19,val10,int20,val11,int21, val12,int22
      END IF

      still_more_levels : DO
         READ ( UNIT = iunit , FMT = meas_format ) &
                p,pf, z,zf, t,tf, td,tdf, spd,spdf, dir,dirf,u,uf, v, vf, rh, rhf, th, thf
         IF ( keep ) THEN
            WRITE( UNIT = iunit+1, FMT = meas_format ) &
                   p,pf, z,zf, t,tf, td,tdf, spd,spdf, dir,dirf,u,uf, v, vf, rh, rhf, th, thf
         END IF
         IF ( ABS(p+777777.0).LT.1) THEN
            EXIT still_more_levels
         END IF
      END DO still_more_levels

      READ ( UNIT = iunit , FMT = end_format ) kx,iflag1,iflag2
      IF ( keep ) THEN
         WRITE ( UNIT = iunit+1 , FMT = end_format )  kx, iflag1,iflag2
      END IF
      
   END DO still_more_obs

END PROGRAM window_obs
