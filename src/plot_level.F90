PROGRAM plot_obs

!  PROGRAM to plot the obs that were used in obsgrid

   USE header1
   USE map_stuff
   USE ncarg
   USE read_station
   USE date_pack

  INCLUDE 'netcdf.inc'

!  IMPLICIT NONE

   
   integer :: funit, ios
   character(len=120) :: flnm, oa_file, cgm_file
   logical :: is_used
   INTEGER :: start_year, start_month, start_day, start_hour, start_minute, start_second
   INTEGER ::   end_year,   end_month,   end_day,   end_hour,   end_minute,   end_second
   INTEGER :: interval, idiff, n_times
   CHARACTER (LEN=19) :: start_date, end_date, rdate
   CHARACTER ( LEN = 132)  :: obs_filename
   LOGICAL :: trim_domain, remove_unverified_data
   INTEGER :: trim_value, grid_id, remove_data_above_qc_flag
   INTEGER :: met_ncid
   INTEGER :: rcode
   LOGICAL :: use_first_guess, f4d, lagtem
   INTEGER :: intf4d


   REAL    , PARAMETER :: pi = 3.14159265, rpd = pi/180.
   INTEGER , PARAMETER :: nsta = 3000

   !  Some character space for the labels on the plots.

   CHARACTER (LEN=84) :: title, info
   CHARACTER (LEN=10) :: date
   CHARACTER (LEN= 8) :: fname,name,idstation,units
   CHARACTER (LEN= 8) :: sid(nsta)
   CHARACTER (LEN= 6) :: lab
   CHARACTER (LEN= 1) :: char1

   !  The values for u, v, T, SLP, rh, and (x,y) for each station.
   REAL , DIMENSION(nsta) :: u , v , t , slp , rh , xc , yc

   !  QC values for each variable, at each station location for a level.
   INTEGER , DIMENSION(nsta) :: wqc , tqc , rhqc , pqc 

   !  Pressure level for this level.
   INTEGER :: opress

   !  Some sort of domain size.
   INTEGER :: imax , jmax

   !  Error return code.
   INTEGER :: ier

   !  Pressure used for testing which level to process.
   INTEGER :: ipress

   !  How many observations at this level.
   INTEGER :: num

   !  Some things for the set call.
   REAL :: px1,px2,py1,py2,xx1,xx2,yy1,yy2
   INTEGER :: lll

   !  Label sizes.
   REAL :: size , size5
   INTEGER :: isz , m1 , m2 , ma

   INTEGER :: i

   !  Wind speed and direction.
   REAL :: skt , dir

   
    !  Input file stuff
    namelist  /record1/ start_year, start_month, start_day, start_hour, start_minute, start_second, &
                          end_year,   end_month,   end_day,   end_hour,   end_minute,   end_second, &
                        interval
    namelist /record2/ obs_filename, remove_data_above_qc_flag, & 
                      trim_domain, trim_value, grid_id, remove_unverified_data

    namelist /record7/ use_first_guess, f4d, intf4d, lagtem

   ! default 
   start_minute = 0
   start_second = 0
   end_minute = 0
   end_second = 0
   grid_id = 1
   
   print*," "
   
   ! Read parameters from Fortran namelist
   DO funit=15,100
      inquire(unit=funit, opened=is_used)
      IF (.not. is_used) EXIT
   END DO
   OPEN(funit,file='namelist.oa',status='old',form='formatted')
   READ(funit,record1)
   READ(funit,record2)
   READ(funit,record7)
   CLOSE(funit)

   ! Which interval should we use
   IF ( f4d ) interval = intf4d
   
   ! Build starting date string
   WRITE(start_date, '(i4.4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)') &
         start_year,'-',start_month,'-',start_day,'_',start_hour,':',start_minute,':',start_second
   
   ! Build ending date string
   WRITE(end_date, '(i4.4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)') &
         end_year,'-',end_month,'-',end_day,'_',end_hour,':',end_minute,':',end_second
 
   ! See which times we are interested in
   CALL geth_idts(end_date, start_date, idiff)
   IF ( idiff < 0 ) idiff = 1
   n_times = idiff / interval
 
   ! Build metoa file
   WRITE(oa_file, '("metoa_em.d",I2.2".",A19,".nc")') grid_id, start_date
   print*," Attempting to open file: ", trim(oa_file)
   rcode = nf_open(oa_file, 0, met_ncid)
   IF ( rcode /= 0 ) THEN
      print*," ERROR opening file: ", trim(oa_file)
      print*,"       Cannot continue without this information"
      STOP
   ENDIF
   print*," "

   !  Need to know the domain size to generate a map.
   CALL read_header ( met_ncid , imax , jmax ) 

   rdate = start_date
   LOOP_TIMES : DO itimes = 1, n_times+1

     print*," "
     ! name the filename
     WRITE(flnm,FMT='("plotobs_out.d",i2.2,".",A19,".0000")') grid_id, rdate

     ! Build graphics output file
     WRITE(cgm_file, '("levels.d",i2.2,".",A19,".cgm")') grid_id, rdate

     !  Initialize the NCAR Graphics stuff.
     CALL start_it(cgm_file)

     ! Date information
     READ(rdate,FMT='(    I4.4)') start_year
     READ(rdate,FMT='( 5X,I2.2)') start_month
     READ(rdate,FMT='( 8X,I2.2)') start_day
     READ(rdate,FMT='(11X,I2.2)') start_hour
     WRITE(date,FMT='(i4.4,i2.2,i2.2,i2.2)') start_year,start_month,start_day,start_hour
     print*,"Working on date: ", date
     print*,'Plots written to: ', trim(cgm_file)
     
     opress = -999
     ier = 0
     ma = 1
  
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     !  Loop over all of the stations that we are going to process for this
     !  selected time period.
     
     process_level : DO
  
        !  Read a single level of this data.
  
        CALL rdsta (u, v, t, rh, slp, ipress, opress, sid, xc, yc, wqc, tqc, rhqc, pqc, &
                    nsta, num, ier, trim(flnm) )
  
        IF ( ma .NE. 1 ) THEN
           PRINT '(A)',title
           CYCLE process_level
        END IF

        !  With the big header info, generate the map background.  Make the
        !  line width the minimum possible.
  
        CALL plotmap ( met_ncid, grid_id )
        CALL setusv('LW',1000)
  
        !  Scale the label size based on the pressure level - more surface obs, so
        !  make the labels smaller.
        
        IF (opress .ne. 1001) THEN
           isz = 9
        ELSE
           isz = 8
        END IF
  
        CALL sflush
  
        CALL pcsetc('FC - FUNCTION CODE CHARACTER', '|')
        CALL getset(px1,px2,py1,py2,xx1,xx2,yy1,yy2,lll)
  
        size=MAX( (px2-px1),(py2-py1) ) / 35.
        size5=size*.2
  
        !  For this level, loop over each of the observations.
  
        each_station : do i = 1, num
  
           !  Set the text color index to color #1.
  
           CALL gstxci(1)
  
           !  If the 5 character station ID is not something stupid, use the name,
           !  ELSE, just stick in a location sign.
  
           IF ( index(sid(i),'-----').ne.0) THEN
              lab = '  +  '
           ELSE
              WRITE(lab,fmt='(a5)') sid(i)
           END IF

           !  Station ID in black, print out the station ID.
  
           CALL setcl (0,sid(i))
           CALL plchhq(xc(i),yc(i),'|F4V-80H110|'//lab,size5,0.0,0.0)
           !  Process the temperature and, if this is the surface, the slp
  
           CALL setcl (tqc(i),sid(i))
           m2 = (isz-4)/2
           m1 = m2 - 1
           lab = ' '
           IF (t(i) .gt. -998.) WRITE(lab,fmt='(f6.1)') t(i)-273.15
           CALL plchhq(xc(i),yc(i),'|F4V60H-120|'//lab,size5,0.0,0.0)
      !    CALL wtstr(xc(i)-m2,yc(i)+m1,lab,isz,0,0)
           CALL sflush
           IF (opress .eq. 1001) THEN
              CALL setcl (pqc(i),sid(i))
              lab = ' '
              IF (slp(i) .gt. -998.) WRITE(lab,fmt='(f6.1)') slp(i)/100.
              CALL plchhq(xc(i),yc(i),'|F4V60H110|'//lab,size5,0.0,0.0)
    !         CALL wtstr(xc(i)+m2,yc(i)+m1,lab,isz,0,0)
           END IF
  
           !  Process the dew point.
  
           CALL setcl (rhqc(i),sid(i))
           lab = ' '
           IF (rh(i) .gt. -998.) WRITE(lab,fmt='(i3)') nint(rh(i))
           CALL plchhq(xc(i),yc(i),'|F4V-60H-70|'//lab,size5,0.0,0.0)
  !        CALL wtstr(xc(i)-m1,yc(i)-m1,lab,isz,0,0)
  
           !  Process the wind speed and direction.
  
           CALL setcl (wqc(i),sid(i))
           IF (u(i) .gt. -999. .and. v(i) .gt. -999.) THEN
              IF (abs(u(i)) .lt. .00001 .and. abs(v(i)) .lt. .00001) THEN
                 skt = 0.
                 dir = 0.
              ELSE
                 skt = SQRT(u(i)*u(i) + v(i)*v(i))
                 dir = ACOS(-v(i)/skt) * (180./3.14159265)
              END IF
              IF (u(i) .gt. 0.) dir = 360. - dir
              skt = skt * 1.9
              CALL barb (dir, skt, xc(i), yc(i) )
           END IF
  
        END DO each_station
  
        !  Put the color key on the plot.
  
        IF ( ma .eq. 1) THEN
           CALL sflush
           CALL gstxci(1)
           CALL gsplci(1)
           WRITE (title, fmt='(" Real-time input observations  ",a10, "   LEVEL = ",i4, "   NO. OF OBS = ",i8)' ) date, opress, num
           CALL wtstr(float(jmax/2),1.04*float(imax),title(1:84),12,0,0)
           CALL set (0.,1.,0.,1.,0.,1.,0.,1.,1)
           CALL sflush
           CALL gstxci(2)
           info = 'Red: Failed errmax test'
           CALL wtstr(.05,.12,info,12,0,-1)
           CALL sflush
           CALL gstxci(3)
           info = 'Orange: Extrapolated from single level'
           CALL wtstr(.05,.10,info,12,0,-1)
           CALL sflush
           CALL gstxci(4)
           info = 'Violet: Vertically interpolated and failed errmax'
           CALL wtstr(.05,.08,info,12,0,-1)
           CALL sflush
           CALL gstxci(5)
           info = 'Green: Vertically interpolated'
           CALL wtstr(.05,.06,info,12,0,-1)
           CALL sflush
           CALL gstxci(7)
           info = 'Gray: Single level interpolated, failed errmax'
           CALL wtstr(.05,.04,info,12,0,-1)
           CALL sflush
           CALL gstxci(6)
           info = 'Cyan: Other errmax failure'
           CALL wtstr(.05,.02,info,12,0,-1)
           CALL sflush
           CALL frame
        END IF
  
        PRINT '(A)',title
  
        !  Have we finished reading all of the different levels?
  
        IF ( ier .NE. 0 ) EXIT process_level
  
     END DO process_level
  
     !  Shut down NCAR Graphics.
     CALL stop_it

     CALL geth_newdate(rdate, trim(start_date), itimes*interval)

  END DO LOOP_TIMES

  print*," "
  print*," All done"

END PROGRAM plot_obs
