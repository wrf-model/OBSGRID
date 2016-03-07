program plot_raobs

  use report_module
  use skewt_module
  use mapinfo_module
  use date_pack

  INCLUDE 'netcdf.inc'

  logical :: lcf

  type(report) :: rpt

  integer :: iunit = 10
  integer :: funit, ios

  character(len=120) :: flnm, oa_file, cgm_file

  real, dimension(200) :: pres, temp, dwpt, spd, dir
  real, allocatable, dimension(:) :: plvls

  real, dimension(8) :: xarr

  logical :: rfile = .FALSE.
  logical :: read_metoa = .FALSE.
  
  logical :: DO_THE_PLOT, is_used

  INTEGER :: start_year, start_month, start_day, start_hour, start_minute, start_second
  INTEGER ::   end_year,   end_month,   end_day,   end_hour,   end_minute,   end_second
  INTEGER :: interval, idiff, n_times
  CHARACTER (LEN=19) :: start_date, end_date, rdate
  CHARACTER ( LEN = 132)  :: obs_filename
  LOGICAL :: trim_domain, remove_unverified_data
  INTEGER :: trim_value, grid_id, remove_data_above_qc_flag
  LOGICAL :: use_first_guess, f4d, lagtem
  INTEGER :: intf4d

  INTEGER, DIMENSION(4)                        :: iend
  INTEGER                                      :: met_ncid, iv, idims
  INTEGER                                      :: var_type, var_ndims, var_natt
  INTEGER, DIMENSION(6)                        :: var_shape
  REAL, ALLOCATABLE, DIMENSION(:,:,:)          :: met_em_dum1
  INTEGER                                      :: rcode, ndims, nvars, ngatts, nunlimdimid
  INTEGER, ALLOCATABLE, DIMENSION(:)           :: dim_values
  CHARACTER (LEN=31),ALLOCATABLE, DIMENSION(:) :: dim_names
  CHARACTER (LEN=80)                           :: cval

  CHARACTER (LEN=80)                           :: file_type

  namelist  /record1/ start_year, start_month, start_day, start_hour, start_minute, start_second, &
                        end_year,   end_month,   end_day,   end_hour,   end_minute,   end_second, &
                      interval
  namelist /record2/ obs_filename, remove_data_above_qc_flag, &
                    trim_domain, trim_value, grid_id, remove_unverified_data
  namelist /record7/ use_first_guess, f4d, intf4d, lagtem
  namelist /plot_sounding/ file_type, read_metoa
 
  ! default
  file_type = 'useful'
  grid_id    = 1
  start_minute = 0
  start_second = 0
  end_minute = 0
  end_second = 0

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
  READ(funit,plot_sounding)
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
  !BPR BEGIN
  !WRITE(oa_file, '("metoa_em.d",I2.2".",A19,".nc")') grid_id, start_date
  WRITE(oa_file, '("metoa_em.d",I2.2,".",A19,".nc")') grid_id, start_date
  !BPR END
  
  IF ( read_metoa ) THEN
     print*," Attempting to open file: ", trim(oa_file)
     rcode = nf_open(oa_file, 0, met_ncid)
     IF ( rcode /= 0 ) THEN
        rfile = .FALSE.
        print*," ERROR opening file: ", trim(oa_file)
        print*,"       Continue without this information"
        print*," "
     ELSE   
        rfile = .TRUE.
        print*," Found the file - will use map info from this file"
        print*," "
        rcode = nf_inq(met_ncid, ndims, nvars, ngatts, nunlimdimid)
        ALLOCATE ( dim_values(ndims) )
        ALLOCATE ( dim_names(ndims) )
        DO idims = 1, ndims
           rcode = nf_inq_dim(met_ncid, idims, dim_names(idims), dim_values(idims))
        ENDDO
  
        iend = 1
        cval = "PRES"
        rcode = nf_inq_varid ( met_ncid, cval, iv )
        rcode = nf_inq_var(met_ncid, iv, cval, var_type, var_ndims, var_shape, var_natt)
        DO idims=1,var_ndims-1
           iend(idims)  = dim_values(var_shape(idims))
        ENDDO
        IF ( ALLOCATED(met_em_dum1) ) DEALLOCATE(met_em_dum1)
        ALLOCATE(plvls(iend(3)))
        ALLOCATE ( met_em_dum1(iend(1),iend(2),iend(3)) )
        rcode = nf_get_var_real(met_ncid, iv, met_em_dum1)
        plvls(:) = met_em_dum1(1,1,:)

        call get_mapinfo( met_ncid )
     ENDIF
  ENDIF 


  rdate = start_date
  LOOP_TIMES : DO itimes = 1, n_times+1

     WRITE(flnm,FMT='("qc_obs_",A,".d",i2.2,".",A19,".0000")') trim(file_type), grid_id, rdate
     open(iunit, file=flnm, form='formatted', status='old', action='read',iostat=ios)
     IF ( ios == 0 ) THEN
        print*," Reading file: ", trim(flnm)
     ELSE
        print*," Looking for file: ", trim(flnm)
        print*,"         not available - must be end of input file"
        print*," "
        exit LOOP_TIMES
     ENDIF

     ! Build metoa file
     WRITE(cgm_file, '("soundings_",A,".d",i2.2,".",A19,".cgm")') trim(file_type), grid_id, rdate
     print *,' Plots written to: ', trim(cgm_file)
     print*," "
     call skewt_opngks(cgm_file)

     call read_r_report(iunit, ierr, rpt)

     do while (ierr == 0)

        if (rfile) then
           call lltoxy_mapinfo(rpt%xlat, rpt%xlon, xloc, yloc)
        else
           xloc = -999999.0
           yloc = -999999.0
        endif
   
        ilw = 2000
        ictl = 6   ! Select Temperature color
        icdl = 16  ! Select Dewpoint color
        lcf = .FALSE.
        woffs = 0.
        badval = -999999.
        if (rpt%num_lvls > 2) then
           ! Build the data arrays for the skewt plotting routines:
           j = 0
           pres = -999999.
           temp = -999999.
           dwpt = -999999.
           spd  = -999999.
           dir  = -999999.
           DO_THE_PLOT = .FALSE.
           ILOOP : do i = 1, rpt%num_lvls
              if ((rpt%pressure(i)%qc .lt. 2**15).and.(rpt%pressure(i)%data > 0.)) then
                 if (rfile) then
                    do k = 1, size(plvls)
                       if (abs(plvls(k)-rpt%pressure(i)%data) .lt. 1) then
                          exit
                       endif
                       if (k.eq.size(plvls)) cycle ILOOP
                    enddo
                 endif
                 j = j + 1
                 DO_THE_PLOT = .TRUE.
                 pres(j) = rpt%pressure(i)%data
                 if (rpt%temperature(i)%qc .lt. 2**15) temp(j) = rpt%temperature(i)%data
                 if (rpt%dewpoint(i)%qc .lt. 2**15) dwpt(j) = rpt%dewpoint(i)%data
                 if (rpt%speed(i)%qc .lt. 2**15) spd(j) = rpt%speed(i)%data
                 if (rpt%direction(i)%qc .lt. 2**15) dir(j) = rpt%direction(i)%data
              endif
           enddo ILOOP
           if (j < 2) DO_THE_PLOT = .FALSE.
           if (DO_THE_PLOT) then
              call skewt(pres, temp, dwpt, spd, dir, j, &
                TRIM(rpt%id_string),&
                rpt%hdate, rpt%name_string, &
                xloc, yloc,&
                rpt%xlat, rpt%xlon, ilw, ictl, icdl, lcf, woffs, badval)
              call frame
           endif
        endif

        call read_r_report(iunit, ierr, rpt)
     enddo

     call skewt_clsgks

     CALL geth_newdate(rdate, trim(start_date), itimes*interval)

  END DO LOOP_TIMES

  print*," All done"

end program plot_raobs

