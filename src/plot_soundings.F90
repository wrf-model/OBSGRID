program plot_raobs

  use report_module
  use skewt_module
  use mapinfo_module

  logical :: lcf

  type(report) :: rpt

  integer :: iunit = 10
  integer :: aunit = 99

  character(len=120) :: flnm, oafile

  real, dimension(200) :: pres, temp, dwpt, spd, dir
  real, allocatable, dimension(:) :: plvls

  real, dimension(8) :: xarr

  logical :: rfile = .FALSE.
  
  integer, dimension(50,20) :: bhi
  real, dimension(20,20) :: bhr
  integer :: ndim
  integer, dimension(4) :: stdim, endim
  real :: xtime
  character(len=4) :: stagger, order
  character(len=24) :: date
  character(len=9) :: name
  character(len=25) :: units
  character(len=46) :: desc
  logical :: DO_THE_PLOT


  open ( unit=aunit , file='LITTLE_R_DOMAIN1' , status='old', &
         form='unformatted' , access='sequential' )

  call skewt_opngks("soundings.cgm")

  numarg = 2

     rfile = .TRUE.

     do
        read(aunit, iostat=ierr) iflag
        if (ierr /= 0) exit
        if (iflag == 0) then
           read(aunit) bhi, bhr
           call get_mapinfo(bhi, bhr)
        else if (iflag == 1) then
           read(aunit) ndim, stdim, endim, xtime, stagger, order, date, &
                name, units, desc
           if (name == "PRESSURE") then
              allocate(plvls(endim(1)))
              read(aunit) plvls
              exit
           else
              read(aunit) xdum
           endif
        else if (iflag == 2) then
           
        else
           stop "ERROR:  PROBLEM READING IFLAG FROM ANALYSIS FILE."
        endif
     enddo
     
     close(aunit)

     if (.not. allocated(plvls)) then
        print*, 'PRESSURE not found.'
        rfile = .FALSE.
     endif

  print *,'Enter raobs input file name (such as "useful_out_1993-03-13_00:00:00.0000"):'
  read(5,fmt='(A)') flnm
  flnm=trim(flnm)
  print *,'flnm=',flnm

  open(iunit, file=flnm, form='formatted', status='old', action='read')

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

  print *,'plots are called soundings.cgm'

end program plot_raobs

subroutine print_usage(flnm)
  implicit none
  character(len=*) :: flnm
  write(*,'(/,"      Usage:  ", A, " obs_file [analysis_file]",/)') trim(flnm)
  write(*,'(7x,"where ''obs_file'' is the name of a little_r-formatted observations file.")')
  write(*,'(/,13x,"''analysis_file'' (optional) is the name of a little_r analysis ")')
  write(*,'(13x,"file (e.g., LITTLE_R_DOMAIN1).",//)')
end subroutine print_usage
