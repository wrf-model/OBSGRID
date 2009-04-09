module report_module
! Module to deal with the observations in the little_r observations format.

  type field
     real :: data
     integer :: qc
  end type field

  type report

     real :: xlat, xlon
     character(len=40) :: id_string
     character(len=40) :: name_string

     character(len=40) :: platform_string
     character(len=40) :: source_string

     real :: elev

     integer :: num_vld_fld, num_error, num_warning
     integer :: seq_num, num_dups

     logical :: is_sound, bogus, discard

     integer :: sut, julian
     character(len=20) :: date_char
     character(len=19) :: hdate

     type(field) :: slp, ref_pres, ground_t, sst, psfc, precip,&
          t_max, t_min, t_min_night, ptend_3, ptend_24, cloud_covr, ceiling

     type(field), dimension(200) :: pressure, height, temperature, dewpoint,&
          speed, direction, u, v, rh, thickness

     integer :: end_num_vld_fld, end_num_error, end_num_warning

     integer :: num_lvls

  end type report

  character(len=84), parameter :: rpt_format =&
       '  ( 2f20.5 , 2a40 , '  &
       // ' 2a40 , 1f20.5 , 5i10 , 3L10 , '  &
       // ' 2i10 , a20 ,  13( f13.5 , i7 ) ) '

  character(len=22), parameter :: meas_format =&
       ' ( 10( f13.5 , i7 ) ) '
  character(len=14), parameter :: end_format =&
       ' ( 3 ( i7 ) ) ' 

contains

  subroutine read_r_report(iunit, ierr, rpt)
! Read a full station report in the little_r observations format.
! Returns structure rpt of type(report), which holds the whole report.
    type(report) :: rpt

    type(field) :: slp, ref_pres, ground_t, sst, psfc, precip,&
         t_max, t_min, t_min_night, ptend_3, ptend_24, cloud_covr, ceiling

    type(field) :: pressure, height, temperature, dewpoint,&
         speed, direction, u, v, rh, thickness

    integer :: iunit
    integer :: ierr

    character(len=40) :: string1, string2, string3, string4

    real :: xlat, xlon, elev
    integer :: julian, sut, num_error, num_warning, num_dups, kval, seq_num
    integer :: end_num_vld_fld, end_num_error, end_num_warning
    logical :: is_sound, bogus, discard
    character(len=20) :: date_char
    integer :: k

    ierr = 0

    read(iunit, FMT=rpt_format, end=1200, err = 1300)&
         xlat, xlon, string1, string2,&
         string3, string4, elev, kval, num_error, num_warning,&
         seq_num, num_dups, is_sound, bogus, discard, &
         sut, julian, date_char,&
         slp, ref_pres, &
         ground_t, sst, psfc, precip, t_max, t_min, t_min_night,&
         ptend_3, ptend_24, cloud_covr, ceiling

    rpt%xlat = xlat
    rpt%xlon = xlon
    rpt%id_string = string1
    rpt%name_string = string2
    rpt%platform_string = string3
    rpt%source_string = string4
    rpt%elev = elev
    rpt%num_vld_fld = kval
    rpt%num_error = num_error
    rpt%num_warning = num_warning
    rpt%seq_num = seq_num
    rpt%num_dups = num_dups
    rpt%is_sound = is_sound
    rpt%bogus = bogus
    rpt%discard = discard
    rpt%sut = sut
    rpt%julian = julian
    rpt%date_char = date_char
    rpt%slp = slp
    rpt%ref_pres = ref_pres
    rpt%ground_t = ground_t
    rpt%sst = sst
    rpt%psfc = psfc
    rpt%precip = precip
    rpt%t_max = t_max
    rpt%t_min = t_min
    rpt%t_min_night = t_min_night
    rpt%ptend_3 = ptend_3
    rpt%ptend_24 = ptend_24
    rpt%cloud_covr = cloud_covr
    rpt%ceiling = ceiling

    date_char = trim(rpt%date_char)
    rpt%hdate = date_char(7:10)//'-'//date_char(11:12)//'-'//date_char(13:14)//&
         '+'//date_char(15:16)//':'//date_char(17:18)//':'//date_char(19:20)


    k = 0
    LVLOOP : do
       read(iunit, FMT=meas_format)&
            pressure, height, temperature, dewpoint, speed, direction, u, v, rh, thickness

       k = k + 1
       rpt%pressure(k) = pressure
       rpt%height(k) = height
       rpt%temperature(k) = temperature
       rpt%dewpoint(k) = dewpoint
       rpt%speed(k) = speed
       rpt%direction(k) = direction
       rpt%u(k) = u
       rpt%v(k) = v
       rpt%rh(k) = rh
       rpt%thickness(k) = thickness

       if (pressure%data == -777777.) exit 

    enddo LVLOOP

    rpt%num_lvls = k

    read(iunit, fmt=end_format) end_num_vld_fld, end_num_error, end_num_warning

    rpt%end_num_vld_fld = end_num_vld_fld
    rpt%end_num_error = end_num_error
    rpt%end_num_warning = end_num_warning
    return

1200 continue
    ierr = 1
    return

1300 continue
    ierr = 2
    return

  end subroutine read_r_report

end module report_module
