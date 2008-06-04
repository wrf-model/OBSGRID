module mapinfo_module
! Module to deal with map projections 


! type MAPINFO_TYPE is a defined type that stores all sorts of the 
! details about the map projection.
  type mapinfo_type
     character(len=2) :: project
     real :: dskm, phic, xlonc
     integer :: ix, jx, ixmoad, jxmoad
     real :: xsouth, xwest, xratio, truelat1, truelat2
  end type mapinfo_type

  type(mapinfo_type) :: mapinfo

contains

  subroutine get_mapinfo( met_ncid )

    INCLUDE 'netcdf.inc'

    real :: rdummy
    integer :: i, idummy
    integer :: met_ncid
    integer :: rcode
    integer :: ndims, nvars, ngatts, nunlimdimid
    integer :: wedim, sndim
    integer :: dim_val
    CHARACTER (LEN=31) :: dim_name

    character(LEN=2), parameter, dimension(3) :: NPROJ = (/'LC','ST','ME'/)

    rcode = nf_inq(met_ncid, ndims, nvars, ngatts, nunlimdimid)
    dims_loop : DO i = 1, ndims
       rcode = nf_inq_dim(met_ncid, i, dim_name, dim_val)
       IF ( dim_name == 'west_east_stag'    ) wedim  = dim_val
       IF ( dim_name == 'south_north_stag'  ) sndim  = dim_val
    ENDDO dims_loop

    rcode = NF_GET_ATT_REAL(met_ncid, nf_global, "CEN_LAT", rdummy )
            mapinfo%phic = rdummy
    rcode = NF_GET_ATT_REAL(met_ncid, nf_global, "STAND_LON", rdummy )
            mapinfo%xlonc = rdummy
    rcode = NF_GET_ATT_REAL(met_ncid, nf_global, "TRUELAT1", rdummy )
            mapinfo%truelat1 = rdummy   
    rcode = NF_GET_ATT_REAL(met_ncid, nf_global, "TRUELAT2", rdummy )
            mapinfo%truelat2 = rdummy   
    rcode = NF_GET_ATT_REAL(met_ncid, nf_global, "DX", rdummy )
            mapinfo%dskm = rdummy*0.001
        
    !! HARDCODE - we are only working with metoa files - 
    !!            so no need to keep track of expanded domain
    mapinfo%ix = sndim
    mapinfo%jx = wedim
    mapinfo%ixmoad = sndim
    mapinfo%jxmoad = wedim

    !! Expanded domains
    !! (13,1) is domain
    !! (8,1) is expanded
    !!if ((bhi(13) == 1).and.(bhi(8) == 1)) then
       !!mapinfo%ix = bhi(9)
       !!mapinfo%jx = bhi(10)
       !!mapinfo%ixmoad = bhi(9)
       !!mapinfo%jxmoad = bhi(10)
    !!endif

    rcode = NF_GET_ATT_INT(met_ncid, nf_global, "j_parent_start", idummy )
            mapinfo%xsouth = real(idummy)
    rcode = NF_GET_ATT_INT(met_ncid, nf_global, "i_parent_start", idummy )
            mapinfo%xwest = real(idummy)
    rcode = NF_GET_ATT_INT(met_ncid, nf_global, "parent_grid_ratio", idummy )
            mapinfo%xratio = float(idummy)
    rcode = NF_GET_ATT_INT(met_ncid, nf_global, "MAP_PROJ", idummy )
            mapinfo%project = NPROJ(idummy)

  end subroutine get_mapinfo

  subroutine lltoxy_mapinfo(xlat, xlon, x, y)
! Given the latitude and longitude, computes the X and Y in the model grid.
! If the MAPINFO structure has been filled, we can use all the information
! stored there and not have to pass everything through the argument list.
! If the MAPINFO structure has not been filled, we're in trouble.  So 
! get_mapinfo must be called sometime before lltoxy_mapinfo is called.

!*****************************************************************************!
!  Notes    - Modeled after (i.e., swiped from) XYTOLL in the plots.o library !
!*****************************************************************************!
    implicit none
    real, intent(in)  ::    xlat    ! latitude of the point of interest.
    real, intent(in)  ::    xlon    ! longitude of the point of interest.
    real, intent(out) ::    x       ! x location of the given (lat,lon) point
    real, intent(out) ::    y       ! y location of the given (lat,lon) point

!  Parameters: ---------------------------------------------------------------!

    real, parameter :: pi = 3.14159265   ! you know!  pi = 180 degrees
!   real, parameter :: re = 6370.949     ! the radius of the earth in km
    real, parameter :: re = 6370.        ! the radius of the earth in km consistent with rest of MM5 system
    real, parameter :: ce = 2.*pi*re     ! the circumference of the earth in km
    real, parameter :: degrad = pi/180.

!  Integer variables: --------------------------------------------------------!

    integer   :: imax              ! maximum vertical grid point    (local)
    integer   :: jmax              ! maximum horizontal grid point  (local)

!  Real variables: -----------------------------------------------------------!

    real          :: flat1
    real          :: confac            ! cone factor
    real          :: rcln              ! center longitude in radians  (local)
    real          :: rclt              ! center latitude in radians   (local)
    real          :: cj                ! center x coord. for grid     (local)
    real          :: ci                ! center y coord. for grid     (local)
    real          :: dj                ! distance from the central
    !  meridian to the point        (local)
    real          :: di                ! distance from pole to point  (local)
    real          :: bm                ! calculation variable         (local)

    real :: rlat, rlon, djovrdi, djsdis

!****************************  subroutine begin  *****************************C

    imax = NINT((float(mapinfo%ixmoad)-(2.*mapinfo%xsouth)+1.) * mapinfo%xratio)+1
    jmax = NINT((float(mapinfo%jxmoad)-(2.*mapinfo%xwest )+1.) * mapinfo%xratio)+1

    rlat =  xlat * degrad
    rlon =  xlon * degrad
    rclt = mapinfo%phic * degrad
    rcln = mapinfo%xlonc * degrad
    flat1 = mapinfo%truelat1 * degrad
    cj = float(jmax + 1) * 0.5
    ci = float(imax + 1) * 0.5

    if (mapinfo%project(1:2) .eq. 'ME') then
       di = re * log(tan (0.5*(rlat + pi * 0.5)))
       dj = re *(rlon - rcln)
       y = ci +(di + re * log(cos(rclt)/(1 + sin(rclt))))/mapinfo%dskm
    else if (mapinfo%project(1:2) .eq. 'CE') then
       di = (ce*0.5)*(rlat-rclt)/pi
       dj = ce * (rlon-rcln)/(pi*2.)
       y = ci + di/mapinfo%dskm
    else if (mapinfo%project(1:2) .eq. 'LC') then
       if ((xlat == -90.) .and. (mapinfo%phic > 0)) then
          x = -1.E25
          y = -1.E25
          return
       endif
       call lccone_mapinfo(confac)
       if (mapinfo%phic .ge. 0.0) then
          bm = tan(-0.5*(rlat - pi * 0.5))
          djovrdi = -tan((rlon - rcln)*confac)
          djsdis = ((bm / tan((pi * 0.5 - abs(flat1))*0.5)) &
               **confac * sin(pi * 0.5 - abs(flat1))/(confac/re))**2
          di = -sqrt(djsdis/(1+djovrdi*djovrdi))
          dj = di * djovrdi
          y = ci + &
               (di + re/confac * sin(pi*0.5 - flat1)* &
               (tan((pi * 0.5 - rclt) * 0.5) / &
               tan((pi * 0.5 - flat1) * 0.5))**confac) &
               /mapinfo%dskm
!print*,"map info ", y
!print*,ci, di, re, confac, rclt, flat1, mapinfo%dskm
       else
          bm = tan( 0.5*(rlat + pi * 0.5))
          djovrdi = tan((rlon - rcln)*confac)
          djsdis = ((bm / tan((pi * 0.5 - abs(flat1))*0.5)) &
               **confac * sin(pi * 0.5 - abs(flat1))/(confac/re))**2
          di = sqrt(djsdis/(1+djovrdi*djovrdi))
          dj = di * djovrdi
          y = ci + (di + re/confac * &
               sin(-pi * 0.5 - flat1) * &
               (tan((-pi * 0.5 - rclt) * 0.5) / &
               tan((-pi*0.5 - flat1) * 0.5))**confac) &
               /mapinfo%dskm
       end if
    else if (mapinfo%project(1:2) .eq. 'ST') then
       if (mapinfo%phic .ge. 0.0) then
          djovrdi = - tan(rlon - rcln)
          bm = tan(-0.5*(rlat - pi * 0.5))
          djsdis = (re * bm * &
               (1.0 + cos( pi * 0.5 - flat1)))**2
          di = -sign(sqrt(djsdis/(1+djovrdi**2)),cos(rlon-rcln))
          dj = di * djovrdi
          y =  ci + (di + re * sin(pi * 0.5 - rclt) * &
               (1.0 + cos(pi * 0.5 - flat1)) / &
               (1.0 + cos(pi * 0.5 - rclt)) )/mapinfo%dskm
       else
          djovrdi =   tan(rlon - rcln)
          bm = tan( 0.5*(rlat + pi * 0.5))
          djsdis = (re * bm * &
               (1.0 + cos(-pi * 0.5 - flat1)))**2
          di = sign(sqrt(djsdis/(1+djovrdi**2)),cos(rlon-rcln))
          dj = di * djovrdi
          y = ci + (di + re * sin(-pi * 0.5 - rclt) * &
               (1.0 + cos(-pi * 0.5 - flat1)) / &
               (1.0 + cos(-pi * 0.5 - rclt)))/mapinfo%dskm
       end if
    endif
    x =  cj + dj/mapinfo%dskm

!*****************************  Subroutine End  ******************************C

  end subroutine lltoxy_mapinfo

  subroutine lccone_mapinfo(confac)
! Computes the cone factor for Lambert Conformal projections.  This routine
! uses information from the MAPINFO structure.
    real, parameter :: conv=0.01745329251994
    integer :: msign
    real :: confac
    msign = int(sign(1.0, mapinfo%phic))
    if (abs(mapinfo%truelat1-mapinfo%truelat2).lt.1.E-2) then
       confac = sin(mapinfo%truelat1*conv)
    else
       confac = log10(cos(mapinfo%truelat1*conv))-log10(cos(mapinfo%truelat2*conv))
       confac = confac/(log10(tan((45.-float(msign)*mapinfo%truelat1/2.)*conv))-&
            log10(tan((45.-float(msign)*mapinfo%truelat2/2.)*conv)))
    endif

  end subroutine lccone_mapinfo

end module mapinfo_module
