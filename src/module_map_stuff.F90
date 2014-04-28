MODULE map_stuff

CONTAINS

!-------------------------------------------
   SUBROUTINE plotmap ( met_ncid, grid_id )

      INCLUDE 'netcdf.inc'

      REAL :: rdummy
      INTEGER :: idummy, rcode
      INTEGER :: met_ncid, grid_id
      INTEGER :: ndims, nvars, ngatts, nunlimdimid
      INTEGER :: wedim, sndim, map_proj
      INTEGER :: dim_val
      CHARACTER (LEN=31) :: dim_name

      CHARACTER (LEN=2) :: project
      CHARACTER (LEN=2) , DIMENSION(3) , PARAMETER :: nproj = (/ 'LC','ST','ME' /)

      REAL :: dx, ds, phic, xlonc, stdlon, xsouth, xwest, xratio, truelat1, truelat2 
      REAL :: pl1 , pl2 , pl3 , pl4 , polat , rot
      REAL :: xleft , xright
      REAL :: vl , vr , vb , vt , wl , wr , wb , wt
      INTEGER :: ixmoad , jxmoad , kproj , ix , jx , jlts , jgrid , idot , iusout , ier , ls
      integer :: nf_tstart2(3), nf_tcount2(3), i
      real, allocatable, dimension(:,:) :: scr2d

      rcode = nf_inq(met_ncid, ndims, nvars, ngatts, nunlimdimid)
      dims_loop : DO i = 1, ndims
         rcode = nf_inq_dim(met_ncid, i, dim_name, dim_val)
         IF ( dim_name == 'west_east_stag'    ) wedim  = dim_val
         IF ( dim_name == 'south_north_stag'  ) sndim  = dim_val
      ENDDO dims_loop

      rcode = NF_GET_ATT_INT(met_ncid, nf_global, "MAP_PROJ", idummy )
      project = nproj(idummy)
      map_proj = idummy

      !  Grid distance in km, center lat/lon.
      rcode = NF_GET_ATT_REAL(met_ncid, nf_global, "DX", dx ) 
      ds = dx/1000.
      rcode = NF_GET_ATT_REAL(met_ncid, nf_global, "CEN_LAT", phic )
      rcode = NF_GET_ATT_REAL(met_ncid, nf_global, "CEN_LON", xlonc )
      rcode = NF_GET_ATT_REAL(met_ncid, nf_global, "STAND_LON", stdlon )

      !  The i and j dimensions, in Mother Of All Domain space.
      ixmoad = sndim
      jxmoad = wedim

      !  The starting locations of this domain, wrt to MOAD.
      rcode = NF_GET_ATT_INT(met_ncid, nf_global, "j_parent_start", idummy )
      xsouth = real(idummy)
      rcode = NF_GET_ATT_INT(met_ncid, nf_global, "i_parent_start", idummy )
      xwest = real(idummy)

      if (grid_id .gt. 1) then
      allocate (scr2d(wedim-1,sndim-1))
      rcode = nf_inq_varid (met_ncid, 'XLAT', idummy)
      if (rcode .ne. 0) then 
        rcode = nf_inq_varid (met_ncid, 'XLAT_M', idummy)
      write(0,*) 'error = ',nf_strerror(rcode)
      endif
      nf_tstart2(1)=1
      nf_tstart2(2)=1
      nf_tstart2(3)=1
      nf_tcount2(1)=wedim-1
      nf_tcount2(2)=sndim-1
      nf_tcount2(3)=1
      rcode = nf_get_vara_real (met_ncid, idummy, nf_tstart2,nf_tcount2, scr2d)
      write(0,*) 'error = ',nf_strerror(rcode)
      pl1 = scr2d(1,1)
      pl3 = scr2d(wedim-1,sndim-1)
      rcode = nf_inq_varid (met_ncid, 'XLONG', idummy)
      if (rcode .ne. 0) then 
        rcode = nf_inq_varid (met_ncid, 'XLONG_M', idummy)
      write(0,*) 'error = ',nf_strerror(rcode)
      endif 
      rcode = nf_get_vara_real (met_ncid, idummy, nf_tstart2,nf_tcount2, scr2d)
      pl2 = scr2d(1,1)
      pl4 = scr2d(wedim-1,sndim-1)
      endif

      !  The ratio of the grid distance of the mother domain to the grid distance
      !  from this domain.
      rcode = NF_GET_ATT_INT(met_ncid, nf_global, "parent_grid_ratio", idummy )
      xratio = REAL(idummy)

      !  The latitude where the projection is true (map factors = 1).
      rcode = NF_GET_ATT_REAL(met_ncid, nf_global, "TRUELAT1", truelat1 )
      rcode = NF_GET_ATT_REAL(met_ncid, nf_global, "TRUELAT2", truelat2 )

      !  Get the lat/lon of the lower left hand corner.

      write(0,*) 'ds = ',ds,' phic = ',phic,' xlonc = ',xlonc
      write(0,*) 'truelat1 = ',truelat1,' truelat2 = ',truelat2
      ix = sndim
      jx = wedim
      if (grid_id .eq. 1) then
      ! CALL xytoll(1.,1.,pl1,pl2,project,ds,phic,xlonc,ixmoad,jxmoad,xsouth,xwest,xratio,truelat1,truelat2)
        CALL ij_to_ll( map_proj, truelat1, truelat2, stdlon, phic, xlonc, 0.0, 0.0,  &
                     wedim/2.0, sndim/2.0, dx, dx, 0.0, 0.0, 1.0, 1.0, pl1, pl2 )

      !  Get the lat/lon of the upper right corner point.

      ! CALL xytoll(REAL(jx),REAL(ix),pl3,pl4,project,ds,phic,xlonc,ixmoad,jxmoad,xsouth,xwest,xratio,truelat1,truelat2)
        CALL ij_to_ll( map_proj, truelat1, truelat2, stdlon, phic, xlonc, 0.0, 0.0,  &
                     wedim/2.0, sndim/2.0, dx, dx, 0.0, 0.0, real(wedim), real(sndim), pl3, pl4 )
      endif

      write(0,*) 'pl1 = ',pl1,' pl2 = ',pl2
      write(0,*) 'xratio = ',xratio,' xsouth = ',xsouth,' xwest = ',xwest
      write(0,*) 'ixmoad = ',ixmoad,' jxmoad = ',jxmoad
      write(0,*) 'jx = ',jx,' ix = ',ix,' project = ',project
      
      ! kproj:  1 = stereographic, 3 = lambert conformal, 9 = mercator
      !  polat, lat of pole?
    
      !  We need to go from MM5 info to NCAR Graphics map definitions.

      IF (project.EQ.'ST') THEN
         kproj = 1
         polat =  90.
         IF (phic .LT. 0.) polat = -90.
         rot = 0.
      ELSE IF (project.EQ.'LC') THEN
         kproj = 3
         rot = truelat1
         polat = truelat2
      ELSEIF (project.EQ.'ME') THEN
         kproj = 9
         polat = 0.
         rot = 0.
      ELSE
         PRINT '(A,I8,A)','Projection is weird, projection # = ',project,'.'
         STOP 'plotmap'
      END IF
      !     jlts = 2:
      !         pl1 = latitude of lower-left corner point
      !         pl2 = longitude of lower-left corner point
      !         pl3 = latitude of upper-right corner point
      !         pl4 = longitude of upper-right corner point
      
      jlts = 2

      IF (pl2 .LT. 0) THEN
         xleft = 360 + pl2
      ELSE
         xleft = pl2
      END IF
      IF (pl4 .LT. 0) THEN
         xright = 360 + pl4
      ELSE
         xright = pl4
      END IF
      
      !  How far apart are the lat/lon lines.

      IF (xright - xleft .GT. 50.) THEN
         jgrid = 10
      ELSE IF (xright - xleft .GT. 30.) THEN
         jgrid = 5            ! every 5 deg
      ELSE
         jgrid = 5            ! every 5 deg
      END IF

      !  Is the map to be dots or a line.

      idot = 0

      !  Would you like to see the US state boundaries?

      iusout = 1

      !  Make the map within the middle 80% of the domain.

      CALL mappos(.1,.9,.1,.9)

      !  This is the "generate a map" call.

      write(0,*) 'kproj = ',kproj,' polat = ',polat,' xlonc = ',xlonc
      write(0,*) 'rot = ',rot,' pl1 = ',pl1,' pl2 = ',pl2,' pl3 = ',pl3,' pl4 = ',pl4
      write(0,*) 'jlts = ',jlts,' jgrid = ',jgrid

      CALL supmap(kproj,polat,stdlon,rot,pl1,pl2,pl3,pl4,jlts,jgrid,iusout,idot,ier)

      !  Was everything OK?

      IF (ier.NE.0) THEN
         PRINT '(A,I8,A)','We wanted to have the map return code be 0, it was ',ier,'.  Oops.'
         STOP 'MAP_hosed'
      END IF

      !  We are going to remember the set call for the station locations.

      CALL getset(vl,vr,vb,vt,wl,wr,wb,wt,ls)
      CALL set(vl,vr,vb,vt,1.,REAL(jx), 1., REAL(ix), ls)

   END SUBROUTINE plotmap

!-----------------------------------------------------------

   SUBROUTINE lltoxy (xlat,xlon,x,y,&
        project,ds,phic,xlonc,ixmoad,jxmoad,xsouth,xwest,&
        xratio,truelat1,truelat2)
   
   !*****************************************************************************C
   !  Notes    - Modeled after XYTOLL in the plots.o library                     C
   !*****************************************************************************C
     implicit none
   
   !  Parameters
   
     REAL, parameter :: pi = 3.14159265   ! you know!  pi = 180 degrees
!    REAL, parameter :: re = 6370.949     ! the radius of the earth in km
     REAL, parameter :: re = 6370.        ! the radius of the earth in km consistent with rest of MM5 system
     REAL, parameter :: ce = 2.*pi*re     ! the circumference of the earth in km
   
   !  Character variables
   
     character(LEN=2) :: project          ! projection indicator    (in)
   
   !  Integer variables
   
     INTEGER   :: imax              ! maximum vertical gird point
     INTEGER   :: jmax              ! maximum horizontal grid point
     INTEGER   :: isign             ! N (1) or S (-1) hemisphere phic
   
   !  Real variables
   
     REAL          :: x                 ! x coord. to be changed   (in)
     REAL          :: y                 ! y coord. to be changed   (in)
     REAL          :: xlat              ! resulting latitude     (out)
     REAL          :: xlon              ! resulting longitude    (out)
     REAL          :: ds                ! grid distance in km
     REAL          :: phic              ! center latitude
     REAL          :: xlonc             ! center longitude
     REAL          :: truelat1
     REAL          :: truelat2
     REAL          :: confac            ! cone factor
     REAL          :: rcln              ! center longitude in radians  (local)
     REAL          :: rclt              ! center latitude in radians   (local)
     REAL          :: cj                ! center x coord. for grid     (local)
     REAL          :: ci                ! center y coord. for grid     (local)
     REAL          :: dj                ! distance from the central
                                        !  meridian to the point        (local)
     REAL          :: di                ! distance from pole to point  (local)
     REAL          :: bm                ! calculation variable         (local)
   
     INTEGER :: ixmoad, jxmoad
     REAL :: xsouth, xratio, rlat, rlon, djovrdi, djsdis, xwest
   
   
   !****************************  SUBROUTINE begin  *****************************C
   
     imax = NINT((REAL(ixmoad)-(2.*xsouth)+1.) * xratio)+1
     jmax = NINT((REAL(jxmoad)-(2.*xwest )+1.) * xratio)+1
   
     rlat =  xlat * pi / 180.0
     rlon =  xlon * pi / 180.0
     rclt = phic * pi / 180.0
     rcln = xlonc * pi / 180.0
     cj = REAL(jmax + 1) * 0.5
     ci = REAL(imax + 1) * 0.5
   
     IF (project(1:2) .EQ. 'ME') THEN
        di = re * log(tan ((rlat + pi * 0.5)/(2.0)))
        dj = re *(rlon - rcln)
        y = ci +(di + re * log(COS(rclt)/(1 + SIN(rclt))))/ds
     ELSE IF (project(1:2) .EQ. 'CE') THEN
        di = (ce*0.5)*(rlat-rclt)/pi
        dj = ce * (rlon-rcln)/(pi*2.)
        y = ci + di/ds
     ELSE IF (project(1:2) .EQ. 'LC') THEN
        IF ((xlat == -90.) .and. (phic > 0)) THEN
           x = -1.E25
           y = -1.E25
           return
        END IF
        IF (phic.GE.0) isign = 1
        IF (phic.LT.0) isign = -1
        call lccone(truelat1,truelat2,isign,confac)
        IF (phic .GE. 0.0) THEN
           bm = TAN(-(rlat - pi * 0.5) /  2.0)
        ELSE
           bm = TAN((rlat + pi * 0.5) /  2.0)
        END IF
        IF (phic .GE. 0.0) THEN
   !           (dj/di)
           djovrdi = -TAN((rlon - rcln)*confac)
        ELSE
           djovrdi = TAN((rlon - rcln)*confac)
        END IF
   !        (dj**2 + di**2)
        djsdis = ((bm / TAN((pi * 0.5 - ABS(truelat1 * pi/180.0))/2.0)) &
             **confac * &
             SIN(pi * 0.5 - ABS(truelat1 * pi/180.0))/(confac/re))**2
        IF (phic.GE.0.0) THEN
           di = -SQRT(djsdis/(1+djovrdi*djovrdi))
        ELSE
           di = SQRT(djsdis/(1+djovrdi*djovrdi))
        END IF
        dj = di * djovrdi
        IF (phic .GE. 0.0) THEN
           y = ci + &
                (di + re/confac * SIN(pi*0.5 - (truelat1 * pi/180.0))* &
                (TAN((pi * 0.5 - rclt) * 0.5) / &
                TAN((pi * 0.5 - (truelat1 * pi/180.0)) * 0.5))**confac) &
                /ds
        ELSE
           y = ci + (di + re/confac * &
                SIN(-pi * 0.5 - (truelat1 * pi/180.0)) * &
                (TAN((-pi * 0.5 - rclt) * 0.5) / &
                TAN((-pi*0.5 - (truelat1 * pi/180.0)) * 0.5))**confac) &
                /ds
        END IF
     ELSE IF (project(1:2) .EQ. 'ST') THEN
        IF (phic .GE. 0.0) THEN
   !           (dj/di)
           djovrdi = - TAN(rlon - rcln)
        ELSE
           djovrdi =   TAN(rlon - rcln)
        END IF
        IF (phic .GE. 0.0) THEN
           bm = TAN((rlat - pi * 0.5)/(-2.0))
        ELSE
           bm = TAN((rlat + pi * 0.5)/( 2.0))
        END IF
        IF (phic .GE. 0.0) THEN
   !           dj**2 + di**2
           djsdis = (re * bm * &
                (1.0 + COS( pi * 0.5 - (truelat1 * pi/180.0))))**2
           di = -SQRT(djsdis/(1+djovrdi*djovrdi))
        ELSE
   !           dj**2 + di**2
           djsdis = (re * bm * &
                (1.0 + COS(-pi * 0.5 - (truelat1 * pi/180.0))))**2
           di = SQRT(djsdis/(1+djovrdi*djovrdi))
        END IF
        dj = di * djovrdi
        IF (phic .GT. 0.0) THEN
           y =  ci + (di + re * SIN(pi * 0.5 - rclt) * &
                (1.0 + COS(pi * 0.5 - (truelat1 * pi/180.0))) / &
                (1.0 + COS(pi * 0.5 - rclt)) )/ds
        ELSE
           y = ci + (di + re * SIN(-pi * 0.5 - rclt) * &
                (1.0 + COS(-pi * 0.5 - (truelat1 * pi/180.0))) / &
                (1.0 + COS(-pi * 0.5 - rclt)))/ds
        END IF
     END IF
     x =  cj + dj/ds
   
   END SUBROUTINE lltoxy

!-----------------------------------------------------------
   
   SUBROUTINE xytoll (x,y,xlat,xlon,&
        project,ds,phic,xlonc,ixmoad,jxmoad,xsouth,xwest,&
        xratio,truelat1,truelat2)
   
   !*****************************************************************************C
   !  xytoll   - This is a MAPDRV routine                                        C
   !  Section  - Labels                                                          C
   !  Purpose  - To  transform  mesoscale gird point coordinates into  latitude, C
   !             longitude coordinates.                                          C
   !                                                                             C
   !  On entry - X  and  Y are an ordered pair representing a grid point in  the C
   !             mesoscale grid.  XYLLON is a common block that contains the in- C
   !             formation necessary for describing the domain.                  C
   !                                                                             C
   !  On exit  - XLAT, XLON contain  the latitude and longitude respectively     C
   !             that resulted from the transformation.                          C
   !                                                                             C
   !  Assume   - Nothing.                                                        C
   !                                                                             C
   !  Notes    - The formula's used in this routine were taken from the  PROGRAM C
   !             TERRAIN DOCUMENTATION AND USER'S GUIDE.                         C
   !                                                                             C
   !  Author   - Jeremy Asbill   Date - September 17, 1990      for the MM4 club C
   !*****************************************************************************C
     implicit none
   
   !  Parameters
   
     REAL, parameter :: pi = 3.14159265   ! you know!  pi = 180 degrees
!    REAL, parameter :: re = 6370.949     ! the radius of the earth in km
     REAL, parameter :: re = 6370.        ! the radius of the earth in km consistent with rest of MM5 system
     REAL, parameter :: ce = 2.*pi*re     ! the circumference of the earth in km
   
   !  Character variables
   
     character(LEN=2) :: project           ! projection indicator            (in)
   
   !  Integer variables
   
     INTEGER ::       imax             ! gridpoints in i(y) direction    (in)
     INTEGER ::       jmax             ! gridpoints in j(x) direction    (in)
     INTEGER ::       isign            ! N (1) or S (-1) hemisphere phic
   
   !  Real variables
   
     REAL ::            x                ! x coord. to be changed          (in)
     REAL ::            y                ! y coord. to be changed          (in)
     REAL ::            xlat             ! resulting latitude             (out)
     REAL ::            xlon             ! resulting longitude            (out)
     REAL ::            ds                ! grid distance in km
     REAL ::            phic              ! center latitude
     REAL ::            xlonc             ! center longitude
     REAL ::            truelat1
     REAL ::            truelat2
     REAL ::            confac           ! cone factor
     REAL ::            rcln             ! center longitude in radians  (local)
     REAL ::            rclt             ! center latitude in radians   (local)
     REAL ::            cj               ! center x coord. for grid     (local)
     REAL ::            ci               ! center y coord. for grid     (local)
     REAL ::            dj               ! distance from the central
   !                                          meridian to the point        (local)
     REAL ::            di               ! distance from pole to point  (local)
     REAL ::            bm               ! calculation variable         (local)
   
   
     INTEGER :: ixmoad, jxmoad
     REAL :: xsouth, xratio, xwest
   !****************************  SUBROUTINE begin  *****************************C
   
     imax = NINT((REAL(ixmoad)-(2.*xsouth)+1.) * xratio)+1
     jmax = NINT((REAL(jxmoad)-(2.*xwest )+1.) * xratio)+1
   
   !  Convert the center latitude and longitude of the domain to radians
   
     rclt = phic * pi/180.0
     rcln = xlonc * pi/180.0
   
   !  Find the center values of the grid in mesoscale grid coordinates
   
     cj = REAL(jmax + 1) * 0.5
     ci = REAL(imax + 1) * 0.5
   
   !  Calculate the distance from the vertical axis to (J,I)
   
     dj = (x - cj) * ds
   
   !  The rest is figured out differently for each type of projection, so ...
   !  If the projection is mercator ('ME') THEN ...
   
     IF (project(1:2) .EQ. 'ME') THEN
   
   !  Calculate the distance the point in question is from the pole
   
        di = -re * log(COS(rclt)/(1 + SIN(rclt))) + &
             (y - ci) * ds
   
   !  Calculate the latitude desired in radians
   
        xlat = 2.0 * aTAN(exp(di/re)) - pi * 0.5
   
   !  Calculate the longitude desired in radians
   
        xlon = rcln + dj/re
   
   !  If the projection is cylindrical equidistant ('CE') THEN ...
   
     ELSE IF (project(1:2) .EQ. 'CE') THEN
   
   !  Calculate the distance from the horizontal axis to (J,I)
   
        di = (y - ci) * ds
   
   !  Determine the shift north-south
   
        xlat = rclt + (pi * di/(ce * 0.5))
   
   !  Determine the shift east-west
   
        xlon = rcln + (2 * pi * dj/ce)
   
   !  If the projection is lambert conic conformal ('LC') THEN ...
   
     ELSE IF (project(1:2) .EQ. 'LC') THEN
        IF (phic.GE.0) isign = 1
        IF (phic.LT.0) isign = -1
        call lccone(truelat1,truelat2,isign,confac)
   
   !  Calculate the distance from the pole to J,I
   
        IF (phic .GE. 0.0) THEN
           di = -re/confac * SIN(pi * 0.5 - (truelat1 * pi/180.0)) * &
                (TAN((pi * 0.5 - rclt) * 0.5) / &
                TAN((pi * 0.5 - (truelat1 * pi/180.0)) * 0.5))**confac + &
                (y - ci) * ds
        ELSE
           di = -re/confac * SIN(-pi * 0.5 - (truelat1 * pi/180.0)) * &
                (TAN((-pi * 0.5 - rclt) * 0.5) / &
                TAN((-pi*0.5 - (truelat1 * pi/180.0)) * 0.5))**confac + &
                (y - ci) * ds
        END IF
   
   !  Calculate out the Big Messy equation refered to as c1 in the document
   !  from which this formula was taken
   
        bm = TAN((pi * 0.5 - ABS(truelat1 * pi/180.0))/2.0) * &
             (confac/re * SQRT(dj**2 + di**2) / &
             SIN(pi * 0.5 - ABS(truelat1 * pi/180.0)))**(1.0/confac)
   
   !  Calculate the desired latitude in radians
   
        IF (phic .GE. 0.0) THEN
           xlat = pi * 0.5 - 2.0 * aTAN(bm)
        ELSE
           xlat = -pi * 0.5 + 2.0 * aTAN(bm)
        END IF
   
   !  Calculate the desired longitude in radians
   
        IF (phic .GE. 0.0) THEN
           xlon = rcln + (1.0/confac) * ATAN2(dj,-di)
        ELSE
           xlon = rcln + (1.0/confac) * ATAN2(dj,di)
        END IF
   
   !  If the projection is polar stereographic ('ST') THEN ...
   
     ELSE IF (project(1:2) .EQ. 'ST') THEN
   
   !  Calculate the distance J,I lies from the "true" point
   
        IF (phic .GT. 0.0) THEN
           di = -re * SIN(pi * 0.5 - rclt) * &
                (1.0 + COS(pi * 0.5 - (truelat1 * pi/180.0))) / &
                (1.0 + COS(pi * 0.5 - rclt)) + &
                (y - ci) * ds 
        ELSE
           di = -re * SIN(-pi * 0.5 - rclt) * &
                (1.0 + COS(-pi * 0.5 - (truelat1 * pi/180.0))) / &
                (1.0 + COS(-pi * 0.5 - rclt)) + &
                (y - ci) * ds
        END IF
   
   !  Calculate the Big Messy quantity as would be done, for lambert conformal
   !  projections.  This quantity is different in value, same in purpose of
   !  BM above
   
        IF (phic .GE. 0.0) THEN
           bm = (1/re) * SQRT(dj**2 + di**2) / &
                (1.0 + COS(pi * 0.5 - (truelat1 * pi/180.0)))
        ELSE
           bm = (1/re) * SQRT(dj**2 + di**2) / &
                (1.0 + COS(-pi * 0.5 - (truelat1 * pi/180.0)))
        END IF
   
   !  Calculate the desired latitude in radians
   
        IF (phic .GE. 0.0) THEN
           xlat = pi * 0.5 - 2.0 * aTAN(bm)
        ELSE
           xlat = -pi * 0.5 + 2.0 * aTAN(bm)
        END IF
   
   !  Calculate the desired longitude in radians
   
        IF (phic .GE. 0.0) THEN
           xlon = rcln + ATAN2(dj,-di)
        ELSE
           xlon = rcln + ATAN2(dj,di)
        END IF
     END IF
   
   !  Convert the calculated lat,lon pair into degrees
   
     xlat = xlat* 180.0/pi
     xlon = xlon * 180.0/pi
   
   !  Make sure no values are greater than 180 degrees and none
   !  are less than -180 degrees
   
     IF (xlon .GT. 180.0)  xlon = xlon - 360.0
     IF (xlon .LT. -180.0) xlon = xlon + 360.0
   
   END SUBROUTINE xytoll

!-----------------------------------------------------------
   
   SUBROUTINE lccone (truelat1,truelat2,sign,confac)
     REAL, parameter :: conv=0.01745329251994
     INTEGER :: sign
     REAL :: truelat1,truelat2,confac
     IF (ABS(truelat1-truelat2).LT.1.E-2) THEN
        confac = SIN(truelat1*conv)
     ELSE
        confac = LOG10(COS(truelat1*conv))-LOG10(COS(truelat2*conv))
        confac = confac/(LOG10(TAN((45.-REAL(sign)*truelat1/2.)*conv))-&
             LOG10(TAN((45.-REAL(sign)*truelat2/2.)*conv)))
     END IF
   END SUBROUTINE lccone

!-----------------------------------------------------------
   
      SUBROUTINE ll_to_ij( map_proj, truelat1, truelat2, stdlon, &
                           lat1, lon1, pole_lat, pole_lon, &
                           knowni, knownj, dx, dy, latinc, loninc,  &
                           lat, lon, i, j )


         ! Converts input lat/lon values to the cartesian (i,j) value
         ! for the given projection. 

         INTEGER map_proj
         REAL truelat1, truelat2, stdlon
         REAL lat1, lon1, pole_lat, pole_lon, knowni, knownj
         REAL dx, dy, latinc, loninc, lat, lon

         REAL clain, dlon, rsw, deltalon, deltalat
         REAL reflon, scale_top, ala1, alo1, ala, alo, rm, polei, polej
         REAL REbydx   ! Earth radius divided by dx
         REAL deltalon1, tl1r, ctl1r, arg, cone, hemi
         REAL i, j
         REAL lat1n, lon1n, olat, olon

         REAL PI, RAD_PER_DEG, DEG_PER_RAD, RE_M

!!!      lat1     ! SW latitude (1,1) in degrees (-90->90N)
!!!      lon1     ! SW longitude (1,1) in degrees (-180->180E)
!!!      dx       ! Grid spacing in meters at truelats
!!!      dlat     ! Lat increment for lat/lon grids
!!!      dlon     ! Lon increment for lat/lon grids
!!!      stdlon   ! Longitude parallel to y-axis (-180->180E)
!!!      truelat1 ! First true latitude (all projections)
!!!      truelat2 ! Second true lat (LC only)
!!!      hemi     ! 1 for NH, -1 for SH
!!!      cone     ! Cone factor for LC projections
!!!      polei    ! Computed i-location of pole point
!!!      polej    ! Computed j-location of pole point
!!!      rsw      ! Computed radius to SW corner
!!!      knowni   ! X-location of known lat/lon
!!!      knownj   ! Y-location of known lat/lon
!!!      RE_M     ! Radius of spherical earth, meters
!!!      REbydx   ! Earth radius divided by dx

         PI = 3.141592653589793
         RAD_PER_DEG = PI/180.
         DEG_PER_RAD = 180./PI
         RE_M     = 6370000.  ! Radius of spherical earth, meters
         REbydx = RE_M / dx

         hemi = 1.0
         if ( truelat1 < 0.0 ) then
           hemi = -1.0
         endif 


         !MERCATOR
         IF ( map_proj .eq. 3 ) THEN

            !  Preliminary variables
            clain = COS(RAD_PER_DEG*truelat1)
            dlon = dx / (RE_M * clain)

            ! Compute distance from equator to origin, and store in the rsw tag.
            rsw = 0.
            IF (lat1 .NE. 0.) THEN
               rsw = (ALOG(TAN(0.5*((lat1+90.)*RAD_PER_DEG))))/dlon
            ENDIF
      
            deltalon = lon - lon1
            IF (deltalon .LT. -180.) deltalon = deltalon + 360.
            IF (deltalon .GT. 180.) deltalon = deltalon - 360.
            i = knowni + (deltalon/(dlon*DEG_PER_RAD))
            j = knownj + (ALOG(TAN(0.5*((lat + 90.) * RAD_PER_DEG)))) /  &
                   dlon - rsw
       
         !PS
         ELSEIF ( map_proj .eq. 2 ) THEN

            reflon = stdlon + 90.

            ! Compute numerator term of map scale factor
            scale_top = 1. + hemi * SIN(truelat1 * RAD_PER_DEG)
      
            ! Compute radius to lower-left (SW) corner
            ala1 = lat1 * RAD_PER_DEG
            rsw = REbydx*COS(ala1)*scale_top/(1.+hemi*SIN(ala1))
      
            ! Find the pole point
            alo1 = (lon1 - reflon) * RAD_PER_DEG
            polei = knowni - rsw * COS(alo1)
            polej = knownj - hemi * rsw * SIN(alo1)
      
            ! Find radius to desired point
            ala = lat * RAD_PER_DEG
            rm = REbydx * COS(ala) * scale_top/(1. + hemi *SIN(ala))
            alo = (lon - reflon) * RAD_PER_DEG
            i = polei + rm * COS(alo)
            j = polej + hemi * rm * SIN(alo)

         !LAMBERT
         ELSEIF ( map_proj .eq. 1 ) THEN

            IF (ABS(truelat2) .GT. 90.) THEN
               truelat2=truelat1
            ENDIF

            IF (ABS(truelat1-truelat2) .GT. 0.1) THEN
               cone=(ALOG(COS(truelat1*RAD_PER_DEG))-            &
                     ALOG(COS(truelat2*RAD_PER_DEG))) /          &
               (ALOG(TAN((90.-ABS(truelat1))*RAD_PER_DEG*0.5 ))- &
                ALOG(TAN((90.-ABS(truelat2))*RAD_PER_DEG*0.5 )) )
            ELSE
               cone = SIN(ABS(truelat1)*RAD_PER_DEG )
            ENDIF
      
            ! Compute longitude differences and ensure we stay out of the
            ! forbidden "cut zone"
            deltalon1 = lon1 - stdlon
            IF (deltalon1 .GT. +180.) deltalon1 = deltalon1 - 360.
            IF (deltalon1 .LT. -180.) deltalon1 = deltalon1 + 360.
      
            ! Convert truelat1 to radian and compute COS for later use
            tl1r = truelat1 * RAD_PER_DEG
            ctl1r = COS(tl1r)
      
            ! Compute the radius to our known lower-left (SW) corner
            rsw = REbydx * ctl1r/cone * &
                  (TAN((90.*hemi-lat1)*RAD_PER_DEG/2.) / &
                   TAN((90.*hemi-truelat1)*RAD_PER_DEG/2.))**cone
      
            ! Find pole point
            arg = cone*(deltalon1*RAD_PER_DEG)
            polei = hemi*knowni - hemi * rsw * SIN(arg)
            polej = hemi*knownj + rsw * COS(arg)
      
            ! Compute deltalon between known longitude and standard lon and ensure
            ! it is not in the cut zone
            deltalon = lon - stdlon
            IF (deltalon .GT. +180.) deltalon = deltalon - 360.
            IF (deltalon .LT. -180.) deltalon = deltalon + 360.
      
            ! Radius to desired point
            rm = REbydx * ctl1r/cone *   &
                 (TAN((90.*hemi-lat)*RAD_PER_DEG/2.) /   &
                  TAN((90.*hemi-truelat1)*RAD_PER_DEG/2.))**cone
      
            arg = cone*(deltalon*RAD_PER_DEG)
            i = polei + hemi * rm * SIN(arg)
            j = polej - rm * COS(arg)
      
            ! Finally, if we are in the southern hemisphere, flip the i/j
            ! values to a coordinate system where (1,1) is the SW corner
            ! (what we assume) which is different than the original NCEP
            ! algorithms which used the NE corner as the origin in the 
            ! southern hemisphere (left-hand vs. right-hand coordinate?)
            i = hemi * i
            j = hemi * j


        !lat-lon
        ELSEIF ( map_proj .eq. 6 ) THEN

          if ( pole_lat /= 90. ) then
            call rotate_coords(lat,lon,olat,olon, &
              pole_lat,pole_lon,stdlon,-1)
            lat = olat
            lon = olon + stdlon
          end if

          ! make sure center lat/lon is good
          if ( pole_lat /= 90. ) then
            call rotate_coords(lat1,lon1,olat,olon, &
              pole_lat,pole_lon,stdlon,-1)
            lat1n = olat
            lon1n = olon + stdlon
            deltalat = lat - lat1n
            deltalon = lon - lon1n
          else
            deltalat = lat - lat1
            deltalon = lon - lon1
          end if

          ! Compute i/j
          i = deltalon/loninc
          j = deltalat/latinc

          i = i + knowni
          j = j + knownj

         ELSE

           print*,"ERROR: Do not know map projection ", map_proj

         ENDIF

      RETURN
      END  SUBROUTINE ll_to_ij


      SUBROUTINE ij_to_ll( map_proj, truelat1, truelat2, stdlon,  &
                           lat1, lon1, pole_lat, pole_lon,  &
                           knowni, knownj, dx, dy, latinc, loninc,  &
                           ai, aj, lat, lon )

        ! Converts input lat/lon values to the cartesian (i,j) value
        ! for the given projection.

        INTEGER map_proj
        REAL truelat1, truelat2, stdlon
        REAL lat1, lon1, pole_lat, pole_lon, knowni, knownj
        REAL dx, dy, latinc, loninc, ai, aj

        REAL clain, dlon, rsw, deltalon, deltalat
        REAL reflon, scale_top, ala1, alo1, ala, alo, rm, polei, polej
        REAL REbydx   ! Earth radius divided by dx
        REAL deltalon1, tl1r, ctl1r, arg, cone, hemi

        REAL PI, RAD_PER_DEG, DEG_PER_RAD, RE_M

        REAL inew, jnew, r, r2
        REAL chi,chi1,chi2
        REAL xx, yy, lat, lon
  
        REAL rlat, rlon, olat, olon, lat1n, lon1n
        REAL phi_np, lam_np, lam_0, dlam
        REAL sinphi, cosphi, coslam, sinlam
        REAL gi2, arccos


!!!     lat1     ! SW latitude (1,1) in degrees (-90->90N)
!!!     lon1     ! SW longitude (1,1) in degrees (-180->180E)
!!!     dx       ! Grid spacing in meters at truelats
!!!     dlat     ! Lat increment for lat/lon grids
!!!     dlon     ! Lon increment for lat/lon grids
!!!     stdlon   ! Longitude parallel to y-axis (-180->180E)
!!!     truelat1 ! First true latitude (all projections)
!!!     truelat2 ! Second true lat (LC only)
!!!     hemi     ! 1 for NH, -1 for SH
!!!     cone     ! Cone factor for LC projections
!!!     polei    ! Computed i-location of pole point
!!!     polej    ! Computed j-location of pole point
!!!     rsw      ! Computed radius to SW corner
!!!     knowni   ! X-location of known lat/lon
!!!     knownj   ! Y-location of known lat/lon
!!!     RE_M     ! Radius of spherical earth, meters
!!!     REbydx   ! Earth radius divided by dx
   
        PI = 3.141592653589793
        RAD_PER_DEG = PI/180.
        DEG_PER_RAD = 180./PI
        RE_M     = 6370000.  ! Radius of spherical earth, meters
        REbydx = RE_M / dx

        hemi = 1.0
        if ( truelat1 < 0.0 ) then
          hemi = -1.0
        endif


        !MERCATOR
        IF ( map_proj .eq. 3 ) THEN     

          !  Preliminary variables
          clain = COS(RAD_PER_DEG*truelat1)
          dlon = dx / (RE_M * clain)

          ! Compute distance from equator to origin, and store in the rsw tag.
          rsw = 0.
          IF (lat1 .NE. 0.) THEN
             rsw = (ALOG(TAN(0.5*((lat1+90.)*RAD_PER_DEG))))/dlon
          ENDIF

          lat = 2.0*ATAN(EXP(dlon*(rsw + aj-knownj)))*deg_per_rad - 90.
          lon = (ai-knowni)*dlon*deg_per_rad + lon1
          IF (lon.GT.180.) lon = lon - 360.
          IF (lon.LT.-180.) lon = lon + 360.


        !PS
        ELSEIF ( map_proj .eq. 2 ) THEN     

          ! Compute the reference longitude by rotating 90 degrees to the east
          ! to find the longitude line parallel to the positive x-axis.
          reflon = stdlon + 90.
    
          ! Compute numerator term of map scale factor
          scale_top = 1. + hemi * SIN(truelat1 * rad_per_deg)

          ! Compute radius to known point
          ala1 = lat1 * RAD_PER_DEG
          rsw = REbydx*COS(ala1)*scale_top/(1.+hemi*SIN(ala1))

          ! Find the pole point
          alo1 = (lon1 - reflon) * RAD_PER_DEG
          polei = knowni - rsw * COS(alo1)
          polej = knownj - hemi * rsw * SIN(alo1)

          ! Compute radius to point of interest
          xx = ai - polei
          yy = (aj - polej) * hemi
          r2 = xx**2 + yy**2
    
          ! Now the magic code
          IF (r2 .EQ. 0.) THEN
             lat = hemi * 90.
             lon = reflon
          ELSE
             gi2 = (rebydx * scale_top)**2.
             lat = deg_per_rad * hemi * ASIN((gi2-r2)/(gi2+r2))
             arccos = ACOS(xx/SQRT(r2))
             IF (yy .GT. 0) THEN
                lon = reflon + deg_per_rad * arccos
             ELSE
                lon = reflon - deg_per_rad * arccos
             ENDIF
          ENDIF

          ! Convert to a -180 -> 180 East convention
          IF (lon .GT. 180.) lon = lon - 360.
          IF (lon .LT. -180.) lon = lon + 360.

        !LAMBERT
        ELSEIF ( map_proj .eq. 1 ) THEN     

          IF (ABS(truelat2) .GT. 90.) THEN
            truelat2=truelat1
          ENDIF

          IF (ABS(truelat1-truelat2) .GT. 0.1) THEN
            cone=(ALOG(COS(truelat1*RAD_PER_DEG))-   &
                  ALOG(COS(truelat2*RAD_PER_DEG))) /          &
                  (ALOG(TAN((90.-ABS(truelat1))*RAD_PER_DEG*0.5 ))-  &
                  ALOG(TAN((90.-ABS(truelat2))*RAD_PER_DEG*0.5 )) )
          ELSE
            cone = SIN(ABS(truelat1)*RAD_PER_DEG )
          ENDIF
   
          ! Compute longitude differences and ensure we stay out of the
          ! forbidden "cut zone"
          deltalon1 = lon1 - stdlon 
          IF (deltalon1 .GT. +180.) deltalon1 = deltalon1 - 360.
          IF (deltalon1 .LT. -180.) deltalon1 = deltalon1 + 360.

          ! Convert truelat1 to radian and compute COS for later use
          tl1r = truelat1 * RAD_PER_DEG
          ctl1r = COS(tl1r)
      
          ! Compute the radius to our known point
          rsw = REbydx * ctl1r/cone *  &
                (TAN((90.*hemi-lat1)*RAD_PER_DEG/2.) /  &
                 TAN((90.*hemi-truelat1)*RAD_PER_DEG/2.))**cone

          ! Find pole point
          alo1 = cone*(deltalon1*RAD_PER_DEG)
          polei = knowni - rsw * SIN(alo1)
          polej = knownj + hemi * rsw * COS(alo1)

          chi1 = (90. - hemi*truelat1)*rad_per_deg
          chi2 = (90. - hemi*truelat2)*rad_per_deg
 
          ! See if we are in the southern hemispere and flip the indices if we are.
          inew = hemi * ai
          jnew = hemi * aj

          ! Compute radius**2 to i/j location
          reflon = stdlon + 90.
          xx = inew - polei
          yy = polej - jnew
          r2 = (xx*xx + yy*yy)
          r = SQRT(r2)/rebydx

          ! Convert to lat/lon
          IF (r2 .EQ. 0.) THEN
            lat = hemi * 90.
            lon = stdlon
          ELSE
            lon = stdlon + deg_per_rad * ATAN2(hemi*xx,yy)/cone
            lon = AMOD(lon+360., 360.)
            IF (chi1 .EQ. chi2) THEN
             chi =2.0*ATAN((r/TAN(chi1))**(1./cone)*TAN(chi1*0.5))
            ELSE
             chi =2.0*ATAN((r*cone/SIN(chi1))**(1./cone)*TAN(chi1*0.5))
            ENDIF
            lat = (90.0-chi*deg_per_rad)*hemi
          ENDIF
 
          IF (lon .GT. +180.) lon = lon - 360.
          IF (lon .LT. -180.) lon = lon + 360.


        !lat-lon
        ELSEIF ( map_proj .eq. 6 ) THEN     

          inew = ai - knowni
          jnew = aj - knownj
    
          if (inew < 0.)        inew = inew + 360./loninc
          if (inew >= 360./dx)  inew = inew - 360./loninc
    
          ! Compute deltalat and deltalon
          deltalat = jnew*latinc
          deltalon = inew*loninc
        
          if ( pole_lat /= 90. ) then
            call rotate_coords(lat1,lon1,olat,olon, &
              pole_lat,pole_lon,stdlon,-1)
            lat1n = olat
            lon1n = olon + stdlon
            lat = deltalat + lat1n
            lon = deltalon + lon1n
          else
            lat = deltalat + lat1
            lon = deltalon + lon1
          end if

   
          if ( pole_lat /= 90. ) then
            lon = lon - stdlon
            call rotate_coords(lat,lon,olat,olon, &
              pole_lat,pole_lon,stdlon,1)
            lat = olat
            lon = olon
          end if

          if (lon < -180.) lon = lon + 360.
          if (lon >  180.) lon = lon - 360.

        ELSE
     
           print*,"ERROR: Do not know map projection ", map_proj

        ENDIF

        RETURN

      END SUBROUTINE ij_to_ll


      SUBROUTINE rotate_coords(ilat,ilon,olat,olon, &
                 lat_np,lon_np,lon_0,direction)
        REAL ilat, ilon
        REAL olat, olon
        REAL lat_np, lon_np, lon_0
        INTEGER direction

        ! >=0, default : computational -> geographical
        ! < 0          : geographical  -> computational
  
        REAL rlat, rlon
        REAL phi_np, lam_np, lam_0, dlam
        REAL sinphi, cosphi, coslam, sinlam
        REAL PI, RAD_PER_DEG, DEG_PER_RAD 

        PI = 3.141592653589793
        RAD_PER_DEG = PI/180.
        DEG_PER_RAD = 180./PI
  
        ! Convert all angles to radians
        phi_np = lat_np * rad_per_deg
        lam_np = lon_np * rad_per_deg
        lam_0  = lon_0  * rad_per_deg
        rlat = ilat * rad_per_deg
        rlon = ilon * rad_per_deg
  
        IF (direction < 0) THEN
           ! The equations are exactly the same except for one small difference
           ! with respect to longitude ...
           dlam = PI - lam_0
        ELSE
           dlam = lam_np
        END IF
        sinphi = COS(phi_np)*COS(rlat)*COS(rlon-dlam) +  &
                 SIN(phi_np)*SIN(rlat)
        cosphi = SQRT(1.-sinphi*sinphi)
        coslam = SIN(phi_np)*COS(rlat)*COS(rlon-dlam) -  &
                 COS(phi_np)*SIN(rlat)
        sinlam = COS(rlat)*SIN(rlon-dlam)
        IF ( cosphi /= 0. ) THEN
           coslam = coslam/cosphi
           sinlam = sinlam/cosphi
        END IF
        olat = deg_per_rad*ASIN(sinphi)
        olon = deg_per_rad*(ATAN2(sinlam,coslam)-dlam-lam_0+lam_np)

      END SUBROUTINE rotate_coords

END MODULE map_stuff
