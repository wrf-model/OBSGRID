MODULE map_stuff

CONTAINS

!-------------------------------------------
   SUBROUTINE plotmap (bhi, bhr)

      INTEGER , DIMENSION(50,20) :: bhi
      REAL    , DIMENSION(20,20) :: bhr
      CHARACTER (LEN=2) :: project
      CHARACTER (LEN=2) , DIMENSION(3) , PARAMETER :: nproj = (/ 'LC','ST','ME' /)

      REAL :: ds , phic , xlonc , xsouth , xwest , xratio , truelat1 , truelat2 
      REAL :: pl1 , pl2 , pl3 , pl4 , polat , rot
      REAL :: xleft , xright
      REAL :: vl , vr , vb , vt , wl , wr , wb , wt
      INTEGER :: ixmoad , jxmoad , kproj , ix , jx , jlts , jgrid , idot , iusout , ier , ls

      !  What projection are we using?

      project = nproj(bhi(7,1))

      !  Grid distance in km, center lat/lon.
      
      ds = bhr(9,1)/1000.
      phic = bhr(2,1)
      xlonc = bhr(3,1)

      !  The i and j dimensions, in Mother Of All Domain space.

      ixmoad = bhi(5,1)
      jxmoad = bhi(6,1)

      !  The starting locations of this domain, wrt to MOAD.

      IF ((bhi(8,1).EQ.0) .and. (bhi(13,1).EQ.1)) THEN 
         xsouth = bhr(10,1)
         xwest = bhr(11,1)
      ELSEIF (bhi(13,1) .GT. 1) THEN
         xsouth = bhr(10,1)
         xwest = bhr(11,1)
      ELSE
         xsouth = bhr(10,1) - bhi(11,1)
         xwest = bhr(11,1) - bhi(12,1)
      END IF

      !  The ratio of the grid distance of the mother domain to the grid distance
      !  from this domain.

      xratio = REAL(bhi(20,1))

      !  The latitude where the projection is true (map factors = 1).

      truelat1 = bhr(5,1)
      truelat2 = bhr(6,1)

      !  Get the lat/lon of the lower left hand corner.

      CALL xytoll(1.,1.,pl1,pl2,project,ds,phic,xlonc,ixmoad,jxmoad,xsouth,xwest,xratio,truelat1,truelat2)

      !  Is there some sort of expanded domain with which to concern ourselves?

      IF ((bhi(8,1) .EQ. 1) .and. (bhi(13,1) .EQ. 1)) THEN
         ix = bhi(9,1)
         jx = bhi(10,1)
      ELSE
         ix = bhi(16,1)
         jx = bhi(17,1)
      END IF

      !  Get the lat/lon of the upper right corner point.

      CALL xytoll(REAL(jx),REAL(ix),pl3,pl4,project,ds,phic,xlonc,ixmoad,jxmoad,xsouth,xwest,xratio,truelat1,truelat2)
      
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

      CALL supmap(kproj,polat,xlonc,rot,pl1,pl2,pl3,pl4,jlts,jgrid,iusout,idot,ier)

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
   
END MODULE map_stuff
