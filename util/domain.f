      program domain
c
c     ... this file is a simple way to decide on the domain
c         to choose in the mm4 system
c
c         1) copy this file (domain.f) to $TMPDIR
c                # cd $TMPDIR
c                # cp ~mesouser/Decks/Terrain/domain.f .
c         2) compile and load the file
c                # ncargf77 domain.f
c         3) run it, and answer the questions
c                # a.out
c         4) look at the plots from a vt100 terminal from shavano
c                # ictrans -d vt100 gmeta
c
      character *80 fname
      character *1  cdum
      common /xyll/ xn,ds,xlatc,xlonc,imax,jmax,kproj
c
c     ... joe gks open jazz
c
      call opngks
      call gsclip (0)
c
c     ... read in the pertinent info from the user
c
      print *,'enter imax,jmax (i refers to y, j refers to x direction)'
      read (5,*) imax,jmax
c
      print *,'enter grid distance in km '
      read(5,*) ds
c
      print *,'enter center latitude, center longitude '
      read(5,*) xlatc,xlonc
c
      print *,'enter the lower left corner i,j '
      read(5,*) plly,pllx
c
      print *,'enter the upper right corner i,j '
      read(5,*) pury,purx
c
      print *,'enter type of projection: (1)lambert (2)polar ',
     *            '(3) mercator'
      read(5,*)kproj
c
c     ... user supplied map info into ezmap choices
c
      if(kproj.eq.1) jproj=3
      if(kproj.eq.2) jproj=1
      if(kproj.eq.3) jproj=9
      if (kproj.eq.3) then
         xn=0.0
         polat=0.
         rot=0.
      else if (kproj.eq.1) then
         xn=0.716
         if (xlatc.ge.0) then
            rot=30.
            polat=60.
          else
            rot=-30
            polat=-60.
          end if
      else if(kproj.eq.2) then
           xn=1.0
           rot=0.
           polat=90.
           if (xlatc.lt.0) polat=-90.
      end if
      polong=xlonc
c
c     ... do not grind them with details, set up some simple things
c
      jgrid=10
      jlts=-2
      iusout=1
      idot=0
c
c     ... get lower left point in lat/lon
c
      call xyll1(1.,1.,xlat,xlon)
      pl1=xlat
      pl2=xlon
c
c     ... get upper right point in lat/lon
c
      call xyll1(purx,float(imax),xlat,xlon)
      pl3=xlat
      pl4=xlon
      print *,'lower left   lat,lon=',pl1,pl2
      print *,'upper right  lat,lon=',pl3,pl4
c
c     ... feed in lower left, upper right, etc to map routine
c
      call supmap(jproj,polat,polong,rot,pl1,pl2,pl3,pl4,jlts,
     1            jgrid,iusout,idot,ier)
      if (ier.ne.0) then
         print *,'error:  ier=',ier
         go to 111
      end if
      call setusv('LW',1000)
      call perim(INT(purx)-1,1,imax-1,1)


      call frame
c
c     ... all back stop
c
      call clsgks
      stop 99999
111   continue
      print *,'hey, this wasn''t so easy after all'
      call frame
c     call gclwk (1) 
      call clsgks
      stop 111
      end

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine xyll1(x,y,xlat,xlon)
c
      common /xyll/ xn,ds,xlatc,xlonc,imax,jmax,kproj
c
c     ... xn:          cone factor 
c                         1.0 for polar stereographic.
c                         0.716 for lambert conformal.
c                         0.0 for mercator.
c     ... ds:          grid spacing in km.
c     ... xlatc,xlonc: central latitude and longitude.
c     ... imax,jmax:   grid number in y & x direction
c     ... kproj:       projection
c                         1 lambert
c                         2 polar
c                         3 mercator
c
      a=6370.
      conv=57.29578
      centri=(float(imax)+1.)/2.
      centrj=(float(jmax)+1.)/2.

      if(kproj.eq.1) then  ! lambert conformal
         phi1=30.0  ! true at 30
         phic=90.0-xlatc
         phi1=phi1/conv
         phic=phic/conv
         cnst=a*sin(phi1)/(xn*tan(phi1/2.)**xn)
         xxc=0.0
         yyc=-1.*cnst*(tan(phic/2.)**xn)
         xx=(x-centrj)*ds+xxc
         yy=(y-centri)*ds+yyc
         r=sqrt(xx*xx+yy*yy)
         yxn=1./xn
         phi=(r*xn/(sin(phi1)*a))**yxn
         phi=phi*tan(phi1/2.)
         phi=2.0*atan(phi)
         phi=phi*conv
         xlat=90.0-phi
         s=asin(xx/r)
         s=s*conv
         s=s/xn
         xlon=xlonc+s
c
      else if(kproj.eq.2) then  ! polar stereo
         pi=3.14159265358
         xlamc=xlonc*pi/180.
         phic=xlatc*pi/180.
c        psi1=pi/4.5  ! true at 40
         psi1=pi/3.   ! true at 60
         xx=(x-float(jmax+1)/2.)* ds
         yy=-a/xn*sin(pi/2.-phic)*
     *      (  (1.+cos(pi/2.-psi1)) / (1.+cos(pi/2.-phic))  )**xn +
     *      (y-float(imax+1)/2.)* ds
         c1=1./a * sqrt(xx**2 + yy**2) *
     *      (  1./(1. + cos(pi/2.-psi1)))
         xlat=pi/2. - 2.*atan(c1)
         xlon=xlamc+ 1./xn * atan2(xx,-yy)
         xlat=xlat*180./pi
         xlon=xlon*180./pi
c
      else if(kproj.eq.3) then  ! mercator
         phi1 = 0.0
         phi1 = phi1/conv
         c2   = a*cos(phi1)
         phictr=xlatc/conv
         xcntr= 0.0
         cell = cos(phictr)/(1.0+sin(phictr))
         ycntr= -c2*alog(cell)
         xx   = (x-centrj)*ds
         yy   = (y-centri)*ds+ycntr
         xlon = xlonc + ((xx - xcntr)/c2) * conv
         cell = exp(yy/c2)
         xlat = 2.0*(conv*atan(cell)) - 90.
c
      endif
c
      return
      end
