module skewt_module
contains

! FUNCTIONS:  MAPPINGS FROM (P,T) TO CM ON SKEWT
  real function skewt_fy(P) result(fy)
    real :: P ! Pressure in mb
    fy = 132.182 - 44.061 * ALOG10(P)
  end function skewt_fy

  real function skewt_fx(T,Y) result(fx)
    real :: T ! T in C
    real :: Y
    fx = 0.54 * T + 0.90692 * Y
  end function skewt_fx

  subroutine skewt_opngks(fname)

    character(len=*) :: fname
    character(len=80) :: cdum, fdum

    COMMON /SKWDRW/ ISKDRW

    ISKDRW = 0    !  Make sure to redraw the skew-t background.

    call gopks(6,idum)

    fdum = trim(fname)

    call gesc(-1391,1,fdum,1,1,cdum)
    call gopwk(1,3,1)
    call gopwk(2,9,3)
    call gacwk(1)
    call gacwk(2)
    call gsclip(0)

    call dfclrs

    return
  end subroutine skewt_opngks

subroutine skewt_clsgks

  call gdawk(1)
  call gclwk(1)
  call gdawk(2)
  call gclwk(2)
  call gclks

  return
end subroutine skewt_clsgks


SUBROUTINE SKEWT (PRES,TEMP,DWPT,SPD,DIR,NN,TITLINE,&
     HDATE,ISTAT,XLOC,YLOC,YLAT,YLON,ILW,ICTL,ICDL,LCF,WOFFS,&
     BADVAL)
!
!     SKEWT- PLOTS SOUNDINGS ON A SKEWT, LOG P THERMODYNAMIC DIAGRAM
!
!     PRES- PRESSURE ARRAY FOR THERMODYNAMIC DATA (Pa)
!     TEMP- TEMPERATURE ARRAY (K)
!     DWPT- DEW POINT ARRAY (CELSIUS)
!     SPD- WIND SPEED ARRAY (M/S)
!     DIR- WIND DIRECTION ARRAY (DEGREES-0 IS NORTH)
!     NN- NUMBER OF DATA LEVELS
!     TITLINE - 80-character title.
!     HDATE - MODEL INITIAL TIME YYMMDDHH
!     ISTAT- STRING CONTAINING STATION NAME (40 CHAR.)
!     XLOC, YLOC - X and Y locations
!     YLAT, YLON - Lat and Lon
!     ICTL - Color index for the temperature line.
!     ICDL - Color index for the temperature line.
!     ILW -  LINE WIDTH ( 2000 recommended for thin, 6000 recommended for thick)
!     LCF - LOGICAL CALL FRAME.  If .TRUE., FRAME is called.
!     BADVAL- VALUE ASSIGNED TO MISSING DATA  --TEMP,DWPT TRACES ARE
!                   TERMINATED AS SOON AS THIS NUMBER IS ENCOUNTERED.
!
!        OUTPUT PARAMETERS.....
!     ALL INPUT PARAMETERS REMAIN UNCHANGED BY THIS ROUTINE.
!
  real, dimension(nn) :: pres,temp,dwpt,spd,dir
  CHARACTER(len=120) :: LAB
  character(len=40)  :: ISTAT
  character(len=*)   :: hdate
!
  CHARACTER(LEN=*) ::    TITLINE
  INTEGER ::          ITSTRT,&
       ITITLEN,&
       LCNTR
  COMMON /SKWDRW/ ISKDRW
!
  LOGICAL LCF
!
!
!  DEGREES TO RADIANS
  PARAMETER (DTR = 0.0174532925)
!  WIND BARB DATA
  PARAMETER (XM = 24.2)
!
!  LINE WIDTH VARIABLES
!
  DATA ISWIDE / 1000 /

  logical :: NFP

!
!  STATEMENT FUNCTIONS:  MAPPINGS FROM (P,T) TO CM ON SKEWT
!
  if (titline.eq."40754") then
     do i = 1, nn
        print*, pres(i), temp(i), dwpt(i)
     enddo
  endif

!
!  Initialization
!

  DO I = 1, 120
     LAB(I:I) = ' '
  ENDDO
!
!  Check Character Strings
!
  DO I = 1, 40
     if ( (ichar(istat(i:i)).lt.32) .or.&
          (ichar(istat(i:i)).gt.127)) then
        print*, 'ISTAT PROBLEM.'
        istat(i:i) = '*'
     endif
  enddo
  DO I = 1, len(titline)
     if ( (ichar(titline(i:i)).lt.32) .or.&
          (ichar(titline(i:i)).gt.127)) then
        print*, 'TITLINE PROBLEM.'
        titline(i:i) = '*'
     endif
  enddo

!
!  TEST TO SEE IF A BACKGROUND HAS BEEN DRAWN, IF NOT CALL SKEWT_BACKGROUND
!
  CALL GSPLCI (4 ) ! Background color
  CALL GSPMCI (1 )
  CALL GSTXCI (1 )

  IF (ISKDRW .EQ. 0) THEN
     CALL SKEWT_BACKGROUND
     ISKDRW = 1
  END IF

!KWM      CALL GSPLCI (9 )
!KWM      CALL GSPMCI (9 )
!KWM      CALL GSTXCI (9 )
  CALL GSPLCI (4 )
  CALL GSPMCI (4 )
  CALL GSTXCI (4 )
!
!  SKEWT BACKGROUND HAS BEEN GENERATED-- PLOT THE SOUNDING
!
  CALL GFLAS3 (1)
  CALL SET(.05,.95,.05,.95,-19.0,27.1,-.9346217,44.061,1)
!

!  PUT ON TITLE

  ITSTRT = INDEX(TITLINE,':')
  IF (ITSTRT .NE. 0) THEN
     DO LCNTR=1,80-ITSTRT
        TITLINE(LCNTR:LCNTR) = TITLINE(LCNTR+ITSTRT:LCNTR+ITSTRT)
     ENDDO
     DO LCNTR=80-ITSTRT+2,80
        TITLINE(LCNTR:LCNTR)=' '
     ENDDO
     ITITLEN = 80
     LOOP_100 : DO LCNTR=80,1,-1
        IF (TITLINE(LCNTR:LCNTR) .NE. ' ') THEN
           ITITLEN = LCNTR
           GOTO 110
        END IF
     enddo LOOP_100
110  CONTINUE
  ELSE
     ITITLEN=LEN(TITLINE)
  END IF

  CALL WTSTR (4.05,-2.4,TITLINE(1:ITITLEN),12,0,0)
  if ((xloc > -900000) .and. (yloc > -900000)) then
     WRITE(LAB,106) ISTAT, HDATE(1:13), XLOC, YLOC, YLAT, YLON
106  FORMAT(1X,A40,2X,A13,'  X/Y = ',F6.2,1x,F6.2,&
         2X,'LAT/LON = ',F7.3,2X,F8.3)
  else
     WRITE(LAB,206) ISTAT, HDATE(1:13), YLAT, YLON
206  FORMAT(1X,A40,2X,A13,28X,'LAT/LON = ',F7.3,2X,F8.3)
  endif
  CALL WTSTR (-19.5,-1.4,trim(LAB),8,0,-1)
  CALL GSPMCI (1 )
  CALL GSTXCI (1 )
  CALL GSPLCI (ICTL)   !  Make wind-barb line the color of the T line.
  IF(NN.GT.0) CALL LINE(XM+woffs,-.9346217,XM+woffs,44.061)
  CALL GSPLCI (1 )
!
!  SOLID DASH PATTERN, INCREASED SPOT SIZE (DEFAULT=8)
!
  CALL DASHDB (65535)
  CALL GETUSV ('LW',ISNORM)
  CALL SETUSV ('LW',ISWIDE)
!
  CALL SETUSV('LW',ilw)
  CALL GSPLCI (ictl)
  CALL GSPMCI (ictl)
  CALL GSTXCI (ictl)

  II = 0
  NFP = .TRUE.
  DO I=1,NN
     IF ((TEMP(I) .gt. -111111.).and.(PRES(I).gt.-111111.)) THEN
        II =  II + 1
        Y=skewt_FY(PRES(I) * 0.01)
        X=skewt_FX(TEMP(I)-273.15,Y)
        IF (II==1) CALL FRSTPT(X,Y)
        CALL VECTOR(X,Y)
     ENDIF


!KWM     IF ((TEMP(I) .gt. -111111.).and.(PRES(I).gt.-111111.)) THEN
!KWM        II =  II + 1
!KWM        Y=skewt_FY(PRES(I) * 0.01)
!KWM        X=skewt_FX(TEMP(I)-273.15,Y)
!KWM        IF (NFP) CALL FRSTPT(X,Y)
!KWM        CALL VECTOR(X,Y)
!KWM        NFP = .FALSE.
!KWM     else
!KWM        if ( II < NN ) then
!KWM           call plotif(0,0,2)
!KWM           NFP = .TRUE.
!KWM        endif
!KWM     ENDIF
  enddo

  CALL SETUSV('LW',ilw)
  CALL GSPLCI (icdl)
  CALL GSPMCI (icdl)
  CALL GSTXCI (icdl)

  II = 0
  LOOP_70 : DO I=1,NN
     IF ((DWPT(I).gt.-111111.).and.(PRES(I).gt.-111111.))THEN
        II = II + 1
        Y=skewt_FY(PRES(I) * 0.01)
        X=skewt_FX(DWPT(I)-273.15,Y)
        IF(II.EQ.1)CALL FRSTPT(X,Y)
        CALL VECTOR(X,Y)
     ENDIF
  enddo LOOP_70

  CALL SETUSV('LW',2000)
  CALL GSPLCI (8 )
  CALL GSPMCI (8 )
  CALL GSTXCI (8 )
!
!         PLOT WIND VECTORS
!
  IF (NN.LE.0) GO TO 76

  CALL SETUSV ('LW',ISWIDE)
  DO I=1,NN
     if ((abs(spd(i)-badval).lt.1).or.(abs(dir(i)-badval).lt.1)) cycle
     IF (DIR(I) > 360.) cycle
     if (pres(i) < -111111.) cycle
     if (spd(i) < -111111.) cycle
     if (dir(i) < -111111.) cycle
     ANG=DIR(I)*DTR
     U = -SPD(I)*SIN(ANG)
     V = -SPD(I)*COS(ANG)
     Y1=skewt_FY(PRES(I) * 0.01)
     CALL SKEWT_WNDBARB (XM+WOFFS,Y1,U,V)
  enddo
76 CONTINUE
  CALL GSPLCI (1 )
  CALL GSPMCI (1 )
  CALL GSTXCI (1 )
!
!  RESET TO NORMAL SPOT SIZE AND EXIT
!
  CALL SETUSV ('LW',ISNORM)
!
  IF (LCF) THEN
     CALL FRAME
  ENDIF
!
  RETURN
END SUBROUTINE SKEWT

SUBROUTINE SKEWT_BACKGROUND
!
!  SKEWT_BACKGROUND:  PLOTS BACKGROUND FOR A SKEWT, LOG P THERMODYNAMIC DIAGRAM
!
  COMMON /SKWDRW/ ISKDRW
!
!  ENCODE BUFFER
!
  CHARACTER(LEN=4) :: ITIT
!
!  SKEWT BORDER
!
  DIMENSION XB(7),YB(7)
  DATA XB/-19.,27.1,27.1,18.6,18.6,-19.,-19./
  DATA YB/-.9346217,-.9346217,9.,17.53,44.061,44.061,-.9346217/
!
!  PRESSURE LINE SPECS
!
  DIMENSION PLV(11),PLN(11,2)
  DATA PLV/100.,200.,300.,400.,500.,600.,700.,800.,900.,&
       1000.,1050./
  DATA PLN/11*-19.,4*18.6,22.83,26.306,5*27.1/
!
!  TEMPERATURE LINE SPECS
!
  DIMENSION TP(15,2)
  DATA TP/8*1050.,855.,625.,459.,337.,247.,181.,132.,730.,580.,500.,&
       430.,342.,251.,185.,135.,7*100./
!
!  MIXING RATIO SPECS
!
  REAL ::        RAT(8)
  CHARACTER(LEN=2) :: LRAT(8)
  DATA RAT/20.,12.,8.,5.,3.,2.,1.,0.4/
  DATA LRAT/'20','12',' 8',' 5',' 3',' 2',' 1','.4'/
!
!  DRY/SATURATED ADAIBAT BUFFERS
!
  DIMENSION SX(162),SY(162),Y45(162)
!
!  DEGREES TO RADIANS, ABSOLUTE ZERO
!
  PARAMETER (ABZ = 273.16)
!
!  MAPPINGS FROM (P,T) TO CM ON SKEWT
!
!
!  DRAW SKEWT BORDER
!
  CALL GFLAS1 (1)
  CALL SET(.05,.95,.05,.95,-19.0,27.1,-.9346217,44.061,1)
  CALL CURVE(XB,YB,7)
!
!  DRAW THE PRESSURE LINES
!
  LOOP_10 : DO K=1,11
     Y1=skewt_FY(PLV(K))
     IF(K.NE.1.AND.K.NE.11) CALL LINE(PLN(K,1),Y1,PLN(K,2),Y1)
     ITS=NINT(PLV(K))
     WRITE (ITIT,101) ITS
101  FORMAT(I4)
!        CALL PWRY(-20.9,Y1,ITIT,4,1.9,0,0)
     CALL WTSTR (-19.2,Y1,ITIT,11,0,1)
  enddo LOOP_10
!
!  DRAW TEMPERATURE LINES
!
  T=40.
  LOOP_20 : DO I=1,15
     Y1=skewt_FY(TP(I,1))
     Y2=skewt_FY(TP(I,2))
     X1=skewt_FX(T,Y1)
     X2=skewt_FX(T,Y2)
     CALL LINE(X1,Y1,X2,Y2)
     ITS=NINT(T)
     IF(ITS.EQ.20) GO TO 19
     X2=X2+0.4
     Y2=Y2+0.441
     WRITE (ITIT,101) ITS
!           CALL PWRY(X2,Y2,ITIT,4,12,.83422,0)
     CALL WTSTR (X2,Y2,ITIT,12,47,-1)
19   T=T-10.
  enddo LOOP_20
!
!  TICK MARKS AT 500 MB
!
  Y1=13.2627
  Y2=13.75
  T=-52.
  LOOP_25 : DO I=1,31
     T=T+2.
     IF(AMOD(T,10.).EQ.0.)cycle LOOP_25 
     X1=skewt_FX(T,Y1)
     X2=skewt_FX(T,Y2)
     CALL LINE(X1,Y1,X2,Y2)
  enddo LOOP_25
!
!  DRAW MIXING RATIO LINES
!
  CALL DASHDB (3855)      ! PATTERN = 0000111100001111
  Y1=skewt_FY(1050.)
  Y2=skewt_FY(700.)
  LOOP_30 : DO I=1,8
     X1=skewt_FX(SKEWT_TMR(RAT(I),1050.)-ABZ,Y1)
     X2=skewt_FX(SKEWT_TMR(RAT(I), 700.)-ABZ,Y2)
     CALL LINED(X1,Y1,X2,Y2)
!        CALL PWRY(X2,Y2+0.6,LRAT(I),2,10,0,1)
     CALL WTSTR (X2,Y2+0.6,LRAT(I),10,0,0)
  enddo LOOP_30
!
!  DRAW SATURATED ADIABATS
!
  CALL DASHDB (31710)     ! PATTERN = 0111101111011110
  TS=32.
  LOOP_40 : DO I=1,7
     P=1060.
     TK=TS+ABZ
     AOS=SKEWT_OS(TK,1000.)
     LOOP_35 : DO J=1,86
        P=P-10.
        ATSA=SKEWT_TSA(AOS,P)-ABZ
        SY(J)=skewt_FY(P)
        SX(J)=skewt_FX(ATSA,SY(J))
     enddo LOOP_35
     CALL CURVED(SX,SY,86)
     ITS=NINT(TS)
     WRITE (ITIT,102) ITS
102  FORMAT(I2)
!        CALL PWRY(SX(86),SY(86)+0.6,ITIT,2,10,0,1)
     CALL WTSTR (SX(86),SY(86)+0.6,ITIT(1:2),10,0,0)
     TS=TS-4.0
  enddo LOOP_40
!
!  DRAW DRY ADIABAT LINES
!
  CALL DASHDB (21845)     ! PATTERN = 0101010101010101
!     CALL DASHD(4444B)
  T=51.
  LOOP_45 : DO I=1,162
     Y45(I)=66.67*(5.7625544-ALOG(T+ABZ))
     T=T-1.0
  enddo LOOP_45
  T=450.
  TD=52.
  LOOP_55 : DO I=1,20
     T=T-10.
     K=0
     YD=66.67*(ALOG(T)-5.7625544)
     LOOP_50 : DO J=1,162
        YPD=Y45(J)+YD
        TX=TD-FLOAT(J)
        IF(YPD.GT.44.061) GO TO 54
        IF(YPD.LT.-.9346217) cycle LOOP_50 
        XPD=skewt_FX(TX,YPD)
        IF(XPD.LT.-19.0)GO TO 54
        IF(XPD.GT.27.1)cycle LOOP_50 
        IF(XPD.GT.18.6.AND.T.GT.350.0)cycle LOOP_50 
        K=K+1
        SX(K)=XPD
        SY(K)=YPD
     enddo LOOP_50
!
54   CALL CURVED(SX,SY,K)
     ITS=NINT(T)
     WRITE (ITIT,103) ITS
103  FORMAT(I3)
!        IF(ITS .GE. 320) THEN
!           CALL PWRY(SX(K-3),43.0,ITIT,3,10,0,1)
!           ELSE
!           CALL PWRY(-18.0,SY(K-3),ITIT,3,10,0,1)
!        END IF
     X=SX(K-3)
     Y=SY(K-3)
     IF(X.LT.-15.0) X = -17.95
     IF(Y.GT.40.0)  Y = 42.9
     CALL WTSTR (X,Y,ITIT(1:3),10,0,0)
  enddo LOOP_55
!
  CALL GFLAS2
  ISKDRW = 1
  RETURN
END SUBROUTINE SKEWT_BACKGROUND

SUBROUTINE SKEWT_WNDBARB (XBASE,YBASE,U,V)

!***********************************************************************
! WNDBARB - FOR THE GRAPH PORTION OF GRIN
!  THIS ROUTINE DRAWS A SINGLE WIND BARB PER CALL.
!
! ON INPUT - FOUR VARIABLES COME IN.  XBASE CONTAINS THE HORIZONTAL COOR
!            DINATE OF THE BASE OF THE WIND BARB.  YBASE CONTAINS THE
!            VERTICAL COORDINATE OF THE BASE OF THE WIND BARB.  U CONTAI
!            THE EAST-WEST WIND COMPONENT IN METERS PER SECOND.  V CONTA
!            THE NORTH-SOUTH WIND COMPONENT IN METERS PER SECOND.
!
! ON OUTPUT - ONE BARB HAS BEEN DRAWN TO UNIT NUMBER 2, WHICH CORRESPOND
!            TO GMETA.CGM, THE GKS OUTPUT META CODE FILE.
!
! ASSUMPTIONS - THIS ROUTINE ASSUMES THAT GKS HAS BEEN OPENED AND A WORK
!            STATION HAS BEEN SET.
!
! REVISED BY - JEREMY ASBILL ON MAY 3, 1990.
!***********************************************************************

!  INPUT VARIABLE DECLARATIONS ...

  REAL ::      XBASE,YBASE,&
       U,V

!  LOCAL PARAMETER ...
!    SC SPECIFIES IN THE NORMALIZED (FRACTIONAL) GRAPHICS COORDINATE
!       SYSTEM HOW LONG THE BARB SHAFT IS TO BE
!    COORDINATE SYSTEMS ARE EXPLAINED IN NCAR GRAPHICS USER'S GUIDE VERS
!    2.00 ON PAGE 46.

  PARAMETER (SC = 0.05493)

!  LOCAL VARIABLE DECLARATIONS ...

  INTEGER ::   LLSV                 ! SAVE VARIABLE, SCALING FOR SET

  LOGICAL DONE                   ! T => SUBROUTINE ENDS, F => LOOP A

  REAL ::      WINDVCT,             &! WIND VECTOR MAGNITUDE
       FLSV,FRSV,FBSV,FTSV, &! SAVE VARIABLES, FRACTIONAL COORDI
       ULSV,URSV,UBSV,UTSV, &! SAVE VARIABLES, INCOMING USER COO
       NEWXBASE,            &! FRACTIONAL X COORD. FOR BARB BASE
       NEWYBASE,            &! FRACTIONAL Y COORD. FOR BARB BASE
       XCOMP,               &! X COMPONENT OF GRAPHICAL VECTOR (
       YCOMP,               &! Y COMPONENT OF GRAPHICAL VECTOR (
       PK,                  &! PLACE KEEPER
       FETHLENX,            &! X COMPONENT OF GRAPHICAL VECT. (F
       FETHLENY             ! Y COMPONENT OF GRAPHICAL VECT. (F

!  LOCAL ARRAY DECLARATIONS ...

  INTEGER ::   IJUNK(5)             ! CALCULATION ARRAY FOR SFSGFA

  REAL ::      POINTX(3),           &! USED TO SPECIFY POINTS TO DRAW BE
       POINTY(3),           &! USED TO SPECIFY POINTS TO DRAW BE
       JUNK(7)              ! CALCULATION ARRAY FOR SFSGFA

!***************************** SUBROUTINE BEGIN ************************

!  INITIALIZE LOOP, BOOLEAN INDICATOR

  DONE = .FALSE.

!  CALCULATE THE WIND VECTOR MAGNITUDE IN KNOTS

  IF ((U .EQ. 0) .AND. (V .EQ. 0)) THEN
     WINDVCT = 1.0
  ELSE
     WINDVCT = SQRT(U**2 + V**2) * 1.94
  END IF

!  SAVE INCOMING USER COORDINATES AND CHANGE BACK TO NORMALIZED COORDINA
!  DOCUMENTATION FOR SET AND GETSET CAN BE FOUND IN NCAR GRAPHICS USER'S
!  GUIDE VERSION 2.00 ON PAGES 49 (GETSET) AND 53 (SET).

  CALL GETSET (FLSV,FRSV,FBSV,FTSV,ULSV,URSV,UBSV,UTSV,LLSV)
  CALL SET    (FLSV,FRSV,FBSV,FTSV, 0.0, 1.0, 0.0, 1.0, 1)

!  DETERMINE WHERE THE BASE OF THE BARB IS IN THE NORMALIZED COORDINATES

  NEWXBASE = (XBASE - ULSV)/(URSV - ULSV)
  NEWYBASE = (YBASE - UBSV)/(UTSV - UBSV)

!  CALCULATE THE X DISTANCE AND Y DISTANCE FROM THE BASE OF THE BARB THA
!  DEFINES THE BARBS TIP (NORMALIZED COORD'S)

  XCOMP = -SC * U * 1.94/WINDVCT
  YCOMP = -SC * V * 1.94/WINDVCT

!  DETERMINE THE ACTUAL LOCATION IN NORMALIZED COORDINATES OF THE BARB'S

  POINTX(1) = NEWXBASE + XCOMP
  POINTY(1) = NEWYBASE + YCOMP

!  DRAW THE BARB SHAFT, DOCUMENTATION FOR THE LINE SUBROUTINE CAN BE FOU
!  IN NCAR GRAPHICS USER'S GUIDE VERSION 2.00 ON PAGE 50

  CALL LINE (NEWXBASE,NEWYBASE,POINTX(1),POINTY(1))

!  DETERMINE THE FEATHER LENGTH

  FETHLENX = 0.3 * YCOMP
  FETHLENY = -0.3 * XCOMP

!  SET THE PLACE KEEPER AND BOOST THE WIND MAGNITUDE

  PK = 0.9
  WINDVCT = WINDVCT + 2.5

!  BEGIN MAKING FEATHERS

10 CONTINUE

!    DRAW A FLAG FOR EVERY 50 KNOTS WIND MAGNITUDE

  IF (WINDVCT .GE. 50.0) THEN

!      DETERMINE THE POSITION OF THE FLAG TIP, POINT_(2)
!      AND DETERMINE POSITION WHERE FLAG BOTTOM MEETS THE SHAFT, POINT_(

     POINTX(2) = POINTX(1) + FETHLENX + 0.0005
     POINTY(2) = POINTY(1) + FETHLENY + 0.0005
     POINTX(3) = PK * XCOMP + NEWXBASE
     POINTY(3) = PK * YCOMP + NEWYBASE

!      DRAW FLAG

     CALL LINE (POINTX(1),POINTY(1),POINTX(2),POINTY(2))
     CALL LINE (POINTX(3),POINTY(3),POINTX(2),POINTY(2))

!      FILL IN FLAG, DOCUMENTATION FOR SFSGFA CAN BE FOUND IN NCAR
!      GRAPHICS GUIDE TO NEW UTILITIES VERSION 3.00 ON PAGE 4-8

     CALL SFSETR ('SP',0.000001)
     CALL SFSGFA (POINTX,POINTY,3,JUNK,5,IJUNK,7,2)

!      REMOVE 50 KNOTS FROM WIND MAGNITUDE (ALREADY DRAWN IN)

     WINDVCT = WINDVCT - 50.0

!      DETERMINE NEW BEGIN POINT FOR NEXT FLAG OR FEATHER

     PK = PK - 0.05
     POINTX(1) = PK * XCOMP + NEWXBASE
     POINTY(1) = PK * YCOMP + NEWYBASE
     PK = PK - 0.1

!    DRAW A FULL FEATHER FOR WIND MAGNITUDE OF EVERY 10 KNOTS

  ELSE IF (WINDVCT .GE. 10.0) THEN

!      CALCULATE POSITION OF FEATHER END

     POINTX(2) = POINTX(1) + FETHLENX + 0.0005
     POINTY(2) = POINTY(1) + FETHLENY + 0.0005

!      DRAW FEATHER

     CALL LINE (POINTX(1),POINTY(1),POINTX(2),POINTY(2))

!      REMOVE 10 KNOTS FROM WIND MAGNITUDE (ALREADY DRAWN IN)

     WINDVCT = WINDVCT - 10.0

!      DETERMINE NEW START POINT FOR NEXT FEATHER OR FLAG

     POINTX(1) = PK * XCOMP + NEWXBASE
     POINTY(1) = PK * YCOMP + NEWYBASE
     PK = PK - 0.1

!    DRAW A HALF FEATHER FOR EVERY 5 KNOTS OF WIND MAGNITUDE

  ELSE IF (WINDVCT .GE. 5.0) THEN

!      CALCULATE POSITION OF TIP OF HALF FEATHER

     POINTX(2) = POINTX(1) + 0.5 * FETHLENX + 0.0005
     POINTY(2) = POINTY(1) + 0.5 * FETHLENY + 0.0005

!      DRAW IN FEATHER

     CALL LINE (POINTX(1),POINTY(1),POINTX(2),POINTY(2))

!      TELL LOOP TO QUIT

     DONE = .TRUE.
  ELSE
     DONE = .TRUE.
  END IF

!  IF THERE IS STILL MORE WIND MAGNITUDE (>= 5 KNOTS) LOOP AGAIN

  IF (.NOT. DONE) GOTO 10

!  RESET USER COORDINATES TO THE INCOMING VALUES

  CALL SET (FLSV,FRSV,FBSV,FTSV,ULSV,URSV,UBSV,UTSV,LLSV)

!****************************** SUBROUTINE END *************************

  RETURN
END SUBROUTINE SKEWT_WNDBARB

FUNCTION  SKEWT_TMR(W, P)
!  TMR(KELVIN),W(GRAMS WATER VAPOR/KILOGRAM DRY AIR),P(MILLIBAR)
  X =  ALOG10 (W * P / (622.+ W))
  SKEWT_TMR = 10.**(.0498646455*X+2.4082965)-7.07475+38.9114*&
       ((10.**( .0915*X ) - 1.2035 )**2 )
  RETURN
END FUNCTION SKEWT_TMR

FUNCTION SKEWT_TSA(OS, P)
!  TSA AND OS(KELVIN),P(MILLIBARS)
!  SIGN(A,B) REPLACES THE ALGEBRAIC SIGN OF A WITH THE SIGN OF B
  A = OS
  TQ = 253.16
  D = 120
  LOOP_1 : DO  I = 1,12
     D = D/2.
!  IF THE TEMPERATURE DIFFERENCE,X, IS SMALL,EXIT THIS LOOP
     X=A*EXP(-2.6518986*SKEWT_W(TQ,P)/TQ)-TQ*((1000./P)**.286)
     IF(ABS(X).LT.0.01)GO TO 2
     TQ = TQ + SIGN(D,X)
  enddo LOOP_1
2 SKEWT_TSA=TQ
  RETURN
END FUNCTION SKEWT_TSA

FUNCTION SKEWT_W(T, P)
!BPR BEGIN
!!  W(GRAMS WATER VAPOR/KILOGRAM DRY AIR ), P(MILLIBAR )
!  SKEWT_W(GRAMS WATER VAPOR/KILOGRAM DRY AIR ), P(MILLIBAR )
!BPR END
  IF (T .GE. 999.) THEN
!BPR BEGIN
!    W = 0.0
     SKEWT_W = 0.0
!BPR END
  ELSE
     X =  SKEWT_ESAT(T)
     SKEWT_W = 621.97 * X / (P - X)
  ENDIF
  RETURN
END FUNCTION SKEWT_W

FUNCTION  SKEWT_ESAT(T)
!  ESAT(MILLIBARS),T(KELVIN)
  PARAMETER (ABZ=273.16)
  TC=T-ABZ
  SKEWT_ESAT=6.1078*EXP((17.2693882*TC)/(TC+237.3))
  RETURN
END FUNCTION SKEWT_ESAT

FUNCTION  SKEWT_OS(T, P)
!  OS AND T (KELVIN) , P (MILLIBARS )
  SKEWT_OS = T*((1000./P)**.286)/(EXP(-2.6518986*SKEWT_W(T,P)/T))
  RETURN
END FUNCTION SKEWT_OS

subroutine skewt_mappings(P, T, X, Y)
!
!  MAPPINGS FROM (P,T) TO CM ON SKEWT
!
  Y = 132.182 - 44.061 * ALOG10(P)
  X = 0.54 * T + 0.90692 * Y

  return
end subroutine skewt_mappings

SUBROUTINE DFCLRS
!
! Define a set of RGB color triples for colors 1 through 15.
! I got this from the Penn State 1-d model plot programs.
  DIMENSION RGBV(3,20)
!
! Define the RGB color triples needed below.
!------------------------------------------------------------------
!                    RED   GREEN   BLUE       #.  Color
!------------------------------------------------------------------
  DATA RGBV / 0.00 , 0.00 , 0.00 ,   &!  1. WHITE
       0.66 , 0.66 , 0.66 ,   &!  2. LIGHT GRAY
       0.40 , 0.40 , 0.40 ,   &!  3. DARK GRAY
       0.00 , 0.00 , 0.00 ,   &!  4. BLACK
       0.20 , 0.56 , 0.80 ,   &!  5. SKY BLUE 
       0.00 , 0.00 , 1.00 ,   &!  6. BLUE
       0.80 , 0.80 , 0.00 ,   &!  7. LIGHT YELLOW
       1.00 , 0.00 , 1.00 ,   &!  8. MAGENTA
       1.00 , 1.00 , 0.00 ,   &!  9. YELLOW
       0.00 , 1.00 , 0.00 ,   &! 10. GREEN
       0.14 , 0.25 , 0.14 ,   &! 11. FOREST GREEN
       0.00 , 1.00 , 1.00 ,   &! 12. CYAN
       0.50 , 0.00 , 1.00 ,   &! 13. PURPLE
       1.00 , 0.00 , 0.38 ,   &! 14. LIGHT RED
       1.00 , 0.50 , 0.00 ,   &! 15. ORANGE
       1.00 , 0.00 , 0.00 ,   &! 16. RED
       0.75 , 0.50 , 1.00 ,   &! 17. MID-BLUE
       0.00 , 0.15 , 0.30 ,   &! 18. DULL MID-BLUE
       0.20 , 0.40 , 0.20 ,   &! 19. BRIGHT FOREST GREEN
       0.60 , 0.30 , 0.00 /    ! 20. DULL ORANGE
!------------------------------------------------------------------
! Define 21 different color indices, for indices 0 through 20.  The
! color corresponding to index 0 is black and the color corresponding
! to index 1 is white.
!
  CALL GSCR (1,0,1.0,1.0,1.0)
!
  LOOP_101 : DO I=1,20
     CALL GSCR (1,I,RGBV(1,I),RGBV(2,I),RGBV(3,I))
  enddo LOOP_101

  RETURN
END SUBROUTINE DFCLRS

end module skewt_module
