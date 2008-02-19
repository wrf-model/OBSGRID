MODULE ncarg

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE stop_it

!     CALL clsgks
      CALL gdawk(1)
      CALL gclwk(1)
      CALL gdawk(2)
      CALL gclwk(2)
      CALL gclks

   END SUBROUTINE stop_it

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE start_it( fname )

      INTEGER :: idum
      CHARACTER(LEN=*) :: fname
      CHARACTER(LEN=80) :: cdum, fdum

      !  Initialize the NCAR Graphics routines.
      fdum = trim(fname)

!     CALL opngks
      CALL gopks(6,idum)
      CALL gesc(-1391,1,fdum,1,1,cdum)
      CALL gopwk(1,3,1)
      CALL gopwk(2,9,3)
      CALL gacwk(1)
      CALL gacwk(2)
   
      !  Stay within the given boundaries.

      CALL gsclip ( 0 )
   
      !  Set up the colors that will be used for the text, line and polymarkers.
   
      CALL gscr ( 1, 0, 1.00, 1.00, 1.00 )  ! white background
      CALL gscr ( 1, 1, 0.00, 0.00, 0.00 )  ! black foreground
      CALL gscr ( 1, 2, 0.80, 0.00, 0.00 )  ! Red
      CALL gscr ( 1, 3, 0.80, 0.40, 0.00 )  ! Orange
      CALL gscr ( 1, 4, 0.25, 0.00, 0.25 )  ! Violet
      CALL gscr ( 1, 5, 0.10, 0.52, 0.10 )  ! dark-Green
      CALL gscr ( 1, 6, 0.00, 0.50, 1.00 )  ! light Blue (Cyan)
      CALL gscr ( 1, 7, 0.60, 0.60, 0.60 )  ! Medium Grey

   END SUBROUTINE start_it

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE setcl (iqc,idstation)

      !  Set the color based on the quality control flag.

      INTEGER , INTENT(IN) :: iqc
      CHARACTER (LEN= 8) :: idstation

      CALL sflush
      CALL gstxci(1)
      CALL gsplci(1)
      IF (iqc .EQ. 256) THEN   ! vertically interpolated
         CALL gstxci(5)
         CALL gsplci(5)
      ELSE IF (iqc .EQ. 32768) THEN   ! failed errmax
         CALL gstxci(2)
         CALL gsplci(2)
      ELSE IF (iqc .EQ. 512) THEN    ! extrapolated from single level
         CALL gstxci(3)
         CALL gsplci(3)
      ELSE IF (iqc .EQ. 33024) THEN
         CALL gstxci(4)
         CALL gsplci(4)
      ELSE IF (iqc .EQ. 33280) THEN
         CALL gstxci(7)
         CALL gsplci(7)
      ELSE IF (iqc .GE. 2**15) THEN
         CALL gstxci(6)
         CALL gsplci(6)
         WRITE(6,*) 'idstation = ',idstation,' qc = ',iqc
      END IF
   END SUBROUTINE setcl

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE barb ( wd, ws, x0, y0)

      REAL , INTENT(IN) :: wd, ws, x0, y0

      REAL :: brblen , flglen , x0f , y0f , blen , a , wd1 , bdir , xcos , ysin , xe1 , ye1
      REAL :: xcos2 , ysin2 , dlen, clen , xs , xe , ys , ye , five
      INTEGER :: iws , ia , len , nten

      !  Given wd in degrees, ws in knots, build a barb.

      REAL , DIMENSION(3) :: ptx, pty

      REAL , EXTERNAL :: cfux , cfuy , cufx , cufy

      !  If the wind speed is slow, just put a circle there.

      IF (ws .EQ. 0.) THEN

         CALL pwritx (x0, y0,"'KGU'Y",6,12,0,0)

      ELSE

         !  The absolute length of the barb and flag.

         brblen = 0.025
         flglen = brblen * 0.5

         !  Get coords in fractional system

         x0f = cufx(x0)
         y0f = cufy(y0)

         !  Where to start the barb.  In the ole timey days, this was called a
         !  pen down.

         CALL frstpt (x0,y0)

         !  This tells us how many pennants (each is 50 kts).         

         iws = IFIX (( ws + 2.5)/50.)
         blen = (1.0 + 0.2*REAL(iws)) * brblen

         !  The direction of the barb.

         a = 0.
         wd1 = wd + (ATAN(a)) * 180./3.14159265
         bdir = (90. - WD1) * 0.017453
         xcos = COS (bdir)
         ysin = SIN (bdir)
         xe1 = blen * xcos + x0f
         ye1 = blen * ysin + y0f

         !  Now draw the line for the barb.

         call vector (cfux(xe1), cfuy(ye1))

         !  90 degrees from the barb - this may be a northern-hemispheric-centrist chauvanism.

         xcos2 = COS (bdir - 1.570796)
         ysin2 = SIN (bdir - 1.570796)

         !  Draw each pennant on the barb.

         DO ia = 1, iws
            clen = (1.0 + 0.2*REAL(ia+1)) * brblen
            dlen = (1.0 + 0.2*REAL(ia) ) * brblen
            ptx(1) = dlen * xcos + x0f
            pty(1) = dlen * ysin + y0f
            ptx(2) = flglen * xcos2 + ptx(1)
            pty(2) = flglen * ysin2 + pty(1)
            ptx(3) = clen * xcos + x0f
            pty(3) = clen * ysin + y0f
            CALL frstpt (cfux(ptx(1)), cfuy(pty(1)))
            CALL vector (cfux(ptx(2)), cfuy(pty(2)))
            CALL vector (cfux(ptx(3)), cfuy(pty(3)))
            CALL vector (cfux(ptx(1)), cfuy(pty(1)))
         END DO

         !  Subtract out the pennants, how many 10 kt flags are left over.
         
         nten = IFIX (( ws - REAL(iws)*50. + 2.5)/10.)

         !  Plop down each of the 10 kt flags.

         DO ia = 1, nten
            xs = brblen * (1.0 - 0.2*REAL(ia-1)) * xcos + x0f
            ys = brblen * (1.0 - 0.2*REAL(ia-1)) * ysin + y0f
            xe = flglen * xcos2 + xs
            ye = flglen * ysin2 + ys
            CALL frstpt (cfux(xs), cfuy(ys))
            CALL vector (cfux(xe), cfuy(ye))
         END DO

         !  And is there a "half" flag, worth 5 kts, remaining?

         five = ws - REAL(iws)*50. - REAL(nten)*10.

         IF (five .GT. 2.5) THEN
            IF (nten .eq. 0 .AND. iws .EQ. 0) nten = 1
            xs = brblen * (1.0 - 0.2*REAL(nten)) * xcos + x0f
            ys = brblen * (1.0 - 0.2*REAL(nten)) * ysin + y0f
            xe = 0.5 * flglen * xcos2 + xs
            ye = 0.5 * flglen * ysin2 + ys
            CALL frstpt (cfux(xs), cfuy(ys))
            CALL vector (cfux(xe), cfuy(ye))
         ENDIF
      ENDIF

   END SUBROUTINE barb

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE ncarg
