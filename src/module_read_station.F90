MODULE read_station

CONTAINS

   SUBROUTINE rdsta (u, v, t, rh, slp, ipress, opress, sid,  xc, yc, wqc, tqc, rhqc, pqc, nsta, num, ier, filename)

      INTEGER , INTENT(IN) :: nsta
      INTEGER , INTENT(OUT) :: num , ier
      REAL , INTENT(OUT) , DIMENSION(nsta) :: u, v, t, slp, rh, xc, yc
      INTEGER , INTENT(OUT) , DIMENSION(nsta) :: wqc, tqc, rhqc, pqc
      INTEGER , INTENT(INOUT) :: opress , ipress
      CHARACTER (LEN= *) , INTENT(IN) :: filename
      CHARACTER (LEN= 8) , DIMENSION(nsta) :: sid
      CHARACTER (LEN=80) , SAVE :: fmt
      CHARACTER (LEN= 8) :: name , idstation
      CHARACTER (LEN= 1) :: char1

      LOGICAL :: got_a_file
      INTEGER :: ok

      INTEGER :: nobs , loop , iobnum , iqc , i
      REAL :: value , diff , x , y 

      !  First time in, or not?  Makes a difference.

      IF (opress .lt. 0.) THEN
         INQUIRE ( FILE = filename , EXIST = got_a_file )
  
         !  And the worst case scenario is, that file ain't around.

         IF ( .NOT. got_a_file ) THEN
            PRINT '(A,A,A)','The file, ',TRIM(filename),' does not exist.'
            STOP 'NO_plotobs_out_file_found'
         END IF

         !  First time in, open the file.

         OPEN ( UNIT    = 9            , &
                FILE    = filename     , &
                FORM    = 'FORMATTED'  , &
                ACCESS  = 'SEQUENTIAL' , &
                STATUS  = 'OLD'          )

         !  Bypass the beginning print in the obs file.

         READ (9,'(a80)') fmt
         !PRINT '(A,A)','format statement is = ',fmt

      ELSE

         opress = ipress

      END IF

      !  Initialize the data and flags to undefined values.

      u    = -999.
      v    = -999.
      t    = -999.
      rh   = -999.
      slp  = -999.
      xc   = 0.
      yc   = 0.
      wqc  = 0
      tqc  = 0
      rhqc = 0
      pqc  = 0
      num  = 0
      ier  = 0

      do_a_level : DO

         !  A bit of header info in front of each variable.

         READ (9,'(A)',IOSTAT=ok) char1

         !  Can we continue?  This is where we would find the end of file.

         IF ( ok .NE. 0 ) THEN
            ier = 1
            EXIT do_a_level
         END IF

         READ (9,'(I8)') nobs
         READ (9,'(A)') char1
         READ (9,'(A)') char1

         !  Do all of the obs for this variable, at this level.

         do_a_variable : DO loop = 1, nobs
            READ(9,FMT=fmt)name,ipress,iobnum,idstation,value,diff,x,y,iqc
            !WRITE (*,FMT=fmt) name,ipress,iobnum,idstation,value,diff,x,y,iqc
            !READ(9,FMT=fmt)name,ipress,iobnum,idstation,value,x,y,iqc
            !WRITE (*,FMT=fmt) name,ipress,iobnum,idstation,value,x,y,iqc
! READ (9,1212) name,ipress,iobnum,idstation,value,diff,x,y,iqc
! WRITE(6,1212) name,ipress,iobnum,idstation,value,diff,x,y,iqc
! 1212  FORMAT( 3X,A8,3X,I6,3X,I5,3X,A8,3X,2(G13.6,3X),2(F7.2,3X),I7 )

            !  Is this the right level?  If not, jump out of the level loop, which means
            !  exit the routine.

            IF (opress .gt. 0. .and. ipress .ne. opress) THEN
               BACKSPACE (9)
               BACKSPACE (9)
               BACKSPACE (9)
               BACKSPACE (9)
               BACKSPACE (9)
               EXIT do_a_level
            END IF

            opress = ipress
            IF ( name .EQ. 'UU      ') THEN
               u(loop) = value
               wqc(loop) = iqc
               xc(loop) = x
               yc(loop) = y
               sid(loop) = idstation
               num = num + 1
            ELSE IF (name .EQ. 'VV      ') THEN
               DO i = 1, num
                  IF (sid(i) .EQ. idstation .and. x .EQ. xc(i) .and. y .EQ. yc(i)) THEN
                     v(i) = value
                     CYCLE do_a_variable
                  END IF
               END DO
               num = num + 1   ! didn't find a matching ob
               v(num) = value
               wqc(num) = iqc
               xc(num) = x
               yc(num) = y
               sid(num) = idstation
            ELSE IF (name .EQ. 'TT      ') THEN
               DO i = 1, num
                  IF (sid(i) .EQ. idstation .and. x .EQ. xc(i) .and. y .EQ. yc(i)) THEN
                     t(i) = value
                     tqc(i) = iqc
   !  WRITE(6,*) 'found match sid = ',sid(i)
                     CYCLE do_a_variable
                  END IF
               END DO
   !  WRITE(6,*) 'didnt find sid = ',sid(loop)
               num = num + 1   ! didn't find a matching ob
               t(num) = value
               tqc(num) = iqc
               xc(num) = x
               yc(num) = y
               sid(num) = idstation
            ELSE IF (name .EQ. 'RH      ') THEN
               DO i = 1, num
                  IF (sid(i) .EQ. idstation .and. x .EQ. xc(i) .and. y .EQ. yc(i)) THEN
                     rh(i) = value
                     rhqc(i) = iqc
                     CYCLE do_a_variable
                  END IF
               END DO
               num = num + 1   ! didn't find a matching ob
   ! WRITE(6,*) 'num = ',num,' value = ',value
               rh(num) = value
               rhqc(num) = iqc
               xc(num) = x
               yc(num) = y
               sid(num) = idstation
            ELSE IF (name .EQ. 'PMSL') THEN
               DO i = 1, num
                  IF (sid(i) .EQ. idstation .and. x .EQ. xc(i) .and. y .EQ. yc(i)) THEN
                     slp(i) = value
                     pqc(i) = iqc
                     CYCLE do_a_variable
                  END IF
               END DO
               num = num + 1   ! didn't find a matching ob
               slp(num) = value
               pqc(num) = iqc
               xc(num) = x
               yc(num) = y
               sid(num) = idstation
            END IF
         END DO do_a_variable

      END DO do_a_level

   END SUBROUTINE rdsta

END MODULE read_station
