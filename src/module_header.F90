MODULE header1

   !  Ye old big header.

   INTEGER            , DIMENSION(50,20) :: bhi
   REAL               , DIMENSION(20,20) :: bhr

   CONTAINS

   SUBROUTINE read_header ( imax , jmax )

      IMPLICIT NONE

      INTEGER , INTENT(OUT) :: imax , jmax
      INTEGER :: ier , iflag

      !  The only MM5-ish file we need is something that has a big header,
      !  and just enough info to get a map generated.
     
      OPEN ( UNIT=99 , FILE='LITTLE_R_DOMAIN1' , STATUS='OLD' , &
             ACCESS='SEQUENTIAL' , FORM='UNFORMATTED' )
   
      READ(99, iostat=ier) iflag
      IF(ier .ne. 0) THEN
         PRINT '(A)','Error READing big header flag, link something to fort.99.'
         STOP 'Link_anything_to_99'
      END IF

      IF (iflag .EQ. 0) THEN
         READ(99) bhi, bhr
         IF(ier .ne. 0) THEN
            PRINT '(A)','Error READing big header, meet me half way and link an MM5 file to fort.99'
            STOP 'Link_MM5_thing_to_99'
         END IF
   
         !  Pull the horizontal domain size from the header info.
   
         imax = bhi(16,1)
         jmax = bhi(17,1)
         IF ((bhi(15,1).EQ.0) .and. (bhi(8,1).eq.1) .and. (bhi(1,1).lt.3)) THEN
            imax = bhi(9,1)
            jmax = bhi(10,1)
         END IF
      ELSE
         PRINT '(A)','Trouble READing V3 file, as in - I don''t think it is one.'
         STOP 'V3_is_doubtful'
      END IF

   END SUBROUTINE read_header

END MODULE header1
