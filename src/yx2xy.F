
!------------------------------------------------------------------------------

SUBROUTINE yx2xy ( yx , ifirst , isecond , xy )

!  This subroutine performs a simple 2D transpose, which is required
!  for the usual (x,y) to (i,j) orientation switches.  The variable
!  "yx" is always the input gridded field.  The array dimensions
!  are the first and second for the input array (ifirst and isecond,
!  respectively).

   INTEGER , INTENT ( IN )                             :: ifirst , & 
                                                          isecond
   REAL , DIMENSION (ifirst ,isecond) , INTENT ( IN  ) :: yx
   REAL , DIMENSION (isecond,ifirst ) , INTENT ( OUT ) :: xy

   INTEGER                                             :: i , j 

   DO i = 1 , isecond
      DO j = 1 , ifirst
         xy(i,j) = yx(j,i)
      END DO
   END DO

END SUBROUTINE yx2xy
