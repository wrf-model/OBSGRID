MODULE printout

CONTAINS

!------------------------------------------------------------------------------
   
   SUBROUTINE sample1 ( field , name , kbu )
   
   !  This routine provides sample values of the input 1D field.
   
      IMPLICIT NONE 

      REAL , DIMENSION ( : )                         :: field
      CHARACTER *(*)                                 :: name
      INTEGER                                        :: kbu 
   
      REAL                                           :: fmax , &
                                                        fmin , &
                                                        fmean
      INTEGER                                        :: k
   
      !  Here is the header for this printout.  List the variable
      !  name and some required dashes.
   
      WRITE ( UNIT = * , FMT = ' (//"    VARAIBLE = ",a8) ' ) name
      WRITE ( UNIT = * , &
      FMT = ' ("            Maximum       Minimum         Mean    ") ' )
      WRITE ( UNIT = * , &
      FMT = ' ("             Value         Value         Value    ") ' )
      WRITE ( UNIT = * , &
      FMT = ' ("--------------------------------------------------") ' )
   
      !  Find and/or compute the extremes and mean of the field.
      !  These values are printed out directly from inside this 
      !  routine, not passed back through the argument list.
   
      !  Set the  max/min to reasonable initial values using
      !  FORTRAN intrinsics.
   
      fmax = -1. * 10.**(20)
      fmin =       10.**(20)
      
      !  The mean of each layer comes from a sum, which we
      !  initialize to zero.
 
      fmean = 0.
   
      k_loop : DO k = 1 , kbu
         fmax = max ( field(k) , fmax )
         fmin = min ( field(k) , fmin )
         fmean = fmean + field(k)
      END DO k_loop
      fmean = fmean / real ( kbu )
      
      !  We now know the max/min/mean of this layer, for this variable.
   
      WRITE ( UNIT = * , FMT = ' (8x,        3(2x,g12.5)) ' ) fmax , fmin , fmean
   
      WRITE ( UNIT = * , &
      FMT = ' ("--------------------------------------------------"/) ' )
   
   END SUBROUTINE sample1

!------------------------------------------------------------------------------
   
   SUBROUTINE sample2 ( field , name , iewd , jnsd , istagu , istagv )
   
   !  This routine provides sample values of the input 2D field.
   
      IMPLICIT NONE 

      REAL , DIMENSION ( : , : )                     :: field
      CHARACTER *(*)                                 :: name
      INTEGER                                        :: iewd , &
                                                        jnsd , &
                                                        istagu , istagv
   
      REAL                                           :: fmax , &
                                                        fmin , &
                                                        fmean
      INTEGER                                        :: i , j , k
   
      !  Here is the header for this printout.  List the variable
      !  name and some required dashes.
   
      WRITE ( UNIT = * , FMT = ' (//"    VARAIBLE = ",a8) ' ) name
      WRITE ( UNIT = * , &
      FMT = ' (" K-Level    Maximum       Minimum         Mean    ") ' )
      WRITE ( UNIT = * , &
      FMT = ' (" Index       Value         Value         Value    ") ' )
      WRITE ( UNIT = * , &
      FMT = ' ("--------------------------------------------------") ' )
   
      !  Find and/or compute the extremes and mean of the field
      !  for this single input level.  These values are printed out
      !  directly from inside this routine, not passed back through
      !  the argument list.
   
      !  Set the  max/min to reasonable initial values using
      !  FORTRAN intrinsics.
   
      fmax = -1. * 10.**(20)
      fmin =       10.**(20)
      
      !  The mean of each layer comes from a sum, which we
      !  initialize to zero.
   
      fmean = 0.
   
      i_loop : DO i = 1 , iewd - istagu
         j_loop : DO j = 1 , jnsd - istagv
            fmax = max ( field(j,i) , fmax )
            fmin = min ( field(j,i) , fmin )
            fmean = fmean + field(j,i)
         END DO j_loop
      END DO i_loop
      fmean = fmean / real ( ( iewd - istagu ) * ( jnsd - istagv ) )
      
      !  We now know the max/min/mean of this layer, for this variable.
   
      WRITE ( UNIT = * , FMT = ' (3x,i3 ,2x, 3(2x,g12.5) ) ' ) 1 , fmax , fmin , fmean
   
      WRITE ( UNIT = * , &
      FMT = ' ("--------------------------------------------------"/) ' )
   
   END SUBROUTINE sample2

!------------------------------------------------------------------------------
   
   SUBROUTINE sample3 ( field , name , iewd , jnsd , kbu , istagu , istagv )
   
      IMPLICIT NONE 

   !  This routine provides sample values of the input 3D field.
   
      REAL , DIMENSION ( : , : , : )                 :: field
      CHARACTER *(*)                                 :: name
      INTEGER                                        :: iewd , &
                                                        jnsd , &
                                                        kbu  , &
                                                        istagu , istagv
   
      REAL                                           :: fmax , &
                                                        fmin , &
                                                        fmean
      INTEGER                                        :: i , j , k
   
      !  Here is the header for this printout.  List the variable
      !  name and some required dashes.
   
      WRITE ( UNIT = * , FMT = ' (//"    VARAIBLE = ",a8) ' ) name
      WRITE ( UNIT = * , &
      FMT = ' (" K-Level    Maximum       Minimum         Mean    ") ' )
      WRITE ( UNIT = * , &
      FMT = ' (" Index       Value         Value         Value    ") ' )
      WRITE ( UNIT = * , &
      FMT = ' ("--------------------------------------------------") ' )
   
      !  Find and/or compute the extremes and mean of the field
      !  for each k-level.  These values are printed out
      !  directly from inside this routine, not passed back through
      !  the argument list.
   
      find_extremes : DO k = 1 , kbu
   
         !  Set the  max/min to reasonable initial values using
         !  FORTRAN intrinsics.
   
         fmax = -1. * 10.**(20)
         fmin =       10.**(20)
      
         !  The mean of each layer comes from a sum, which we
         !  initialize to zero.
   
         fmean = 0.
   
         i_loop : DO i = 1 , iewd - istagu
            j_loop : DO j = 1 , jnsd - istagv
               fmax = max ( field(j,i,k) , fmax )
               fmin = min ( field(j,i,k) , fmin )
               fmean = fmean + field(j,i,k)
            END DO j_loop
         END DO i_loop
         fmean = fmean / real ( ( iewd - istagu ) * ( jnsd - istagv ) )
      
         !  We now know the max/min/mean of this layer, for this variable.
   
         WRITE ( UNIT = * , FMT = ' (3x,i3 ,2x, 3(2x,g12.5) ) ' ) k , fmax , fmin , fmean
   
      END DO find_extremes
   
      WRITE ( UNIT = * , &
      FMT = ' ("--------------------------------------------------"/) ' )
   
   END SUBROUTINE sample3

!------------------------------------------------------------------------------

END MODULE printout
