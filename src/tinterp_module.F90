MODULE tinterp

   USE small_header_data
   USE input_data

CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE store_fa ( t , u , v , rh , slp_x , &
                      iew_alloc , jns_alloc , kbu_alloc , num3d , num2d , icount ) 

!  Store the separate arrays back into the all_3d and all_2d arrays.  This is to allow the
!  final analysis to be saved for use in the interpolation schemes for sfc FDDA.

   IMPLICIT NONE

   INTEGER , INTENT(IN) ::  iew_alloc , jns_alloc , kbu_alloc , num3d , num2d , icount

   REAL , INTENT(IN) , DIMENSION(jns_alloc,iew_alloc,kbu_alloc) :: t , u , v , rh
   REAL , INTENT(IN) , DIMENSION(jns_alloc,iew_alloc)           :: slp_x

   INTEGER :: loop_count , tt

   IF ( icount .EQ. 1 ) THEN
      tt = first_time
   ELSE
      tt = second_time
   END IF

   the_3d_search : DO loop_count = 1 , num3d
      IF      ( all_3d(loop_count,tt)%small_header%name(1:8) .EQ. 'TT      ' ) THEN
         all_3d(loop_count,tt)%array = t 
      ELSE IF ( all_3d(loop_count,tt)%small_header%name(1:8) .EQ. 'UU      ' ) THEN
         all_3d(loop_count,tt)%array = u 
      ELSE IF ( all_3d(loop_count,tt)%small_header%name(1:8) .EQ. 'VV      ' ) THEN
         all_3d(loop_count,tt)%array = v 
      ELSE IF ( all_3d(loop_count,tt)%small_header%name(1:8) .EQ. 'RH      ' ) THEN
         all_3d(loop_count,tt)%array = rh
      END IF
   END DO the_3d_search

   the_2d_search : DO loop_count = 1 , num2d
      IF      ( all_2d(loop_count,tt)%small_header%name(1:8) .EQ. 'PMSL    ' ) THEN
         all_2d(loop_count,tt)%array = slp_x * 100.
      END IF
   END DO the_2d_search

END SUBROUTINE store_fa

!------------------------------------------------------------------------------

SUBROUTINE temporal_interp ( t , u , v , uA , vA , uC , vC , h , rh , pres , slp_x , slp_C , snow , pressure , & 
iew_alloc , jns_alloc , kbu_alloc , num3d , num2d , &
icount_fdda , icount_1 , icount_2 ) 

!  This routine temporally interpolates date to the FDDA time from the surrounding
!  traditional analysis time periods.  Only the passed variables (3d and 2d) are
!  interpolated, since those are currrently the only ones on which an objective
!  analysis is performed.  

   IMPLICIT NONE

   INTEGER , INTENT(IN) ::  iew_alloc , jns_alloc , kbu_alloc , num3d , num2d , icount_fdda , icount_1 , icount_2
   REAL , INTENT(IN)  , DIMENSION(kbu_alloc)                     :: pressure

   REAL , INTENT(OUT) , DIMENSION(jns_alloc,iew_alloc,kbu_alloc) :: t , u , v , h , rh
   REAL , INTENT(OUT) , DIMENSION(jns_alloc,iew_alloc,kbu_alloc) :: uA , vA , uC , vC , pres
   REAL , INTENT(OUT) , DIMENSION(jns_alloc,iew_alloc)           :: slp_x , slp_C , snow

   INTEGER :: loop_count , k , kk
   INTEGER , DIMENSION(4) :: k_want = 0
   REAL    , DIMENSION(4) :: p_want = (/ 1001. , 850., 700., 500. /)
   !BPR BEGIN
   INTEGER  :: first_time_weight_int, second_time_weight_int, weight_divide_int
   REAL  :: first_time_weight_real, second_time_weight_real, weight_divide_real
   !BPR END

   !  Get the levels that we are interested in for the SFC FDDA.  We only need the surface for
   !  the objective analysis, but we also need the 850, 700, and 500 hPa levels for diagnosing
   !  the surface pressure.

   k_want(1) = 1
   wanted : DO kk = 2 , 4
      search_press : DO k = 2 , kbu_alloc
         IF ( ABS( pressure(k) - p_want(kk) ) .LT. 1. ) THEN
            k_want(kk) = k
            EXIT search_press
         END IF
      END DO search_press
   END DO wanted

   IF ( ( k_want(2) .EQ. 0 ) .OR. &
        ( k_want(3) .EQ. 0 ) .OR. &
        ( k_want(4) .EQ. 0 ) ) THEN
      PRINT '(A)','Error finding the 850, 700, 500 hPa levels.'
print *,'pressure=',pressure
      STOP 'Toasted_on_levels'
   END IF

   !  Linearly interpolate the 3d data.

   !Weighting of first time (how far away temporally is the second time from the
   !desired time)
   first_time_weight_int = icount_2-icount_fdda
   !Weighting of second time (how far away temporally is the first time from the
   !desired time)
   second_time_weight_int = icount_fdda-icount_1
   !Now check for special case of both times being exactly at the desired time
   !We used integers so that this comparison could be made  
   IF((first_time_weight_int.eq.0).and.(second_time_weight_int.eq.0)) THEN
     !Just use the first time since both times are the same
     first_time_weight_int = 1
   ENDIF
   weight_divide_int = first_time_weight_int + second_time_weight_int
   IF(weight_divide_int.eq.0) THEN
    PRINT *,'ERROR: temporal_interp: About to attempt divide by zero'
    STOP
   ENDIF
   !Now convert from integer to real
   first_time_weight_real = REAL(first_time_weight_int)
   second_time_weight_real = REAL(second_time_weight_int)
   weight_divide_real = REAL(weight_divide_int)
   

   the_3d_search : DO loop_count = 1 , num3d
      IF      ( all_3d(loop_count,second_time)%small_header%name(1:8) .EQ. 'TT      ' ) THEN
         t (1:jns_alloc-1,1:iew_alloc-1,:) = &
         ( first_time_weight_real * all_3d(loop_count,first_time )%array(1:jns_alloc-1,1:iew_alloc-1,:) + &
           second_time_weight_real * all_3d(loop_count,second_time)%array(1:jns_alloc-1,1:iew_alloc-1,:) ) / &
           weight_divide_real 
         t ( jns_alloc,:iew_alloc-1,:) = t ( jns_alloc-1,:iew_alloc-1,:)
         t (:jns_alloc, iew_alloc  ,:) = t (:jns_alloc  , iew_alloc-1,:)
      ELSE IF ( all_3d(loop_count,second_time)%small_header%name(1:8) .EQ. 'UU      ' ) THEN
         u (1:jns_alloc  ,1:iew_alloc  ,:) = &
         ( first_time_weight_real * all_3d(loop_count,first_time )%array(1:jns_alloc  ,1:iew_alloc  ,:) + &
           second_time_weight_real * all_3d(loop_count,second_time)%array(1:jns_alloc  ,1:iew_alloc  ,:) ) / &
           weight_divide_real 
      ELSE IF ( all_3d(loop_count,second_time)%small_header%name(1:8) .EQ. 'VV      ' ) THEN
         v (1:jns_alloc  ,1:iew_alloc  ,:) = &
         ( first_time_weight_real * all_3d(loop_count,first_time )%array(1:jns_alloc  ,1:iew_alloc  ,:) + &
           second_time_weight_real * all_3d(loop_count,second_time)%array(1:jns_alloc  ,1:iew_alloc  ,:) ) / &
           weight_divide_real 
      ELSE IF ( all_3d(loop_count,second_time)%small_header%name(1:8) .EQ. 'UA      ' ) THEN
         uA (1:jns_alloc  ,1:iew_alloc  ,:) = &
         ( first_time_weight_real * all_3d(loop_count,first_time )%array(1:jns_alloc  ,1:iew_alloc  ,:) + &
           second_time_weight_real * all_3d(loop_count,second_time)%array(1:jns_alloc  ,1:iew_alloc  ,:) ) / &
           weight_divide_real 
      ELSE IF ( all_3d(loop_count,second_time)%small_header%name(1:8) .EQ. 'VA      ' ) THEN
         vA (1:jns_alloc  ,1:iew_alloc  ,:) = &
         ( first_time_weight_real * all_3d(loop_count,first_time )%array(1:jns_alloc  ,1:iew_alloc  ,:) + &
           second_time_weight_real * all_3d(loop_count,second_time)%array(1:jns_alloc  ,1:iew_alloc  ,:) ) / &
           weight_divide_real 
      ELSE IF ( all_3d(loop_count,second_time)%small_header%name(1:8) .EQ. 'UC      ' ) THEN
         uC (1:jns_alloc  ,1:iew_alloc  ,:) = &
         ( first_time_weight_real * all_3d(loop_count,first_time )%array(1:jns_alloc  ,1:iew_alloc  ,:) + &
           second_time_weight_real * all_3d(loop_count,second_time)%array(1:jns_alloc  ,1:iew_alloc  ,:) ) / &
           weight_divide_real 
      ELSE IF ( all_3d(loop_count,second_time)%small_header%name(1:8) .EQ. 'VC      ' ) THEN
         vC (1:jns_alloc  ,1:iew_alloc  ,:) = &
         ( first_time_weight_real * all_3d(loop_count,first_time )%array(1:jns_alloc  ,1:iew_alloc  ,:) + &
           second_time_weight_real * all_3d(loop_count,second_time)%array(1:jns_alloc  ,1:iew_alloc  ,:) ) / &
           weight_divide_real 
      ELSE IF ( all_3d(loop_count,second_time)%small_header%name(1:8) .EQ. 'GHT     ' ) THEN
         h (1:jns_alloc  ,1:iew_alloc  ,:) = &
         ( first_time_weight_real * all_3d(loop_count,first_time )%array(1:jns_alloc  ,1:iew_alloc  ,:) + &
           second_time_weight_real * all_3d(loop_count,second_time)%array(1:jns_alloc  ,1:iew_alloc  ,:) ) / &
           weight_divide_real 
         h( jns_alloc,:iew_alloc-1,:) = h( jns_alloc-1,:iew_alloc-1,:)
         h(:jns_alloc, iew_alloc  ,:) = h(:jns_alloc  , iew_alloc-1,:)
      ELSE IF ( all_3d(loop_count,second_time)%small_header%name(1:8) .EQ. 'RH      ' ) THEN
         rh(1:jns_alloc-1,1:iew_alloc-1,:) = &
         ( first_time_weight_real * all_3d(loop_count,first_time )%array(1:jns_alloc-1,1:iew_alloc-1,:) + &
           second_time_weight_real * all_3d(loop_count,second_time)%array(1:jns_alloc-1,1:iew_alloc-1,:) ) / &
           weight_divide_real 
         rh( jns_alloc,:iew_alloc-1,:) = rh( jns_alloc-1,:iew_alloc-1,:)
         rh(:jns_alloc, iew_alloc  ,:) = rh(:jns_alloc  , iew_alloc-1,:)
      ELSE IF ( all_3d(loop_count,second_time)%small_header%name(1:8) .EQ. 'PRES    ' ) THEN
         pres(1:jns_alloc-1,1:iew_alloc-1,:) = &
         ( first_time_weight_real * all_3d(loop_count,first_time )%array(1:jns_alloc-1,1:iew_alloc-1,:) + &
           second_time_weight_real * all_3d(loop_count,second_time)%array(1:jns_alloc-1,1:iew_alloc-1,:) ) / &
           weight_divide_real 
         pres( jns_alloc,:iew_alloc-1,:) = pres( jns_alloc-1,:iew_alloc-1,:)
         pres(:jns_alloc, iew_alloc  ,:) = pres(:jns_alloc  , iew_alloc-1,:)
      END IF
   END DO the_3d_search

   !  Linearly interpolate the 2d data.

   the_2d_search : DO loop_count = 1 , num2d
      IF      ( all_2d(loop_count,second_time)%small_header%name(1:8) .EQ. 'PMSL    ' ) THEN
         slp_x(1:jns_alloc-1,1:iew_alloc-1) = &
         ( first_time_weight_real * all_2d(loop_count,first_time )%array(1:jns_alloc-1,1:iew_alloc-1) + &
           second_time_weight_real * all_2d(loop_count,second_time)%array(1:jns_alloc-1,1:iew_alloc-1) ) / &
           weight_divide_real 
         slp_x( jns_alloc,:iew_alloc-1) = slp_x( jns_alloc-1,:iew_alloc-1)
         slp_x(:jns_alloc, iew_alloc  ) = slp_x(:jns_alloc  , iew_alloc-1)
      ELSE IF ( all_2d(loop_count,second_time)%small_header%name(1:8) .EQ. 'PMSL_C  ' ) THEN
         slp_C(1:jns_alloc-1,1:iew_alloc-1) = &
         ( first_time_weight_real * all_2d(loop_count,first_time )%array(1:jns_alloc-1,1:iew_alloc-1) + &
           second_time_weight_real * all_2d(loop_count,second_time)%array(1:jns_alloc-1,1:iew_alloc-1) ) / &
           weight_divide_real 
         slp_C( jns_alloc,:iew_alloc-1) = slp_C( jns_alloc-1,:iew_alloc-1)
         slp_C(:jns_alloc, iew_alloc  ) = slp_C(:jns_alloc  , iew_alloc-1)
      ELSE IF ( all_2d(loop_count,second_time)%small_header%name(1:8) .EQ. 'SNOW    ' ) THEN
         snow(1:jns_alloc-1,1:iew_alloc-1) = &
         ( first_time_weight_real * all_2d(loop_count,first_time )%array(1:jns_alloc-1,1:iew_alloc-1) + &
           second_time_weight_real * all_2d(loop_count,second_time)%array(1:jns_alloc-1,1:iew_alloc-1) ) / &
           weight_divide_real 
         snow( jns_alloc,:iew_alloc-1) = snow( jns_alloc-1,:iew_alloc-1)
         snow(:jns_alloc, iew_alloc  ) = snow(:jns_alloc  , iew_alloc-1)
      END IF
   END DO the_2d_search

END SUBROUTINE temporal_interp

!------------------------------------------------------------------------------

SUBROUTINE lagtem_assign ( t , u , v , uA , vA , uC , vC , h , rh , pres , slp_x , slp_C , & 
snow , iew_alloc , jns_alloc , kbu_alloc , num3d , num2d , &
icount_fdda , icount_1 , icount_2 ) 

!  This routine assigns the previous time period's data to be this time period's
!  first guess.  Is this cutting edge research or what?

   IMPLICIT NONE

   INTEGER , INTENT(IN) ::  iew_alloc , jns_alloc , kbu_alloc , num3d , num2d , icount_fdda , icount_1 , icount_2

   REAL , INTENT(OUT) , DIMENSION(jns_alloc,iew_alloc,kbu_alloc) :: t , u , v , h , rh
   REAL , INTENT(OUT) , DIMENSION(jns_alloc,iew_alloc,kbu_alloc) :: uA , vA , uC , vC , pres
   REAL , INTENT(OUT) , DIMENSION(jns_alloc,iew_alloc)           :: slp_x , slp_C , snow

   INTEGER :: loop_count , tt

   !  Which time level to assign?  If the intermediate FDDA time is the first one that we would
   !  look at in here (i.e., the data is 00 Z and 12 Z, and we are starting with 03 Z in this
   !  routine call), then we assign the data from the "first" analysis time level (00 Z).  If we 
   !  are processing subsequent times upto but not including the final time between these two
   !  analysis times (i.e., same 00 Z and 12 Z, so we do work with 06 Z and 09 Z, but not 12 Z),
   !  the first guess is the previous lagtem'ed datas final analysis.  If we are at the final
   !  time between these two analysis periods (same 00 Z and 12 Z, and we are doing 12 Z sfc
   !  FDDA), then the first guess we use is the final analysis from the second analysis period.

   !  This is the first FDDA loop between these two periods.  Simply assign the data from the
   !  final analysis ffrom the "first" analysis time period.  This is the 03 Z example between 00 Z and 12 Z.

   IF ( icount_fdda - icount_1 .EQ. 1 ) THEN

      !  Search the 3d data.
   
      the_3d_search_a : DO loop_count = 1 , num3d
         IF      ( all_3d(loop_count,first_time)%small_header%name(1:8) .EQ. 'TT      ' ) THEN
            t = all_3d(loop_count,first_time )%array
            t ( jns_alloc,:iew_alloc-1,:) = t ( jns_alloc-1,:iew_alloc-1,:)
            t (:jns_alloc, iew_alloc  ,:) = t (:jns_alloc  , iew_alloc-1,:)
         ELSE IF ( all_3d(loop_count,first_time)%small_header%name(1:8) .EQ. 'UU      ' ) THEN
            u = all_3d(loop_count,first_time )%array
         ELSE IF ( all_3d(loop_count,first_time)%small_header%name(1:8) .EQ. 'VV      ' ) THEN
            v = all_3d(loop_count,first_time )%array
         ELSE IF ( all_3d(loop_count,first_time)%small_header%name(1:8) .EQ. 'UA      ' ) THEN
            uA = all_3d(loop_count,first_time )%array
         ELSE IF ( all_3d(loop_count,first_time)%small_header%name(1:8) .EQ. 'VA      ' ) THEN
            vA = all_3d(loop_count,first_time )%array
         ELSE IF ( all_3d(loop_count,first_time)%small_header%name(1:8) .EQ. 'UC      ' ) THEN
            uC = all_3d(loop_count,first_time )%array
         ELSE IF ( all_3d(loop_count,first_time)%small_header%name(1:8) .EQ. 'VC      ' ) THEN
            vC = all_3d(loop_count,first_time )%array
         ELSE IF ( all_3d(loop_count,first_time)%small_header%name(1:8) .EQ. 'GHT     ' ) THEN
            h = all_3d(loop_count,first_time )%array
            h( jns_alloc,:iew_alloc-1,:) = h( jns_alloc-1,:iew_alloc-1,:)
            h(:jns_alloc, iew_alloc  ,:) = h(:jns_alloc  , iew_alloc-1,:)
         ELSE IF ( all_3d(loop_count,first_time)%small_header%name(1:8) .EQ. 'RH      ' ) THEN
            rh= all_3d(loop_count,first_time )%array
            rh( jns_alloc,:iew_alloc-1,:) = rh( jns_alloc-1,:iew_alloc-1,:)
            rh(:jns_alloc, iew_alloc  ,:) = rh(:jns_alloc  , iew_alloc-1,:)
         ELSE IF ( all_3d(loop_count,first_time)%small_header%name(1:8) .EQ. 'PRES    ' ) THEN
            pres= all_3d(loop_count,first_time )%array
            pres( jns_alloc,:iew_alloc-1,:) = pres( jns_alloc-1,:iew_alloc-1,:)
            pres(:jns_alloc, iew_alloc  ,:) = pres(:jns_alloc  , iew_alloc-1,:)
         END IF
      END DO the_3d_search_a
   
      !  Search the 2d data.
   
      the_2d_search_a : DO loop_count = 1 , num2d
         IF      ( all_2d(loop_count,first_time)%small_header%name(1:8) .EQ. 'PMSL    ' ) THEN
            slp_x = all_2d(loop_count,first_time )%array
            slp_x( jns_alloc,:iew_alloc-1) = slp_x( jns_alloc-1,:iew_alloc-1)
            slp_x(:jns_alloc, iew_alloc  ) = slp_x(:jns_alloc  , iew_alloc-1)
         ELSE IF ( all_2d(loop_count,first_time)%small_header%name(1:8) .EQ. 'PMSL_C  ' ) THEN
            slp_C = all_2d(loop_count,first_time )%array
            slp_C( jns_alloc,:iew_alloc-1) = slp_C( jns_alloc-1,:iew_alloc-1)
            slp_C(:jns_alloc, iew_alloc  ) = slp_C(:jns_alloc  , iew_alloc-1)
         ELSE IF ( all_2d(loop_count,first_time)%small_header%name(1:8) .EQ. 'SNOW    ' ) THEN
            snow = all_2d(loop_count,first_time )%array
            snow( jns_alloc,:iew_alloc-1) = snow( jns_alloc-1,:iew_alloc-1)
            snow(:jns_alloc, iew_alloc  ) = snow(:jns_alloc  , iew_alloc-1)
         END IF
      END DO the_2d_search_a

   !  This is a tweener time!  The first guess is the final analysis of the last time's first
   !  guess, and guess what, the data is already loaded in the right place.  No way?  WAY!

   ELSE IF ( ( icount_fdda .GT. icount_1 + 1 ) .AND. ( icount_fdda .LT. icount_2 ) ) THEN
!print *,'We''re in the tweener section, man, and we ain''t doin'' nuthin'

   !  This is the last FDDA loop between these two periods.  Simply assign the data from the
   !  final analysis from the "second" analysis time period.  This is the 12 Z example between 00 Z and 12 Z.

   ELSE IF ( icount_fdda .EQ. icount_2 ) THEN

      !  Search the 3d data.
   
      the_3d_search_b : DO loop_count = 1 , num3d
         IF      ( all_3d(loop_count,second_time)%small_header%name(1:8) .EQ. 'TT      ' ) THEN
            t = all_3d(loop_count,second_time )%array
            t ( jns_alloc,:iew_alloc-1,:) = t ( jns_alloc-1,:iew_alloc-1,:)
            t (:jns_alloc, iew_alloc  ,:) = t (:jns_alloc  , iew_alloc-1,:)
         ELSE IF ( all_3d(loop_count,second_time)%small_header%name(1:8) .EQ. 'UU      ' ) THEN
            u = all_3d(loop_count,second_time )%array
         ELSE IF ( all_3d(loop_count,second_time)%small_header%name(1:8) .EQ. 'VV      ' ) THEN
            v = all_3d(loop_count,second_time )%array
         ELSE IF ( all_3d(loop_count,second_time)%small_header%name(1:8) .EQ. 'UA      ' ) THEN
            uA = all_3d(loop_count,second_time )%array
         ELSE IF ( all_3d(loop_count,second_time)%small_header%name(1:8) .EQ. 'VA      ' ) THEN
            vA = all_3d(loop_count,second_time )%array
         ELSE IF ( all_3d(loop_count,second_time)%small_header%name(1:8) .EQ. 'UC      ' ) THEN
            uC = all_3d(loop_count,second_time )%array
         ELSE IF ( all_3d(loop_count,second_time)%small_header%name(1:8) .EQ. 'VC      ' ) THEN
            vC = all_3d(loop_count,second_time )%array
         ELSE IF ( all_3d(loop_count,second_time)%small_header%name(1:8) .EQ. 'GHT     ' ) THEN
            h = all_3d(loop_count,second_time )%array
            h( jns_alloc,:iew_alloc-1,:) = h( jns_alloc-1,:iew_alloc-1,:)
            h(:jns_alloc, iew_alloc  ,:) = h(:jns_alloc  , iew_alloc-1,:)
         ELSE IF ( all_3d(loop_count,second_time)%small_header%name(1:8) .EQ. 'RH      ' ) THEN
            rh= all_3d(loop_count,second_time )%array
            rh( jns_alloc,:iew_alloc-1,:) = rh( jns_alloc-1,:iew_alloc-1,:)
            rh(:jns_alloc, iew_alloc  ,:) = rh(:jns_alloc  , iew_alloc-1,:)
         ELSE IF ( all_3d(loop_count,second_time)%small_header%name(1:8) .EQ. 'PRES    ' ) THEN
            pres= all_3d(loop_count,second_time )%array
            pres( jns_alloc,:iew_alloc-1,:) = pres( jns_alloc-1,:iew_alloc-1,:)
            pres(:jns_alloc, iew_alloc  ,:) = pres(:jns_alloc  , iew_alloc-1,:)
         END IF
      END DO the_3d_search_b
   
      !  Search the 2d data.
   
      the_2d_search_b : DO loop_count = 1 , num2d
         IF      ( all_2d(loop_count,second_time)%small_header%name(1:8) .EQ. 'PMSL    ' ) THEN
            slp_x = all_2d(loop_count,second_time )%array
            slp_x( jns_alloc,:iew_alloc-1) = slp_x( jns_alloc-1,:iew_alloc-1)
            slp_x(:jns_alloc, iew_alloc  ) = slp_x(:jns_alloc  , iew_alloc-1)
         ELSE IF ( all_2d(loop_count,second_time)%small_header%name(1:8) .EQ. 'PMSL_C  ' ) THEN
            slp_C = all_2d(loop_count,second_time )%array
            slp_C( jns_alloc,:iew_alloc-1) = slp_C( jns_alloc-1,:iew_alloc-1)
            slp_C(:jns_alloc, iew_alloc  ) = slp_C(:jns_alloc  , iew_alloc-1)
         ELSE IF ( all_2d(loop_count,second_time)%small_header%name(1:8) .EQ. 'SNOW    ' ) THEN
            snow = all_2d(loop_count,second_time )%array
            snow( jns_alloc,:iew_alloc-1) = snow( jns_alloc-1,:iew_alloc-1)
            snow(:jns_alloc, iew_alloc  ) = snow(:jns_alloc  , iew_alloc-1)
         END IF
      END DO the_2d_search_b

   END IF

END SUBROUTINE lagtem_assign

!------------------------------------------------------------------------------

END MODULE tinterp

