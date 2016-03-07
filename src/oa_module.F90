MODULE obj_analysis

   REAL , PARAMETER , PRIVATE :: pi = 3.14159265358
   REAL , PARAMETER , PRIVATE :: twopi = 2. * 3.14159265358

CONTAINS

!
!---------------------------------------------------------------------------

SUBROUTINE clean_rh ( rh , iew , jns , rh_min , rh_max )

   INTEGER                        :: iew , jns
   REAL , DIMENSION ( iew , jns ) :: rh
   REAL                           :: rh_min , rh_max

   WHERE ( rh .GT. rh_max ) rh = rh_max
   WHERE ( rh .LT. rh_min ) rh = rh_min

END SUBROUTINE clean_rh

!
!---------------------------------------------------------------------------

SUBROUTINE cressman ( obs_value , obs , xob , yob , numobs , &
gridded , iew , jns , &
crsdot , name , radius_influence , &
dxd , u_banana , v_banana , pressure , passes , smooth_type , lat_center , &
!BPR BEGIN
!use_first_guess ) 
use_first_guess , fg_value, scale_cressman_rh_decreases ) 
!BPR END

   ! arguments:
   !    obs_value: difference between 1st guess and obs
   !               ( need any duplicated obs removed from this array )
!BPR BEGIN
   !               ACTUALLY obs_value appears to be a value NOT a difference
!BPR END
   !    obs:       difference of obs_value and the interpolated first_guess value
   !    numobs:    number of station obs
   !    xob,yob:   x, y locations of station obs
   !    gridded:   background/first field as input, final analysis as output
   !    iew,ins:   1st and 2nd dimensions of array 'gridded'
   !    name:      variable name, T, RH, p, U, or V
   !    dxd        grid distance of this domain (m)
   !    u_banana   first guess field of U for banana scheme
   !    v_banana   first guess field of V for banana scheme
   !    pressure   pressure of this level (Pa, 1001 signifies surface)
   !    use_first_guess T/F to include first-guess in the objective analysis
!BPR BEGIN
   !    fg_value: Model value at the ob location
!BPR END
    
   IMPLICIT NONE

   INTEGER, INTENT ( IN )                     :: numobs
   INTEGER, INTENT ( IN )                     :: iew, jns
   REAL,    INTENT ( IN ), DIMENSION ( numobs )      :: obs_value , obs, xob, yob
!BPR BEGIN
   REAL,    INTENT ( IN ), DIMENSION ( numobs )      :: fg_value 
!BPR END
   REAL,    INTENT ( INOUT ), DIMENSION( iew , jns ) :: gridded
   INTEGER, INTENT ( IN )                     :: radius_influence
   CHARACTER (8), INTENT ( IN )               :: name
   REAL ,   INTENT ( IN ) , DIMENSION( iew , jns ) :: u_banana , v_banana
   REAL ,   INTENT ( IN )                     :: dxd , pressure
   INTEGER                                    :: crsdot
   INTEGER , INTENT( IN )                     :: passes , smooth_type
   REAL ,   INTENT ( IN )                     :: lat_center
   LOGICAL , INTENT ( IN )                    :: use_first_guess

   REAL, DIMENSION ( iew,jns )                :: denominator,    &
                                                 numerator ,     &
                                                 perturbation
   REAL                                       :: dxob, dyob, distance , &
                                                 weight
   INTEGER                                    :: iob, job
   INTEGER                                    :: numpts,         &
                                                 i,  j , n ,     &
                                                 iew_min  ,      &
                                                 jns_min  ,      &
                                                 iew_max  ,      &
                                                 jns_max
   !BPR BEGIN
   !CHARACTER (LEN=8)                          :: choice
   REAL                                       :: scale_factor
   REAL                                       :: real_crsdot_divide_2, real_i, real_j
   REAL                                       :: distance_squared
   INTEGER                                    :: radius_influence_squared
   LOGICAL                                    :: scale_cressman_rh_decreases
   !BPR END


   INCLUDE 'error.inc'

   INTERFACE 
      INCLUDE 'error.int'
   END INTERFACE

   !  Initialize the Cressman weights (numerator and denominator) and
   !  the total perturbation field to zero.

   numerator    = 0
   denominator  = 0
   perturbation = 0

   !BPR BEGIN
   !Calculate this here so we do not have to do it many, many times below
   real_crsdot_divide_2 = REAL(crsdot)/2.
   radius_influence_squared = radius_influence * radius_influence
   !BPR END

   !  Loop over each observation.

   obs_loop : DO n = 1 , numobs

      !  Every observation has a first guess wind field location that will be specified
      !  as an objective analysis that is circular, elliptical or banana-shaped.

      !BPR BEGIN
      !choice = 'undefine'
      !BPR END

      !  To use the observation, it must be far enough inside our domain
      !  (using the cross/dot configuration for the staggered variables).  If
      !  it is not, we simply zoom back to the top of this loop.

      !BPR BEGIN
      !IF ( ( xob(n) .LE.   1 + REAL(crsdot)/2. ) .OR. & 
      !     ( yob(n) .LE.   1 + REAL(crsdot)/2. ) .OR. & 
      !     ( xob(n) .GE. iew - REAL(crsdot)/2. ) .OR. & 
      !     ( yob(n) .GE. jns - REAL(crsdot)/2. ) ) THEN
      IF ( ( xob(n) .LE.   1 + real_crsdot_divide_2 ) .OR. & 
           ( yob(n) .LE.   1 + real_crsdot_divide_2 ) .OR. & 
           ( xob(n) .GE. iew - real_crsdot_divide_2 ) .OR. & 
           ( yob(n) .GE. jns - real_crsdot_divide_2 ) ) THEN
      !BPR END
         CYCLE obs_loop
      END IF

      !  Compute the window in which this observation has influence.  For
      !  some points inside this window the observation will have zero 
      !  influence, but this allows a direct method rather than a search.

      !BPR BEGIN
      !iew_min = MAX ( NINT ( (xob(n)-REAL(crsdot)/2.) ) - 4*radius_influence -1 ,            1 )
      !jns_min = MAX ( NINT ( (yob(n)-REAL(crsdot)/2.) ) - 4*radius_influence -1 ,            1 )
      !iew_max = MIN ( NINT ( (xob(n)-REAL(crsdot)/2.) ) + 4*radius_influence +1 , iew - crsdot ) 
      !jns_max = MIN ( NINT ( (yob(n)-REAL(crsdot)/2.) ) + 4*radius_influence +1 , jns - crsdot )
      iew_min = MAX ( NINT ( (xob(n)-real_crsdot_divide_2) ) - 4*radius_influence -1 ,            1 )
      jns_min = MAX ( NINT ( (yob(n)-real_crsdot_divide_2) ) - 4*radius_influence -1 ,            1 )
      iew_max = MIN ( NINT ( (xob(n)-real_crsdot_divide_2) ) + 4*radius_influence +1 , iew - crsdot ) 
      jns_max = MIN ( NINT ( (yob(n)-real_crsdot_divide_2) ) + 4*radius_influence +1 , jns - crsdot )
      !BPR END

      !  For the grid points in the window, compute this observation's
      !  influence.

      !BPR BEGIN
      !Significantly edited code for compuatational efficiency.  The exact same
      !algorithms are used, they are just restructured to speed things up since
      !we want to minimize the code called on the innermost loop where we are
      !looping over the x and y direction for every observation.

      !  Distance from each of the windowed grid points to the observation.  The distance is in 
      !  grid point units.  If this is a horizontal wind field or a moisture field, we can use the
      !  banana scheme, which elongates the acceptable distance from an observation and also is
      !  able to handle first-guess wind curvature.

      IF ( ( name(1:8) .EQ. 'UU      ' ) .OR. &
           ( name(1:8) .EQ. 'VV      ' ) .OR. &
           ( name(1:8) .EQ. 'RH      ' )  ) THEN
       CALL dist_uu_vv_rh ( iew_min, iew_max, jns_min, jns_max, iew , jns , real_crsdot_divide_2, name(1:8) , &
                   obs(n), xob(n) , yob(n) , &
                   dxd , pressure , u_banana , v_banana , radius_influence , &
                   lat_center , gridded, fg_value(n), use_first_guess, numerator , denominator, & !BPR BEGIN
                   scale_cressman_rh_decreases )
!                  lat_center , gridded, fg_value(n), use_first_guess, numerator , denominator  ) !BPR END
      ELSE

       CALL dist_not_uu_vv_rh ( iew_min, iew_max, jns_min, jns_max, real_crsdot_divide_2,  name(1:8) , &
                  obs(n), xob(n) , yob(n) ,  radius_influence , numerator, denominator )

      END IF
      !BPR END


   END DO obs_loop

   !  We have looped over all of the observations, and have computed the required
   !  sums for the objective analysis for each grid point.  Compute the grid
   !  point perturbations from these values.

   WHERE ( denominator .GT. 0 ) 
      perturbation = numerator / denominator
   ENDWHERE

   !  Which smoother will we use? Smooth the final analysis with a 5-point stencil 
   !  (1-2-1 on the lateral boundaries), or the traditional smoother-desmoother.
   !  Both accept a 2D field with the same arguments.

   IF      ( smooth_type .EQ. 1 ) THEN
      CALL smooth_5            ( perturbation , iew, jns, passes , crsdot ) 
   ELSE IF ( smooth_type .EQ. 2 ) THEN
      CALL smoother_desmoother ( perturbation , iew, jns, passes , crsdot )
   END IF

   !  Add the first guess to the perturbations from the observations at
   !  each (i,j).

   IF ( use_first_guess ) THEN
      gridded = gridded + perturbation
   ELSE
      gridded = perturbation
   END IF

END SUBROUTINE cressman

!
!---------------------------------------------------------------------------

!  This routine computes the distance from an observation to a grid point.
!  The Euclidean distance is modified to allow for elongation due to
!  strong stream-wise flow, as well as for curvature.
!  The distance is then used to update the denominator and numerator used to 
!  determine the grid point perturbations for the objective analysis

SUBROUTINE dist_uu_vv_rh ( iew_min, iew_max, jns_min, jns_max, iew , jns , real_crsdot_divide_2,  name , &
                  obs_n, xob , yob , dxd , &
                  pressure , u_banana , v_banana , radius_influence , &
                  lat_center , gridded, fg_value_n, use_first_guess, numerator_pert, denominator_pert, & !BPR BEGIN
                  scale_cressman_rh_decreases )
!                 lat_center , gridded, fg_value_n, use_first_guess, numerator_pert, denominator_pert ) !BPR END
  
   IMPLICIT NONE
 
   !  Input variables.

   INTEGER , INTENT(IN)                      :: iew_min  , &
                                                iew_max  , &
                                                jns_min  , &
                                                jns_max 
   INTEGER , INTENT(IN)                      :: iew      , &
                                                jns    
   REAL    , INTENT(IN)                      :: real_crsdot_divide_2
   CHARACTER (LEN=8), INTENT(IN)             :: name
   REAL    , INTENT(IN)                      :: obs_n    , &
                                                xob      , &
                                                yob      , &
                                                dxd      , &
                                                pressure
   REAL    , INTENT(IN) , DIMENSION(:,:)     :: u_banana , &
                                                v_banana
   INTEGER    , INTENT(IN)                   :: radius_influence
   REAL    , INTENT(IN)                      :: lat_center
   REAL,    INTENT (IN), DIMENSION(:,:)      :: gridded
   REAL,    INTENT (IN)                      :: fg_value_n
   LOGICAL , INTENT (IN)                     :: use_first_guess
   LOGICAL , INTENT (IN)                     :: scale_cressman_rh_decreases !BPR BEGIN/END
   !  Input/Output variables.
   REAL    , INTENT(INOUT), DIMENSION (:,:)  :: numerator_pert, &
                                                denominator_pert

   !  Local variables.
   REAL :: real_i , real_j  
   REAL :: radius_influence_squared

   REAL :: numerator_curvature, denominator_curvature
   REAL :: radius_of_curvature , k
   REAL :: u_ob , dudx , speed_ob , dydx
   REAL :: v_ob , dvdx
   REAL :: x1 , x2 , y1 , y2 , d12
   INTEGER :: i1 , i2 , j1 , j2
   INTEGER :: ia, ja, ib, jb, ic, jc, id, jd
   REAL :: ua, va, ub, vb, uc, vc, ud, vd, vor
   INTEGER :: ii,jj
   REAL :: distance_squared
   

   INTEGER :: i_center , j_center

   REAL :: theta_ob , theta_ij
   
   REAL :: rij
   
   REAL :: vcrit , elongation

   REAL :: x , y , theta_diff

   LOGICAL :: USE_CIRCLE, USE_BANANA, USE_ELLIPSE
   LOGICAL :: do_rh_scaling
 
   REAL :: one_divided_by_pi, one_divided_by_twopi
   REAL :: weight
   REAL :: scale_factor

   USE_CIRCLE = .FALSE.
   USE_BANANA = .FALSE.
   USE_ELLIPSE = .FALSE.
   radius_influence_squared = REAL(radius_influence * radius_influence)

   !  The critical speed to decide whether or not to do a simple
   !  circular Cressman or one of the anisotropic schemes.  If the 
   !  observation speed is < the critical speed, do a circular
   !  Cressman distance computation.  At 1000 hPa, the critical speed
   !  is 5 m/s; this goes linearly to 15 m/s at 500 hPa.

   IF ( pressure < 500. ) THEN
      vcrit = 15.
   ELSE
      vcrit = 25. - 0.02 * pressure
   END IF

   !  Compute the radius of curvature for the streamline passing through
   !  this observation location.  This is the information available from 
   !  eqn 4.16 in Manning and Haagenson, 1994 (MH94).  Note that MH94 eqn 4.16
   !  is different from that computation given below.

   !  dy/dx = slope of streamline through grid point, which is v/u
   !  d^2y/dx^2 = rate of change of slope of streamline through grid point, which
   !              is d(v/u)/dx if taken along stream line trajectory, which is
   !              u*dv/dx - v*du/dx
   !              -----------------
   !                    u^2
   !  k = curvature = abs(dy/dx) / ( 1 + (d^2y/dx^2)^2 ) ^(3/2)
   !  radius of curvature = 1 / k

   !  Values of u and v nearest the observation point in the first-guess data and the
   !  respective speed at this point.

   u_ob = u_banana(NINT(xob),NINT(yob))
   v_ob = v_banana(NINT(xob),NINT(yob))
   speed_ob = SQRT ( u_ob**2 + v_ob**2 ) 

   !  We can see if we need to do any fancy computations by just looking at the speed
   !  of the observation.  If the first-guess speed is less than the cutoff, we do a 
   !  simple circular Cressman and vamoose.

   IF ( speed_ob .GE. vcrit ) THEN
   
      !  The dy/dx has to be defined, so u can't be = zero.  If it is, then we end up
      !  doing an ellipse, not a banana, anyways.

      IF ( ABS(u_ob) .GE. 1.E-5 ) THEN
         dydx = v_ob / u_ob
      ELSE
         dydx = 1000.
      END IF

      !  Differences are computed from approximately 3 grid points on either side of the 
      !  observation location, normal and tanget to the flow at the observation location.
 
      x1 = xob + 3. * ( u_ob/speed_ob )
      y1 = yob + 3. * ( v_ob/speed_ob )
      x2 = xob - 3. * ( u_ob/speed_ob )
      y2 = yob - 3. * ( v_ob/speed_ob )

      !  These values need to be inside the domain since we will use them as indices
      !  for the u and v first-guess winds.

      i1 = MAX ( MIN ( NINT(x1) , iew ) , 1 )
      j1 = MAX ( MIN ( NINT(y1) , jns ) , 1 )
      i2 = MAX ( MIN ( NINT(x2) , iew ) , 1 )
      j2 = MAX ( MIN ( NINT(y2) , jns ) , 1 )

      !  And the distance between these two points (in grid point units) is 

      d12 = SQRT( (REAL(i1) - REAL(i2))**2 + (REAL(j1) - REAL(j2))**2 )

      !  So we can now compute the du/dx and the dv/dx, where we are assumming that
      !  we are on a streamline through this observation point.  We extend the 
      !  difference calculation in approximately 3 grid points to either side of
      !  the observation point on the tangent to the streamline.

      dudx = ( u_banana(i1,j1) - u_banana(i2,j2) ) / ( d12 * dxd ) 
      dvdx = ( v_banana(i1,j1) - v_banana(i2,j2) ) / ( d12 * dxd ) 
   
      !  Build the numerator and denominator of the curvature equation (k from above).
      !  Since the radius of curvature is the inverse of the curvature, there is no
      !  real need to do a temporary computation and then an inverse, so we just do it
      !  directly (denominator / numerator).  The radius of curvature is now in units
      !  of (m).
   
      numerator_curvature = ABS ( u_ob * dvdx - v_ob * dudx ) / u_ob**2
      denominator_curvature = ( SQRT ( 1. + dydx**2 ) ) ** 3

      IF ( ( ABS ( u_ob) .LT. 1.E-5 ) .OR. ( ABS ( v_ob) .LT. 1.E-5 ) .OR. &
           ( ABS ( numerator_curvature ) .LT. 1.E-5 ) ) THEN
         radius_of_curvature = 1000.
      ELSE
!        k = numerator / denominator
!        radius_of_curvature = 1. / k
         radius_of_curvature = denominator_curvature / numerator_curvature
      END IF

   END IF

   !  Circular Cressman implies a Euclidean distance without elongation.  This is used
   !  when the speed of the observation is less than the critical value, for this 
   !  pressure level.

   IF ( speed_ob .LT. vcrit ) THEN

      USE_CIRCLE = .true.

   !  If the wind speed is large enough and the radius of curvature is small enough, then
   !  the banana shaped ellipse is used.

   ELSE IF ( ABS(radius_of_curvature) .LT. ( 3.* radius_influence ) * dxd ) THEN

      USE_BANANA = .true.


      one_divided_by_twopi = 1. / twopi
      one_divided_by_pi = 1. / pi

      !  The radius of curvature needs to have a sign associated with it so that
      !  from the given point on the radius, we can find the center.  From the 
      !  observation location, find the four surrounding points that are three
      !  grid units away.  However, as before, they need to be in the domain.

      ia = MAX ( NINT ( xob ) - 3 ,   1 ) 
      ja =              yob             
      ib =              xob              
      jb = MAX ( NINT ( yob ) - 3 ,   1 ) 
      ic = MIN ( NINT ( xob ) + 3 , iew ) 
      jc =              yob             
      id =              xob              
      jd = MIN ( NINT ( yob ) + 3 , jns ) 

      !  With these four grid points, get the first-guess u and v at these points.

      ua = u_banana(ia,ja)
      va = v_banana(ia,ja)
      ub = u_banana(ib,jb)
      vb = v_banana(ib,jb)
      uc = u_banana(ic,jc)
      vc = v_banana(ic,jc)
      ud = u_banana(id,jd)
      vd = v_banana(id,jd)

      !  Get the SIGN of the relative vorticity at the observation point with
      !  these surrounding values.
   
      vor = ( vc - va ) - ( ud - ub )

      !  Modify the SIGN of the radius_of_curvature with the SIGN of the relative
      !  vorticity.

      radius_of_curvature = radius_of_curvature * SIGN ( 1.,vor )

      !  Compute the location of the center of curvature for the streamline
      !  passing through this observation location.  This is the same as eqn 4.20a
      !  and 4.20b in MH94, except that i0, ic, j0, jc are not in the reversed
      !  model directions (i is east west, j is north south).

      i_center = NINT(xob) - radius_of_curvature * v_ob / speed_ob
      j_center = NINT(yob) + radius_of_curvature * u_ob / speed_ob
      

      !  Compute the angle from the positve x-direction to the line connecting the
      !  center of the curvature to the observation.
      !  The variable theta_ob is the same as the output from eqn 4.21 MH94.

      theta_ob = ATAN2 ( REAL(j_center)-yob     , REAL(i_center)-xob     )

      elongation = 1. + 0.7778 / vcrit * speed_ob


   !  We have the fast enough speed (the observation speed is larger than the critical
   !  speed), but too large radius of curvature (radius_curvature > 3X radius of influence).
   !  This is, therefore, the ellipse scheme.

   ELSE

      USE_ELLIPSE = .true.

      elongation = 1. + 0.7778 / vcrit * speed_ob

   END IF


   !Determine if RH scaling should be done
   IF ( ( scale_cressman_rh_decreases ) .AND. &
        ( ( name(1:8) .EQ. 'RH      ' ) .AND. ( use_first_guess ) .AND. ( obs_n .lt. 0.0) ) ) THEN
    do_rh_scaling = .TRUE.
   ELSE 
    do_rh_scaling = .FALSE.
   ENDIF

   !Now that we have done everything we can do outside the loop, loop over all
   !relevant x, y points for whatever shape was chosen
   JLOOP: DO jj = jns_min , jns_max
    real_j = REAL(jj)
    ILOOP: DO ii = iew_min , iew_max
     real_i = REAL(ii)

     IF( USE_CIRCLE ) THEN

      distance_squared = &
       ( ( xob - real_crsdot_divide_2 - real_i ) * ( xob - real_crsdot_divide_2 - real_i ) ) +  &
       ( ( yob - real_crsdot_divide_2 - real_j ) * ( yob - real_crsdot_divide_2 - real_j ) )

     ELSEIF( USE_BANANA ) THEN

      !  Compute the angle from the positve x-direction to the line connecting the
      !  center of the curvature to the (i,j) locations in the window.
      !  theta_ij is eqn 4.22 from MH94.

      theta_ij = ATAN2 ( REAL(j_center)-real_j , REAL(i_center)-real_i )

      !  Compute the distance from the center of the curvature to the grid point
      !  (i,j).  This corresponds to the definition given for rij on the top of page 63
      !  in MH94.

      rij = SQRT ( (REAL(i_center)-real_i)**2 + (REAL(j_center)-real_j)**2 )

      !  All of the required pieces are assembled for use in MH94 eqns 4.15a,
      !  4.14b and 4.14.  First, make sure that the values of theta_diff are bounded
      !  by +- 2 pi.  Then, bound the values between +- pi.

      theta_diff = theta_ob - theta_ij
      IF ( theta_diff .GT.  twopi ) theta_diff = theta_diff -  &
       REAL ( INT ( theta_diff*one_divided_by_twopi ) ) * twopi
      IF ( theta_diff .LT. -twopi ) theta_diff = theta_diff + &
        REAL ( INT ( theta_diff*one_divided_by_twopi ) ) * twopi
      IF ( theta_diff .GT.     pi ) theta_diff = theta_diff - twopi
      IF ( theta_diff .LT.    -pi ) theta_diff = theta_diff + twopi

!     x = radius_of_curvature * theta_diff
      x = radius_of_curvature * theta_diff * one_divided_by_pi  ! this dividing by pi is home grown
      y = ABS(radius_of_curvature) - rij

      distance_squared = x*x/elongation + y*y

     ELSEIF( USE_ELLIPSE ) THEN

      x = ( u_ob * ( real_i - xob ) + v_ob * ( real_j - yob ) ) / speed_ob
      y = ( v_ob * ( real_i - xob ) - u_ob * ( real_j - yob ) ) / speed_ob

      distance_squared = x*x/elongation + y*y

     ELSE
      PRINT *,'ERROR: dist_uu_vv_rh: It appears that no shape was chosen for radius of influence'
      STOP
     ENDIF


     !  If the observation is within a grid point's radius of influence,
     !  sum its effect.

     IF ( distance_squared .LT. radius_influence_squared ) THEN
      scale_factor = 1.0
      IF ( do_rh_scaling ) THEN
        !do_rh_scaling is TRUE if ALL of the following are TRUE:
        !1) Field is RH, 2) use_first_guess is true, 3) the observed RH was
        !lower than the first-guess RH
        !Calculate a scale_factor that is:
        ! first-guess RH at current location (i.e., location where ob is being applied
        !     DIVIDED BY
        ! First-guess RH at the location of the observation
        scale_factor = min( 1.0, gridded(ii,jj)/fg_value_n)
      END IF
      !weight = ( radius_influence**2 - distance**2 ) / ( radius_influence**2 + distance**2 )
      weight = ( radius_influence_squared - distance_squared ) / ( radius_influence_squared + distance_squared )
      numerator_pert(ii,jj)   = numerator_pert(ii,jj)   + weight * weight * obs_n * scale_factor
      denominator_pert(ii,jj) = denominator_pert(ii,jj) + weight
     END IF

    END DO ILOOP
   END DO JLOOP
 
END SUBROUTINE dist_uu_vv_rh


SUBROUTINE dist_not_uu_vv_rh ( iew_min, iew_max, jns_min, jns_max, real_crsdot_divide_2,  name , &
                  obs_n, xob , yob ,  radius_influence , numerator, denominator )
  
   IMPLICIT NONE
 
   !  Input variables.

   INTEGER , INTENT(IN)                      :: iew_min  
   INTEGER , INTENT(IN)                      :: iew_max 
   INTEGER , INTENT(IN)                      :: jns_min
   INTEGER , INTENT(IN)                      :: jns_max 
   REAL    , INTENT(IN)                      :: real_crsdot_divide_2
   CHARACTER (LEN=8), INTENT(IN)             :: name
   REAL    , INTENT(IN)                      :: obs_n 
   REAL    , INTENT(IN)                      :: xob     
   REAL    , INTENT(IN)                      :: yob
   INTEGER , INTENT(IN)                      :: radius_influence
   !  Input/Output variables.
   REAL    , INTENT(INOUT), DIMENSION (:,:)  :: denominator,    &
                                                numerator

   !  Local variables.
   REAL    :: real_i , real_j  
   REAL    :: distance_squared, radius_influence_squared
   INTEGER :: ii,jj
   REAL    :: weight

   radius_influence_squared = REAL(radius_influence * radius_influence)


   JLOOP: DO jj = jns_min , jns_max
    real_j = REAL(jj)
    ILOOP: DO ii = iew_min , iew_max
     real_i = REAL(ii)

     distance_squared = &
      ( ( xob - real_crsdot_divide_2 - real_i ) * ( xob - real_crsdot_divide_2 - real_i ) ) +  &
      ( ( yob - real_crsdot_divide_2 - real_j ) * ( yob - real_crsdot_divide_2 - real_j ) )

     IF ( distance_squared .LT. radius_influence_squared ) THEN

      weight = ( radius_influence_squared - distance_squared ) / ( radius_influence_squared + distance_squared )
      numerator(ii,jj)   = numerator(ii,jj)   + weight * weight * obs_n
      denominator(ii,jj) = denominator(ii,jj) + weight

     END IF
   
    END DO ILOOP
   END DO JLOOP

END SUBROUTINE dist_not_uu_vv_rh

!
!---------------------------------------------------------------------------

!BPR BEGIN
!SUBROUTINE get_background_for_oa ( t , u , v , rh , slp_x , &
SUBROUTINE get_background_for_oa ( t , u , v , rh , slp_x , pres , &
!BPR END
dum2d , kp , name , &
iew_alloc , jns_alloc , kbu_alloc )

   INCLUDE 'first_guess_size.inc'
   INCLUDE 'first_guess.inc'

   REAL , DIMENSION ( iew_alloc , jns_alloc ) :: dum2d
   INTEGER                                    :: kp
   CHARACTER ( LEN = 8 )                      :: name


   which_var_to_get : SELECT CASE ( name ) 

      CASE ( 'UU      ' ) which_var_to_get
         CALL yx2xy (     u(1,1,kp) , jns_alloc , iew_alloc , dum2d )

      CASE ( 'VV      ' ) which_var_to_get
         CALL yx2xy (     v(1,1,kp) , jns_alloc , iew_alloc , dum2d )

      CASE ( 'TT      ' ) which_var_to_get
         CALL yx2xy (     t(1,1,kp) , jns_alloc , iew_alloc , dum2d )

      CASE ( 'RH      ' ) which_var_to_get
         CALL yx2xy (    rh(1,1,kp) , jns_alloc , iew_alloc , dum2d )

      CASE ( 'PMSL    ' ) which_var_to_get
         CALL yx2xy ( slp_x         , jns_alloc , iew_alloc , dum2d )

!BPR BEGIN
      CASE ( 'PSFC    ' ) which_var_to_get
         IF(kp.ne.1) THEN
          PRINT *,'ERROR: get_background_for_oa: When get_background_for_oa is called'
          PRINT *,' for PSFC the vertical level must be set to 1 but kp = ',kp
          STOP 
         END IF
         CALL yx2xy ( pres(1,1,kp)  , jns_alloc , iew_alloc , dum2d )
!BPR END

   END SELECT which_var_to_get

END SUBROUTINE get_background_for_oa

!
!---------------------------------------------------------------------------

SUBROUTINE mqd ( obs , xob , yob , numobs , &
gridded , iew , jns , &
crsdot , name , passes , smooth_type , use_first_guess ) 

!  Multiquadric interpolation based on Nuss and Titley (1994 MWR).
!  The free parameters are lambda, the smoothing factor and c, the
!  multiquadric parameter. If the procedure bombs, the most likely
!  problem is that the value for c is incorrect. For long, narrow
!  domains, it might be a problem. When there are lots of obs, this
!  routine is a memory hog.

   ! arguments:
   !    obs:       difference between 1st guess and obs
   !               ( need any duplicated obs removed from this array )
   !    numobs:    number of station obs
   !    xob,yob:   x, y locations of station obs
   !    gridded:   background/first field as input, final analysis as output
   !    iew,ins:   1st and 2nd dimensions of array 'gridded'
   !    name:      variable name, T, RH, p, U, or V
   !    use_first_guess: T/F use the first-guess fields in the objective analysis
    
   ! other variables:
   !    qi:        distance function for obs
   !    qgi:       distance function for gridded data
   !    dummy:     used for matrix calculation
   !    workarray: to hold the input first guess array
   !    errm:      error account for obs

   IMPLICIT NONE

   INTEGER, INTENT ( IN )                     :: numobs
   INTEGER, INTENT ( IN )                     :: iew, jns
   REAL,    INTENT ( IN ), DIMENSION ( numobs )      :: obs, xob, yob
   REAL,    INTENT ( INOUT ), DIMENSION( iew , jns ) :: gridded
   CHARACTER (8), INTENT ( IN )               :: name
   INTEGER                                    :: crsdot
   INTEGER , INTENT ( IN )                    :: passes , smooth_type 
   LOGICAL , INTENT ( IN )                    :: use_first_guess

   REAL, DIMENSION ( numobs, numobs )         :: qi
   REAL                                       :: errm,           &
                                                 c
   INTEGER                                    :: numpts,         &
                                                 i,  j,          &
                                                 ig, jg
   REAL, DIMENSION ( iew,jns )                :: workarray2
   REAL, PARAMETER                            :: lambda = 0.0025
   INTEGER , DIMENSION ( numobs ) :: pvt
   INTEGER , PARAMETER            :: job = 01
   REAL    , DIMENSION ( numobs ) :: z , & 
                                     work
   REAL                           :: condition , det(2)
   INTEGER                        :: imult , jmult , kmult


   REAL, DIMENSION ( iew, numobs )            :: qgi
   REAL, DIMENSION (numobs)                   :: dummy

   INCLUDE 'error.inc'
   INTERFACE 
      INCLUDE 'error.int'
   END INTERFACE

   c      = 0.01 * ( max ( iew, jns) )
   numpts = iew * jns

   !  Set errm, the mean error value for the variable being analyzed. 
   !  The results aren't too sensitive to this parameter, but the 
   !  values must be sane.

   SELECT CASE ( name )
      CASE ( 'UU      ' ) ; errm = 2.
      CASE ( 'VV      ' ) ; errm = 2.
      CASE ( 'PMSL    ' ) ; errm = 1.
      !BPR BEGIN
      CASE ( 'PSFC    ' ) ; errm = 1.
      !BPR END
      CASE ( 'TT      ' ) ; errm = 0.5
      CASE ( 'RH      ' ) ; errm = 10.
   END SELECT

   !  Compute the radial basis function (how far is each observation
   !  away from every other observation).  The diagonal entries are
   !  identically zero, and are treated below.

   do j = 1, numobs
      do i = 1, numobs
         qi(i,j) = -1.* SQRT(((xob(j)-xob(i))**2 +               &
                              (yob(j)-yob(i))**2)/(c*c)+1.)
      end do
   end do

   !  Account for observational uncertainty (Nuss and Titley 1994),
   !  which are the diagonal entries.

   do j = 1, numobs
      qi(j,j) = qi(j,j) + numobs * lambda * errm
   end do

   !  Find the inverse of Qi (radial basis function).  These are single
   !  precision LINPAK routines.  Double precision routines are available.

   CALL sgeco ( qi , numobs , numobs , pvt , condition , z )
   CALL sgedi ( qi , numobs , numobs , pvt , det , work , job )

   !  Compute the radial basis function for gridded points (similar to
   !  how far away is each observation away from every grid point).

   dummy = matmul(Qi, obs)

   DO jg = 1, jns
      DO ig = 1, iew
         DO i = 1, numobs
            qgi(ig,i) = -1.* SQRT(((REAL(ig)-xob(i))**2 +          &
                                  (REAL(jg)-yob(i))**2)/(c*c)+1.)
         END DO
      END DO
      workarray2(:,jg) = matmul(Qgi,dummy)
   END DO


   !  Which smoother will we use? Smooth the final analysis with a 5-point stencil 
   !  (1-2-1 on the lateral boundaries), or the traditional smoother-desmoother.
   !  Both accept a 2D field with the same arguments.

   IF      ( smooth_type .EQ. 1 ) THEN
      CALL smooth_5            ( workarray2 , iew, jns, passes , crsdot ) 
   ELSE IF ( smooth_type .EQ. 2 ) THEN
      CALL smoother_desmoother ( workarray2 , iew, jns, passes , crsdot )
   END IF

   !  Add the first guess to the perturbations from the observations at
   !  each (i,j).

   IF ( use_first_guess ) THEN
      gridded = gridded + workarray2
   ELSE
      gridded = workarray2
   END IF

END SUBROUTINE mqd

!
!---------------------------------------------------------------------------

!BPR BEGIN
!SUBROUTINE put_background_from_oa ( t , u , v , rh , slp_x , &
SUBROUTINE put_background_from_oa ( t , u , v , rh , slp_x , pres , &
!BPR END
dum2d , kp , name , &
iew_alloc , jns_alloc , kbu_alloc )

   INCLUDE 'first_guess_size.inc'
   INCLUDE 'first_guess.inc'

   REAL , DIMENSION ( iew_alloc , jns_alloc ) :: dum2d
   INTEGER                                    :: kp
   CHARACTER ( LEN = 8 )                      :: name


   which_var_to_get : SELECT CASE ( name ) 

      CASE ( 'UU      ' ) which_var_to_get
         CALL yx2xy ( dum2d , iew_alloc , jns_alloc ,     u(1,1,kp) ) 

      CASE ( 'VV      ' ) which_var_to_get
         CALL yx2xy ( dum2d , iew_alloc , jns_alloc ,     v(1,1,kp) ) 

      CASE ( 'TT      ' ) which_var_to_get
         CALL yx2xy ( dum2d , iew_alloc , jns_alloc ,     t(1,1,kp) ) 

      CASE ( 'RH      ' ) which_var_to_get
         CALL yx2xy ( dum2d , iew_alloc , jns_alloc ,    rh(1,1,kp) ) 

      CASE ( 'PMSL    ' ) which_var_to_get
         CALL yx2xy ( dum2d , iew_alloc , jns_alloc , slp_x         ) 

      !BPR BEGIN
      CASE ( 'PSFC    ' ) which_var_to_get
         IF(kp.ne.1) THEN
          PRINT *,'ERROR: pot_background_from_oa: When put_background_from_oa is called'
          PRINT *,' for PSFC the vertical level must be set to 1 but kp = ',kp
          STOP 
         ENDIF
         CALL yx2xy ( dum2d , iew_alloc , jns_alloc , pres(1,1,kp)  ) 
      !BPR END

   END SELECT which_var_to_get

END SUBROUTINE put_background_from_oa

!
!---------------------------------------------------------------------------

SUBROUTINE smoother_desmoother ( slab , imx , jmx , passes , crsdot )

   IMPLICIT NONE

   INTEGER                        :: imx , jmx , passes , crsdot
   REAL , DIMENSION ( imx , jmx ) :: slab , & 
                                     slabnew

   REAL , DIMENSION ( 2 )         :: xnu
   INTEGER                        :: i , j , loop , n 

   xnu  =  (/ 0.50 , -0.52 /)

   !  The odd number passes of this are the "smoother", the even
   !  number passes are the "de-smoother" (note the differnt signs on xnu).

   smoothing_passes : DO loop = 1 , passes * 2

      n  =  2 - MOD ( loop , 2 )
 
      DO i = 2 , imx - 1 - crsdot
         DO j = 2 , jmx - 1 - crsdot
            slabnew(i,j) = slab(i,j) + xnu(n) *  & 
            ((slab(i,j+1) + slab(i,j-1)) * 0.5-slab(i,j))
         END DO
      END DO
 
      DO i = 2 , imx - 1 - crsdot
         DO j = 2 , jmx - 1 - crsdot
            slab(i,j) = slabnew(i,j)
         END DO
      END DO
 
      DO j = 2 , jmx - 1 - crsdot
         DO i = 2 , imx - 1 - crsdot
            slabnew(i,j) = slab(i,j) + xnu(n) *  &
            ((slab(i+1,j) + slab(i-1,j)) * 0.5-slab(i,j))
         END DO
      END DO
 
      DO i = 2 , imx - 1 - crsdot
         DO j = 2 , jmx - 1 - crsdot
            slab(i,j) = slabnew(i,j)
         END DO
      END DO
 
   END DO smoothing_passes

END SUBROUTINE smoother_desmoother

!
!---------------------------------------------------------------------------

SUBROUTINE smooth_5 ( field , iew , jns , passes , crsdot )

   INTEGER                        :: iew , jns , &
                                     passes    , &
                                     crsdot
   REAL , DIMENSION ( iew , jns ) :: field

   REAL , DIMENSION ( iew , jns ) :: temp
   INTEGER                        :: i , j , num_passes

   !  How may passes of this smoother are we using.

   smoothing_passes : DO num_passes = 1 , passes

      !  Apply 5-point stencil smoother on interior of the domain.
   
      DO j = 2 , jns - 1 - crsdot
         DO i = 2 , iew - 1 - crsdot
            temp(i,j) = ( field(i  ,j  ) * 4. +  & 
                          field(i+1,j  )      +  & 
                          field(i-1,j  )      +  & 
                          field(i  ,j+1)      +  & 
                          field(i  ,j-1)      )  * 0.125
         END DO
      END DO

      !  Apply 3-point stencil smoother on the boundaries.
   
      i = 1
      DO j = 2 , jns - 1 - crsdot
         temp(i,j) = ( field(i  ,j  ) * 2. +  & 
                       field(i  ,j+1)      +  & 
                       field(i  ,j-1)      )  * 0.25
      END DO

      i = iew - crsdot
      DO j = 2 , jns - 1 - crsdot
         temp(i,j) = ( field(i  ,j  ) * 2. +  & 
                       field(i  ,j+1)      +  & 
                       field(i  ,j-1)      )  * 0.25
      END DO
   
      j = 1
      DO i = 2 , iew - 1 - crsdot
         temp(i,j) = ( field(i  ,j  ) * 2. +  & 
                       field(i+1,j  )      +  & 
                       field(i-1,j  )      ) * 0.25
      END DO
   
      j = jns - crsdot
      DO i = 2 , iew - 1 - crsdot
         temp(i,j) = ( field(i  ,j  ) * 2. +  & 
                       field(i+1,j  )      +  & 
                       field(i-1,j  )      ) * 0.25
      END DO
   
      !  Store smoothed field back into original array.
   
      DO j = 2 , jns - 1 - crsdot
         DO i = 2 , iew - 1 - crsdot
            field(i,j) = temp(i,j)
         END DO
      END DO
   
      !  Store smoothed boundary field back into original array.
   
      DO j = 2 , jns - 1 - crsdot
         field(1         ,j) = temp(1         ,j)
         field(iew-crsdot,j) = temp(iew-crsdot,j)
      END DO
   
      DO i = 2 , iew - 1 - crsdot
         field(i,1         ) = temp(i,1         )
         field(i,jns-crsdot) = temp(i,jns-crsdot)
      END DO

   END DO smoothing_passes

END SUBROUTINE smooth_5

!---------------------------------------------------------------------------
!NOTE THAT THIS IS AN ELEMENTAL FUNCTION SO WE CAN PASS IT SCALARS OR ARRAYS
!AND THE RESULT WILL BE CONFORMABLE TO THE INPUT ARRAYS
!Find water vapor mixing ratio given relative humidity, temperature, and
!pressure
!Algorithm from subroutine mixing_ratio in final_analysis_module.F90
!If error then returns -99999, since we cannot "STOP" because this is an
!elemental function
ELEMENTAL FUNCTION qv_from_rh(rh,tk,pressure)

 IMPLICIT NONE

 REAL, INTENT(IN) :: rh  !Relative humidity (%)
 REAL, INTENT(IN) :: tk  !Temperature (Kelvin)
 REAL, INTENT(IN) :: pressure !Pressure (hPa)
 REAL :: qv_from_rh      !Water vapor mixing ratio (kg/kg)

 REAL :: es, qs
 REAL :: rh_use

 REAL,         PARAMETER     :: EPS         = 0.622
 REAL,         PARAMETER     :: SVP1        = 0.6112
 REAL,         PARAMETER     :: SVP2        = 17.67
 REAL,         PARAMETER     :: SVP3        = 29.65
 REAL,         PARAMETER     :: SVPT0       = 273.15

 REAL,         PARAMETER     :: small_diff  = 0.00001
 LOGICAL                     :: error_found

 error_found = .false.
 
 !Do some error checking on the inputs

 !RH
 rh_use = rh
 IF((rh.lt.0.0).or.(rh.gt.100.0)) THEN
  !If its out of range by more than a little bit then error out
  IF((rh.lt.(0.0-small_diff)).OR.(rh.gt.(100.0+small_diff))) THEN
   !PRINT *,'ERROR: qv_from_rh: RH outside expected ranges: ',rh
   error_found = .true. 
   qv_from_rh = -99999.0
  ELSE 
  !If its just slightly outside of range then set it to edge of range
   rh_use=max(min(rh,100.0),0.0)
  ENDIF
 ENDIF

 !Temperature
 IF(tk<50.0) THEN
  !PRINT *,'ERROR: qv_from_rh: Temperature of less than 50 K passed in: ',tk
  error_found = .true. 
   qv_from_rh = -99998.0
 ENDIF

 !Pressure
 IF((pressure<0.0).or.(pressure>1200.0)) THEN
  !PRINT *,'ERROR: qv_from_rh: Pressure outside expected range passed in: ',pressure,' (hPa aka mb)'
  error_found = .true. 
  qv_from_rh = -99997.0
 ENDIF
  
 IF(error_found) THEN
  !qv_from_rh = -99999.0
 ELSE
  es = svp1 * 10. * EXP(svp2 * (tk - svpt0) / (tk - svp3))
  qs = eps * es / (pressure - es) 
  qv_from_rh = MAX(0.01 * rh_use * qs,0.0)
 ENDIF
 
END FUNCTION qv_from_rh

!---------------------------------------------------------------------------
!NOTE THAT THIS IS AN ELEMENTAL FUNCTION SO WE CAN PASS IT SCALARS OR ARRAYS
!AND THE RESULT WILL BE CONFORMABLE TO THE INPUT ARRAYS
!Find relative humidity given water vapor mixing ratio, temperature, and
!pressure
!Based on algorithm in qv_from_rh which is in turned based on algorithm from 
!subroutine mixing_ratio in final_analysis_module.F90
!If error then returns -99999, since we cannot "STOP" because this is an
!elemental function
ELEMENTAL FUNCTION rh_from_qv(qv,tk,pressure)

 IMPLICIT NONE

 REAL, INTENT(IN) :: qv  !Water vapor mixing ratio (kg/kg)
 REAL, INTENT(IN) :: tk  !Temperature (Kelvin)
 REAL, INTENT(IN) :: pressure !Pressure (hPa)
 REAL :: rh_from_qv      !Relative humidity (%)

 REAL :: es, qs
 REAL :: qv_use

 REAL,         PARAMETER     :: EPS         = 0.622
 REAL,         PARAMETER     :: SVP1        = 0.6112
 REAL,         PARAMETER     :: SVP2        = 17.67
 REAL,         PARAMETER     :: SVP3        = 29.65
 REAL,         PARAMETER     :: SVPT0       = 273.15

 REAL,         PARAMETER     :: small_diff  = 0.00001
 LOGICAL                     :: error_found

 error_found = .false.
 
 qv_use = qv
 !Do some error checking on the inputs

 !Water vapor mixing ratio
 IF(qv.lt.0.0) THEN
  !If its out of range by more than a little bit then error out
  IF(qv.lt.(0-small_diff)) THEN
   !PRINT *,'ERROR: rh_from_qv: QV outside expected ranges: ',qv
   error_found = .true. 
  ELSE 
  !If its just slightly outside of range then set it to edge of range
   qv_use=0.0
  ENDIF
 ENDIF

 !Temperature
 IF(tk<50.0) THEN
  !PRINT *,'ERROR: rh_from_qv: Temperature of less than 50 K passed in: ',tk
  error_found = .true. 
 ENDIF

 !Pressure
 IF((pressure<0.0).or.(pressure>1200.0)) THEN
  !PRINT *,'ERROR: rh_from_qv: Pressure outside expected range passed in: ',pressure,' (hPa aka mb)'
  error_found = .true. 
 ENDIF
  
 IF(error_found) THEN
  rh_from_qv = -99999.0
 ELSE
  es = svp1 * 10. * EXP(svp2 * (tk - svpt0) / (tk - svp3))
  qs = eps * es / (pressure - es) 
  rh_from_qv = MIN ( MAX( ( 100.0 * qv_use ) / qs,0.0), 100.0 )
 ENDIF
 
END FUNCTION rh_from_qv




END MODULE obj_analysis
