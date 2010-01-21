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
use_first_guess ) 

   ! arguments:
   !    obs_value: difference between 1st guess and obs
   !               ( need any duplicated obs removed from this array )
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
    
   IMPLICIT NONE

   INTEGER, INTENT ( IN )                     :: numobs
   INTEGER, INTENT ( IN )                     :: iew, jns
   REAL,    INTENT ( IN ), DIMENSION ( numobs )      :: obs_value , obs, xob, yob
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
   CHARACTER (LEN=8)                          :: choice

   INCLUDE 'error.inc'

   INTERFACE 
      INCLUDE 'error.int'
   END INTERFACE

   !  Initialize the Cressman weights (numerator and denominator) and
   !  the total perturbation field to zero.

   numerator    = 0
   denominator  = 0
   perturbation = 0

   !  Loop over each observation.

   obs_loop : DO n = 1 , numobs

      !  Every observation has a first guess wind field location that will be specified
      !  as an objective analysis that is circular, elliptical or banana-shaped.

      choice = 'undefine'

      !  To use the observation, it must be far enough inside our domain
      !  (using the cross/dot configuration for the staggered variables).  If
      !  it is not, we simply zoom back to the top of this loop.

      IF ( ( xob(n) .LE.   1 + REAL(crsdot)/2. ) .OR. & 
           ( yob(n) .LE.   1 + REAL(crsdot)/2. ) .OR. & 
           ( xob(n) .GE. iew - REAL(crsdot)/2. ) .OR. & 
           ( yob(n) .GE. jns - REAL(crsdot)/2. ) ) THEN
         CYCLE obs_loop
      END IF

      !  Compute the window in which this observation has influence.  For
      !  some points inside this window the observation will have zero 
      !  influence, but this allows a direct method rather than a search.

      iew_min = MAX ( NINT ( (xob(n)-REAL(crsdot)/2.) ) - 4*radius_influence -1 ,            1 )
      jns_min = MAX ( NINT ( (yob(n)-REAL(crsdot)/2.) ) - 4*radius_influence -1 ,            1 )
      iew_max = MIN ( NINT ( (xob(n)-REAL(crsdot)/2.) ) + 4*radius_influence +1 , iew - crsdot ) 
      jns_max = MIN ( NINT ( (yob(n)-REAL(crsdot)/2.) ) + 4*radius_influence +1 , jns - crsdot )

      !  For the grid points in the window, compute this observation's
      !  influence.

      DO j = jns_min , jns_max
         DO i = iew_min , iew_max

            !  Distance from each of the windowed grid points to the observation.  The distance is in 
            !  grid point units.  If this is a horizontal wind field or a moisture field, we can use the
            !  banana scheme, which elongates the acceptable distance from an observation and also is
            !  able to handle first-guess wind curvature.

            IF ( ( name(1:8) .EQ. 'UU      ' ) .OR. &
                 ( name(1:8) .EQ. 'VV      ' ) .OR. &
                 ( name(1:8) .EQ. 'RH      ' )  ) THEN
               CALL dist ( i , j , iew , jns , crsdot , name(1:8) , obs_value(n) , xob(n) , yob(n) , &
                           dxd , pressure , u_banana , v_banana , radius_influence , &
                           lat_center , choice , distance )
   
            ELSE 
               distance = SQRT ( ( ( xob(n) - REAL(crsdot)/2. - REAL(i) ) * ( xob(n) - REAL(crsdot)/2. - REAL(i) ) ) +  &
                                 ( ( yob(n) - REAL(crsdot)/2. - REAL(j) ) * ( yob(n) - REAL(crsdot)/2. - REAL(j) ) ) ) 

            END IF

            !  If the observation is within a grid point's radius of influence,
            !  sum its effect.

            IF ( distance .LT. radius_influence ) THEN
               weight = ( radius_influence**2 - distance**2 ) / ( radius_influence**2 + distance**2 )
               numerator(i,j)   = numerator(i,j)   + weight * weight * obs(n)
               denominator(i,j) = denominator(i,j) + weight
            END IF

         END DO
      END DO

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

SUBROUTINE dist ( i , j , iew , jns , crsdot , name , obs_value , xob , yob , dxd , &
                  pressure , u_banana , v_banana , radius_influence , &
                  lat_center , choice , distance )
  
   IMPLICIT NONE
 
   !  Input variables.

   INTEGER , INTENT(IN)                      :: i        , &
                                                j        , &
                                                iew      , &
                                                jns      , &
                                                crsdot   , &
                                                radius_influence
   REAL    , INTENT(IN)                      :: obs_value , &
                                                xob      , &
                                                yob      , &
                                                dxd      , &
                                                pressure , &
                                                lat_center
   REAL    , INTENT(IN) , DIMENSION(iew,jns) :: u_banana , &
                                                v_banana
   CHARACTER (LEN=8), INTENT(IN)             :: name

   !  Output variables.

   REAL    , INTENT(OUT)                     :: distance
   CHARACTER (LEN=8) , INTENT(INOUT)         :: choice

   !  Local variables.

   REAL :: radius_of_curvature , k
   REAL :: u_ob , dudx , speed_ob , dydx
   REAL :: v_ob , dvdx
   REAL :: numerator , denominator
   REAL :: x1 , x2 , y1 , y2 , d12
   INTEGER :: i1 , i2 , j1 , j2
   INTEGER :: ia, ja, ib, jb, ic, jc, id, jd
   REAL :: ua, va, ub, vb, uc, vc, ud, vd, vor
   

   INTEGER :: i_center , j_center

   REAL :: theta_ob , theta_ij
   
   REAL :: rij
   
   REAL :: vcrit , elongation

   REAL :: x , y , theta_diff

   SAVE

   IF ( choice .EQ. 'undefine' ) THEN

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
      
         numerator = ABS ( u_ob * dvdx - v_ob * dudx ) / u_ob**2
         denominator = ( SQRT ( 1. + dydx**2 ) ) ** 3
   
         IF ( ( ABS ( u_ob) .LT. 1.E-5 ) .OR. ( ABS ( v_ob) .LT. 1.E-5 ) .OR. &
              ( ABS ( numerator ) .LT. 1.E-5 ) ) THEN
            radius_of_curvature = 1000.
         ELSE
   !        k = numerator / denominator
   !        radius_of_curvature = 1. / k
            radius_of_curvature = denominator / numerator
         END IF
   
      END IF

   END IF

   !  Circular Cressman implies a Euclidean distance without elongation.  This is used
   !  when the speed of the observation is less than the critical value, for this 
   !  pressure level.

   IF ( speed_ob .LT. vcrit ) THEN

      distance = SQRT ( ( ( xob - REAL(crsdot)/2. - REAL(i) ) * ( xob - REAL(crsdot)/2. - REAL(i) ) ) +  &
                        ( ( yob - REAL(crsdot)/2. - REAL(j) ) * ( yob - REAL(crsdot)/2. - REAL(j) ) ) ) 
   
      choice = 'circle  '

   !  If the wind speed is large enough and the radius of curvature is small enough, then
   !  the banana shaped ellipse is used.

   ELSE IF ( ABS(radius_of_curvature) .LT. ( 3.* radius_influence ) * dxd ) THEN

      IF ( choice .EQ. 'undefine' ) THEN

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
      
      END IF

      !  Compute the angle from the positve x-direction to the line connecting the
      !  center of the curvature to a couple of different points - the observation and
      !  the (i,j) locations in the window.  The variable theta_ob is the same as the
      !  output from eqn 4.21 MH94, and theta_ij is eqn 4.22 from MH94.

      theta_ob = ATAN2 ( REAL(j_center)-yob     , REAL(i_center)-xob     )
      theta_ij = ATAN2 ( REAL(j_center)-REAL(j) , REAL(i_center)-REAL(i) )

      !  Compute the distance from the center of the curvature to the grid point
      !  (i,j).  This corresponds to the definition given for rij on the top of page 63
      !  in MH94.

      rij = SQRT ( (REAL(i_center)-REAL(i))**2 + (REAL(j_center)-REAL(j))**2 )

      !  All of the required pieces are assembled for use in MH94 eqns 4.15a, 
      !  4.14b and 4.14.  First, make sure that the values of theta_diff are bounded
      !  by +- 2 pi.  Then, bound the values between +- pi.

      theta_diff = theta_ob - theta_ij
      IF ( theta_diff .GT.  twopi ) theta_diff = theta_diff - REAL ( INT ( theta_diff/twopi ) ) * twopi
      IF ( theta_diff .LT. -twopi ) theta_diff = theta_diff + REAL ( INT ( theta_diff/twopi ) ) * twopi
      IF ( theta_diff .GT.     pi ) theta_diff = theta_diff - twopi
      IF ( theta_diff .LT.    -pi ) theta_diff = theta_diff + twopi

!     x = radius_of_curvature * theta_diff
      x = radius_of_curvature * theta_diff / pi  ! this dividing by pi is home grown
      y = ABS(radius_of_curvature) - rij

      elongation = 1. + 0.7778 / vcrit * speed_ob
      distance = SQRT ( x*x/elongation + y*y )

      choice = 'banana  '

   !  We have the fast enough speed (the observation speed is larger than the critical
   !  speed), but too large radius of curvature (radius_curvature > 3X radius of influence).
   !  This is, therefore, the ellipse scheme.

   ELSE

      x = ( u_ob * ( REAL(i) - xob ) + v_ob * ( REAL(j) - yob ) ) / speed_ob
      y = ( v_ob * ( REAL(i) - xob ) - u_ob * ( REAL(j) - yob ) ) / speed_ob

      elongation = 1. + 0.7778 / vcrit * speed_ob
      distance = SQRT ( x*x/elongation + y*y )

      choice = 'ellipse '

   END IF
 
END SUBROUTINE dist

!
!---------------------------------------------------------------------------

SUBROUTINE get_background_for_oa ( t , u , v , rh , slp_x , &
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

SUBROUTINE put_background_from_oa ( t , u , v , rh , slp_x , &
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

END MODULE obj_analysis
