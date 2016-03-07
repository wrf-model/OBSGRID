!BPR BEGIN
MODULE get_fg_at_point

CONTAINS

!Retrieve first guess temperature at chosen point
SUBROUTINE GET_FG_T_AT_POINT(fg_3d_t,fg_3d_h,latitude,longitude,height_msl,fg_t_at_point)

   USE map_utils
   USE map_utils_helper

   IMPLICIT NONE

   REAL, DIMENSION(:,:,:),INTENT(IN)  :: fg_3d_t       !First guess 3D temperature
   REAL, DIMENSION(:,:,:),INTENT(IN)  :: fg_3d_h       !First guess 3D height
   REAL,                  INTENT(IN)  :: latitude      !Latitude where we want temperature
   REAL,                  INTENT(IN)  :: longitude     !Longitude where we want temperature
   REAL,                  INTENT(IN)  :: height_msl    !Height (m MSL) where want want temperature
   REAL,                  INTENT(OUT) :: fg_t_at_point !First guess temperature at chosen point

   REAL                    :: x_location, y_location
   INTEGER                 :: jns_alloc,iew_alloc,kbu_alloc
   INTEGER, DIMENSION(4)   :: nearest_x_locations, nearest_y_locations
   INTEGER                 :: cur_point
   LOGICAL                 :: found
   REAL                    :: h_below, h_above
   REAL                    :: t_below, t_above
   REAL                    :: vinterp_weight
   REAL, DIMENSION(4)      :: fg_t_at_near_point
   INTEGER                 :: k
   REAL                    :: x_location_nonint, y_location_nonint

   INCLUDE 'missing.inc'
   INCLUDE 'constants.inc'

   !Get location of ob in domain grid cell space
   CALL latlon_to_ij(projd , latitude  , longitude , x_location , y_location )

   jns_alloc = size(fg_3d_t,1)
   iew_alloc = size(fg_3d_t,2)
   kbu_alloc = size(fg_3d_t,3)

   !If ob is outside domain then mark as missing and return
   IF( (x_location .LT. 1) .OR. (x_location .GT. iew_alloc-1 ) .OR. &
       (y_location .LT. 1) .OR. (y_location .GT. jns_alloc-1 ) ) THEN
    fg_t_at_point = missing_r
    RETURN
   ENDIF

   !Create four surrounding points
   nearest_x_locations(1)=floor(x_location)
   nearest_y_locations(1)=floor(y_location)
   nearest_x_locations(2)=floor(x_location)
   nearest_y_locations(2)=ceiling(y_location)
   nearest_x_locations(3)=ceiling(x_location)
   nearest_y_locations(3)=floor(y_location)
   nearest_x_locations(4)=ceiling(x_location)
   nearest_y_locations(4)=ceiling(y_location)
   

   !Vertically interpolate temperature to desired height at four surrounding grid points
   !If desired point is below bottom of first guess temperature, extrapolate using standard
   !atmosphere
   DO cur_point = 1,4
    found = .FALSE.
    find_trapping : DO k = 2 , kbu_alloc-1    
     IF ( ( height_msl .GE. fg_3d_h(nearest_y_locations(cur_point), &
                            nearest_x_locations(cur_point),k ) ) .AND. &
          ( height_msl .LE. fg_3d_h(nearest_y_locations(cur_point), &
                            nearest_x_locations(cur_point),k+1 ) ) ) THEN
      found = .TRUE.
      h_below = fg_3d_h(nearest_y_locations(cur_point),nearest_x_locations(cur_point),k )
      h_above = fg_3d_h(nearest_y_locations(cur_point),nearest_x_locations(cur_point),k+1 )
      t_below = fg_3d_t(nearest_y_locations(cur_point),nearest_x_locations(cur_point),k )
      t_above = fg_3d_t(nearest_y_locations(cur_point),nearest_x_locations(cur_point),k+1 )
      EXIT find_trapping
     END IF
    END DO find_trapping
    !If we have first guess levels sandwiching the desired height
    IF( found ) THEN 
     !Vertically interpolate temperature to desired height
     !height_msl
     vinterp_weight = (h_above - height_msl ) / ( h_above - h_below ) 
     fg_t_at_near_point(cur_point) =  vinterp_weight * t_below + (1 - vinterp_weight) * t_above 

    !The desired height is below the lowest non-surface first-guess level
    ELSEIF( height_msl .LT. fg_3d_h(nearest_y_locations(cur_point),nearest_x_locations(cur_point),2 ) ) THEN
     fg_t_at_near_point(cur_point) = fg_3d_t(nearest_y_locations(cur_point),nearest_x_locations(cur_point),2 ) - &
                     ( dTdz_std_atm_lt_11000m *  &
                      ( fg_3d_h(nearest_y_locations(cur_point),nearest_x_locations(cur_point),2 ) - height_msl ) )

    !The desired height is above the top first-guess level
    ELSEIF( height_msl .GT. &
            fg_3d_h(nearest_y_locations(cur_point),nearest_x_locations(cur_point),kbu_alloc ) ) THEN
     fg_t_at_point = missing_r
     RETURN
    !We should never get to this branch so error out
    ELSE
     PRINT *,'ERROR: GET_FG_T_AT_POINT: Entered branch that should never be entered: height_msl = ',height_msl
     STOP

    ENDIF

   ENDDO
   
   !Horizontally interpolate to desired location
   x_location_nonint = x_location - floor(x_location)
   y_location_nonint = y_location - floor(y_location)

   fg_t_at_point =  ( 1.-x_location_nonint ) * ( ( 1.-y_location_nonint ) * fg_t_at_near_point(1)   &
                    +                                 y_location_nonint   * fg_t_at_near_point(2) ) &
                    +    x_location_nonint *   ( ( 1.-y_location_nonint ) * fg_t_at_near_point(3)   &
                    +                                 y_location_nonint   * fg_t_at_near_point(4) )


END SUBROUTINE GET_FG_T_AT_POINT

!------------------------------------------------------------------------------
END MODULE get_fg_at_point
!BPR END
