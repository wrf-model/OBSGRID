!------------------------------------------------------------------------------

SUBROUTINE proc_get_info_header ( print_header , &
iewd , jnsd , kbu , grid_id , map_projection , expanded , iewe , jnse , &
dxc , lat_center , lon_center , cone_factor , true_lat1 , true_lat2 , pole , &
dxd , ptop , help )
   
!  This routine returns a value to the user, from a long list of available
!  options.  This routine is intended to be user callable.  If the keyword
!  help is passed (regardless of whether help is assigned T or F), then 
!  a list of available parameters is output.  This is a utility routine 
!  may be called from anywhere, though the record header must be 
!  initialized by a call to proc_header prior to this call (from anywhere
!  else in the program).
   
   IMPLICIT NONE

   LOGICAL            :: print_header
   
   !  INTEGER information 
   
   INTEGER , OPTIONAL :: iewd              , & ! domain i size
                         jnsd              , & ! domain j size
                         kbu               , & ! num vert levs, bottom up
                         grid_id           , & ! domain ID
                         map_projection    , & ! 1=LC, 2=PS, 3=MC
                         expanded          , & ! 1=expanded; 0=not expanded
                         iewe              , & ! expanded grid i size
                         jnse                  ! expanded grid j size
                                                        
   !  REAL information 
   
   REAL ,    OPTIONAL :: dxc               , & ! coarse grid distance
                         lat_center        , & ! center latitude
                         lon_center        , & ! central longitude
                         cone_factor       , & ! cone factor
                         true_lat1         , & ! lat of true projection
                         true_lat2         , & ! lat of true projection
                         pole              , & ! location of pole in latitude
                         dxd                   ! domain grid distance
   
   !  INTEGER information from input program.
   

   !  REAL information from REGRID.
   
   REAL ,    OPTIONAL :: ptop                  ! top of analysis, lid
                                                        
   !  This is to allow folks to see what is currently available
   !  from inside here, names and such.

   LOGICAL , OPTIONAL :: help                  ! T/F: print stuff only

   !  For a loop index for some print out, we need this guy defined.

   INTEGER            :: loop_count

   !  The other end of this common is in proc_header.  This is a way
   !  to allow calls to retrieve data, without requiring the 
   !  record header to be included in the argument list.
   
   INCLUDE 'header.common'

   !  Is this a request for simple printout of the available data, or
   !  a valid data request?

   present_help : IF ( PRESENT ( help ) ) THEN

      !  Display all of the available information from this routine, and
      !  return to the calling routine.  Note that even though help is
      !  a LOGICAL variable, it does not need to be initialized, only 
      !  PRESENT.  Even a FALSEly initialized help will activate this
      !  block of the IF statement.

      WRITE ( UNIT = * , FMT = * ) '-------------------------------'
      WRITE ( UNIT = * , FMT = * ) 'INTEGER information              '
      WRITE ( UNIT = * , FMT = * ) 'iewd              domain i size'
      WRITE ( UNIT = * , FMT = * ) 'jnsd              domain j size'
      WRITE ( UNIT = * , FMT = * ) 'kbu               num vert levels'
      WRITE ( UNIT = * , FMT = * ) 'grid_id           grid ID'
      WRITE ( UNIT = * , FMT = * ) 'map_projection    1=LC, 2=PS, 3=MC'
      WRITE ( UNIT = * , FMT = * ) 'expanded          1=expanded; 0=not expanded'
      WRITE ( UNIT = * , FMT = * ) 'iewe              expanded grid i size'
      WRITE ( UNIT = * , FMT = * ) 'jnse              expanded grid j size'
      WRITE ( UNIT = * , FMT = * ) '-------------------------------'
      WRITE ( UNIT = * , FMT = * ) 'REAL information              '
      WRITE ( UNIT = * , FMT = * ) 'dxc               coarse grid distance'
      WRITE ( UNIT = * , FMT = * ) 'lat_center        center latitude'
      WRITE ( UNIT = * , FMT = * ) 'lon_center        central longitude'
      WRITE ( UNIT = * , FMT = * ) 'cone_factor       cone factor'
      WRITE ( UNIT = * , FMT = * ) 'true_lat1         lat of true projection'
      WRITE ( UNIT = * , FMT = * ) 'true_lat2         lat of true projection'
      WRITE ( UNIT = * , FMT = * ) 'pole              location of pole in latitude'
      WRITE ( UNIT = * , FMT = * ) 'dxd               domain grid distance'
      WRITE ( UNIT = * , FMT = * ) 'ptop              top of analysis, lid'
      WRITE ( UNIT = * , FMT = * ) '-------------------------------'

   ELSE present_help

      !  The values that are requested should be printed out in here
      !  since we can assume the value of the argument is interesting.

      IF ( print_header ) WRITE ( UNIT = * , &
      FMT = '(/"  Data values available from header:")')

      !  These may be needed for some of the other requests, so they have
      !  to be done prior to the other assignments.

      IF ( PRESENT ( grid_id       ) ) THEN
         grid_id       = bhi( 1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'grid_id = ',  grid_id 
      ENDIF
      
      IF ( PRESENT ( map_projection   ) ) THEN
         map_projection   = bhi( 2)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'map_projection = ',  map_projection 
      ENDIF
      
      kbu     = bhi(3)
      IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'kbu = ',  kbu 

      IF ( PRESENT ( iewe             ) ) THEN
         iewe             = bhi( 4)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'iewe = ',  iewe 
      ENDIF

      IF ( PRESENT ( jnse             ) ) THEN
         jnse             = bhi( 5)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'jnse = ',  jnse 
      ENDIF

      IF ( PRESENT ( expanded         ) ) THEN
         expanded         = bhi( 6)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'expanded = ',  expanded
      ENDIF
      

      IF ( PRESENT ( iewd             ) ) THEN
         iewd             = bhi( 7)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'iewd = ',  iewd 
      ENDIF

      IF ( PRESENT ( jnsd             ) ) THEN
         jnsd             = bhi( 8)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'jnsd = ',  jnsd 
      ENDIF



      IF ( PRESENT ( dxc              ) ) THEN
         dxc              = bhr( 1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'dxc = ',  dxc 
      ENDIF

      IF ( PRESENT ( lat_center       ) ) THEN
         lat_center       = bhr( 2)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'lat_center = ',  lat_center 
      ENDIF

      IF ( PRESENT ( lon_center       ) ) THEN
         lon_center       = bhr( 3)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'lon_center = ',  lon_center 
      ENDIF

      IF ( PRESENT ( cone_factor      ) ) THEN
         cone_factor      = bhr( 4)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'cone_factor = ',  cone_factor 
      ENDIF

      IF ( PRESENT ( true_lat1        ) ) THEN
         true_lat1        = bhr( 5)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'true_lat1 = ',  true_lat1 
      ENDIF

      IF ( PRESENT ( true_lat2        ) ) THEN
         true_lat2        = bhr( 6)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'true_lat2 = ',  true_lat2 
      ENDIF

      IF ( PRESENT ( pole             ) ) THEN
         pole             = bhr( 7)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'pole = ',  pole
      ENDIF

      IF ( PRESENT ( dxd              ) ) THEN
         dxd              = bhr( 8)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'dxd = ',  dxd 
      ENDIF

      IF ( PRESENT ( ptop             ) ) THEN
         ptop             = bhr( 9)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'ptop = ',  ptop 
      ENDIF

      IF ( print_header ) WRITE ( UNIT = * , FMT = '(//)' )

   END IF present_help
   
END SUBROUTINE proc_get_info_header
   
