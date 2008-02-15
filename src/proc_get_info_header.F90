!------------------------------------------------------------------------------

SUBROUTINE proc_get_info_header ( print_header , program , &
iewc , jnsc , map_projection , expanded , iewe , jnse , &
dxc , lat_center , lon_center , cone_factor , true_lat1 , true_lat2 , pole , &
domain_id , mother_id , nest_level , &
iewd , jnsd , iew_startm , jns_startm , &
ratio_wrt_coarse , ratio_wrt_mother , &
dxd , xew_startc , yns_startc , &
kbu , &
ptop , help )
   
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
   
   INTEGER , OPTIONAL :: iewc              , & ! coarse grid i size 
                         jnsc              , & ! coarse grid j size
                         map_projection    , & ! 1=LC, 2=PS, 3=MC
                         expanded          , & ! 1=expanded; 0=not expanded
                         iewe              , & ! expanded grid i size
                         jnse              , & ! expanded grid j size
                         domain_id         , & ! domain ID
                         mother_id         , & ! mother ID
                         nest_level        , & ! nest level (0=coarse)
                         iewd              , & ! domain i size
                         jnsd              , & ! domain j size
                         iew_startm        , & ! LL location in mother i
                         jns_startm        , & ! LL location in mother j
                         ratio_wrt_coarse  , & ! grid ratio wrt coarse
                         ratio_wrt_mother      ! grid ratio wrt mother
                                                        
   !  REAL information 
   
   REAL ,    OPTIONAL :: dxc               , & ! coarse grid distance
                         lat_center        , & ! center latitude
                         lon_center        , & ! central longitude
                         cone_factor       , & ! cone factor
                         true_lat1         , & ! lat of true projection
                         true_lat2         , & ! lat of true projection
                         pole              , & ! location of pole in latitude
                         dxd               , & ! domain grid distance
                         xew_startc        , & ! LL location in coarse x
                         yns_startc            ! LL location in coarse y 
   
   !  INTEGER information from input program.
   
   INTEGER , OPTIONAL :: program           , & ! 2=REGRID, 3=little_r, 8=interpolated model output
                         kbu                   ! num vert levs, bottom up

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
      WRITE ( UNIT = * , FMT = * ) 'iewc              coarse grid i size '
      WRITE ( UNIT = * , FMT = * ) 'jnsc              coarse grid j size'
      WRITE ( UNIT = * , FMT = * ) 'map_projection    1=LC, 2=PS, 3=MC'
      WRITE ( UNIT = * , FMT = * ) 'expanded          1=expanded; 0=not expanded'
      WRITE ( UNIT = * , FMT = * ) 'iewe              expanded grid i size'
      WRITE ( UNIT = * , FMT = * ) 'jnse              expanded grid j size'
      WRITE ( UNIT = * , FMT = * ) 'domain_id         domain ID'
      WRITE ( UNIT = * , FMT = * ) 'mother_id         mother ID'
      WRITE ( UNIT = * , FMT = * ) 'nest_level        nest level (0=coarse)'
      WRITE ( UNIT = * , FMT = * ) 'iewd              domain i size'
      WRITE ( UNIT = * , FMT = * ) 'jnsd              domain j size'
      WRITE ( UNIT = * , FMT = * ) 'iew_startm        LL location in mother i'
      WRITE ( UNIT = * , FMT = * ) 'jns_startm        LL location in mother j'
      WRITE ( UNIT = * , FMT = * ) 'ratio_wrt_coarse  grid ratio wrt coarse'
      WRITE ( UNIT = * , FMT = * ) 'ratio_wrt_mother  grid ratio wrt mother'
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
      WRITE ( UNIT = * , FMT = * ) 'xew_startc        LL location in coarse x'
      WRITE ( UNIT = * , FMT = * ) 'yns_startc        LL location in coarse y '
      WRITE ( UNIT = * , FMT = * ) '-------------------------------'
      WRITE ( UNIT = * , FMT = * ) 'INTEGER information from REGRID, little_r or interpolated MM5.'
      WRITE ( UNIT = * , FMT = * ) 'program           2=REGRID, 3=little_r, 8=interpolated MM5'
      WRITE ( UNIT = * , FMT = * ) 'kbu               num vert levels'
      WRITE ( UNIT = * , FMT = * ) '-------------------------------'
      WRITE ( UNIT = * , FMT = * ) 'REAL information from REGRID.'
      WRITE ( UNIT = * , FMT = * ) 'ptop              top of analysis, lid'
      WRITE ( UNIT = * , FMT = * ) '-------------------------------'

   ELSE present_help

      !  The values that are requested should be printed out in here
      !  since we can assume the value of the argument is interesting.

      IF ( print_header ) WRITE ( UNIT = * , &
      FMT = '(/"  Data values available from header:")')

      !  These may be needed for some of the other requests, so they have
      !  to be done prior to the other assignments.
      
      program = bhi(  1,1)
      IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'program = ',  program 
      kbu     = bhi(12,program)
      IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'kbu = ',  kbu 
      
      !  The available choices, pulled from the record header.
      
      IF ( PRESENT ( iewc             ) ) THEN
         iewc             = bhi(  6,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'iewc = ',  iewc 
      ENDIF

      IF ( PRESENT ( jnsc             ) ) THEN
         jnsc             = bhi(  5,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'jnsc = ',  jnsc 
      ENDIF

      IF ( PRESENT ( map_projection   ) ) THEN
         map_projection   = bhi(  7,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'map_projection = ',  map_projection 
      ENDIF

      IF ( PRESENT ( expanded         ) ) THEN
         expanded         = bhi(  8,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'expanded = ',  expanded
      ENDIF

      IF ( PRESENT ( iewe             ) ) THEN
         iewe             = bhi( 10,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'iewe = ',  iewe 
      ENDIF

      IF ( PRESENT ( jnse             ) ) THEN
         jnse             = bhi(  9,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'jnse = ',  jnse 
      ENDIF

      IF ( PRESENT ( domain_id        ) ) THEN
         domain_id        = bhi( 13,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'domain_id = ',  domain_id 
      ENDIF

      IF ( PRESENT ( mother_id        ) ) THEN
         mother_id        = bhi( 14,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'mother_id = ',  mother_id 
      ENDIF

      IF ( PRESENT ( nest_level       ) ) THEN
         nest_level       = bhi( 15,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'nest_level = ',  nest_level 
      ENDIF

      IF ( PRESENT ( iewd             ) ) THEN
         iewd             = bhi( 17,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'iewd = ',  iewd 
      ENDIF

      IF ( PRESENT ( jnsd             ) ) THEN
         jnsd             = bhi( 16,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'jnsd = ',  jnsd 
      ENDIF

      IF ( PRESENT ( iew_startm       ) ) THEN
         iew_startm       = bhi( 19,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'iew_startm = ',  iew_startm 
      ENDIF

      IF ( PRESENT ( jns_startm       ) ) THEN
         jns_startm       = bhi( 18,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'jns_startm = ',  jns_startm 
      ENDIF

      IF ( PRESENT ( ratio_wrt_coarse ) ) THEN
         ratio_wrt_coarse = bhi( 20,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'ratio_wrt_coarse = ',  ratio_wrt_coarse 
      ENDIF

      IF ( PRESENT ( ratio_wrt_mother ) ) THEN
         ratio_wrt_mother = bhi( 21,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'ratio_wrt_mother = ',  ratio_wrt_mother 
      ENDIF

      IF ( PRESENT ( dxc              ) ) THEN
         dxc              = bhr(  1,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'dxc = ',  dxc 
      ENDIF

      IF ( PRESENT ( lat_center       ) ) THEN
         lat_center       = bhr(  2,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'lat_center = ',  lat_center 
      ENDIF

      IF ( PRESENT ( lon_center       ) ) THEN
         lon_center       = bhr(  3,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'lon_center = ',  lon_center 
      ENDIF

      IF ( PRESENT ( cone_factor      ) ) THEN
         cone_factor      = bhr(  4,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'cone_factor = ',  cone_factor 
      ENDIF

      IF ( PRESENT ( true_lat1        ) ) THEN
         true_lat1        = bhr(  5,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'true_lat1 = ',  true_lat1 
      ENDIF

      IF ( PRESENT ( true_lat2        ) ) THEN
         true_lat2        = bhr(  6,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'true_lat2 = ',  true_lat2 
      ENDIF

      IF ( PRESENT ( pole             ) ) THEN
         pole             = bhr(  7,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'pole = ',  pole
      ENDIF

      IF ( PRESENT ( dxd              ) ) THEN
         dxd              = bhr(  9,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'dxd = ',  dxd 
      ENDIF

      IF ( PRESENT ( xew_startc       ) ) THEN
         xew_startc       = bhr( 11,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'xew_startc = ',  xew_startc 
      ENDIF

      IF ( PRESENT ( yns_startc       ) ) THEN
         yns_startc       = bhr( 10,1)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'yns_startc = ',  yns_startc 
      ENDIF

      IF ( PRESENT ( ptop             ) ) THEN
         ptop             = bhr(  2,2)
         IF ( print_header ) WRITE ( UNIT = * , FMT = * ) 'ptop = ',  ptop 
      ENDIF

      IF ( print_header ) WRITE ( UNIT = * , FMT = '(//)' )

   END IF present_help
   
END SUBROUTINE proc_get_info_header
   
