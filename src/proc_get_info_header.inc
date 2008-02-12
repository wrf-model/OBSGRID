   INTEGER            :: iewc              , & ! coarse grid i size 
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
                                                        
   REAL               :: dxc               , & ! coarse grid distance
                         lat_center        , & ! center latitude
                         lon_center        , & ! central longitude
                         cone_factor       , & ! cone factor
                         true_lat1         , & ! lat of true projection
                         true_lat2         , & ! lat of true projection
                         pole              , & ! location of pole in latitude
                         dxd               , & ! domain grid distance
                         xew_startc        , & ! LL location in coarse x
                         yns_startc            ! LL location in coarse y 
   
   INTEGER            :: program           , & ! 2=DATAGRID, 3=RAWINS
                         kbu               , & ! num vert levs, bottom up
                         num3d             , & ! num 3D arrays
                         num2d             , & ! num 2D arrays
                         num1d                 ! num 1D arrays
   REAL               , DIMENSION ( 1000 )   &
                      :: pressure              ! array of pres (bottom up)
                                                        
   REAL               :: ptop                  ! top of analysis, lid
                                                        
   LOGICAL            :: help                  ! T/F: print stuff only
