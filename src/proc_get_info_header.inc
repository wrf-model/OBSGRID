   INTEGER            :: iewd              , & ! domain i size
                         jnsd              , & ! domain j size
                         kbu               , & ! num vert levs, bottom up
                         grid_id           , & ! grid ID
                         map_projection    , & ! 1=LC, 2=PS, 3=MC
                         expanded          , & ! 1=expanded; 0=not expanded
                         iewe              , & ! expanded grid i size
                         jnse                  ! expanded grid j size
                                                        
   REAL               :: dxc               , & ! coarse grid distance
                         lat_center        , & ! center latitude
                         lon_center        , & ! central longitude
                         cone_factor       , & ! cone factor
                         true_lat1         , & ! lat of true projection
                         true_lat2         , & ! lat of true projection
                         pole              , & ! location of pole in latitude
                         dxd                   ! domain grid distance

   REAL               :: ptop                  ! top of analysis, lid
   
   INTEGER            :: num3d             , & ! num 3D arrays
                         num2d             , & ! num 2D arrays
                         num1d                 ! num 1D arrays

   REAL               , DIMENSION ( 1000 )   &
                      :: pressure              ! array of pres (bottom up)
                                                        
                                                        
   LOGICAL            :: help                  ! T/F: print stuff only
