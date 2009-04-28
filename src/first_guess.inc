   REAL , DIMENSION(jns_alloc,iew_alloc,kbu_alloc) &
                                   :: t  , &
                                      u  , &
                                      v  , &
                                      uA , &
                                      vA , &
                                      uC , &
                                      vC , &
                                      h  , &
                                      rh , &
                                      pres

   REAL , DIMENSION(jns_alloc,iew_alloc) &
                                   :: terrain      , &
                                      latitude_x   , &
                                      longitude_x  , &
                                      latitude_d   , &
                                      longitude_d  , &
                                      slp_x        , &
                                      slp_C        , &
                                      snow         , &
                                      sst
