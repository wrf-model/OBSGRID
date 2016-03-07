   !  Record 1 NAMELIST variables.

   INTEGER                 ::  start_year                   , &
                               start_month                  , &
                               start_day                    , &
                               start_hour                   , &
                               start_minute                 , &
                               start_second                 , &
                                 end_year                   , &
                                 end_month                  , &
                                 end_day                    , &
                                 end_hour                   , &
                                 end_minute                 , &
                                 end_second               
   INTEGER                 :: interval    

   !  Record 2 NAMELIST variables.

   CHARACTER ( LEN = 132 ) :: fg_filename     
   CHARACTER ( LEN = 132 ) ::  obs_filename    
   LOGICAL :: trim_domain, remove_unverified_data
   INTEGER :: trim_value, grid_id, remove_data_above_qc_flag

   !  Record 3 NAMELIST variables.

   INTEGER                 :: max_number_of_obs    

   LOGICAL                 :: fatal_if_exceed_max_obs    

   !  Record 4 NAMELIST variables.

   LOGICAL                 :: qc_test_error_max            , &
                              qc_test_buddy                , &
                              qc_test_vert_consistency     , &
   !  BPR BEGIN
   !                          qc_test_convective_adj       
                              qc_test_convective_adj       , &
                              qc_psfc
   !  BPR END

   REAL                    :: max_error_t                  , &
                              max_error_uv                 , &
                              max_error_z                  , &
                              max_error_rh                 , &
   ! BPR BEGIN
                              max_error_dewpoint           , &
   ! BPR END  
                              max_error_p                  , &
                              max_buddy_t                  , &
                              max_buddy_uv                 , &
                              max_buddy_z                  , &
                              max_buddy_rh                 , &
   ! BPR BEGIN
                              max_buddy_dewpoint           , &
   ! BPR END
                              max_buddy_p                  , &
                              buddy_weight                 , &
                              max_p_extend_t               , &
                              max_p_extend_w    
!BPR BEGIN
     LOGICAL               :: use_p_tolerance_one_lev   
     INTEGER               :: max_p_tolerance_one_lev_qc
!BPR END


   !  Record 5 NAMELIST variables.
 
   LOGICAL                 :: print_found_obs              , & 
                              print_header                 , & 
                              print_analysis               , &                
                              print_qc_vert                , &
                              print_qc_dry                 , & 
                              print_obs_files              , &
                              print_error_max              , &
                              print_buddy                  , &
                              print_oa    

   !  Record 7 NAMELIST variables.

   LOGICAL                 :: use_first_guess              , &
                              f4d                          , &
                              lagtem

   INTEGER                 :: intf4d

   !  Record 8 NAMELIST variables.

   INTEGER                 :: smooth_type                  , &
                              smooth_sfc_wind              , &
                              smooth_sfc_temp              , &
                              smooth_sfc_rh                , &
                              smooth_sfc_slp               , &
                              smooth_upper_wind            , &
                              smooth_upper_temp            , &
                              smooth_upper_rh      

   !  Record 9 NAMELIST variables.

   CHARACTER ( LEN = 132 ) :: oa_type    , oa_3D_type

   INTEGER                 :: mqd_minimum_num_obs          , &
                              mqd_maximum_num_obs    
   INTEGER , DIMENSION(10) :: radius_influence    

   LOGICAL                 :: oa_min_switch                , &
                              oa_max_switch
   INTEGER                 :: oa_3D_option    
   !BPR BEGIN
   LOGICAL                 :: oa_psfc
   LOGICAL                 :: scale_cressman_rh_decreases
   REAL                    :: radius_influence_sfc_mult
   INTEGER                 :: max_p_tolerance_one_lev_oa
   !BPR END 
