   REAL , PARAMETER  :: g    =    9.81  , & !  acceleration of gravity ( m/s **2)
                        cp   = 1004.    , & !  specific heat at constant pressure
                        gasr =  287.95  , & !  gas constant for dry air
                        rcp  = gasr / cp       
!BPR BEGIN
   REAL , PARAMETER  :: L_over_Rv = 5418.12
   REAL , PARAMETER  :: Rv_over_L = 1. / L_over_Rv
   REAL , PARAMETER  :: t_at_sea_level_std_atm = 288.15
   REAL , PARAMETER  :: dTdz_std_atm_lt_11000m = -0.0065
   REAL , PARAMETER  :: very_small_rh = 0.0001
!BPR END
