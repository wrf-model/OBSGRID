   !  The following parameters are used as flags to indicate either a
   !  MISSING data value (real number) in the "measurement" TYPE or an
   !  "end of data" flag (two fields with end_data  flag mean end of
   !  measurements for one observation).  The undefined values are when
   !  a GTS observation is read but determined to be in error (not according
   !  to the code).  Geopotential height is referenced as a positive value
   !  undefined, which explains the duplicate.

   INTEGER , PARAMETER                            ::  undefined1    =  999999
   REAL    , PARAMETER                            ::  undefined1_r  =  999999.
   INTEGER , PARAMETER                            ::  undefined2    = -999999
   REAL    , PARAMETER                            ::  undefined2_r  = -999999.
   INTEGER , PARAMETER                            ::  missing       = -888888
   REAL    , PARAMETER                            ::  missing_r     = -888888.
   INTEGER , PARAMETER                            ::  end_data      = -777777
   REAL    , PARAMETER                            ::  end_data_r    = -777777.

   !  Define various ways for bad data to be flagged.  Each of these is
   !  2^n, so that the errors may be summed.

   INTEGER , PARAMETER                            ::  p_from_h_first_guess           = 2 **  1
   INTEGER , PARAMETER                            ::  p_std_atm_and_height           = 2 **  2
   INTEGER , PARAMETER                            ::  p_std_atm_and_temp             = 2 **  3

   INTEGER , PARAMETER                            ::  zero_t_td                      = 2 **  4
   INTEGER , PARAMETER                            ::  zero_spd                       = 2 **  5
   INTEGER , PARAMETER                            ::  negative_spd                   = 2 **  6
   INTEGER , PARAMETER                            ::  wrong_direction                = 2 **  7

   INTEGER , PARAMETER                            ::  vert_interpolated              = 2 **  8
   INTEGER , PARAMETER                            ::  extend_influence               = 2 **  9

   INTEGER,  PARAMETER                            ::  wrong_t_sign                   = 2 ** 10
   INTEGER,  PARAMETER                            ::  t_fail_supa_inver              = 2 ** 11
   INTEGER,  PARAMETER                            ::  wrong_wind_data                = 2 ** 12
   INTEGER,  PARAMETER                            ::  convective_adjustment          = 2 ** 13
   INTEGER , PARAMETER                            ::  no_buddies                     = 2 ** 14

   INTEGER , PARAMETER                            ::  no_qc_possible                 = 2 ** 15

   INTEGER , PARAMETER                            ::  fails_error_max                = 2 ** 16
   INTEGER , PARAMETER                            ::  fails_buddy_check              = 2 ** 17
   INTEGER , PARAMETER                            ::  outside_of_domain              = 2 ** 18
