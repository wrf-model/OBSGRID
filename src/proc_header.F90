!------------------------------------------------------------------------------
SUBROUTINE proc_header ( filename , bhid , bhrd , nml )

!  This routine processes the modeling system record header on input.

   USE namelist

   IMPLICIT NONE
   INCLUDE 'netcdf.inc'

   CHARACTER ( LEN = 132 ) ,  INTENT ( IN )      :: filename
   TYPE ( all_nml ) , INTENT ( IN )              :: nml
   INTEGER , DIMENSION (50)                      :: bhid
   REAL    , DIMENSION (20)                      :: bhrd

   INTEGER                                       :: i, idummy
   REAL                                          :: rdummy
   INTEGER                                       :: ntimes, wedim, sndim, btdim
   INTEGER                                       :: rcode
   INTEGER                                       :: met_ncid
   INTEGER                                       :: ndims, nvars, ngatts, nunlimdimid
   INTEGER                                       :: dim_val
   CHARACTER (LEN=31)                            :: dim_name
   CHARACTER (LEN=24), ALLOCATABLE, DIMENSION(:) :: times
   REAL, ALLOCATABLE, DIMENSION(:,:,:,:)         :: ptop

   INTEGER                :: map_proj
   REAL                   :: truelat1, truelat2, cone
   REAL, PARAMETER        :: PI = 3.141592653589793
   REAL, PARAMETER        :: RAD_PER_DEG = PI/180.


   !  The record header is initialized with the call to this
   !  routine.  That data is assigned to these arrays that
   !  are accessible from a COMMON block in the proc_get_info_header
   !  routine.  

   INCLUDE 'header.common'

   rcode = nf_open(filename, 0, met_ncid)

   rcode = nf_inq(met_ncid, ndims, nvars, ngatts, nunlimdimid)
   dims_loop : DO i = 1, ndims
      rcode = nf_inq_dim(met_ncid, i, dim_name, dim_val)
      IF ( dim_name == 'Time'              ) ntimes = dim_val
      IF ( dim_name == 'west_east_stag'    ) wedim  = dim_val
      IF ( dim_name == 'south_north_stag'  ) sndim  = dim_val
      IF ( dim_name == 'num_metgrid_levels') btdim  = dim_val
   ENDDO dims_loop

   ALLOCATE ( times (ntimes) )
   rcode = nf_inq_varid    ( met_ncid, "Times", i )
   rcode = nf_get_var_text ( met_ncid, i, times )
   ALLOCATE ( ptop (wedim-1, sndim-1, btdim, ntimes) )
   rcode = nf_inq_varid    ( met_ncid, "PRES", i )
   rcode = nf_get_var_real ( met_ncid, i, ptop )
   rcode = NF_GET_ATT_INT(met_ncid, nf_global, "MAP_PROJ", map_proj )


   bhid( 1) = nml%record_2%grid_id      !!! grid ID
   bhid( 2) = map_proj
   bhid( 3) = btdim                     !!! number of vertical levels in 3d data
   bhid( 4) = wedim                     !!! full incoming we dim (staggered)
   bhid( 5) = sndim                     !!! full incoming sn dim (staggered)

   if ( nml%record_2%trim_domain ) then
      bhid( 6) = 1                         !!! domain expanded
      bhid( 7) = wedim - 2 * nml%record_2%trim_value 
      bhid( 8) = sndim - 2 * nml%record_2%trim_value
   else
      bhid( 6) = 0                         !!! domain not expanded
      bhid( 7) = wedim                  
      bhid( 8) = sndim                 
   endif


   rcode = NF_GET_ATT_REAL(met_ncid, nf_global, "DX", rdummy )
        bhrd( 1) = rdummy               !!! dx
        bhrd( 8) = rdummy               !!! dx
   rcode = NF_GET_ATT_REAL(met_ncid, nf_global, "CEN_LAT", rdummy )
        bhrd( 2) = rdummy               !!! center lat
   rcode = NF_GET_ATT_REAL(met_ncid, nf_global, "STAND_LON", rdummy )
        bhrd( 3) = rdummy               !!! stand lon
   rcode = NF_GET_ATT_REAL(met_ncid, nf_global, "TRUELAT1", truelat1 )
        bhrd( 5) = truelat1             !!! truelat1
   rcode = NF_GET_ATT_REAL(met_ncid, nf_global, "TRUELAT2", truelat2 )
        bhrd( 6) = truelat2             !!! truelat2

        bhrd( 9) =  ptop(1,1,btdim,1)   !!! Top pressure used in analysis, pressure defining model lid (Pa)


    cone = 1.
    if ( map_proj .eq. 1 ) then
      IF (ABS(truelat1-truelat2) .GT. 0.1) THEN
         cone=(ALOG(COS(truelat1*RAD_PER_DEG))-            &
               ALOG(COS(truelat2*RAD_PER_DEG))) /          &
         (ALOG(TAN((90.-ABS(truelat1))*RAD_PER_DEG*0.5 ))- &
          ALOG(TAN((90.-ABS(truelat2))*RAD_PER_DEG*0.5 )) )
      ELSE
         cone = SIN(ABS(truelat1)*RAD_PER_DEG )
      ENDIF
    endif
        bhrd( 4) = cone                 !!! cone factor

   rcode = nf_close(met_ncid) 

   !  The reason to put the header data in COMMON is if it happens to be filled!

   bhi  = bhid
   bhr  = bhrd

END SUBROUTINE proc_header
