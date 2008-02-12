!------------------------------------------------------------------------------
SUBROUTINE proc_header ( filename , bhid , bhrd , nml )

!  This routine processes the modeling system record header on input.

   USE namelist

   IMPLICIT NONE
   INCLUDE 'netcdf.inc'

   CHARACTER ( LEN = 132 ) ,INTENT ( IN )        :: filename
   INTEGER                , INTENT (OUT ) , &
                            DIMENSION (50,20)    :: bhid
   REAL                   , INTENT (OUT ) , &
                            DIMENSION (20,20)    :: bhrd
   TYPE ( all_nml ) , INTENT ( IN )              :: nml

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


        bhid( 1,1) = 2
        bhid( 5,1) = sndim                     !!! sn dim (unstaggered)
        bhid( 6,1) = wedim                     !!! we dim (unstaggered)
   rcode = NF_GET_ATT_INT(met_ncid, nf_global, "MAP_PROJ", map_proj )
        bhid( 7,1) = map_proj
        bhid( 9,1) = bhid( 5,1)
        bhid(10,1) = bhid( 6,1)
   rcode = NF_GET_ATT_INT(met_ncid, nf_global, "grid_id", idummy )
        bhid(13,1) = idummy                    !!! domain ID
        bhid(15,1) = 0                         !!! nest level
        IF ( idummy > 1 ) bhid(15,1) = 1
   rcode = NF_GET_ATT_INT(met_ncid, nf_global, "parent_id", idummy )
        bhid(14,1) = idummy                    !!! mother domain ID
        bhid(16,1) = bhid( 5,1)
        bhid(17,1) = bhid( 6,1)
   rcode = NF_GET_ATT_INT(met_ncid, nf_global, "j_parent_start", idummy )
        bhid(18,1) = idummy                    !!! location in mother domain
        bhrd(10,1) = real(idummy) 
   rcode = NF_GET_ATT_INT(met_ncid, nf_global, "i_parent_start", idummy )
        bhid(19,1) = idummy                    !!! location in mother domain
        bhrd(11,1) = real(idummy) 
   rcode = NF_GET_ATT_INT(met_ncid, nf_global, "parent_grid_ratio", idummy )
        bhid(20,1) = idummy                    !!! grid ratio wrt coarse domain
        bhid(21,1) = idummy                    !!! grid ratio wrt mother domain

   READ(times(1)(1:4),'(i4)')idummy
        bhid( 5,2) = idummy               !!! year of the start time
   READ(times(1)(6:7),'(i2)')idummy
        bhid( 6,2) = idummy               !!! month of the start time
   READ(times(1)(9:10),'(i2)')idummy
        bhid( 7,2) = idummy               !!! day of the start time
   READ(times(1)(12:13),'(i2)')idummy
        bhid( 8,2) = idummy               !!! hour of the start time
   READ(times(1)(15:16),'(i2)')idummy
        bhid( 9,2) = idummy               !!! minute of the start time
   READ(times(1)(18:19),'(i2)')idummy
        bhid(10,2) = idummy               !!! second of the start time

        bhid(12,2) = btdim                !!! Anticipated number of vertical levels in 3d data

   rcode = NF_GET_ATT_REAL(met_ncid, nf_global, "DX", rdummy )
        bhrd( 1,1) = rdummy               !!! dx
        bhrd( 9,1) = rdummy               !!! dx
   rcode = NF_GET_ATT_REAL(met_ncid, nf_global, "CEN_LAT", rdummy )
        bhrd( 2,1) = rdummy               !!! center lat
   rcode = NF_GET_ATT_REAL(met_ncid, nf_global, "STAND_LON", rdummy )
        bhrd( 3,1) = rdummy               !!! stand lon
   rcode = NF_GET_ATT_REAL(met_ncid, nf_global, "TRUELAT1", truelat1 )
        bhrd( 5,1) = truelat1             !!! truelat1
   rcode = NF_GET_ATT_REAL(met_ncid, nf_global, "TRUELAT2", truelat2 )
        bhrd( 6,1) = truelat2             !!! truelat2

        bhrd( 2,2) =  ptop(1,1,btdim,1)   !!! Top pressure used in analysis, pressure defining model lid (Pa)


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
        bhrd( 4,1) = cone                 !!! cone factor

   rcode = nf_close(met_ncid) 

   !  The reason to put the header data in COMMON is if it happens to be filled!

   bhi  = bhid
   bhr  = bhrd

END SUBROUTINE proc_header
