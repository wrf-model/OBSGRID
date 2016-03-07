MODULE final_analysis

   USE first_guess

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE write_analysis ( met_ncid, oa_ncid , &
t , u , v , uA , vA , uC , vC , h , rh , pres , terrain , &
latitude_x , longitude_x , latitude_d , longitude_d , &
slp_x , slp_C , sst , tobbox , odis , &
iew_alloc , jns_alloc , kbu_alloc , iewd , jnsd , date_char, print_analysis , &
!BPR BEGIN
!oa_type , oa_3D_type , oa_3D_option , rad_influence )
oa_type , oa_3D_type , oa_3D_option , rad_influence , oa_psfc )
!BPR END

!  This routine assembles the correct data fields and outputs them with the
!  appropriate small header.

   USE input_data

   IMPLICIT NONE 
   INCLUDE 'netcdf.inc'

   INCLUDE 'big_header.inc'
   INCLUDE 'first_guess_size.inc'
   INCLUDE 'first_guess.inc'
   INCLUDE 'error.inc'
 
   INTEGER                                                :: met_ncid, oa_ncid
   REAL , DIMENSION ( jns_alloc , iew_alloc )             :: tobbox , odis
   INTEGER                                                :: iewd , jnsd
   REAL , DIMENSION ( jns_alloc , iew_alloc , kbu_alloc ) :: pertubation
   REAL , DIMENSION ( jnsd , iewd , kbu_alloc )           :: dum3d
   REAL , DIMENSION ( iewd , jnsd , kbu_alloc )           :: met_em_dum3d
   REAL , DIMENSION ( jnsd , iewd )                       :: dum2d
   REAL , DIMENSION ( iewd , jnsd )                       :: met_em_dum2d
   REAL , ALLOCATABLE , DIMENSION(:,:,:)                  :: met_em_3d
   REAL , ALLOCATABLE , DIMENSION(:,:)                    :: met_em_2d
   CHARACTER (LEN=19)                                     :: date_char

   CHARACTER *(*)                                 :: oa_type , oa_3D_type
   INTEGER , DIMENSION(10)                        :: rad_influence
   INTEGER                                        :: oa_3D_option
   !BPR BEGIN
   LOGICAL                                        :: oa_psfc
   !BPR END

   INTERFACE 
      INCLUDE 'error.int'
      INCLUDE 'sample.int'
   END INTERFACE
   
   INTEGER :: loop_count , &
              num3d , &
              num2d , &
              num1d , &
              num0d
   REAL,PARAMETER    :: dummy_value = -99999.
   LOGICAL :: print_analysis

   INTEGER                                           :: ivar_number
   INTEGER                                           :: kp, imet, jmet
   INTEGER                                           :: rcode, i, na, ii
   INTEGER                                           :: ndims, nvars, ngatts, nunlimdimid
   CHARACTER (LEN=31),ALLOCATABLE, DIMENSION(:)      :: dname
   INTEGER,           ALLOCATABLE, DIMENSION(:)      :: dval, dval2
   CHARACTER (LEN=31)                                :: cname
   CHARACTER (LEN=50)                                :: cval
   INTEGER                                           :: ilen, itype, ival, idm, natt
   INTEGER, DIMENSION(6)                             :: ishape
   REAL                                              :: rval
   REAL, DIMENSION(16)                               :: rval_corners, corner_lats, corner_lons
   CHARACTER (LEN=19), ALLOCATABLE, DIMENSION(:,:,:) :: text
   REAL, ALLOCATABLE, DIMENSION(:,:,:)               :: d3_data1, d3_data2
   REAL, ALLOCATABLE, DIMENSION(:,:)                 :: d2_data1, d2_data2
   INTEGER, DIMENSION(4)                             :: start_dims, dims_in, dims_out
   CHARACTER (LEN=16)                                :: rads
   INTEGER                                           :: i_rad


   !  We need to keep track of where we are sticking the data for the FDDA option.

   IF ( initial_time ) THEN
      tt = first_time
   ELSE
      tt = second_time
   END IF

   !! First get all the dimension, global variables and unchanged fields from the input file

   ! find out some details about the output file
   rcode = nf_inq(met_ncid, ndims, nvars, ngatts, nunlimdimid)

   ! sort out the dimensions
   IF (ALLOCATED(dname)) DEALLOCATE(dname)
       ALLOCATE (dname(ndims  ))
   IF (ALLOCATED(dval)) DEALLOCATE(dval)
       ALLOCATE  (dval(ndims  ))
   IF (ALLOCATED(dval2)) DEALLOCATE(dval2)
       ALLOCATE (dval2(ndims))
   DO i = 1, ndims
       rcode = nf_inq_dim(met_ncid, i, dname(i), dval(i))
       dval2(i)  = dval(i)
           IF (dname(i) == 'west_east') THEN
               dval2(i) = iewd-1
       ELSEIF (dname(i) == 'west_east_stag') THEN
               dval2(i) = iewd
       ELSEIF (dname(i) == 'south_north') THEN
               dval2(i) = jnsd-1
       ELSEIF (dname(i) == 'south_north_stag') THEN
               dval2(i) = jnsd
       ENDIF

       if ( dname(i) == "Time" ) then
         rcode = nf_def_dim(oa_ncid, dname(i), NF_UNLIMITED, i)
       else
         rcode = nf_def_dim(oa_ncid, dname(i), dval2(i), i)
       end if
   ENDDO

   ! deal with the global variables
   DO i = 1, ngatts
          rcode = nf_inq_attname(met_ncid, NF_GLOBAL, i,    cname)
          rcode = nf_inq_atttype(met_ncid, NF_GLOBAL, cname, itype)
          rcode = nf_inq_attlen (met_ncid, NF_GLOBAL, cname, ilen)

          IF ( itype == 2 ) THEN        ! characters
            rcode = nf_get_att_text (met_ncid, NF_GLOBAL, cname, cval)
            IF(cname(1:5) .eq. 'TITLE') then
               cval = "OUTPUT FROM OBSGRID"
               ilen = len_trim(cval)
            ENDIF 
            rcode = nf_put_att_text(oa_ncid, NF_GLOBAL, cname, ilen, cval(1:ilen))

          ELSEIF ( itype == 4 ) THEN     ! integers
            rcode = nf_get_att_int (met_ncid, NF_GLOBAL, cname, ival)
                IF(cname == 'WEST-EAST_GRID_DIMENSION' .OR. &
                   cname == 'WEST-EAST_PATCH_END_STAG' .OR. &
                   cname == 'i_parent_end'             ) THEN
                   ival = iewd
            ELSEIF(cname == 'SOUTH-NORTH_GRID_DIMENSION' .OR. &
                   cname == 'SOUTH-NORTH_PATCH_END_STAG' .OR. &
                   cname == 'j_parent_end'             ) THEN
                   ival = jnsd
            ELSEIF(cname == 'WEST-EAST_PATCH_END_UNSTAG'    ) THEN
                   ival = iewd-1
            ELSEIF(cname == 'SOUTH-NORTH_PATCH_END_UNSTAG'   ) THEN
                   ival = jnsd-1
            ENDIF
            rcode = nf_put_att_int(oa_ncid, NF_GLOBAL, cname, itype, ilen, ival)

          ELSEIF ( itype == 5 ) THEN    ! real
            IF ( cname == 'corner_lats' .OR. cname == 'corner_lons' ) THEN
               rcode = nf_get_att_real (met_ncid, NF_GLOBAL, cname, rval_corners)
               rcode = nf_put_att_real(oa_ncid, NF_GLOBAL, cname, itype, ilen, rval_corners)
            ELSE
               rcode = nf_get_att_real (met_ncid, NF_GLOBAL, cname, rval)
               rcode = nf_put_att_real(oa_ncid, NF_GLOBAL, cname, itype, ilen, rval)
            ENDIF 
          ENDIF 
   ENDDO     
   ilen = len_trim(oa_type)
   rcode = nf_put_att_text(oa_ncid, NF_GLOBAL, "OA_SCHEME", ilen, oa_type(1:ilen))
   ilen = len_trim(oa_3D_type)
   rcode = nf_put_att_text(oa_ncid, NF_GLOBAL, "OA_3D_SCHEME", ilen, oa_3D_type(1:ilen))
   rcode = nf_put_att_int(oa_ncid, NF_GLOBAL, "OA_3D_option", 4, 1, oa_3D_option)
   DO i_rad = 1,10
     IF ( rad_influence(i_rad) > 0 ) THEN
       WRITE(rads,'("RAD_INFLUENCE_",i2.2)')i_rad
       rcode = nf_put_att_int(oa_ncid, NF_GLOBAL, rads, 4, 1, rad_influence(i_rad))
     END IF
   END DO

   ! train output file for fields to come
   DO i = 1, nvars
          rcode = nf_inq_var(met_ncid, i, cval, itype, idm, ishape, natt)
          rcode = nf_def_var(oa_ncid, cval, itype, idm, ishape, i)
          DO na = 1, natt
             rcode = nf_inq_attname(met_ncid, i, na, cname)
             rcode = nf_copy_att   (met_ncid, i, cname, oa_ncid, i)
          ENDDO
   ENDDO
   rcode = nf_enddef(oa_ncid)

   ! copy all variables to the new file
   start_dims = 1
   DO i = 1, nvars
          rcode = nf_inq_var(met_ncid, i, cval, itype, idm, ishape, natt)
          dims_in  = 1
          dims_out = 1
          DO ii = 1,idm-1
            dims_in(ii)  = dval(ishape(ii))
            dims_out(ii) = dval2(ishape(ii))
          ENDDO

          IF     (itype == 2) THEN        
            allocate (text(dims_in(1), dims_in(2), dims_in(3)))
            rcode = nf_get_var_text (met_ncid, i, text)
            text = date_char
            rcode = nf_put_vara_text(oa_ncid, i, start_dims, dims_in, text)
            deallocate (text)

          ELSEIF (itype == 5 .AND. idm == 4) THEN 
            allocate (d3_data1(dims_in(1), dims_in(2), dims_in(3) ))
            allocate (d3_data2(dims_out(1),dims_out(2),dims_out(3)))
            rcode = nf_get_var_real(met_ncid, i, d3_data1)
            CALL unexpand3 ( d3_data1, dims_in(2), dims_in(1), dims_in(3), d3_data2, dims_out(2), dims_out(1) ) 
            rcode = nf_put_vara_real (oa_ncid, i, start_dims, dims_out, d3_data2)
            deallocate(d3_data1)
            deallocate(d3_data2)

          ELSEIF (itype == 5 .AND. idm == 3) THEN 
            allocate (d2_data1(dims_in(1), dims_in(2) ))
            allocate (d2_data2(dims_out(1),dims_out(2)))
            rcode = nf_get_var_real(met_ncid, i, d2_data1)
            CALL unexpand2 ( d2_data1, dims_in(2), dims_in(1), d2_data2, dims_out(2), dims_out(1) ) 
            rcode = nf_put_vara_real (oa_ncid, i, start_dims, dims_out, d2_data2)
                IF ( trim(cval) == 'XLAT_M' ) THEN
                   corner_lats( 1) = d2_data2(1          ,1          )
                   corner_lats( 2) = d2_data2(1          ,dims_out(2))
                   corner_lats( 3) = d2_data2(dims_out(1),dims_out(2))
                   corner_lats( 4) = d2_data2(dims_out(1),1          )
            ELSEIF ( trim(cval) == 'XLAT_U' ) THEN
                   corner_lats( 5) = d2_data2(1          ,1          )
                   corner_lats( 6) = d2_data2(1          ,dims_out(2))
                   corner_lats( 7) = d2_data2(dims_out(1),dims_out(2))
                   corner_lats( 8) = d2_data2(dims_out(1),1          )
                   corner_lats(13) = corner_lats( 5) - (d2_data2(1          ,2          )-d2_data2(1          ,1            ))/2.0
                   corner_lats(14) = corner_lats( 6) + (d2_data2(1          ,dims_out(2))-d2_data2(1          ,dims_out(2)-1))/2.0
                   corner_lats(15) = corner_lats( 7) + (d2_data2(dims_out(1),dims_out(2))-d2_data2(dims_out(1),dims_out(2)-1))/2.0
                   corner_lats(16) = corner_lats( 8) - (d2_data2(dims_out(1),2          )-d2_data2(dims_out(1),1            ))/2.0
            ELSEIF ( trim(cval) == 'XLAT_V' ) THEN
                   corner_lats( 9) = d2_data2(1          ,1          )
                   corner_lats(10) = d2_data2(1          ,dims_out(2))
                   corner_lats(11) = d2_data2(dims_out(1),dims_out(2))
                   corner_lats(12) = d2_data2(dims_out(1),1          )
            ELSEIF ( trim(cval) == 'XLONG_M' ) THEN
                   corner_lons( 1) = d2_data2(1          ,1          )
                   corner_lons( 2) = d2_data2(1          ,dims_out(2))
                   corner_lons( 3) = d2_data2(dims_out(1),dims_out(2))
                   corner_lons( 4) = d2_data2(dims_out(1),1          )
            ELSEIF ( trim(cval) == 'XLONG_U' ) THEN
                   corner_lons( 5) = d2_data2(1          ,1          )
                   corner_lons( 6) = d2_data2(1          ,dims_out(2))
                   corner_lons( 7) = d2_data2(dims_out(1),dims_out(2))
                   corner_lons( 8) = d2_data2(dims_out(1),1          )
            ELSEIF ( trim(cval) == 'XLONG_V' ) THEN
                   corner_lons( 9) = d2_data2(1          ,1          )
                   corner_lons(10) = d2_data2(1          ,dims_out(2))
                   corner_lons(11) = d2_data2(dims_out(1),dims_out(2))
                   corner_lons(12) = d2_data2(dims_out(1),1          )
                   corner_lons(13) = corner_lons( 9) - (d2_data2(2          ,1          )-d2_data2(1            ,1          ))/2.0 
                   corner_lons(14) = corner_lons(10) - (d2_data2(2          ,dims_out(2))-d2_data2(1            ,dims_out(2)))/2.0
                   corner_lons(15) = corner_lons(11) + (d2_data2(dims_out(1),dims_out(2))-d2_data2(dims_out(1)-1,dims_out(2)))/2.0
                   corner_lons(16) = corner_lons(12) + (d2_data2(dims_out(1),1          )-d2_data2(dims_out(1)-1,1          ))/2.0
            ENDIF
            deallocate(d2_data1)
            deallocate(d2_data2)
         ENDIF
   ENDDO

   !! Replace the corner meta data
   rcode = nf_put_att_real(oa_ncid, NF_GLOBAL, 'corner_lats', 5, 16, corner_lats)
   rcode = nf_put_att_real(oa_ncid, NF_GLOBAL, 'corner_lons', 5, 16, corner_lons)


   !! NOW replace the changed fields in the output file

   ALLOCATE ( met_em_2d (iewd-1, jnsd-1 ) )
   CALL unexpand2 ( slp_x, iew_alloc, jns_alloc, dum2d, iewd, jnsd ) 
   CALL yx2xy ( dum2d, jnsd, iewd, met_em_dum2d )
   DO jmet = 1,jnsd-1
   DO imet = 1,iewd-1
      met_em_2d(imet,jmet) = met_em_dum2d(imet,jmet)
   ENDDO
   ENDDO
   rcode = nf_inq_varid    ( oa_ncid, "PMSL", ivar_number )
   rcode = nf_put_var_real ( oa_ncid, ivar_number, met_em_2d )
   DEALLOCATE (met_em_2d)

   ALLOCATE ( met_em_3d (iewd-1, jnsd-1, kbu_alloc) )
   CALL unexpand3 ( t, iew_alloc, jns_alloc, kbu_alloc, dum3d, iewd, jnsd ) 
   DO kp = 1,kbu_alloc
      CALL yx2xy ( dum3d(1,1,kp), jnsd, iewd, met_em_dum3d(1,1,kp) )
      DO jmet = 1,jnsd-1
      DO imet = 1,iewd-1
         met_em_3d(imet,jmet,kp) = met_em_dum3d(imet,jmet,kp)
      ENDDO
      ENDDO
   ENDDO
   rcode = nf_inq_varid    ( oa_ncid, "TT", ivar_number )
   rcode = nf_put_var_real ( oa_ncid, ivar_number, met_em_3d )
   DEALLOCATE (met_em_3d)

   ALLOCATE ( met_em_3d (iewd-1, jnsd-1, kbu_alloc) )
   CALL unexpand3 ( rh, iew_alloc, jns_alloc, kbu_alloc, dum3d, iewd, jnsd ) 
   DO kp = 1,kbu_alloc
      CALL yx2xy ( dum3d(1,1,kp), jnsd, iewd, met_em_dum3d(1,1,kp) )
      DO jmet = 1,jnsd-1
      DO imet = 1,iewd-1
         met_em_3d(imet,jmet,kp) = met_em_dum3d(imet,jmet,kp)
      ENDDO
      ENDDO
   ENDDO
   rcode = nf_inq_varid    ( oa_ncid, "RH", ivar_number )
   rcode = nf_put_var_real ( oa_ncid, ivar_number, met_em_3d )
   DEALLOCATE (met_em_3d)

 
   !! Get the diff between the unchanged A grid and modified A grid
   !! This pertubation will be interpolated to the C grid before it is added to the 
   !! incoming C grid field 
   ALLOCATE ( met_em_3d (iewd, jnsd-1, kbu_alloc) )
   pertubation = u - uA
   CALL crs2dot_pertU ( pertubation, jns_alloc, iew_alloc, kbu_alloc )
   u = uC + pertubation
   CALL unexpand3 ( u, iew_alloc, jns_alloc, kbu_alloc, dum3d, iewd, jnsd ) 
   DO kp = 1,kbu_alloc
      CALL yx2xy ( dum3d(1,1,kp), jnsd, iewd, met_em_dum3d(1,1,kp) )
      DO jmet = 1,jnsd-1
      DO imet = 1,iewd
         met_em_3d(imet,jmet,kp) = met_em_dum3d(imet,jmet,kp)
      ENDDO
      ENDDO
   ENDDO
   rcode = nf_inq_varid    ( oa_ncid, "UU", ivar_number )
   rcode = nf_put_var_real ( oa_ncid, ivar_number, met_em_3d )
   DEALLOCATE (met_em_3d)

   !! Get the diff between the unchanged A grid and modified A grid
   !! This pertubation will be interpolated to the C grid before it is added to the 
   !! incoming C grid field 
   ALLOCATE ( met_em_3d (iewd-1, jnsd, kbu_alloc) )
   pertubation = v - vA
   CALL crs2dot_pertV ( pertubation, jns_alloc, iew_alloc, kbu_alloc )
   v = vC + pertubation
   CALL unexpand3 ( v, iew_alloc, jns_alloc, kbu_alloc, dum3d, iewd, jnsd ) 
   do kp = 1,kbu_alloc
      CALL yx2xy ( dum3d(1,1,kp), jnsd, iewd, met_em_dum3d(1,1,kp) )
      DO jmet = 1,jnsd
      DO imet = 1,iewd-1
         met_em_3d(imet,jmet,kp) = met_em_dum3d(imet,jmet,kp)
      ENDDO
      ENDDO
   enddo
   rcode = nf_inq_varid    ( oa_ncid, "VV", ivar_number )
   rcode = nf_put_var_real ( oa_ncid, ivar_number, met_em_3d )
   DEALLOCATE (met_em_3d)


   !! Recalculate the surface pressure
   
   ALLOCATE ( d3_data1 (jns_alloc, iew_alloc, kbu_alloc) )
   d3_data1 = pres
   !BPR BEGIN
   !The surface level of the pressure field can be adjusted based on changes
   !made to the sea level pressure by the objective analysis.
   !However, only do this if you did not create a surface pressure objective
   !analysis.
   IF(.NOT.oa_psfc) THEN
   !BPR END
    DO imet = 1,iew_alloc
    DO jmet = 1,jns_alloc
      d3_data1(jmet,imet,1) =    pres(jmet,imet,1) + &
                               ( pres(jmet,imet,1)/slp_C(jmet,imet) ) * &
                               ( slp_x(jmet,imet)-slp_C(jmet,imet) )   
    ENDDO
    ENDDO
   !BPR BEGIN
   ENDIF
   !BPR END

   ALLOCATE ( met_em_3d (iewd-1, jnsd-1, kbu_alloc) )
   ALLOCATE ( met_em_2d (iewd-1, jnsd-1 ) )
   CALL unexpand3 ( d3_data1, iew_alloc, jns_alloc, kbu_alloc, dum3d, iewd, jnsd ) 
   DEALLOCATE(d3_data1)
   CALL yx2xy ( dum3d(1,1,1), jnsd, iewd, met_em_dum3d(1,1,1) )
   DO jmet = 1,jnsd-1
   DO imet = 1,iewd-1
      met_em_3d(imet,jmet,1) = met_em_dum3d(imet,jmet,1)
      met_em_2d(imet,jmet  ) = met_em_dum3d(imet,jmet,1)
   ENDDO
   ENDDO
   DO kp = 2,kbu_alloc
      CALL yx2xy ( dum3d(1,1,kp), jnsd, iewd, met_em_dum3d(1,1,kp) )
      DO jmet = 1,jnsd-1
      DO imet = 1,iewd-1
         met_em_3d(imet,jmet,kp) = met_em_dum3d(imet,jmet,kp)
      ENDDO
      ENDDO
   ENDDO

   rcode = nf_inq_varid    ( oa_ncid, "PSFC", ivar_number )
   IF ( rcode == 0 ) rcode = nf_put_var_real ( oa_ncid, ivar_number, met_em_2d )
   DEALLOCATE (met_em_2d)

   rcode = nf_inq_varid    ( oa_ncid, "PRES", ivar_number )
   rcode = nf_put_var_real ( oa_ncid, ivar_number, met_em_3d )
   DEALLOCATE (met_em_3d)

   rcode = nf_close(met_ncid)
   rcode = nf_close(oa_ncid)


END SUBROUTINE write_analysis

!------------------------------------------------------------------------------

SUBROUTINE write_analysis_fdda ( oa_met , sfc_ncid , total_count, icount_fdda , &
t , u , v , uA , vA , uC , vC , h , rh , pres , &
terrain , slp_x , slp_C , snow , tobbox , odis , pressure , &
iew_alloc , jns_alloc , kbu_alloc , iewd , jnsd , &
fdda_date_24 , ptop , intf4d , lagtem , oa_type , rad_influence )

!  This routine assembles the correct data fields and outputs them with the
!  appropriate small header in the format expected for the surface FDDA file.

   USE input_data

   IMPLICIT NONE 
   INCLUDE 'netcdf.inc'

   INCLUDE 'big_header.inc'
 
   INTEGER                                                :: iewd , jnsd , iew_alloc , jns_alloc , kbu_alloc
   REAL , DIMENSION ( jns_alloc , iew_alloc , kbu_alloc ) :: t , u , v , h , rh , qv
   REAL , DIMENSION ( jns_alloc , iew_alloc , kbu_alloc ) :: uA , vA , uC , vC , pertubation, pres
   REAL , DIMENSION ( jns_alloc , iew_alloc )             :: terrain , slp_x , slp_C , tobbox , odis , psfc , snow
   REAL , DIMENSION ( jnsd , iewd )                       :: dum2d
   REAL , DIMENSION ( kbu_alloc )                         :: pressure

   CHARACTER (LEN=24)                                     :: fdda_date_24
   REAL                                                   :: ptop
   INTEGER                                                :: oa_met , sfc_ncid , icount_fdda , total_count

   INTEGER :: loop_count , &
              num3d , &
              num2d , &
              num1d , &
              num0d

   INTEGER                                           :: rcode, i, na, ii, imet, jmet, ivar_number
   INTEGER                                           :: itime, iwe, isn, iwe_stag, isn_stag
   INTEGER                                           :: ndims, nvars, ngatts, nunlimdimid
   CHARACTER (LEN=31),ALLOCATABLE, DIMENSION(:)      :: dname
   INTEGER,           ALLOCATABLE, DIMENSION(:)      :: dval, dval2
   CHARACTER (LEN=31)                                :: cname
   CHARACTER (LEN=50)                                :: cval
   INTEGER                                           :: ilen, itype, ival, idm, natt
   INTEGER, DIMENSION(6)                             :: ishape
   REAL                                              :: rval
   REAL, DIMENSION(16)                               :: rval_corners
   CHARACTER (LEN=19), ALLOCATABLE, DIMENSION(:,:,:) :: text
   REAL, ALLOCATABLE, DIMENSION(:,:,:)               :: d3_data1, d3_data2
   REAL, ALLOCATABLE, DIMENSION(:,:)                 :: d2_data1, d2_data2
   CHARACTER (LEN=20)                                :: units
   CHARACTER (LEN=45)                                :: description
   REAL, DIMENSION (iewd , jnsd)                     :: met_em_dum2d
   REAL, ALLOCATABLE , DIMENSION(:,:)                :: met_em_2d, th_sfc
   INTEGER, DIMENSION(3)                             :: start_dims, end_dims, start_dims_new
   INTEGER                                           :: is_there
   CHARACTER (LEN=1)                                 :: stagger
   REAL, PARAMETER                                   :: Rd = 287.
   REAL, PARAMETER                                   :: Cp = 7.*Rd/2.

   INTEGER                                           :: intf4d
   LOGICAL                                           :: lagtem
   CHARACTER *(*)                                    :: oa_type
   INTEGER , DIMENSION(10)                           :: rad_influence
   CHARACTER (LEN=16)                                :: rads
   INTEGER                                           :: i_rad

   INTEGER :: unit=12
#ifdef NCARG
real,dimension(iewd,jnsd)::ddd
integer :: i , j
#endif

!BPR BEGIN
integer :: varid
!BPR END

   !  We need to keep track of where we are sticking the data for the FDDA option.

   IF ( initial_time ) THEN
      tt = first_time
   ELSE
      tt = second_time
   END IF

   IF ( icount_fdda == 1 ) THEN
      !! First get all the dimension, global variables and unchanged fields from the input file
      ! find out some details about the output file
      rcode = nf_inq(oa_met, ndims, nvars, ngatts, nunlimdimid)
   
      ! sort out the dimensions
      IF (ALLOCATED(dname)) DEALLOCATE(dname)
          ALLOCATE (dname(ndims))
      IF (ALLOCATED(dval)) DEALLOCATE(dval)
          ALLOCATE  (dval(ndims))
      ii = 0
      DO i = 1, ndims
          rcode = nf_inq_dim(oa_met, i, dname(i), dval(i))
              IF (dname(i) == 'Time') THEN
                  ii = ii + 1
                  rcode = nf_def_dim(sfc_ncid, dname(i), NF_UNLIMITED, ii)
                  itime = ii
          ELSEIF (dname(i) == 'DateStrLen') THEN
                  ii = ii + 1
                  rcode = nf_def_dim(sfc_ncid, dname(i), dval(i), ii)
          ELSEIF (dname(i) == 'west_east') THEN
                  ii = ii + 1
                  rcode = nf_def_dim(sfc_ncid, dname(i), dval(i), ii)
                  iwe = ii
          ELSEIF (dname(i) == 'south_north') THEN
                  ii = ii + 1
                  rcode = nf_def_dim(sfc_ncid, dname(i), dval(i), ii)
                  isn = ii
          ELSEIF (dname(i) == 'west_east_stag') THEN
                  ii = ii + 1
                  rcode = nf_def_dim(sfc_ncid, dname(i), dval(i), ii)
                  iwe_stag = ii
          ELSEIF (dname(i) == 'south_north_stag') THEN
                  ii = ii + 1
                  rcode = nf_def_dim(sfc_ncid, dname(i), dval(i), ii)
                  isn_stag = ii
          ENDIF
   
      ENDDO 
   
      ! deal with the global variables
      loop_global_att : DO i = 1, ngatts
             rcode = nf_inq_attname(oa_met, NF_GLOBAL, i,    cname)
             rcode = nf_inq_atttype(oa_met, NF_GLOBAL, cname, itype)
             rcode = nf_inq_attlen (oa_met, NF_GLOBAL, cname, ilen)

             is_there = 0
             is_there = INDEX(cname,"PATCH") 
                IF ( is_there > 0 ) CYCLE loop_global_att
             is_there = INDEX(cname,"FLAG")  
                IF ( is_there > 0 ) CYCLE loop_global_att
   
             IF ( itype == 2 ) THEN        ! characters
               rcode = nf_get_att_text (oa_met, NF_GLOBAL, cname, cval)
               IF(cname(1:5) .eq. 'TITLE') then
                  cval = "SFCFDDA FILE: SURFACE ANALYSES FOR NUDGING"
                  ilen = len_trim(cval)
               ENDIF 
               rcode = nf_put_att_text(sfc_ncid, NF_GLOBAL, cname, ilen, cval(1:ilen))
   
             ELSEIF ( itype == 4 ) THEN     ! integers
               IF ( cname == 'BOTTOM-TOP_GRID_DIMENSION' ) THEN
                  rcode = nf_put_att_int(sfc_ncid, NF_GLOBAL, cname, itype, ilen, 1)
               ELSE
                  rcode = nf_get_att_int (oa_met, NF_GLOBAL, cname, ival)
                  rcode = nf_put_att_int(sfc_ncid, NF_GLOBAL, cname, itype, ilen, ival)
               ENDIF
   
             ELSEIF ( itype == 5 ) THEN    ! real
               IF ( cname == 'corner_lats' .OR. cname == 'corner_lons' ) THEN
                  rcode = nf_get_att_real (oa_met, NF_GLOBAL, cname, rval_corners)
                  rcode = nf_put_att_real(sfc_ncid, NF_GLOBAL, cname, itype, ilen, rval_corners)
               ELSE
                  rcode = nf_get_att_real (oa_met, NF_GLOBAL, cname, rval)
                  rcode = nf_put_att_real(sfc_ncid, NF_GLOBAL, cname, itype, ilen, rval)
               ENDIF 
             ENDIF 
      ENDDO loop_global_att    
      ilen = len_trim(oa_type)
      rcode = nf_put_att_text(sfc_ncid, NF_GLOBAL, "OA_SCHEME", ilen, oa_type(1:ilen))
      DO i_rad = 1,10
        IF ( rad_influence(i_rad) > 0 ) THEN
          WRITE(rads,'("RAD_INFLUENCE_",i2.2)')i_rad
          rcode = nf_put_att_int(sfc_ncid, NF_GLOBAL, rads, 4, 1, rad_influence(i_rad))
        END IF
      END DO
      IF ( lagtem ) THEN
        rcode = nf_put_att_int(sfc_ncid, NF_GLOBAL, "FLAG_LAGTEM", 4, 1, 1)
      ELSE
        rcode = nf_put_att_int(sfc_ncid, NF_GLOBAL, "FLAG_LAGTEM", 4, 1, 0)
      END IF

      ! train output file for fields to come
      DO i = 1, 23
             IF ( i == 1 ) THEN
                rcode = nf_inq_var(oa_met, i, cval, itype, idm, ishape, natt)
                !BPR BEGIN
                if(rcode.ne.0) then
                 PRINT *, 'ERROR: Attempt to read variable from metoa file failed: ',NF_STRERROR(rcode)
                end if
                !The last argument of nf_def_var is an OUTPUT of nf_def_var
                !If this variable is our loop variable than this call can change the
                !value of our loop value
                !rcode = nf_def_var(sfc_ncid, cval, itype, idm, ishape, i)
                rcode = nf_def_var(sfc_ncid, cval, itype, idm, ishape, varid)
                if(rcode.ne.0) then
                 PRINT *, 'ERROR: Attempt to define a variable in the surface analysis nudging file failed: ',NF_STRERROR(rcode)
                end if
                !BPR END
             ELSE
                ishape(1) = iwe
                ishape(2) = isn
                ishape(3) = itime
                stagger   = ' ' 
                IF ( i == 2 ) THEN
                     cval = 'T2_NDG_OLD'
                     units = 'K'
                     description = 'Temperature'
                ELSEIF ( i == 3 ) THEN
                     cval = 'T2_NDG_NEW'
                     units = 'K'
                     description = 'Temperature'
                ELSEIF ( i == 4 ) THEN
                     cval = 'TH2_NDG_OLD'
                     units = 'K'
                     description = 'Potential Temperature'
                ELSEIF ( i == 5 ) THEN
                     cval = 'TH2_NDG_NEW'
                     units = 'K'
                     description = 'Potential Temperature'
                ELSEIF ( i == 6 ) THEN
                     ishape(1) = iwe_stag
                     cval = 'U10_NDG_OLD'
                     units = 'm/s'
                     description = 'u-wind component'
                     stagger   = 'X' 
                ELSEIF ( i == 7 ) THEN
                     ishape(1) = iwe_stag
                     cval = 'U10_NDG_NEW'
                     units = 'm/s'
                     description = 'u-wind component'
                     stagger   = 'X' 
                ELSEIF ( i == 8 ) THEN
                     ishape(2) = isn_stag
                     cval = 'V10_NDG_OLD'
                     units = 'm/s'
                     description = 'v-wind component'
                     stagger   = 'Y' 
                ELSEIF ( i == 9 ) THEN
                     ishape(2) = isn_stag
                     cval = 'V10_NDG_NEW'
                     units = 'm/s'
                     description = 'v-wind component'
                     stagger   = 'Y' 
                ELSEIF ( i == 10 ) THEN
                     cval = 'RH_NDG_OLD'
                     units = '%'
                     description = 'SURFACE RELATIVE HUMIDITY'
                ELSEIF ( i == 11 ) THEN
                     cval = 'RH_NDG_NEW'
                     units = '%'
                     description = 'SURFACE RELATIVE HUMIDITY'
                ELSEIF ( i == 12 ) THEN
                     cval = 'Q2_NDG_OLD'
                     units = 'kg/kg'
                     description = 'SURFACE MIXING RATIO'
                ELSEIF ( i == 13 ) THEN
                     cval = 'Q2_NDG_NEW'
                     units = 'kg/kg'
                     description = 'SURFACE MIXING RATIO'
                ELSEIF ( i == 14 ) THEN
                     cval = 'PS_NDG_OLD'
                     units = 'Pa'
                     description = 'SURFACE PRESSURE'
                ELSEIF ( i == 15 ) THEN
                     cval = 'PS_NDG_NEW'
                     units = 'Pa'
                     description = 'SURFACE PRESSURE'
                ELSEIF ( i == 16 ) THEN
                     cval = 'PSL_NDG_OLD'
                     units = 'Pa'
                     description = 'SEA-LEVEL PRESSURE'
                ELSEIF ( i == 17 ) THEN
                     cval = 'PSL_NDG_NEW'
                     units = 'Pa'
                     description = 'SEA-LEVEL PRESSURE'
                ELSEIF ( i == 18 ) THEN
                     cval = 'TOB_NDG_OLD'
                     units = 'OBSERVATIONS'      
                     description = 'OBSERVATION DENSITY'
                ELSEIF ( i == 19 ) THEN
                     cval = 'TOB_NDG_NEW'
                     units = 'OBSERVATIONS'      
                     description = 'OBSERVATION DENSITY'
                ELSEIF ( i == 20 ) THEN
                     cval = 'ODIS_NDG_OLD'
                     units = 'km'      
                     description = 'DISTANCE TO NEAREST OBSERVATION'
                ELSEIF ( i == 21 ) THEN
                     cval = 'ODIS_NDG_NEW'
                     units = 'km'      
                     description = 'DISTANCE TO NEAREST OBSERVATION'
                ELSEIF ( i == 22 ) THEN
                     cval = 'SN_NDG_OLD'
                     units = 'kg m-2'      
                     description = 'WATER EQUIVALENT SNOW DEPTH'    
                ELSEIF ( i == 23 ) THEN
                     cval = 'SN_NDG_NEW'
                     units = 'kg m-2'      
                     description = 'WATER EQUIVALENT SNOW DEPTH'    
                ENDIF
                !BPR BEGIN
                !The final argument to nf_dev_var is an OUTPUT and gives the variable id of the
                !variable you are inquiring about
                !If this is your loop variable you can redefine your loop variable and get stuck
                !in an infinite loop!
                !Therefore, instead of using the loop variable "i", we now use "varid"
                !rcode = nf_def_var(sfc_ncid, cval, NF_FLOAT, 3, ishape, i)
                !   rcode = nf_put_att_int(sfc_ncid, i, 'FieldType', NF_INT, 1, 104 )
                !   rcode = nf_put_att_text(sfc_ncid, i, 'MemoryOrder', 3, 'XY ' )
                !      ilen = len_trim(units)
                !   rcode = nf_put_att_text(sfc_ncid, i, 'units', ilen, trim(units) )
                !      ilen = len_trim(description)
                !   rcode = nf_put_att_text(sfc_ncid, i, 'description', ilen, trim(description) )
                !   rcode = nf_put_att_text(sfc_ncid, i, 'stagger', 1, stagger )
                rcode = nf_def_var(sfc_ncid, cval, NF_FLOAT, 3, ishape, varid)
                   rcode = nf_put_att_int(sfc_ncid, varid, 'FieldType', NF_INT, 1, 104 )
                   rcode = nf_put_att_text(sfc_ncid, varid, 'MemoryOrder', 3, 'XY ' )
                      ilen = len_trim(units)
                   rcode = nf_put_att_text(sfc_ncid, varid, 'units', ilen, trim(units) )
                      ilen = len_trim(description)
                   rcode = nf_put_att_text(sfc_ncid, varid, 'description', ilen, trim(description) )
                   rcode = nf_put_att_text(sfc_ncid, varid, 'stagger', 1, stagger )
                !BPR END
             ENDIF
      ENDDO
      rcode = nf_enddef(sfc_ncid)
   ENDIF

   ! write out the new variables to our file

   IF ( icount_fdda < total_count ) THEN
      start_dims    = 1
      end_dims      = 1
      start_dims(2) = icount_fdda
      end_dims(1)   = 19
      rcode = nf_inq_varid     ( sfc_ncid, "Times", ivar_number )
      rcode = nf_put_vara_text ( sfc_ncid, ivar_number, start_dims, end_dims, fdda_date_24 )
   ENDIF

   start_dims        = 1
   start_dims_new    = 1
   end_dims          = 1
   start_dims(3)     = icount_fdda
   start_dims_new(3) = icount_fdda-1

   !! WRITING INFO FOR U10_NDG
   end_dims(1) = iewd
   end_dims(2) = jnsd-1
   pertubation = u - uA
   CALL crs2dot_pertU ( pertubation , jns_alloc, iew_alloc , kbu_alloc )
   u = uC + pertubation
   CALL unexpand3 ( u  , iew_alloc , jns_alloc , 1 , dum2d , iewd , jnsd ) 
   CALL yx2xy ( dum2d, jnsd, iewd, met_em_dum2d )
   ALLOCATE ( met_em_2d (iewd, jnsd-1) )
   DO jmet = 1,jnsd-1
   DO imet = 1,iewd
      met_em_2d(imet,jmet) = met_em_dum2d(imet,jmet)
   ENDDO
   ENDDO
   IF ( icount_fdda < total_count ) THEN
     rcode = nf_inq_varid     ( sfc_ncid, "U10_NDG_OLD", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims, end_dims, met_em_2d )
   ENDIF
   IF ( icount_fdda > 1 ) THEN
     rcode = nf_inq_varid     ( sfc_ncid, "U10_NDG_NEW", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims_new, end_dims, met_em_2d )
   ENDIF
   DEALLOCATE (met_em_2d)

   !! WRITING INFO FOR V10_NDG
   end_dims(1) = iewd-1
   end_dims(2) = jnsd
   pertubation = v - vA
   CALL crs2dot_pertV ( pertubation , jns_alloc, iew_alloc , kbu_alloc )
   v = vC + pertubation
   CALL unexpand3 ( v  , iew_alloc , jns_alloc , 1 , dum2d , iewd , jnsd ) 
   CALL yx2xy ( dum2d, jnsd, iewd, met_em_dum2d )
   ALLOCATE ( met_em_2d (iewd-1, jnsd) )
   DO jmet = 1,jnsd
   DO imet = 1,iewd-1
      met_em_2d(imet,jmet) = met_em_dum2d(imet,jmet)
   ENDDO
   ENDDO
   IF ( icount_fdda < total_count ) THEN
     rcode = nf_inq_varid     ( sfc_ncid, "V10_NDG_OLD", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims, end_dims, met_em_2d )
   ENDIF
   IF ( icount_fdda > 1 ) THEN
     rcode = nf_inq_varid     ( sfc_ncid, "V10_NDG_NEW", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims_new, end_dims, met_em_2d )
   ENDIF
   DEALLOCATE (met_em_2d)
 
   !! WRITING INFO FOR T2_NDG
   end_dims(1) = iewd-1
   end_dims(2) = jnsd-1
   CALL unexpand3 ( t  , iew_alloc , jns_alloc , 1 , dum2d , iewd , jnsd ) 
   CALL yx2xy ( dum2d, jnsd, iewd, met_em_dum2d )
   ALLOCATE ( met_em_2d (iewd-1, jnsd-1) )
   ALLOCATE ( th_sfc (iewd-1, jnsd-1) )
   DO jmet = 1,jnsd-1
   DO imet = 1,iewd-1
      met_em_2d(imet,jmet) = met_em_dum2d(imet,jmet)
   ENDDO
   ENDDO
   IF ( icount_fdda < total_count ) THEN
     rcode = nf_inq_varid     ( sfc_ncid, "T2_NDG_OLD", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims, end_dims, met_em_2d )
   ENDIF
   IF ( icount_fdda > 1 ) THEN
     rcode = nf_inq_varid     ( sfc_ncid, "T2_NDG_NEW", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims_new, end_dims, met_em_2d )
   ENDIF
   th_sfc = met_em_2d
   DEALLOCATE (met_em_2d)

   !! WRITING INFO FOR RH_NDG
   CALL unexpand3 ( rh , iew_alloc , jns_alloc , 1 , dum2d , iewd , jnsd ) 
   CALL yx2xy ( dum2d, jnsd, iewd, met_em_dum2d )
   ALLOCATE ( met_em_2d (iewd-1, jnsd-1) )
   DO jmet = 1,jnsd-1
   DO imet = 1,iewd-1
      met_em_2d(imet,jmet) = met_em_dum2d(imet,jmet)
   ENDDO
   ENDDO
   IF ( icount_fdda < total_count ) THEN
     rcode = nf_inq_varid     ( sfc_ncid, "RH_NDG_OLD", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims, end_dims, met_em_2d )
   ENDIF
   IF ( icount_fdda > 1 ) THEN
     rcode = nf_inq_varid     ( sfc_ncid, "RH_NDG_NEW", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims_new, end_dims, met_em_2d )
   ENDIF
   DEALLOCATE (met_em_2d)

   !! WRITING INFO FOR Q2_NDG
   CALL mixing_ratio ( rh , t , h , terrain , slp_x , pressure , iew_alloc , jns_alloc , kbu_alloc , qv , psfc )
   CALL unexpand3 ( qv , iew_alloc , jns_alloc , 1 , dum2d , iewd , jnsd )
   CALL yx2xy ( dum2d, jnsd, iewd, met_em_dum2d )
   ALLOCATE ( met_em_2d (iewd-1, jnsd-1) )
   DO jmet = 1,jnsd-1
   DO imet = 1,iewd-1
      met_em_2d(imet,jmet) = met_em_dum2d(imet,jmet)
   ENDDO
   ENDDO
   IF ( icount_fdda < total_count ) THEN
     rcode = nf_inq_varid     ( sfc_ncid, "Q2_NDG_OLD", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims, end_dims, met_em_2d )
   ENDIF
   IF ( icount_fdda > 1 ) THEN
     rcode = nf_inq_varid     ( sfc_ncid, "Q2_NDG_NEW", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims_new, end_dims, met_em_2d )
   ENDIF
   DEALLOCATE (met_em_2d)

   !! WRITING INFO FOR PS_NDG & TH2_NDG
   !! We are not using the psfc calculated by the model, but a smmoth one we recalculate

   ALLOCATE ( d2_data1 (jns_alloc, iew_alloc) )
   DO imet = 1,iew_alloc
   DO jmet = 1,jns_alloc
      d2_data1(jmet,imet) =    pres(jmet,imet,1) + &
                             ( pres(jmet,imet,1)/slp_C(jmet,imet) ) * &
                             ( slp_x(jmet,imet)-slp_C(jmet,imet) )
   ENDDO
   ENDDO
   CALL unexpand2 ( d2_data1 , iew_alloc , jns_alloc , dum2d , iewd , jnsd )
   DEALLOCATE(d2_data1)
   CALL yx2xy ( dum2d, jnsd, iewd, met_em_dum2d )
   ALLOCATE ( met_em_2d (iewd-1, jnsd-1) )
   DO jmet = 1,jnsd-1
   DO imet = 1,iewd-1
      met_em_2d(imet,jmet) = met_em_dum2d(imet,jmet)
      th_sfc(imet,jmet) = th_sfc(imet,jmet) * ( 100000./met_em_2d(imet,jmet) ) ** (Rd / Cp)
   ENDDO
   ENDDO
   IF ( icount_fdda < total_count ) THEN
     rcode = nf_inq_varid     ( sfc_ncid, "PS_NDG_OLD", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims, end_dims, met_em_2d )
     rcode = nf_inq_varid     ( sfc_ncid, "TH2_NDG_OLD", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims, end_dims, th_sfc )
   ENDIF
   IF ( icount_fdda > 1 ) THEN
     rcode = nf_inq_varid     ( sfc_ncid, "PS_NDG_NEW", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims_new, end_dims, met_em_2d )
     rcode = nf_inq_varid     ( sfc_ncid, "TH2_NDG_NEW", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims_new, end_dims, th_sfc )
   ENDIF
   DEALLOCATE (met_em_2d)
   DEALLOCATE (th_sfc)

   !! WRITING INFO FOR PSL_NDG
   CALL unexpand2 ( slp_x , iew_alloc , jns_alloc , dum2d , iewd , jnsd ) 
   dum2d(:,iewd)=dum2d(:,iewd-1)
   dum2d(jnsd,:)=dum2d(jnsd-1,:)
   CALL yx2xy ( dum2d, jnsd, iewd, met_em_dum2d )
   ALLOCATE ( met_em_2d (iewd-1, jnsd-1) )
   DO jmet = 1,jnsd-1
   DO imet = 1,iewd-1
      met_em_2d(imet,jmet) = met_em_dum2d(imet,jmet)
   ENDDO
   ENDDO
   IF ( icount_fdda < total_count ) THEN
     rcode = nf_inq_varid     ( sfc_ncid, "PSL_NDG_OLD", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims, end_dims, met_em_2d )
   ENDIF
   IF ( icount_fdda > 1 ) THEN
     rcode = nf_inq_varid     ( sfc_ncid, "PSL_NDG_NEW", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims_new, end_dims, met_em_2d )
   ENDIF
   DEALLOCATE (met_em_2d)

   !! WRITING INFO FOR TOB_NDG
   CALL unexpand2 ( tobbox , iew_alloc , jns_alloc , dum2d , iewd , jnsd ) 
   CALL yx2xy ( dum2d, jnsd, iewd, met_em_dum2d )
   ALLOCATE ( met_em_2d (iewd-1, jnsd-1) )
   DO jmet = 1,jnsd-1
   DO imet = 1,iewd-1
      met_em_2d(imet,jmet) = met_em_dum2d(imet,jmet)
   ENDDO
   ENDDO
   IF ( icount_fdda < total_count ) THEN
     rcode = nf_inq_varid     ( sfc_ncid, "TOB_NDG_OLD", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims, end_dims, met_em_2d )
   ENDIF
   IF ( icount_fdda > 1 ) THEN
     rcode = nf_inq_varid     ( sfc_ncid, "TOB_NDG_NEW", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims_new, end_dims, met_em_2d )
   ENDIF
   DEALLOCATE (met_em_2d)

   !! WRITING INFO FOR ODIS_NDG
   CALL unexpand2 ( odis , iew_alloc , jns_alloc , dum2d , iewd , jnsd ) 
   CALL yx2xy ( dum2d, jnsd, iewd, met_em_dum2d )
   ALLOCATE ( met_em_2d (iewd-1, jnsd-1) )
   DO jmet = 1,jnsd-1
   DO imet = 1,iewd-1
      met_em_2d(imet,jmet) = met_em_dum2d(imet,jmet)
   ENDDO
   ENDDO
   IF ( icount_fdda < total_count ) THEN
     rcode = nf_inq_varid     ( sfc_ncid, "ODIS_NDG_OLD", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims, end_dims, met_em_2d )
   ENDIF
   IF ( icount_fdda > 1 ) THEN
     rcode = nf_inq_varid     ( sfc_ncid, "ODIS_NDG_NEW", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims_new, end_dims, met_em_2d )
   ENDIF
   DEALLOCATE (met_em_2d)

   !! WRITING INFO FOR SN_NDG
   CALL unexpand2 ( snow , iew_alloc , jns_alloc , dum2d , iewd , jnsd ) 
   CALL yx2xy ( dum2d, jnsd, iewd, met_em_dum2d )
   ALLOCATE ( met_em_2d (iewd-1, jnsd-1) )
   DO jmet = 1,jnsd-1
   DO imet = 1,iewd-1
      met_em_2d(imet,jmet) = met_em_dum2d(imet,jmet)
   ENDDO
   ENDDO
   IF ( icount_fdda < total_count ) THEN
     rcode = nf_inq_varid     ( sfc_ncid, "SN_NDG_OLD", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims, end_dims, met_em_2d )
   ENDIF
   IF ( icount_fdda > 1 ) THEN
     rcode = nf_inq_varid     ( sfc_ncid, "SN_NDG_NEW", ivar_number )
     rcode = nf_put_vara_real ( sfc_ncid, ivar_number, start_dims_new, end_dims, met_em_2d )
   ENDIF
   DEALLOCATE (met_em_2d)


   !  When we are done with all the time periods - close the SFCFDDA file
   IF ( total_count == icount_fdda ) rcode = nf_close(sfc_ncid)


END SUBROUTINE write_analysis_fdda

!-------------------------------------------------------------------------------

SUBROUTINE unexpand3 ( orig , iew_big , jns_big , kx , new , iew_small , jns_small ) 

   INTEGER                                         :: iew_big   , jns_big   , &
                                                      kx , &
                                                      iew_small , jns_small
   REAL , DIMENSION ( jns_big   , iew_big   , kx ) :: orig
   REAL , DIMENSION ( jns_small , iew_small , kx ) :: new

   INTEGER                                         :: i , j , k

   DO k = 1 , kx
      DO i = 1 , iew_small
         DO j = 1 , jns_small
            new(j,i,k) = orig(j+(jns_big-jns_small)/2,i+(iew_big-iew_small)/2,k)
         END DO
      END DO
   END DO

END SUBROUTINE unexpand3


!-------------------------------------------------------------------------------

SUBROUTINE unexpand2 ( orig , iew_big , jns_big , new , iew_small , jns_small ) 

   INTEGER                                         :: iew_big   , jns_big   , &
                                                      kx , &
                                                      iew_small , jns_small
   REAL , DIMENSION ( jns_big   , iew_big   )      :: orig
   REAL , DIMENSION ( jns_small , iew_small )      :: new

   INTEGER                                         :: i , j , k

   DO i = 1 , iew_small
      DO j = 1 , jns_small
         new(j,i) = orig(j+(jns_big-jns_small)/2,i+(iew_big-iew_small)/2)
      END DO
   END DO

END SUBROUTINE unexpand2

!-------------------------------------------------------------------------------

SUBROUTINE mixing_ratio     ( rh , t , h , &
                              terrain , slp_x , &
                              pressure , &
                              iew_alloc , jns_alloc , kbu_alloc , &
                              qv , psfc )
                                      
   IMPLICIT NONE 

   INTEGER , INTENT(IN) :: iew_alloc , jns_alloc , kbu_alloc
  
   REAL , INTENT(IN) , DIMENSION(kbu_alloc) :: pressure
   REAL , INTENT(IN) , DIMENSION(jns_alloc,iew_alloc) :: terrain , slp_x
   REAL , INTENT(INOUT) , DIMENSION(jns_alloc,iew_alloc,kbu_alloc) :: rh
   REAL , INTENT(IN) , DIMENSION(jns_alloc,iew_alloc,kbu_alloc) :: t , h

   REAL , INTENT(OUT) , DIMENSION(jns_alloc,iew_alloc,kbu_alloc) :: qv
   REAL , INTENT(OUT) , DIMENSION(jns_alloc,iew_alloc) :: psfc

   INTEGER :: i , j , k 
   REAL                        :: ES
   REAL                        :: QS
   REAL,         PARAMETER     :: EPS         = 0.622
   REAL,         PARAMETER     :: SVP1        = 0.6112
   REAL,         PARAMETER     :: SVP2        = 17.67
   REAL,         PARAMETER     :: SVP3        = 29.65
   REAL,         PARAMETER     :: SVPT0       = 273.15

   !  Compute Qv at all of the mandatory levels - not the surface.

   rh(1:jns_alloc-1,1:iew_alloc-1,:) = MIN ( MAX ( rh(1:jns_alloc-1,1:iew_alloc-1,:) ,  1. ) , 100. )

   DO k = 2 , kbu_alloc
      DO j = 1, jns_alloc-1
         DO i = 1, iew_alloc-1
            es = svp1 * 10. * EXP(svp2 * (t(j,i,k) - svpt0) / (t(j,i,k) - svp3))
            qs = eps * es / (pressure(k) - es)
            qv(j,i,k) = MAX(0.01 * rh(j,i,k) * qs,0.0)
         END DO
      END DO
   END DO

   !  Now, with a 3d Qv, we can do the surface pressure.

   CALL sfcprs ( t , qv , h , slp_x , terrain , pressure , &
   iew_alloc , jns_alloc , kbu_alloc , psfc )

   !  Now, with a surface pressure, we can do the Qv at the surface, which was what
   !  we wanted any who.

   DO j = 1, jns_alloc-1
      DO i = 1, iew_alloc-1
         es = svp1 * 10. * EXP(svp2 * (t(j,i,1) - svpt0) / (t(j,i,1) - svp3))
         qs = eps * es / (psfc(j,i)/100. - es)
         qv(j,i,1) = MAX(0.01 * rh(j,i,1) * qs,0.0)
      END DO
   END DO

END SUBROUTINE mixing_ratio

!-------------------------------------------------------------------------------

SUBROUTINE sfcprs ( t , q , h , slp_x , terrain , pressure , &
iew_alloc , jns_alloc , kbu_alloc , psfc )

   IMPLICIT NONE

   INTEGER , INTENT(IN) ::  iew_alloc , jns_alloc , kbu_alloc
   REAL , INTENT(IN) , DIMENSION(jns_alloc,iew_alloc,kbu_alloc) ::  t , q , h
   REAL , INTENT(IN) , DIMENSION(jns_alloc,iew_alloc) ::  slp_x , terrain
   REAL , INTENT(IN) , DIMENSION(kbu_alloc) :: pressure

   REAL , INTENT(OUT) , DIMENSION(jns_alloc,iew_alloc) ::  psfc

   !  Local variables, counters, etc.

   INTEGER :: k850=0 , k700=0 , k500=0
   INTEGER :: i , j , k    
   REAL , PARAMETER :: R=287. , g = 9.806 , gamma = 6.5e-3
   REAL , PARAMETER :: gammarg = gamma*R/g 
   REAL :: p78 , p57 , gamma78 , gamma57
   REAL , DIMENSION(jns_alloc,iew_alloc) ::  ht , p1 , tslv , tsfc
   REAL :: t850 , t700 , t500 , t1

   !  Find the levels with which we wish to play.

   DO k = 2 , kbu_alloc
      IF      ( NINT(pressure(k) ) .EQ. 850 ) THEN
         k850 = k
      ELSE IF ( NINT(pressure(k) ) .EQ. 700 ) THEN
         k700 = k
      ELSE IF ( NINT(pressure(k) ) .EQ. 500 ) THEN
         k500 = k
      END IF
   END DO

   IF ( ( k850 .EQ. 0 ) .OR. & 
        ( k700 .EQ. 0 ) .OR. & 
        ( k500 .EQ. 0 ) ) THEN
      PRINT '(A)','Problems with find the 850, 700, or 500 mb level.'
      PRINT '(3(A,I3),A)','k850=',k850,', k700=',k700,', k500=',k500,'.'
      STOP 'No_find_levels'
   END IF

   !  First guess at surface pressure to choose the correct lapse rates.

   DO j = 1 , jns_alloc-1
      DO i = 1 , iew_alloc-1
         ht(j,i)= - terrain(j,i) / h(j,i,k850)
         psfc(j,i) = slp_x(j,i) * ( slp_x(j,i)/85000. ) ** ht(j,i)
      END DO
   END DO

   !  Pressure about 100 hPa above the surface.

   DO j = 1 , jns_alloc-1
      DO i = 1 , iew_alloc-1
         IF ( psfc(j,i) .GT. 95000 ) THEN
            p1(j,i) = 85000.
         ELSE IF ( psfc(j,i) .GT. 70000 ) THEN
            p1(j,i) = psfc(j,i) - 10000.
         ELSE 
            p1(j,i) = 50000.
         END IF
      END DO
   END DO

   !  Compute the necessary sea level temp and the surface temp.

   p78 = 1./LOG(850./700.) 
   p57 = 1./LOG(700./500.) 

   DO j = 1 , jns_alloc-1
      DO i = 1 , iew_alloc-1

         !  Virtual temperatures.

         t850 = t(j,i,k850) * ( 1. + 0.608 * q(j,i,k850) ) 
         t700 = t(j,i,k700) * ( 1. + 0.608 * q(j,i,k700) ) 
         t500 = t(j,i,k500) * ( 1. + 0.608 * q(j,i,k500) ) 
  
         !  Lapse rates.

         gamma78 = LOG(t850/t700) * p78
         gamma57 = LOG(t700/t500) * p57

         !  Approximate temperature, about 100 hPa above the surface.

         IF ( psfc(j,i) .GT. 95000 ) THEN
            t1 = t850
         ELSE IF ( psfc(j,i) .GT. 85000 ) THEN
            t1 = t700 * (p1(j,i)/70000.)**gamma78
         ELSE IF ( psfc(j,i) .GT. 70000 ) THEN
            t1 = t500 * (p1(j,i)/50000.)**gamma57
         ELSE 
            t1 = t500
         END IF
     
         !  Extrapolate to get the sea level temperature.

         tslv(j,i) = t1 * ( slp_x(j,i)/p1(j,i) ) ** gammarg

         !  Compute the surface temperature from the sea level temp.

         tsfc(j,i) = tslv(j,i) - gamma * terrain(j,i)
      END DO
   END DO
   
   !  Ta Da - we now do the sea level pressure computation.

   DO j = 1 , jns_alloc-1
      DO i = 1 , iew_alloc-1
         psfc(j,i) = slp_x(j,i) * EXP ( (-terrain(j,i)*g) / (R/2. * (tsfc(j,i) + tslv(j,i))) )
      END DO
   END DO

END SUBROUTINE sfcprs

!-------------------------------------------------------------------------------

END MODULE final_analysis
