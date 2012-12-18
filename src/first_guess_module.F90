MODULE first_guess

   USE small_header_data
   USE input_data
   USE date_pack

   REAL , DIMENSION(:,:,:) , ALLOCATABLE              :: extra_2d_data
   CHARACTER(LEN=80) , DIMENSION(:) , ALLOCATABLE     :: extra_2d_ch
   TYPE(sh) , DIMENSION(:) , ALLOCATABLE              :: extra_2d_sh
   INTEGER                                            :: extra_2d_count

CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE alloc_lots_of_2d_data ( jns_alloc , iew_alloc , num2d_extra )

   IMPLICIT NONE

   INTEGER :: iew_alloc , jns_alloc , num2d_extra

   extra_2d_count = num2d_extra

   IF ( extra_2d_count .GT. 0 ) THEN
      ALLOCATE (extra_2d_data(jns_alloc,iew_alloc,num2d_extra))
      ALLOCATE (extra_2d_ch(num2d_extra))
      ALLOCATE (extra_2d_sh(num2d_extra))
   END IF

END SUBROUTINE alloc_lots_of_2d_data

!------------------------------------------------------------------------------

SUBROUTINE crs2dot(field,dim1,dim2)
  
   IMPLICIT NONE

   INTEGER :: dim1 , dim2
   REAL , DIMENSION(dim1,dim2) :: field,dummy
   INTEGER :: i , j

   dummy(2:dim1-1,2:dim2-1)           = ( field(1:dim1-2,1:dim2-2) + &
                                          field(1:dim1-2,2:dim2-1) + &
                                          field(2:dim1-1,1:dim2-2) + &
                                          field(2:dim1-1,2:dim2-1) ) * 0.25
  
   dummy(2:dim1-1,1:dim2:dim2-1)      = ( field(1:dim1-2,1:dim2-1:dim2-2) + &
                                          field(2:dim1-1,1:dim2-1:dim2-2) ) * 0.5
  
   dummy(1:dim1:dim1-1,2:dim2-1)      = ( field(1:dim1-1:dim1-2,1:dim2-2) + &
                                          field(1:dim1-1:dim1-2,2:dim2-1) ) * 0.5
  
   dummy(1:dim1:dim1-1,1:dim2:dim2-1) =   field(1:dim1-1:dim1-2,1:dim2-1:dim2-2)
  
   field                              =   dummy
  
END SUBROUTINE crs2dot

SUBROUTINE crs2dot_pertU(field,dim1,dim2,dim3)
  
   IMPLICIT NONE

   INTEGER :: dim1 , dim2 , dim3
   REAL , DIMENSION(dim1,dim2,dim3) :: field,dummy

   dummy = 0.0
   dummy(:,2:dim2-1,:)         = ( field(:,1:dim2-2,:) + &
                                   field(:,2:dim2-1,:) ) * 0.5
   dummy(:,1,:)                = field(:,1,:)
   dummy(:,dim2,:)             = field(:,dim2-1,:)
  
   field                       =   dummy
  
END SUBROUTINE crs2dot_pertU

SUBROUTINE crs2dot_pertV(field,dim1,dim2,dim3)
  
   IMPLICIT NONE

   INTEGER :: dim1 , dim2 , dim3
   REAL , DIMENSION(dim1,dim2,dim3) :: field,dummy

   dummy = 0.0
   dummy(2:dim1-1,:,:)         = ( field(1:dim1-2,:,:) + &
                                   field(2:dim1-1,:,:) ) * 0.5
   dummy(1,:,:)                = field(1,:,:)
   dummy(dim1,:,:)             = field(dim1-1,:,:)
  
   field                       =   dummy
  
END SUBROUTINE crs2dot_pertV

!------------------------------------------------------------------------------

SUBROUTINE read_first_guess ( met_ncid , &
bhi , bhr , num3d , num2d , num1d , &
t , u , v , uA , vA , uC , vC , h , rh , pres , terrain , &
latitude_x , longitude_x , latitude_d , longitude_d , &
slp_x , slp_C , sst , snow , pressure , &
iew_alloc , jns_alloc , kbu_alloc , &
current_date_8 , current_time_6 , date_char , icount , print_analysis ) 

!  This routine ingests this time period of the first guess analysis.  The
!  data is passed back through the argument list after a successful
!  sequence of READs.

   USE obj_analysis

   IMPLICIT NONE
   INCLUDE 'netcdf.inc'

   INTEGER                                      :: met_ncid
   INTEGER                                      :: num3d , &
                                                   num2d , &
                                                   num1d 

   INCLUDE 'first_guess_size.inc'
   INCLUDE 'first_guess.inc'

   INTEGER                                      :: current_date_8 , & 
                                                   current_time_6

   INTEGER                                      :: loop_count     , & 
                                                   header_date_8  , &
                                                   header_time_6  , &
                                                   icount
   REAL , DIMENSION(jns_alloc,iew_alloc)        :: dummy
   INTEGER                                      :: extra_counter
   CHARACTER(LEN=8)                             :: dum_char
   INTEGER                                      :: yy, mm, dd, hh, mi, ss
   LOGICAL                                      :: print_analysis , &
                                                   first_loop
   REAL ,    DIMENSION(kbu_alloc)               :: pressure

   REAL, ALLOCATABLE, DIMENSION(:,:,:)          :: met_em_dum1, met_em_dum2
   INTEGER                                      :: kp, i

   INTEGER                                      :: rcode, ndims, nvars, ngatts, nunlimdimid
   INTEGER                                      :: iv, idims
   CHARACTER (LEN=10)                           :: var_name
   INTEGER, DIMENSION(6)                        :: var_shape
   INTEGER                                      :: var_type, var_ndims, var_natt 
   INTEGER, DIMENSION(4)                        :: istart, iend
   INTEGER, ALLOCATABLE, DIMENSION(:)           :: dim_values
   CHARACTER (LEN=31),ALLOCATABLE, DIMENSION(:) :: dim_names
   CHARACTER (LEN=19)                           :: date_char
   REAL, DIMENSION(16)                          :: rval
   CHARACTER (LEN=19)                           :: times_in_file

   INCLUDE 'error.inc'
   INCLUDE 'big_header.inc'

   INTERFACE
      INCLUDE 'error.int'
   END INTERFACE

   !  We need to keep track of where we are sticking the data for the FDDA option.

   IF ( initial_time ) THEN
      tt = first_time
   ELSE
      tt = second_time
   END IF 


   !  Make a copy of the met_em input file, so that we do not overwrite the 
   !  incoming data.

   !!!!! ALSO open the SFC file 
   !!!!! HOW do we know which domain we are dealing with?
   !!!!! Do we have the unexpanded domain information at this point?


   !  This loop allows the routine to READ multiple time periods, looking
   !  for the correct time (it will match current_date_8 and 
   !  current_time_6).  The first time through this loop we allocate the
   !  extra space for the additional 2d fields that may be present.  A test
   !  is required so that if we are not starting at the intitial time period
   !  in the first-guess, we don't try to re-allocate space which has not
   !  been de-allocated.

   first_loop = .TRUE.
   time_loop : DO

      loop3 = 0
      loop2 = 0
      loop1 = 0
      istart = 1
      small_header%staggering = 'M  '
      small_header%current_date = date_char

      rcode = nf_inq(met_ncid, ndims, nvars, ngatts, nunlimdimid)
      ALLOCATE ( dim_values(ndims) ) 
      ALLOCATE ( dim_names(ndims) ) 
      DO idims = 1, ndims
         rcode = nf_inq_dim(met_ncid, idims, dim_names(idims), dim_values(idims))
      ENDDO

      !! Test to see which date we have in the met_em files
      rcode = nf_inq_varid ( met_ncid, "Times", iv )
      rcode = nf_get_var_text(met_ncid, iv, times_in_file)
      IF ( trim(times_in_file) /= trim(date_char) ) THEN
         print*,"   WARNING: Mismatch between datastamp in file and date on file. "
         print*,"            Date on file:                  ", trim(date_char)
         print*,"            Datestamp in file:             ", trim(times_in_file)
         print*,"            OUTPUT file will contain date: ", trim(date_char)
         print*," "
      ENDIF

      variable_loop : DO iv = 1, nvars

         iend = 1

         rcode = nf_inq_var(met_ncid, iv, var_name, var_type, var_ndims, var_shape, var_natt)           
         IF ( var_ndims < 3 ) CYCLE variable_loop
         IF ( var_type /= 5 ) CYCLE variable_loop

         IF ( trim(var_name) /= "TT"       .AND. trim(var_name) /= "UU"       .AND. &
              trim(var_name) /= "VV"       .AND. trim(var_name) /= "GHT"      .AND. &
              trim(var_name) /= "RH"       .AND. trim(var_name) /= "PRES"     .AND. &
              trim(var_name) /= "XLAT_M"   .AND. trim(var_name) /= "XLONG_M"  .AND. &
              trim(var_name) /= "XLAT_V"   .AND. trim(var_name) /= "XLONG_V"  .AND. &
              trim(var_name) /= "SKINTEMP" .AND. trim(var_name) /= "HGT_M"    .AND. &
              trim(var_name) /= "PMSL"     .AND. trim(var_name) /= "SNOW" ) CYCLE variable_loop

         small_header%name = trim(var_name) 

         IF ( print_analysis ) print*,"Reading in: ", trim(var_name)

         ! Get the dimensions of this field
         DO idims=1,var_ndims-1
            iend(idims)  = dim_values(var_shape(idims))
         ENDDO
         IF ( ALLOCATED(met_em_dum1) ) DEALLOCATE(met_em_dum1)
         ALLOCATE ( met_em_dum1(iend(1),iend(2),iend(3)) )
         rcode = nf_get_var_real(met_ncid, iv, met_em_dum1)

   
         IF     ( var_ndims .EQ. 4 ) THEN
            loop3=loop3+1
            IF ( ALLOCATED(met_em_dum2) ) DEALLOCATE(met_em_dum2)
            ALLOCATE ( met_em_dum2(iew_alloc,jns_alloc,kbu_alloc) )
            met_em_dum2(:iend(1),:iend(2),:) = met_em_dum1(:iend(1),:iend(2),:)
            IF ( iend(1) < iew_alloc ) &
                 met_em_dum2(iew_alloc,:iend(2),:) = met_em_dum1(iend(1),:iend(2),:)
            IF ( iend(2) < jns_alloc ) &
                 met_em_dum2(:iend(1),jns_alloc,:) = met_em_dum1(:iend(1),iend(2),:)
            IF ( iend(1) < iew_alloc .AND. iend(2) < jns_alloc ) &
                 met_em_dum2(iew_alloc,jns_alloc,:) = met_em_dum1(iend(1),iend(2),:)

            ALLOCATE ( all_3d(loop3,tt)%array(jns_alloc,iew_alloc,kbu_alloc) ) 
            DO kp = 1,kbu_alloc
               CALL yx2xy ( met_em_dum2(1,1,kp) , iew_alloc , jns_alloc , all_3d(loop3,tt)%array(1,1,kp) )
            ENDDO
            all_3d(loop3,tt)%small_header = small_header

            !! We need to carry 3 arrays for U and V each as we need to keep track of the 
            !! initial C and A grids.
            IF (trim(var_name) == "UU"  ) THEN
               small_header%staggering = 'X  '
               small_header%name = 'UC'
               all_3d(loop3,tt)%small_header = small_header

               uC = all_3d(loop3,tt)%array
               uA = uC
               !BPR BEGIN
               !Appears that uA and uC are both dimensioned as the maximum number
               !of gridpoints in the x and y direction for any of T, U, or V:
               ! (max_x_mass_grid+1,max_y_mass_grid+1)
               !uC is presumably u on the WRF-native C-grid
               ! This should be dimensioned (max_x_mass_grid+1,max_y_mass_grid)
               !uA is presumably u on an unstaggered grid (presumably the WRF mass grid)
               ! This should be dimensioned (max_x_mass_grid,max_y_mass_grid)
               !We need to add "jns_alloc-1" to uC so that there are the same
               !number of grid points on both sides of the equation
               !uA(1:jns_alloc-1,1:iew_alloc-1,:) = ( uC(:,1:iew_alloc-1,:)+ &
               !                                      uC(:,2:iew_alloc  ,:) ) *.5
               uA(1:jns_alloc-1,1:iew_alloc-1,:) = ( uC(1:jns_alloc-1,1:iew_alloc-1,:)+ &
                                                     uC(1:jns_alloc-1,2:iew_alloc  ,:) ) *.5
               !BPR END
               small_header%name = 'UA'
               loop3=loop3+1
               ALLOCATE ( all_3d(loop3,tt)%array(jns_alloc,iew_alloc,kbu_alloc) ) 
               all_3d(loop3,tt)%array = uA
               all_3d(loop3,tt)%small_header = small_header
               small_header%name = 'UU'
               loop3=loop3+1
               ALLOCATE ( all_3d(loop3,tt)%array(jns_alloc,iew_alloc,kbu_alloc) ) 
               all_3d(loop3,tt)%array = uA
               all_3d(loop3,tt)%small_header = small_header
            ENDIF

            IF (trim(var_name) == "VV"  ) THEN               
               small_header%staggering = 'Y  '
               small_header%name = 'VC'
               all_3d(loop3,tt)%small_header = small_header

               vC = all_3d(loop3,tt)%array
               vA = vC
               !BPR BEGIN
               !See note on similar modification to uA/uC above
               !vA(1:jns_alloc-1,1:iew_alloc-1,:) = ( vC(1:jns_alloc-1,:,:)+ &
               !                                      vC(2:jns_alloc  ,:,:) ) *.5
               vA(1:jns_alloc-1,1:iew_alloc-1,:) = ( vC(1:jns_alloc-1,1:iew_alloc-1,:)+ &
                                                     vC(2:jns_alloc  ,1:iew_alloc-1,:) ) *.5
               !BPR END

               small_header%name = 'VA'
               loop3=loop3+1
               ALLOCATE ( all_3d(loop3,tt)%array(jns_alloc,iew_alloc,kbu_alloc) ) 
               all_3d(loop3,tt)%array = vA
               all_3d(loop3,tt)%small_header = small_header
               small_header%name = 'VV'
               loop3=loop3+1
               ALLOCATE ( all_3d(loop3,tt)%array(jns_alloc,iew_alloc,kbu_alloc) ) 
               all_3d(loop3,tt)%array = vA
               all_3d(loop3,tt)%small_header = small_header
            ENDIF

            ! We are keeping track of the 3D pressure fields so that we can get our hands on the surface pressure
            ! But we also need a 1D PRESSURE array for calculations
            IF (trim(var_name) == "PRES"  ) THEN               
               loop1=loop1+1
               small_header%name = 'PRESSURE'
               ALLOCATE ( all_1d(loop1)%array(kbu_alloc) ) 
               all_1d(loop1)%array(1)           = 100100.0
               all_1d(loop1)%array(2:kbu_alloc) = met_em_dum1(1,1,2:)
               all_1d(loop1)%small_header = small_header
            ENDIF

         ELSE IF( var_ndims .EQ. 3 ) THEN
            loop2=loop2+1

            IF ( ALLOCATED(met_em_dum2) ) DEALLOCATE(met_em_dum2)
            ALLOCATE ( met_em_dum2(iew_alloc,jns_alloc,1) )

                IF ( trim(var_name) == "XLAT_V" ) THEN
                   small_header%name = 'XLAT_D'
                   !BPR BEGIN
                   !Appears to be converting from latitude on V points to
                   !latitude on dot points which is presumably exactly in the
                   !center of 4 mass points?
                   !Appears that xlat_d grid point 1,1 is at u=1,v=1
                   !Except for the furthest west and furthest east dot points,
                   !do an average of the v-point latitudes in the west-east
                   !direction
                   !met_em_dum2(2:iend(1)-1,:,1) = ( met_em_dum1(1:iend(1)-1,:,1) + met_em_dum1(2:iend(1),:,1) ) *.5
                   met_em_dum2(2:iend(1),:,1) = ( met_em_dum1(1:iend(1)-1,:,1) + met_em_dum1(2:iend(1),:,1) ) *.5
                   !BPR END
                   deallocate (met_em_dum1)
                   allocate ( met_em_dum1(iew_alloc, jns_alloc-1,1) )
                   rcode = nf_inq_varid ( met_ncid, "XLAT_U", i )
                   rcode = nf_get_var_real ( met_ncid, i, met_em_dum1 )
                   !Now for the furthest west and furthest east dot points, do
                   !an average of the u-point latitudes in the south-north
                   !direction
                   met_em_dum2(1:iew_alloc:iew_alloc,2:jns_alloc-1,1) = ( met_em_dum1(1:iew_alloc:iew_alloc,1:jns_alloc-2,1) &
                                                                        + met_em_dum1(1:iew_alloc:iew_alloc,2:jns_alloc-1,1) &
                                                                        ) *.5
                   !This leaves the corners unfilled, so fill them
                   rcode = nf_get_att_real (met_ncid, nf_global, "corner_lats", rval) 
                   met_em_dum2(1,1,1)                 = rval(13)
                   met_em_dum2(iew_alloc,1,1)         = rval(14)
                   met_em_dum2(iew_alloc,jns_alloc,1) = rval(15)
                   met_em_dum2(1,jns_alloc,1)         = rval(16)
            ELSEIF ( trim(var_name) == "XLONG_V" ) THEN
                   small_header%name = 'XLONG_D'
                   !BPR BEGIN
                   !met_em_dum2(2:iend(1)-1,:,1) = ( met_em_dum1(1:iend(1)-1,:,1) + met_em_dum1(2:iend(1),:,1) ) *.5
                   met_em_dum2(2:iend(1),:,1) = ( met_em_dum1(1:iend(1)-1,:,1) + met_em_dum1(2:iend(1),:,1) ) *.5
                   !BPR END
                   deallocate (met_em_dum1)
                   allocate ( met_em_dum1(iew_alloc, jns_alloc-1,1) )
                   rcode = nf_inq_varid ( met_ncid, "XLONG_U", i )
                   rcode = nf_get_var_real ( met_ncid, i, met_em_dum1 )
                   met_em_dum2(1:iew_alloc:iew_alloc,2:jns_alloc-1,1) = ( met_em_dum1(1:iew_alloc:iew_alloc,1:jns_alloc-2,1) &
                                                                        + met_em_dum1(1:iew_alloc:iew_alloc,2:jns_alloc-1,1) &
                                                                        ) *.5
                   rcode = nf_get_att_real (met_ncid, nf_global, "corner_lons", rval) 
                   met_em_dum2(1,1,1)                 = rval(13)
                   met_em_dum2(iew_alloc,1,1)         = rval(14)
                   met_em_dum2(iew_alloc,jns_alloc,1) = rval(15)
                   met_em_dum2(1,jns_alloc,1)         = rval(16)
            ELSE
                   met_em_dum2(:iend(1),:iend(2),1) = met_em_dum1(:iend(1),:iend(2),1)
                   IF ( iend(1) < iew_alloc ) &
                      met_em_dum2(iew_alloc,:iend(2),1) = met_em_dum1(iend(1),:iend(2),1)
                   IF ( iend(2) < jns_alloc ) &
                      met_em_dum2(:iend(1),jns_alloc,1) = met_em_dum1(:iend(1),iend(2),1)
                   IF ( iend(1) < iew_alloc .AND. iend(2) < jns_alloc ) &
                      met_em_dum2(iew_alloc,jns_alloc,1) = met_em_dum1(iend(1),iend(2),1)
            ENDIF

            ALLOCATE ( all_2d(loop2,tt)%array(jns_alloc,iew_alloc) ) 
            CALL yx2xy ( met_em_dum2(1,1,1) , iew_alloc , jns_alloc , all_2d(loop2,tt)%array )
            all_2d(loop2,tt)%small_header = small_header

            IF (trim(var_name) == "PMSL"  ) THEN
               small_header%name = 'PMSL_C'
               loop2=loop2+1
               ALLOCATE ( all_2d(loop2,tt)%array(jns_alloc,iew_alloc) ) 
               all_2d(loop2,tt)%array(:,:) = all_2d(loop2-1,tt)%array(:,:)
               all_2d(loop2,tt)%small_header = small_header
            END IF

         END IF

      END DO variable_loop

      rcode = nf_close(met_ncid)


      !  We have read in all of the data for this time period.  Was it the right time?

      CALL split_date_char ( small_header%current_date , yy , mm , dd , hh , mi , ss )
      header_date_8 = yy * 10000 + mm * 100 + dd
      header_time_6 = hh * 10000 + mi * 100 + ss

      IF ( ( header_date_8 .NE. current_date_8 ) .OR. &
           ( header_time_6 .NE. current_time_6 ) ) THEN
         PRINT '(A,A,A)','Found wrong time, skipping past ',small_header%current_date,'.'
         CYCLE time_loop
     
      ELSE

         !  We have a valid time, so save the input data and exit out of here.

         num3d = loop3
         num2d = loop2
         num1d = loop1
   
         !  If this is the first time in this routine, we need to ALLOCATE some
         !  space for the extra 2d fields that will be coming this way.  There
         !  are usually only 13 2d arrays.  If there are more than that, we gotta
         !  have room for the data and associated header information.
   
         IF      ( ( first_loop) .AND. ( icount .EQ. 1 ) .AND. ( bhi(1) .EQ. 8 ) ) THEN
            CALL alloc_lots_of_2d_data ( jns_alloc , iew_alloc , num2d-13+1 )
         ELSE IF ( ( first_loop) .AND. ( icount .EQ. 1 ) ) THEN
            CALL alloc_lots_of_2d_data ( jns_alloc , iew_alloc , num2d-13   )
         END IF
         first_loop = .FALSE.
      
         !  Put the 3d data that we know about into storage arrays.
      
         three_D : DO loop_count = 1 , num3d
            IF      ( all_3d(loop_count,tt)%small_header%name(1:8) .EQ. 'TT      ' ) THEN
               t = all_3d(loop_count,tt)%array
               t ( jns_alloc,:iew_alloc-1,kbu_alloc) = t ( jns_alloc-1,:iew_alloc-1,kbu_alloc)
               t (:jns_alloc, iew_alloc  ,kbu_alloc) = t (:jns_alloc  , iew_alloc-1,kbu_alloc)
            ELSE IF ( all_3d(loop_count,tt)%small_header%name(1:8) .EQ. 'UU      ' ) THEN
               u = all_3d(loop_count,tt)%array
            ELSE IF ( all_3d(loop_count,tt)%small_header%name(1:8) .EQ. 'VV      ' ) THEN
               v = all_3d(loop_count,tt)%array
            ELSE IF ( all_3d(loop_count,tt)%small_header%name(1:8) .EQ. 'UA      ' ) THEN
               uA = all_3d(loop_count,tt)%array
            ELSE IF ( all_3d(loop_count,tt)%small_header%name(1:8) .EQ. 'VA      ' ) THEN
               vA = all_3d(loop_count,tt)%array
            ELSE IF ( all_3d(loop_count,tt)%small_header%name(1:8) .EQ. 'UC      ' ) THEN
               uC = all_3d(loop_count,tt)%array
            ELSE IF ( all_3d(loop_count,tt)%small_header%name(1:8) .EQ. 'VC      ' ) THEN
               vC = all_3d(loop_count,tt)%array
            ELSE IF ( all_3d(loop_count,tt)%small_header%name(1:8) .EQ. 'GHT     ' ) THEN
               h = all_3d(loop_count,tt)%array
               h( jns_alloc,:iew_alloc-1,kbu_alloc) = h( jns_alloc-1,:iew_alloc-1,kbu_alloc)
               h(:jns_alloc,iew_alloc   ,kbu_alloc) = h(:jns_alloc  , iew_alloc-1,kbu_alloc)
            ELSE IF ( all_3d(loop_count,tt)%small_header%name(1:8) .EQ. 'RH      ' ) THEN
               DO kp = 1, kbu_alloc
                 CALL clean_rh ( all_3d(loop_count,tt)%array(:,:,kp) , iew_alloc , jns_alloc , 0. , 100. )
               END DO 
               rh= all_3d(loop_count,tt)%array
               rh( jns_alloc,:iew_alloc-1,kbu_alloc) = rh( jns_alloc-1,:iew_alloc-1,kbu_alloc)
               rh(:jns_alloc,iew_alloc   ,kbu_alloc) = rh(:jns_alloc  , iew_alloc-1,kbu_alloc)
            ELSE IF ( all_3d(loop_count,tt)%small_header%name(1:8) .EQ. 'PRES    ' ) THEN
               pres= all_3d(loop_count,tt)%array
               pres( jns_alloc,:iew_alloc-1,kbu_alloc) = pres( jns_alloc-1,:iew_alloc-1,kbu_alloc)
               pres(:jns_alloc,iew_alloc   ,kbu_alloc) = pres(:jns_alloc  , iew_alloc-1,kbu_alloc)
            END IF
            
         END DO three_D
      
         !  Process the 2d data.  Even the fields we don't recognize get saved to be
         !  passed through to the rest of the modeling system.
      
         extra_counter = 0
         two_D : DO loop_count = 1 , num2d
            IF      ( all_2d(loop_count,tt)%small_header%name(1:8) .EQ. 'HGT_M   ' ) THEN
               terrain = all_2d(loop_count,tt)%array
               terrain     ( jns_alloc,:iew_alloc-1) = terrain     ( jns_alloc-1,:iew_alloc-1)
               terrain     (:jns_alloc,iew_alloc   ) = terrain     (:jns_alloc  , iew_alloc-1)
            ELSE IF ( all_2d(loop_count,tt)%small_header%name(1:8) .EQ. 'XLAT_M  ' ) THEN
               latitude_x = all_2d(loop_count,tt)%array
               latitude_x  ( jns_alloc,:iew_alloc-1) = latitude_x  ( jns_alloc-1,:iew_alloc-1)
               latitude_x  (:jns_alloc,iew_alloc   ) = latitude_x  (:jns_alloc  , iew_alloc-1)
            ELSE IF ( all_2d(loop_count,tt)%small_header%name(1:8) .EQ. 'XLONG_M ' ) THEN
               longitude_x = all_2d(loop_count,tt)%array
               longitude_x ( jns_alloc,:iew_alloc-1) = longitude_x ( jns_alloc-1,:iew_alloc-1)
               longitude_x (:jns_alloc,iew_alloc   ) = longitude_x (:jns_alloc  , iew_alloc-1)
            ELSE IF ( all_2d(loop_count,tt)%small_header%name(1:8) .EQ. 'XLAT_D  ' ) THEN
               latitude_d = all_2d(loop_count,tt)%array
            ELSE IF ( all_2d(loop_count,tt)%small_header%name(1:8) .EQ. 'XLONG_D ' ) THEN
               longitude_d = all_2d(loop_count,tt)%array
            ELSE IF ( all_2d(loop_count,tt)%small_header%name(1:8) .EQ. 'PMSL    ' ) THEN
               slp_x = all_2d(loop_count,tt)%array
               slp_x       ( jns_alloc,:iew_alloc-1) = slp_x       ( jns_alloc-1,:iew_alloc-1)
               slp_x       (:jns_alloc,iew_alloc   ) = slp_x       (:jns_alloc  , iew_alloc-1)
            ELSE IF ( all_2d(loop_count,tt)%small_header%name(1:8) .EQ. 'PMSL_C  ' ) THEN
               slp_C = all_2d(loop_count,tt)%array
               slp_C       ( jns_alloc,:iew_alloc-1) = slp_C       ( jns_alloc-1,:iew_alloc-1)
               slp_C       (:jns_alloc,iew_alloc   ) = slp_C       (:jns_alloc  , iew_alloc-1)
            ELSE IF ( all_2d(loop_count,tt)%small_header%name(1:8) .EQ. 'SKINTEMP' ) THEN
               sst = all_2d(loop_count,tt)%array
               sst         ( jns_alloc,:iew_alloc-1) = sst         ( jns_alloc-1,:iew_alloc-1)
               sst         (:jns_alloc,iew_alloc   ) = sst         (:jns_alloc  , iew_alloc-1)
            ELSE IF ( all_2d(loop_count,tt)%small_header%name(1:8) .EQ. 'SNOW' ) THEN
               snow = all_2d(loop_count,tt)%array
               snow         ( jns_alloc,:iew_alloc-1) = snow         ( jns_alloc-1,:iew_alloc-1)
               snow         (:jns_alloc,iew_alloc   ) = snow         (:jns_alloc  , iew_alloc-1)
            END IF
      
         END DO two_D
   
         !  We need to pick up the only important 1d array - pressure.
   
         one_D : DO loop_count = 1 , num1d
            IF      ( all_1d(loop_count)%small_header%name(1:8) .EQ. 'PRESSURE' ) THEN
               pressure = all_1d(loop_count)%array
            END IF
         END DO one_D
      

         !  If we are in this part of the IF test, we found the right time.  We can
         !  therefore exit the time loop.

         EXIT time_loop
  
      END IF
   
   END DO time_loop

END SUBROUTINE read_first_guess

!------------------------------------------------------------------------------

END MODULE first_guess

