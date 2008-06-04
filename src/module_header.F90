MODULE header1

   CONTAINS

   SUBROUTINE read_header ( met_ncid , sndim , wedim )

      IMPLICIT NONE

      INCLUDE 'netcdf.inc'

      INTEGER :: met_ncid
      INTEGER :: rcode
      INTEGER :: i, ndims, nvars, ngatts, nunlimdimid
      INTEGER :: wedim, sndim
      INTEGER :: dim_val
      CHARACTER (LEN=31) :: dim_name


      !  Read the horizontal domain size from the input file

      rcode = nf_inq(met_ncid, ndims, nvars, ngatts, nunlimdimid)
      dims_loop : DO i = 1, ndims
         rcode = nf_inq_dim(met_ncid, i, dim_name, dim_val)
         IF ( dim_name == 'west_east_stag'    ) wedim  = dim_val
         IF ( dim_name == 'south_north_stag'  ) sndim  = dim_val
      ENDDO dims_loop


   END SUBROUTINE read_header

END MODULE header1
