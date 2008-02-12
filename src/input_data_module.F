MODULE input_data

   USE small_header_data

   TYPE input_fields_3d
      REAL , POINTER , DIMENSION(:,:,:) :: array
      TYPE(sh)                          :: small_header
   END TYPE input_fields_3d

   TYPE input_fields_2d
      REAL , POINTER , DIMENSION(:,:)   :: array
      TYPE(sh)                          :: small_header
   END TYPE input_fields_2d

   TYPE input_fields_1d
      REAL , POINTER , DIMENSION(:)     :: array
      TYPE(sh)                          :: small_header
   END TYPE input_fields_1d

   TYPE(input_fields_3d) , ALLOCATABLE , DIMENSION(:,:) :: all_3d
   TYPE(input_fields_2d) , ALLOCATABLE , DIMENSION(:,:) :: all_2d
   TYPE(input_fields_1d) , ALLOCATABLE , DIMENSION(:)   :: all_1d

   INTEGER :: loop3 , loop2 , loop1

   INTEGER :: first_time , second_time , tt
   
   LOGICAL :: initial_time

END MODULE input_data
