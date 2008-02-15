MODULE small_header_data

   TYPE sh
      CHARACTER (LEN= 4)     :: staggering
      CHARACTER (LEN=24)     :: current_date
      CHARACTER (LEN= 9)     :: name
   END TYPE sh

   TYPE(sh) :: small_header

END MODULE small_header_data
