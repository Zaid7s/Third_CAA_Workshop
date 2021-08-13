MODULE Precision_Def
  IMPLICIT NONE
   
  PRIVATE
  
  PUBLIC:: i_def,      r_def,      & ! precision for calculation
            i_def_sp,   r_def_sp,   & ! single precision (4 byte)
            i_def_dp,   r_def_dp      ! double precision (8 byte)
  
  INTEGER, PARAMETER:: i_def    = SELECTED_INT_KIND(9),  &
                        r_def    = SELECTED_REAL_KIND(12), &
                        i_def_sp = SELECTED_INT_KIND(9),   &
                        r_def_sp = SELECTED_REAL_KIND(6),  &
                        i_def_dp = SELECTED_INT_KIND(16),  &
                        r_def_dp = SELECTED_REAL_KIND(12)

END MODULE Precision_Def 
