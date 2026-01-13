module time

! use hdf5_utilities

!-----------------------------------------------------------  
public :: get_time_in_scs
!-----------------------------------------------------------

contains


!---------------------------------------------------
function get_time_in_scs( )

integer, parameter :: p_double = 8
real(kind=p_double) :: get_time_in_scs

! number of ticks since a fixed time in the past
real(kind=p_double), external :: timer_cpu_seconds

! get_time_in_scs = timer_cpu_seconds()
get_time_in_scs = 0.0

end function get_time_in_scs
!---------------------------------------------------

end module time
