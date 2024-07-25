! ****************************************************************************
! Copyright © 2014 Instituto Superior Técnico. All rights reserved. This 
! software is the copyrighted work of Instituto Superior Técnico. 
! Reproduction, in whole or in part, on the Internet, on CD-ROM or any 
! other medium, without the prior written consent of Instituto Superior 
! Técnico is prohibited.
! ****************************************************************************
module time

use hdf5_utilities

!-----------------------------------------------------------------------------  
public :: get_time_in_scs
!-----------------------------------------------------------------------------

contains

!---------------------------------------------------
function get_time_in_scs( )

real(kind=p_double) :: get_time_in_scs

! number of ticks since a fixed time in the past
real(kind=p_double), external :: timer_cpu_seconds

get_time_in_scs = timer_cpu_seconds()

end function get_time_in_scs
!---------------------------------------------------

end module time