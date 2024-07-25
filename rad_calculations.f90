! ****************************************************************************
! Copyright © 2014 Instituto Superior Técnico. All rights reserved. This 
! software is the copyrighted work of Instituto Superior Técnico. 
! Reproduction, in whole or in part, on the Internet, on CD-ROM or any 
! other medium, without the prior written consent of Instituto Superior 
! Técnico is prohibited.
! ****************************************************************************

module rad_calculations

use type_definitions
use track_utilities
use spec_calculations
use ene_calculations
use hdf5_utilities
use time

implicit none

private
  
!-----------------------------------------------------------------------------  
public :: calc_radiation    
!-----------------------------------------------------------------------------

contains


! ****************************************************************************
! ----------------------------------------------------------------------------
subroutine calc_radiation( input, track, detector_ene, detector_spec, myid, comm ) 
! ----------------------------------------------------------------------------
type(t_input), intent(in) :: input
type(t_track), intent(inout) :: track
type(t_detector_ene), intent(inout), optional :: detector_ene
type(t_detector_spec), intent(inout), optional :: detector_spec
integer, intent(in) :: myid, comm

! local variables
real(kind=p_double) :: time1, time2, timesum
logical :: invalidTrack = .false.
integer :: numInvalidTracks = 0
integer :: p
logical :: last_track

time1 = 0.0d0
time2 = 0.0d0
timesum = 0.0d0

last_track = .false.

!###############################################################  
!---------------- beginning of particle loop -------------------
do p = 1, input%npart  
  
  if (p .eq. input%npart) last_track = .true.
    
  call setup_track( input%trackfile, track, input%ndimtrack, invalidTrack, comm )  
  
  if (invalidTrack) then
    numInvalidTracks = numInvalidTracks + 1
    if (p .lt. input%npart) then
      cycle
    else
      exit
    endif    
  endif
  
  call read_track( input%trackfile, track, p, comm )
  
  call track_selection( track, input%track_select_type, input%nbegin, input%nend, &
                        input%nrange, input%x1min, input%x1max, input%x1range, &
                        input%tmin, input%tmax, input%trange, input%enemin, comm )
  
  if (input%coherent .or. .not.(input%m_weight)) track%charge = 1.0d0
    
  call calc_beta( track )
  
  call calc_betaDot( track )
  
  select case (input%diag_type)
  
    case ("standard","farfield","farfieldEndPoints")
  
      !time1 = get_time_in_scs()
      
      call calc_spec( track, detector_spec, input, last_track )  
      
      time2 = get_time_in_scs()
      !print *, "(time2-time1) = ", (time2-time1)
      
    case ("energy")
    
      call calc_ene( track, detector_ene, input%diag_type, myid )     

    case default
    
      if (myid .eq. 0) print *, "diag_type not available !!"
      
  end select    

  call cleanup_track( track )

end do
!---------------- end of particle loop ------------------------

!print *, "myid: ", myid, " time spec = ", timesum

if ((myid .eq. 0) .and. (numInvalidTracks .gt. 0)) then
  print *, "* warning * Found ", &
           numInvalidTracks, "invalid tracks and ignored them"
  print *, " "
endif

!###############################################################

end subroutine calc_radiation
! ----------------------------------------------------------------------------
! ****************************************************************************



end module rad_calculations
