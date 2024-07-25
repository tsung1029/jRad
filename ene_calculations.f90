! ****************************************************************************
! Copyright © 2014 Instituto Superior Técnico. All rights reserved. This 
! software is the copyrighted work of Instituto Superior Técnico. 
! Reproduction, in whole or in part, on the Internet, on CD-ROM or any 
! other medium, without the prior written consent of Instituto Superior 
! Técnico is prohibited.
! ****************************************************************************

module ene_calculations

use hdf5_utilities
use type_definitions

implicit none

private
  
!-----------------------------------------------------------------------------  
public :: calc_ene    
!-----------------------------------------------------------------------------

contains


! ****************************************************************************
! ----------------------------------------------------------------------------
subroutine calc_ene( track, detector, diag_type, myid ) 
! ----------------------------------------------------------------------------

type(t_track), intent(in) :: track
type(t_detector_ene), intent(inout) :: detector
integer, intent(in) :: myid
character(len=*), intent(in) :: diag_type
    
select case (diag_type)

  case('energy') 
    
    call calc_radiated_energy( track, detector, myid  ) 
  
  case('power') 
    	  
    print *, "ERROR: power diagnostic not available" 
    stop
    
  case('angular_power_dist')  
    	  
    print *, "ERROR: angular power distribution not available" 
    stop
          
  case default 
    	  
    print *, "ERROR: diag_type not available" 
    stop
      
end select
  
  
end subroutine calc_ene
! ----------------------------------------------------------------------------
! ****************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine calc_radiated_energy( track, detector, myid  ) 
! ----------------------------------------------------------------------------------------

type(t_track), intent(in) :: track
type(t_detector_ene), intent(inout) :: detector
integer, intent(in) :: myid
! local variables
integer :: indcell1, indcell2, tracksize, ncells1, ncells2
real(kind=p_double) :: xi, xj, xk, x0, y0, z0, dx1det, dx2det
real(kind=p_double) :: x1detmin, x1detmax, x2detmin, x2detmax
character(len=15) :: detector_axis
type(t_auxvecs_pow) :: vecs
integer :: ndimtrack = 0
real(kind=p_double) :: time1, time2, timesum
integer :: mynx

!time1 = 0.0d0
!time2 = 0.0d0
!timesum = 0.0d0

mynx = detector%xaxis%mynx

ndimtrack = track%ndimtrack   
tracksize = track%tracksize 
ncells1 = detector%ncells1
ncells2 = detector%ncells2
x1detmin = detector%x1detmin
x1detmax = detector%x1detmax
x2detmin = detector%x2detmin
x2detmax = detector%x2detmax
detector_axis = trim(detector%detector_axis)
  
!----- allocation of auxiliary arrays for calculations --------
allocate(vecs%a(tracksize), vecs%b(tracksize), vecs%c(tracksize),vecs%d(tracksize))

allocate(vecs%n1MinusBeta1xbetaDot2(tracksize), vecs%n1MinusBeta1xbetaDot3(tracksize))
allocate(vecs%n2MinusBeta2xbetaDot1(tracksize), vecs%n2MinusBeta2xbetaDot3(tracksize))
allocate(vecs%n3MinusBeta3xbetaDot1(tracksize), vecs%n3MinusBeta3xbetaDot2(tracksize))

allocate(vecs%n1(tracksize),vecs%n2(tracksize),vecs%n3(tracksize))

allocate(vecs%r(tracksize),vecs%oneOverR(tracksize))

allocate(vecs%temp(tracksize))

x0 = detector%x0
y0 = detector%y0
z0 = detector%z0 

    
! separate calculations in 1D, or 2D
select case (detector%ndim)
    
  case(1)

    print *, "Not implemented yet in 1D"
	    
  case(2)
    
    ! Detector axis dx
    dx1det = (x1detmax-x1detmin)/(ncells1)
    dx2det = (x2detmax-x2detmin)/(ncells2)
      
    ! ################ cycle in detector line  cells ###################
    !time1 = mpi_wtime()
    do indcell2 = 1, ncells2
      do indcell1 = 1, mynx
	   
        ! coordinates of center of detector cell
        select case (trim(detector_axis))
          case ("x1x2")
            xi = detector%xaxis%xaxislocal(indcell1)
            xj = x2detmin + (indcell2-0.5)*dx2det
            xk = z0
          case ("x1x3")
            xi = detector%xaxis%xaxislocal(indcell1)
            xj = y0
            xk = x2detmin + (indcell2-0.5)*dx2det
          case ("x2x3")
            xi = x0
            xj = detector%xaxis%xaxislocal(indcell1)
            xk = x2detmin + (indcell2-0.5)*dx2det
        end select
	    
        select case (ndimtrack)
          case (1)
            vecs%r = sqrt((track%x1-xi)**2 + xj**2 + xk**2)
            vecs%n1 = (xi-track%x1)/vecs%r
            vecs%n2 = xj/vecs%r
            vecs%n3 = xk/vecs%r
          case (2)
            vecs%r = sqrt((track%x1-xi)**2 + (track%x2-xj)**2 + xk**2)
            vecs%n1 = (xi-track%x1)/vecs%r
            vecs%n2 = (xj-track%x2)/vecs%r
            vecs%n3 = xk/vecs%r 
          case (3)
            vecs%r = sqrt((track%x1-xi)**2 + (track%x2-xj)**2 + (track%x3-xk)**2)
            vecs%n1 = (xi-track%x1)/vecs%r
            vecs%n2 = (xj-track%x2)/vecs%r
            vecs%n3 = (xk-track%x3)/vecs%r
          case default
            print *, "Error: ndimtrack must be 1, 2, or 3"  
        end select
	  	  
        call calc_radiated_energy_cell( track, detector, vecs, indcell1, indcell2 ) 
	  
      enddo 
    enddo
    ! ############## end of cycles on detector cells ###################
    !time2 = mpi_wtime()
    !timesum = timesum + (time2-time1)
  
  case default

    print *, "Error: not implemented"
  	      
end select

!print *, "time in cell loop ", timesum
  
!----- deallocation of auxiliary arrays for calculations --------
deallocate(vecs%a,vecs%b,vecs%c,vecs%d)

deallocate(vecs%n1MinusBeta1xbetaDot2, vecs%n1MinusBeta1xbetaDot3)
deallocate(vecs%n2MinusBeta2xbetaDot1, vecs%n2MinusBeta2xbetaDot3)
deallocate(vecs%n3MinusBeta3xbetaDot1, vecs%n3MinusBeta3xbetaDot2)

deallocate(vecs%n1,vecs%n2,vecs%n3)

deallocate(vecs%r,vecs%oneOverR)

deallocate(vecs%temp)
  
    

end subroutine calc_radiated_energy
! ----------------------------------------------------------------------------
! ****************************************************************************

! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine calc_radiated_energy_cell( track, detector, vecs, indcell1, indcell2  ) 
! ----------------------------------------------------------------------------------------

type(t_track), intent(in) :: track
type(t_detector_ene), intent(inout) :: detector
type(t_auxvecs_pow), intent(inout) :: vecs
integer, intent(in) :: indcell1, indcell2

! local variables
integer :: tracksize
real (kind=p_double) :: dx1det, dx2det
! ----------------------------------------------------------------------------------------

dx1det = (detector%x1detmax-detector%x1detmin)/(detector%ncells1)
dx2det = (detector%x2detmax-detector%x2detmin)/(detector%ncells2)

tracksize = track%tracksize
		
vecs%n1MinusBeta1xbetaDot2 = (vecs%n1-track%beta1) * track%betaDot2
vecs%n1MinusBeta1xbetaDot3 = (vecs%n1-track%beta1) * track%betaDot3
vecs%n2MinusBeta2xbetaDot1 = (vecs%n2-track%beta2) * track%betaDot1
vecs%n2MinusBeta2xbetaDot3 = (vecs%n2-track%beta2) * track%betaDot3
vecs%n3MinusBeta3xbetaDot1 = (vecs%n3-track%beta3) * track%betaDot1
vecs%n3MinusBeta3xbetaDot2 = (vecs%n3-track%beta3) * track%betaDot2
						 
! Calculation of dP/dOmega (power radiated per unit of solid angle)
vecs%a =  vecs%n2*(vecs%n1MinusBeta1xbetaDot2 - vecs%n2MinusBeta2xbetaDot1) + &
          vecs%n3*(vecs%n1MinusBeta1xbetaDot3 - vecs%n3MinusBeta3xbetaDot1)
vecs%b = -vecs%n1*(vecs%n1MinusBeta1xbetaDot2 - vecs%n2MinusBeta2xbetaDot1) + &
          vecs%n3*(vecs%n2MinusBeta2xbetaDot3 - vecs%n3MinusBeta3xbetaDot2)
vecs%c = -vecs%n1*(vecs%n1MinusBeta1xbetaDot3 - vecs%n3MinusBeta3xbetaDot1) - &
          vecs%n2*(vecs%n2MinusBeta2xbetaDot3 - vecs%n3MinusBeta3xbetaDot2)
vecs%d =  1.0d0 - (vecs%n1*track%beta1 + vecs%n2*track%beta2 + vecs%n3*track%beta3)
		
detector%pow2d(indcell1,indcell2) = detector%pow2d(indcell1,indcell2) + &
           abs(track%charge)*real(sum((vecs%a*vecs%a + vecs%b*vecs%b + vecs%c*vecs%c)* &
           track%dt(1)*(dx1det*dx2det/((vecs%r**2)*(vecs%d**5)))))

! CHANGE PART OF DT		
 
end subroutine calc_radiated_energy_cell
! ----------------------------------------------------------------------------
! ****************************************************************************

end module ene_calculations