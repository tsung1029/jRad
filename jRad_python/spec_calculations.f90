! ****************************************************************************
! Copyright © 2014 Instituto Superior Técnico. All rights reserved. This 
! software is the copyrighted work of Instituto Superior Técnico. 
! Reproduction, in whole or in part, on the Internet, on CD-ROM or any 
! other medium, without the prior written consent of Instituto Superior 
! Técnico is prohibited.
! ****************************************************************************

module spec_calculations

use hdf5_utilities
use type_definitions

implicit none

! parameters  
real (kind=p_double), parameter :: facQC = 2.42469d-6

private
  
!-----------------------------------------------------------------------------  
public :: calc_spec    
!-----------------------------------------------------------------------------

contains


! ****************************************************************************
! ----------------------------------------------------------------------------
subroutine calc_spec( track, detector, input, last_track ) 
! ----------------------------------------------------------------------------

type(t_track), intent(in) :: track
type(t_detector_spec), intent(inout) :: detector
type(t_input), intent(in) :: input
logical, intent(in) :: last_track
  
select case (input%diag_type)

  case('standard') 
    
    call calc_spec_standard( track, detector, input, last_track  ) 
    
  case('farfield','farfieldEndPoints') 
    
    call calc_spec_farfield( track, detector, input, last_track )  
    
  case default 
    	  
    print *, "ERROR: diag_type not available" 
    stop
      
end select
  
  
end subroutine calc_spec
! ----------------------------------------------------------------------------
! ****************************************************************************


! ****************************************************************************
! ----------------------------------------------------------------------------
subroutine calc_spec_standard( track, detector, input, last_track ) 
! ----------------------------------------------------------------------------

type(t_track), intent(in) :: track
type(t_detector_spec), intent(inout) :: detector
type(t_input), intent(in) :: input
logical, intent(in) :: last_track

! local variables
integer :: indcell1, indcell2, tracksize, ncells1, ncells2
real(kind=p_double) :: xi, xj, xk, x0, y0, z0, dx1det, dx2det
real(kind=p_double) :: x1detmin, x1detmax, x2detmin, x2detmax
character(len=4) :: detector_axis
type(t_auxvecs) :: vecs
integer :: ndimtrack = 0
real(kind=p_double) :: time1, time2, timesum

time1 = 0.0d0
time2 = 0.0d0
timesum = 0.0d0

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
allocate(vecs%ndotbeta(tracksize),vecs%r(tracksize))

allocate(vecs%temp1cRe(tracksize),vecs%temp1cIm(tracksize))
allocate(vecs%temp2cRe(tracksize),vecs%temp2cIm(tracksize))
allocate(vecs%temp3cRe(tracksize),vecs%temp3cIm(tracksize))

allocate(vecs%exppartRe(tracksize),vecs%exppartIm(tracksize))
allocate(vecs%commonpart1Re(tracksize))
allocate(vecs%commonpart2Re(tracksize),vecs%commonpart2Im(tracksize))

allocate(vecs%temp(tracksize))

call zeroVector( vecs%exppartRe, tracksize  )
call zeroVector( vecs%exppartIm, tracksize  )
call zeroVector( vecs%temp1cRe, tracksize  )
call zeroVector( vecs%temp1cIm, tracksize  )
call zeroVector( vecs%temp2cRe, tracksize  )
call zeroVector( vecs%temp2cIm, tracksize  )
call zeroVector( vecs%temp3cRe, tracksize  )
call zeroVector( vecs%temp3cIm, tracksize  )
call zeroVector( vecs%commonpart1Re, tracksize  )
call zeroVector( vecs%commonpart2Re, tracksize  )
call zeroVector( vecs%commonpart2Im, tracksize  )

x0 = detector%x0
y0 = detector%y0
z0 = detector%z0 
    
! separate calculations in 1D, or 2 or 3D
select case (detector%ndim)
    
  case(1)

    xi = x0
    xj = y0
    xk = z0
    
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
                     
    call calc_spec_standard_cell( track, detector, input, vecs, indcell1, indcell2, &
                                  last_track ) 
    
  case(2)
    
        
    ! Detector axis dx
    dx1det = (x1detmax-x1detmin)/(ncells1)
    
    
    ! ################ cycle in detector line  cells ###################
    do indcell1 = 1, ncells1
    
      select case (trim(detector_axis))
        case ("x1")
          xi = x1detmin + (indcell1-0.5)*dx1det
          xj = y0
          xk = z0
        case ("x2")
          xi = x0
          xj = x1detmin + (indcell1-0.5)*dx1det
          xk = z0
        case ("x3")
          xi = x0
          xj = y0
          xk = x1detmin + (indcell1-0.5)*dx1det
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
              
      call calc_spec_standard_cell( track, detector, input, vecs, indcell1, indcell2, &
                                    last_track ) 
	  
    enddo 
    ! ############## end of cycles on detector cells ################### 
      
  case(3)
      
    ! Detector axis dx
    dx1det = (x1detmax-x1detmin)/(ncells1)
    dx2det = (x2detmax-x2detmin)/(ncells2)
      
    ! ################ cycle in detector line  cells ###################
    do indcell2 = 1, ncells2
      do indcell1 = 1, ncells1
	   
        ! coordinates of center of detector cell
        select case (trim(detector_axis))
          case ("x1x2")
            xi = x1detmin + (indcell1-0.5)*dx1det
            xj = x2detmin + (indcell2-0.5)*dx2det
            xk = z0
          case ("x1x3")
            xi = x1detmin + (indcell1-0.5)*dx1det
            xj = y0
            xk = x2detmin + (indcell2-0.5)*dx2det
          case ("x2x3")
            xi = x0
            xj = x1detmin + (indcell1-0.5)*dx1det
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
	  	  
        call calc_spec_standard_cell( track, detector, input, vecs, indcell1, indcell2, &
                                      last_track ) 
	  
      enddo
    enddo
    ! ############## end of cycles on detector cells ###################
	      
end select

!print *, "time in cell loop ", timesum
  
!----- deallocation of auxiliary arrays for calculations --------
deallocate(vecs%a, vecs%b, vecs%c, vecs%d, vecs%n1MinusBeta1xbetaDot2)
deallocate(vecs%n1MinusBeta1xbetaDot3, vecs%n2MinusBeta2xbetaDot1)
deallocate(vecs%n2MinusBeta2xbetaDot3, vecs%n3MinusBeta3xbetaDot1)
deallocate(vecs%n3MinusBeta3xbetaDot2, vecs%n1,vecs%n2,vecs%n3)
deallocate(vecs%ndotbeta,vecs%r)

deallocate(vecs%temp1cRe,vecs%temp1cIm)
deallocate(vecs%temp2cRe,vecs%temp2cIm)
deallocate(vecs%temp3cRe,vecs%temp3cIm)

deallocate(vecs%exppartRe,vecs%exppartIm)
deallocate(vecs%commonpart1Re)
deallocate(vecs%commonpart2Re,vecs%commonpart2Im)
deallocate(vecs%temp)
  
end subroutine calc_spec_standard
! ----------------------------------------------------------------------------------------
! ****************************************************************************************



! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine calc_spec_standard_cell( track, detector, input, vecs, indcell1, indcell2, &
                                    last_track  ) 
! ----------------------------------------------------------------------------------------

type(t_track), intent(in) :: track
type(t_detector_spec), intent(inout) :: detector
type(t_input), intent(in) :: input
type(t_auxvecs), intent(inout) :: vecs
integer, intent(in) :: indcell1, indcell2
logical, intent(in) :: last_track

! local variables
integer :: ww, wwBegin
real (kind=p_double) :: w
integer :: tracksize
real (kind=p_double) :: temp1csumRe, temp1csumIm, temp2csumRe, temp2csumIm
real (kind=p_double) :: temp3csumRe, temp3csumIm
logical :: leave

leave = .false.

tracksize = track%tracksize


vecs%n1MinusBeta1xbetaDot2 = aMinusBResTimesC( vecs%n1, track%beta1, track%betaDot2, tracksize  )
vecs%n1MinusBeta1xbetaDot3 = aMinusBResTimesC( vecs%n1, track%beta1, track%betaDot3, tracksize )
vecs%n2MinusBeta2xbetaDot1 = aMinusBResTimesC( vecs%n2, track%beta2, track%betaDot1, tracksize )
vecs%n2MinusBeta2xbetaDot3 = aMinusBResTimesC( vecs%n2, track%beta2, track%betaDot3, tracksize )
vecs%n3MinusBeta3xbetaDot1 = aMinusBResTimesC( vecs%n3, track%beta3, track%betaDot1, tracksize )
vecs%n3MinusBeta3xbetaDot2 = aMinusBResTimesC( vecs%n3, track%beta3, track%betaDot2, tracksize )

vecs%ndotbeta = vecs%n1*track%beta1 + vecs%n2*track%beta2 + vecs%n3*track%beta3
  		  
vecs%a = &
    vecs%n2*(vecs%n1MinusBeta1xbetaDot2 - vecs%n2MinusBeta2xbetaDot1) + &
    vecs%n3*(vecs%n1MinusBeta1xbetaDot3 - vecs%n3MinusBeta3xbetaDot1)
vecs%b = &
   -vecs%n1*(vecs%n1MinusBeta1xbetaDot2 - vecs%n2MinusBeta2xbetaDot1) + &
    vecs%n3*(vecs%n2MinusBeta2xbetaDot3 - vecs%n3MinusBeta3xbetaDot2)
vecs%c = &
   -vecs%n1*(vecs%n1MinusBeta1xbetaDot3 - vecs%n3MinusBeta3xbetaDot1) - &
    vecs%n2*(vecs%n2MinusBeta2xbetaDot3 - vecs%n3MinusBeta3xbetaDot2)
vecs%d = 1.0d0 - vecs%ndotbeta
		    				
call calcCommonpart1RePart1( vecs%commonpart1Re, track%dt, vecs%d, tracksize  )
	  
! ################ cycle in frequency ###################	

! test if first cycle is with w=0
w = detector%waxis%waxislocal(1)
wwBegin = 1
if (w .eq. 0.0d0) wwBegin = 2

  

do ww = wwBegin, detector%waxis%mynw

  w = detector%waxis%waxislocal(ww)
      
  vecs%temp = w*(track%t+vecs%r)
  
  vecs%exppartRe = cos(vecs%temp)
  vecs%exppartIm = sin(vecs%temp)
  
  vecs%commonpart2Re = vecs%commonpart1Re*vecs%exppartRe
  vecs%commonpart2Im = vecs%commonpart1Re*vecs%exppartIm
     	  			  		
  vecs%temp1cRe = vecs%commonpart2Re*vecs%a
  vecs%temp1cRe(1) = 0.5d0*vecs%temp1cRe(1)
  vecs%temp1cRe(tracksize) = 0.5d0*vecs%temp1cRe(tracksize)
  
  vecs%temp1cIm = vecs%commonpart2Im*vecs%a
  vecs%temp1cIm(1) = 0.5d0*vecs%temp1cIm(1)
  vecs%temp1cIm(tracksize) = 0.5d0*vecs%temp1cIm(tracksize)
  
  vecs%temp2cRe = vecs%commonpart2Re*vecs%b
  vecs%temp2cRe(1) = 0.5d0*vecs%temp2cRe(1)
  vecs%temp2cRe(tracksize) = 0.5d0*vecs%temp2cRe(tracksize)
  
  vecs%temp2cIm = vecs%commonpart2Im*vecs%b
  vecs%temp2cIm(1) = 0.5d0*vecs%temp2cIm(1)
  vecs%temp2cIm(tracksize) = 0.5d0*vecs%temp2cIm(tracksize)
	
  vecs%temp3cRe = vecs%commonpart2Re*vecs%c
  vecs%temp3cRe(1) = 0.5d0*vecs%temp3cRe(1)
  vecs%temp3cRe(tracksize) = 0.5d0*vecs%temp3cRe(tracksize)
  
  vecs%temp3cIm = vecs%commonpart2Im*vecs%c
  vecs%temp3cIm(1) = 0.5d0*vecs%temp3cIm(1)
  vecs%temp3cIm(tracksize) = 0.5d0*vecs%temp3cIm(tracksize)
  
  ! Calculation of d2I/(dw dOmega) (energy radiated per unit of 
  ! solid angle per frequency unit)
  ! 1st method - integral of all integrand using trapezoidal rule

  if (input%emissivity .eq. "d2W/dwdS") then
  
    vecs%temp1cRe = vecs%temp1cRe / vecs%r
    vecs%temp1cIm = vecs%temp1cIm / vecs%r
    vecs%temp2cRe = vecs%temp2cRe / vecs%r
    vecs%temp2cIm = vecs%temp2cIm / vecs%r
    vecs%temp3cRe = vecs%temp3cRe / vecs%r
    vecs%temp3cIm = vecs%temp3cIm / vecs%r
    
  endif  
    
  temp1csumRe = sum(vecs%temp1cRe)
  temp1csumIm = sum(vecs%temp1cIm)
  temp2csumRe = sum(vecs%temp2cRe)
  temp2csumIm = sum(vecs%temp2cIm)
  temp3csumRe = sum(vecs%temp3cRe)
  temp3csumIm = sum(vecs%temp3cIm)
    
  
  ! Spectrum calculations
  if (input%coherent) then
  
    select case (detector%ndim)

      case(1)
      
        ! vector component 1
        detector%spec_1d_c1_re(ww) = detector%spec_1d_c1_re(ww) + temp1csumRe
        detector%spec_1d_c1_im(ww) = detector%spec_1d_c1_im(ww) + temp1csumIm
        ! vector component 2
        detector%spec_1d_c2_re(ww) = detector%spec_1d_c2_re(ww) + temp2csumRe
        detector%spec_1d_c2_im(ww) = detector%spec_1d_c2_im(ww) + temp2csumIm
        ! vector component 3
        detector%spec_1d_c3_re(ww) = detector%spec_1d_c3_re(ww) + temp3csumRe
        detector%spec_1d_c3_im(ww) = detector%spec_1d_c3_im(ww) + temp3csumIm
        
        if (last_track) then
          detector%spec1d(ww) = detector%spec1d(ww)+ &
               ( (detector%spec_1d_c1_re(ww)**2 + detector%spec_1d_c1_im(ww)**2) + &
                 (detector%spec_1d_c2_re(ww)**2 + detector%spec_1d_c2_im(ww)**2) + &
                 (detector%spec_1d_c3_re(ww)**2 + detector%spec_1d_c3_im(ww)**2) )/4.0d0
        endif
      
      case(2)
    
        ! vector component 1
        detector%spec_2d_c1_re(ww,indcell1) = detector%spec_2d_c1_re(ww,indcell1) + &
                                                temp1csumRe
        detector%spec_2d_c1_im(ww,indcell1) = detector%spec_2d_c1_im(ww,indcell1) + &
                                                temp1csumIm 
        ! vector component 2
        detector%spec_2d_c2_re(ww,indcell1) = detector%spec_2d_c2_re(ww,indcell1) + &
                                                temp2csumRe
        detector%spec_2d_c2_im(ww,indcell1) = detector%spec_2d_c2_im(ww,indcell1) + &
                                                temp2csumIm
        ! vector component 3                                                                                                                        
        detector%spec_2d_c3_re(ww,indcell1) = detector%spec_2d_c3_re(ww,indcell1) + &
                                                temp3csumRe
        detector%spec_2d_c3_im(ww,indcell1) = detector%spec_2d_c3_im(ww,indcell1) + &
                                                temp3csumIm 
                                                
        if (last_track) then
          detector%spec2d(ww,indcell1) = detector%spec2d(ww,indcell1)+ &
               ( (detector%spec_2d_c1_re(ww,indcell1)**2 + &
                  detector%spec_2d_c1_im(ww,indcell1)**2) + &
                 (detector%spec_2d_c2_re(ww,indcell1)**2 + &
                  detector%spec_2d_c2_im(ww,indcell1)**2) + &
                 (detector%spec_2d_c3_re(ww,indcell1)**2 + &
                  detector%spec_2d_c3_im(ww,indcell1)**2) )/4.0d0
        endif                                        
                                                                            	
      case(3)
        
        ! vector component 1                       
        detector%spec_3d_c1_re(ww,indcell1,indcell2) = &
                     detector%spec_3d_c1_re(ww,indcell1,indcell2) + temp1csumRe
        detector%spec_3d_c1_im(ww,indcell1,indcell2) = &
                     detector%spec_3d_c1_im(ww,indcell1,indcell2) + temp1csumIm 
        ! vector component 2                       
        detector%spec_3d_c2_re(ww,indcell1,indcell2) = &
                     detector%spec_3d_c2_re(ww,indcell1,indcell2) + temp2csumRe
        detector%spec_3d_c2_im(ww,indcell1,indcell2) = &
                     detector%spec_3d_c2_im(ww,indcell1,indcell2) + temp2csumIm                          
        ! vector component 3                       
        detector%spec_3d_c3_re(ww,indcell1,indcell2) = &
                     detector%spec_3d_c3_re(ww,indcell1,indcell2) + temp3csumRe
        detector%spec_3d_c3_im(ww,indcell1,indcell2) = &
                     detector%spec_3d_c3_im(ww,indcell1,indcell2) + temp3csumIm                 

        if (last_track) then             
          detector%spec3d(ww,indcell1,indcell2) = detector%spec3d(ww,indcell1,indcell2)+ &
               ( (detector%spec_3d_c1_re(ww,indcell1,indcell2)**2 + &
                  detector%spec_3d_c1_im(ww,indcell1,indcell2)**2) + &
                 (detector%spec_3d_c2_re(ww,indcell1,indcell2)**2 + &
                  detector%spec_3d_c2_im(ww,indcell1,indcell2)**2) + &
                 (detector%spec_3d_c3_re(ww,indcell1,indcell2)**2 + &
                  detector%spec_3d_c3_im(ww,indcell1,indcell2)**2) )/4.0d0
        endif
        
    end select  
  
  else
  
    select case (detector%ndim)

      case(1)
      
        detector%spec1d(ww) = detector%spec1d(ww)+ &
                              abs(track%charge)*( (temp1csumRe**2 + temp1csumIm**2) + &
                                (temp2csumRe**2 + temp2csumIm**2) + &
                                (temp3csumRe**2 + temp3csumIm**2) )/4.0d0
      
      case(2)
    
        detector%spec2d(ww,indcell1) = detector%spec2d(ww,indcell1)+ &
                              abs(track%charge)*( (temp1csumRe**2 + temp1csumIm**2) + &
                                (temp2csumRe**2 + temp2csumIm**2) + &
                                (temp3csumRe**2 + temp3csumIm**2) )/4.0d0
		                            	
      case(3)
                               
        detector%spec3d(ww,indcell1,indcell2) = detector%spec3d(ww,indcell1,indcell2)+ &
                              abs(track%charge)*( (temp1csumRe**2 + temp1csumIm**2) + &
                                (temp2csumRe**2 + temp2csumIm**2) + &
                                (temp3csumRe**2 + temp3csumIm**2) )/4.0d0 
        
    end select
    
  endif        
        
end do  ! ########## end of cycle on frequency ###############


end subroutine calc_spec_standard_cell
! ----------------------------------------------------------------------------
! ****************************************************************************


! ****************************************************************************
! ----------------------------------------------------------------------------
subroutine calc_spec_farfield( track, detector, input, last_track ) 
! ----------------------------------------------------------------------------

type(t_track), intent(in) :: track
type(t_detector_spec), intent(inout) :: detector
type(t_input), intent(in) :: input
logical, intent(in) :: last_track

! local variables
integer :: indcell1, indcell2, tracksize, ncells1, ncells2
real(kind=p_double) :: xi, xj, xk, x0, y0, z0, dx1det, dx2det
real(kind=p_double) :: x1detmin, x1detmax, x2detmin, x2detmax
character(len=4) :: detector_axis
type(t_auxvecs) :: vecs
integer :: ndimtrack = 0
real(kind=p_double) :: time1, time2, timesum

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
allocate(vecs%a(tracksize), vecs%b(tracksize), vecs%c(tracksize))
allocate(vecs%n1(tracksize),vecs%n2(tracksize),vecs%n3(tracksize))
allocate(vecs%ndotr(tracksize),vecs%ndotbeta(tracksize),vecs%r(tracksize))
!allocate(vecs%eta(tracksize),vecs%faccorr(tracksize))

allocate(vecs%temp1cRe(tracksize),vecs%temp1cIm(tracksize))
allocate(vecs%temp2cRe(tracksize),vecs%temp2cIm(tracksize))
allocate(vecs%temp3cRe(tracksize),vecs%temp3cIm(tracksize))

allocate(vecs%exppartRe(tracksize),vecs%exppartIm(tracksize))
allocate(vecs%commonpart1Re(tracksize))
allocate(vecs%commonpart2Re(tracksize),vecs%commonpart2Im(tracksize))

allocate(vecs%temp(tracksize))

call zeroVector( vecs%exppartRe, tracksize  )
call zeroVector( vecs%exppartIm, tracksize  )
call zeroVector( vecs%temp1cRe, tracksize  )
call zeroVector( vecs%temp1cIm, tracksize  )
call zeroVector( vecs%temp2cRe, tracksize  )
call zeroVector( vecs%temp2cIm, tracksize  )
call zeroVector( vecs%temp3cRe, tracksize  )
call zeroVector( vecs%temp3cIm, tracksize  )
call zeroVector( vecs%commonpart1Re, tracksize  )
call zeroVector( vecs%commonpart2Re, tracksize  )
call zeroVector( vecs%commonpart2Im, tracksize  )

x0 = detector%x0
y0 = detector%y0
z0 = detector%z0 

vecs%faccorr = (facQC/track%g)
    
! separate calculations in 1D, or 2 or 3D
select case (detector%ndim)
    
  case(1)

    xi = x0
    xj = y0
    xk = z0
    
    vecs%r = sqrt(xi**2 + xj**2 + xk**2)
      
    vecs%n1 = xi/vecs%r
    vecs%n2 = xj/vecs%r
    vecs%n3 = xk/vecs%r
    
    select case (ndimtrack)
      case (1)
        vecs%ndotr = vecs%n1*track%x1
      case (2)
        vecs%ndotr = vecs%n1*track%x1 + vecs%n2*track%x2
      case (3)
        vecs%ndotr = vecs%n1*track%x1 + vecs%n2*track%x2 + vecs%n3*track%x3
      case default
        print *, "Error: ndimtrack must be 1, 2, or 3"  
    end select
                 
    call calc_spec_farfield_cell( track, detector, input, vecs, indcell1, indcell2, &
                                  last_track ) 
    
	    
  case(2)
    
    ! Detector axis dx
    dx1det = (x1detmax-x1detmin)/(ncells1)
    
    ! ################ cycle in detector line  cells ###################
    do indcell1 = 1, ncells1
     
    
      select case (trim(detector_axis))
        case ("x1")
          xi = x1detmin + (indcell1-0.5)*dx1det
          xj = y0
          xk = z0
        case ("x2")
          xi = x0
          xj = x1detmin + (indcell1-0.5)*dx1det
          xk = z0
        case ("x3")
          xi = x0
          xj = y0
          xk = x1detmin + (indcell1-0.5)*dx1det
      end select
            
      
      vecs%r = sqrt(xi**2 + xj**2 + xk**2)
        
      vecs%n1 = xi/vecs%r
      vecs%n2 = xj/vecs%r
      vecs%n3 = xk/vecs%r
      
      select case (ndimtrack)
        case (1)
          vecs%ndotr = vecs%n1*track%x1
        case (2)
          vecs%ndotr = vecs%n1*track%x1 + vecs%n2*track%x2
        case (3)
          vecs%ndotr = vecs%n1*track%x1 + vecs%n2*track%x2 + vecs%n3*track%x3
        case default
          print *, "Error: ndimtrack must be 1, 2, or 3"  
      end select
            
      call calc_spec_farfield_cell( track, detector, input, vecs, indcell1, indcell2, &
                                    last_track ) 
	  
	    
    enddo 
    ! ############## end of cycles on detector cells ################### 
      
  case(3)
      
    ! Detector axis dx
    dx1det = (x1detmax-x1detmin)/(ncells1)
    dx2det = (x2detmax-x2detmin)/(ncells2)
      
    ! ################ cycle in detector line  cells ###################
    do indcell2 = 1, ncells2
      do indcell1 = 1, ncells1
	   
	    
        ! coordinates of center of detector cell
        select case (trim(detector_axis))
          case ("x1x2")
            xi = x1detmin + (indcell1-0.5)*dx1det
            xj = x2detmin + (indcell2-0.5)*dx2det
            xk = z0
          case ("x1x3")
            xi = x1detmin + (indcell1-0.5)*dx1det
            xj = y0
            xk = x2detmin + (indcell2-0.5)*dx2det
          case ("x2x3")
            xi = x0
            xj = x1detmin + (indcell1-0.5)*dx1det
            xk = x2detmin + (indcell2-0.5)*dx2det
        end select
	    
	  
        vecs%r = sqrt(xi**2 + xj**2 + xk**2)
	  
        vecs%n1 = xi/vecs%r
        vecs%n2 = xj/vecs%r
        vecs%n3 = xk/vecs%r
	  	
        select case (ndimtrack)
          case (1)
            vecs%ndotr = vecs%n1*track%x1
          case (2)
            vecs%ndotr = vecs%n1*track%x1 + vecs%n2*track%x2
          case (3)
            vecs%ndotr = vecs%n1*track%x1 + vecs%n2*track%x2 + vecs%n3*track%x3
          case default
            print *, "Error: ndimtrack must be 1, 2, or 3"  
        end select
	  	  
        call calc_spec_farfield_cell( track, detector, input, vecs, indcell1, indcell2, &
                                      last_track ) 

	  
      enddo
    enddo
    ! ############## end of cycles on detector cells ###################
      
end select

!print *, "time in cell loop ", timesum
  
!----- deallocation of auxiliary arrays for calculations --------
deallocate(vecs%a, vecs%b, vecs%c)
deallocate(vecs%n1,vecs%n2,vecs%n3)
deallocate(vecs%ndotbeta,vecs%ndotr,vecs%r)
!deallocate(vecs%eta,vecs%faccorr)           

deallocate(vecs%temp1cRe,vecs%temp1cIm)
deallocate(vecs%temp2cRe,vecs%temp2cIm)
deallocate(vecs%temp3cRe,vecs%temp3cIm)

deallocate(vecs%exppartRe,vecs%exppartIm)
deallocate(vecs%commonpart1Re)
deallocate(vecs%commonpart2Re,vecs%commonpart2Im)
deallocate(vecs%temp)
  
end subroutine calc_spec_farfield
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine calc_spec_farfield_cell( track, detector, input, vecs, indcell1, indcell2, &
                                    last_track ) 
! ----------------------------------------------------------------------------------------

type(t_track), intent(in) :: track
type(t_detector_spec), intent(inout) :: detector
type(t_input), intent(in) :: input
type(t_auxvecs), intent(inout) :: vecs
integer, intent(in) :: indcell1, indcell2
logical, intent(in) :: last_track

! local variables
integer :: ww, wwBegin
real (kind=p_double) :: w
integer :: tracksize
real (kind=p_double) :: temp1csumRe, temp1csumIm, temp2csumRe, temp2csumIm
real (kind=p_double) :: temp3csumRe, temp3csumIm
real (kind=p_double) :: endpoint1a, endpoint2a, endpoint1b, endpoint2b
real (kind=p_double) :: endpoint1c, endpoint2c
logical :: leave
integer :: myid, mpierr

call MPI_COMM_RANK(MPI_COMM_WORLD, myid, mpierr)

leave = .false.

tracksize = track%tracksize


vecs%ndotbeta = vecs%n1*track%p1 + vecs%n2*track%p2 + &
                vecs%n3*track%p3
                
vecs%ndotbeta = vecs%ndotbeta / track%g                
  		  
vecs%a = &
    vecs%n2*(vecs%n1*track%p2 - vecs%n2*track%p1) - &
    vecs%n3*(vecs%n3*track%p1 - vecs%n1*track%p3)
vecs%b = &
    vecs%n3*(vecs%n2*track%p3 - vecs%n3*track%p2) - &
    vecs%n1*(vecs%n1*track%p2 - vecs%n2*track%p1)
vecs%c = &
    vecs%n1*(vecs%n3*track%p1 - vecs%n1*track%p3) - &
    vecs%n2*(vecs%n2*track%p3 - vecs%n3*track%p2)

vecs%a = vecs%a / track%g
vecs%b = vecs%b / track%g
vecs%c = vecs%c / track%g

if (input%endpoints) then
  endpoint1a = vecs%a(1)/(1.0d0 - vecs%ndotbeta(1))
  endpoint2a = vecs%a(tracksize)/(1.0d0 - vecs%ndotbeta(tracksize))
  endpoint1b = vecs%b(1)/(1.0d0 - vecs%ndotbeta(1)) 
  endpoint2b = vecs%b(tracksize)/(1.0d0 - vecs%ndotbeta(tracksize))
  endpoint1c = vecs%c(1)/(1.0d0 - vecs%ndotbeta(1))
  endpoint2c = vecs%c(tracksize)/(1.0d0 - vecs%ndotbeta(tracksize))
endif
		    				
vecs%commonpart1Re = track%dt
	  
	  
! ################ cycle in frequency ###################	

! test if first cycle is with w=0
w = detector%waxis%waxislocal(1)
wwBegin = 1
if (w .eq. 0.0d0) wwBegin = 2

do ww = wwBegin, detector%waxis%mynw

     
  w = detector%waxis%waxislocal(ww)
  
  vecs%temp = w*(track%t-vecs%ndotr)
  
  vecs%exppartRe = cos(vecs%temp)
  vecs%exppartIm = sin(vecs%temp)
  
  vecs%commonpart2Re = vecs%commonpart1Re*vecs%exppartRe
  vecs%commonpart2Im = vecs%commonpart1Re*vecs%exppartIm
     	  			  		
  vecs%temp1cRe = vecs%commonpart2Re*vecs%a
  vecs%temp1cRe(1) = 0.5d0*vecs%temp1cRe(1)
  vecs%temp1cRe(tracksize) = 0.5d0*vecs%temp1cRe(tracksize)
  
  vecs%temp1cIm = vecs%commonpart2Im*vecs%a
  vecs%temp1cIm(1) = 0.5d0*vecs%temp1cIm(1)
  vecs%temp1cIm(tracksize) = 0.5d0*vecs%temp1cIm(tracksize)
  
  vecs%temp2cRe = vecs%commonpart2Re*vecs%b
  vecs%temp2cRe(1) = 0.5d0*vecs%temp2cRe(1)
  vecs%temp2cRe(tracksize) = 0.5d0*vecs%temp2cRe(tracksize)
  
  vecs%temp2cIm = vecs%commonpart2Im*vecs%b
  vecs%temp2cIm(1) = 0.5d0*vecs%temp2cIm(1)
  vecs%temp2cIm(tracksize) = 0.5d0*vecs%temp2cIm(tracksize)
	
  vecs%temp3cRe = vecs%commonpart2Re*vecs%c
  vecs%temp3cRe(1) = 0.5d0*vecs%temp3cRe(1)
  vecs%temp3cRe(tracksize) = 0.5d0*vecs%temp3cRe(tracksize)
  
  vecs%temp3cIm = vecs%commonpart2Im*vecs%c
  vecs%temp3cIm(1) = 0.5d0*vecs%temp3cIm(1)
  vecs%temp3cIm(tracksize) = 0.5d0*vecs%temp3cIm(tracksize)
  
  ! Calculation of d2I/(dw dOmega) (energy radiated per unit of 
  ! solid angle per frequency unit)
  ! 1st method - integral of all integrand using trapezoidal rule

  if (input%emissivity .eq. "d2W/dwdS") then
  
    vecs%temp1cRe = vecs%temp1cRe / vecs%r
    vecs%temp1cIm = vecs%temp1cIm / vecs%r
    vecs%temp2cRe = vecs%temp2cRe / vecs%r
    vecs%temp2cIm = vecs%temp2cIm / vecs%r
    vecs%temp3cRe = vecs%temp3cRe / vecs%r
    vecs%temp3cIm = vecs%temp3cIm / vecs%r
    
    endpoint1a = endpoint1a / vecs%r(1)
    endpoint2a = endpoint2a / vecs%r(tracksize)
    endpoint1b = endpoint1b / vecs%r(1)
    endpoint2b = endpoint2b / vecs%r(tracksize)
    endpoint1c = endpoint1c / vecs%r(1)
    endpoint2c = endpoint2c / vecs%r(tracksize)
  
  endif
   
  select case (input%endpoints)
  
    case (.false.)
    
      temp1csumRe = sum(vecs%temp1cRe)
      temp1csumIm = sum(vecs%temp1cIm)
      temp2csumRe = sum(vecs%temp2cRe)
      temp2csumIm = sum(vecs%temp2cIm)
      temp3csumRe = sum(vecs%temp3cRe)
      temp3csumIm = sum(vecs%temp3cIm)
  
    case (.true.)
  
      temp1csumRe = sqrt(abs(track%charge))*( ( endpoint2a*vecs%exppartRe(tracksize)/w - &
                      endpoint1a*vecs%exppartRe(1)/w ) + sum(vecs%temp1cIm) )
      temp1csumIm = sqrt(abs(track%charge))*( ( endpoint2a*vecs%exppartIm(tracksize)/w - &
                      endpoint1a*vecs%exppartIm(1)/w ) - sum(vecs%temp1cRe) )
      temp2csumRe = sqrt(abs(track%charge))*( ( endpoint2b*vecs%exppartRe(tracksize)/w - & 
                      endpoint1b*vecs%exppartRe(1)/w ) + sum(vecs%temp2cIm) )
      temp2csumIm = sqrt(abs(track%charge))*( ( endpoint2b*vecs%exppartIm(tracksize)/w - &
                      endpoint1b*vecs%exppartIm(1)/w ) - sum(vecs%temp2cRe) )
      temp3csumRe = ( ( endpoint2c*vecs%exppartRe(tracksize)/w - &
                      endpoint1c*vecs%exppartRe(1)/w ) + sum(vecs%temp3cIm) )
      temp3csumIm = sqrt(abs(track%charge))*( ( endpoint2c*vecs%exppartIm(tracksize)/w - &
                      endpoint1c*vecs%exppartIm(1)/w ) - sum(vecs%temp3cRe) )

  end select 
  

  ! Spectrum calculations
  if (input%coherent) then
  
    select case (detector%ndim)

      case(1)
      
        ! vector component 1
        detector%spec_1d_c1_re(ww) = detector%spec_1d_c1_re(ww) + temp1csumRe
        detector%spec_1d_c1_im(ww) = detector%spec_1d_c1_im(ww) + temp1csumIm
        ! vector component 2
        detector%spec_1d_c2_re(ww) = detector%spec_1d_c2_re(ww) + temp2csumRe
        detector%spec_1d_c2_im(ww) = detector%spec_1d_c2_im(ww) + temp2csumIm
        ! vector component 3
        detector%spec_1d_c3_re(ww) = detector%spec_1d_c3_re(ww) + temp3csumRe
        detector%spec_1d_c3_im(ww) = detector%spec_1d_c3_im(ww) + temp3csumIm
        
        if (last_track) then
          detector%spec1d(ww) = detector%spec1d(ww)+ (w**2)* &
               ( (detector%spec_1d_c1_re(ww)**2 + detector%spec_1d_c1_im(ww)**2) + &
                 (detector%spec_1d_c2_re(ww)**2 + detector%spec_1d_c2_im(ww)**2) + &
                 (detector%spec_1d_c3_re(ww)**2 + detector%spec_1d_c3_im(ww)**2) )/4.0d0
        endif
      
      case(2)
    
        ! vector component 1
        detector%spec_2d_c1_re(ww,indcell1) = detector%spec_2d_c1_re(ww,indcell1) + &
                                                temp1csumRe
        detector%spec_2d_c1_im(ww,indcell1) = detector%spec_2d_c1_im(ww,indcell1) + &
                                                temp1csumIm 
        ! vector component 2
        detector%spec_2d_c2_re(ww,indcell1) = detector%spec_2d_c2_re(ww,indcell1) + &
                                                temp2csumRe
        detector%spec_2d_c2_im(ww,indcell1) = detector%spec_2d_c2_im(ww,indcell1) + &
                                                temp2csumIm
        ! vector component 3                                                                                                                        
        detector%spec_2d_c3_re(ww,indcell1) = detector%spec_2d_c3_re(ww,indcell1) + &
                                                temp3csumRe
        detector%spec_2d_c3_im(ww,indcell1) = detector%spec_2d_c3_im(ww,indcell1) + &
                                                temp3csumIm 
                                                
        if (last_track) then
          detector%spec2d(ww,indcell1) = detector%spec2d(ww,indcell1) +  &
               (w**2)*( (detector%spec_2d_c1_re(ww,indcell1)**2 + &
                         detector%spec_2d_c1_im(ww,indcell1)**2) + &
                        (detector%spec_2d_c2_re(ww,indcell1)**2 + &
                         detector%spec_2d_c2_im(ww,indcell1)**2) + &
                        (detector%spec_2d_c3_re(ww,indcell1)**2 + &
                         detector%spec_2d_c3_im(ww,indcell1)**2) )/4.0d0
        endif                                        
                                                                            	
      case(3)
        
        ! vector component 1                       
        detector%spec_3d_c1_re(ww,indcell1,indcell2) = &
                     detector%spec_3d_c1_re(ww,indcell1,indcell2) + temp1csumRe
        detector%spec_3d_c1_im(ww,indcell1,indcell2) = &
                     detector%spec_3d_c1_im(ww,indcell1,indcell2) + temp1csumIm 
        ! vector component 2                       
        detector%spec_3d_c2_re(ww,indcell1,indcell2) = &
                     detector%spec_3d_c2_re(ww,indcell1,indcell2) + temp2csumRe
        detector%spec_3d_c2_im(ww,indcell1,indcell2) = &
                     detector%spec_3d_c2_im(ww,indcell1,indcell2) + temp2csumIm                          
        ! vector component 3                       
        detector%spec_3d_c3_re(ww,indcell1,indcell2) = &
                     detector%spec_3d_c3_re(ww,indcell1,indcell2) + temp3csumRe
        detector%spec_3d_c3_im(ww,indcell1,indcell2) = &
                     detector%spec_3d_c3_im(ww,indcell1,indcell2) + temp3csumIm                 

        if (last_track) then             
          detector%spec3d(ww,indcell1,indcell2) = detector%spec3d(ww,indcell1,indcell2)+ &
               (w**2)*( (detector%spec_3d_c1_re(ww,indcell1,indcell2)**2 + &
                         detector%spec_3d_c1_im(ww,indcell1,indcell2)**2) + &
                        (detector%spec_3d_c2_re(ww,indcell1,indcell2)**2 + &
                         detector%spec_3d_c2_im(ww,indcell1,indcell2)**2) + &
                        (detector%spec_3d_c3_re(ww,indcell1,indcell2)**2 + &
                         detector%spec_3d_c3_im(ww,indcell1,indcell2)**2) )/4.0d0
        endif
        
    end select
  
  else ! incoherent calculation
    
    select case (detector%ndim)

      case(1)
      
        detector%spec1d(ww) = detector%spec1d(ww)+ &
                              (w**2)*abs(track%charge)*( (temp1csumRe**2 + temp1csumIm**2) + &
                                       (temp2csumRe**2 + temp2csumIm**2) + &
                                       (temp3csumRe**2 + temp3csumIm**2) )/4.0d0
      
      case(2)
      
        detector%spec2d(ww,indcell1) = detector%spec2d(ww,indcell1)+ &
                              (w**2)*abs(track%charge)*( (temp1csumRe**2 + temp1csumIm**2) + &
                                (temp2csumRe**2 + temp2csumIm**2) + &
                                (temp3csumRe**2 + temp3csumIm**2) )/4.0d0
		                            	
      case(3)
		                               
        detector%spec3d(ww,indcell1,indcell2) = detector%spec3d(ww,indcell1,indcell2)+ &
                             (w**2)*abs(track%charge)*( (temp1csumRe**2 + temp1csumIm**2) + &
                                (temp2csumRe**2 + temp2csumIm**2) + &
                                (temp3csumRe**2 + temp3csumIm**2) )/4.0d0 
        
    end select 
    
  endif ! end of spectrum calculations
        
end do  ! ########## end of cycle on frequency ###############


end subroutine calc_spec_farfield_cell
! ----------------------------------------------------------------------------
! ****************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
function aMinusBResTimesC( vec1, vec2, vec3, tracksize  ) result(vec4)
! ----------------------------------------------------------------------------------------

integer, intent(in) :: tracksize
real (kind=p_double), dimension(tracksize), intent(in) :: vec1, vec2, vec3
real (kind=p_double), dimension(tracksize) :: vec4
! local variables
integer :: k

do k = 1, tracksize 
  vec4(k) = (vec1(k) - vec2(k))*vec3(k)
enddo
  
end function aMinusBResTimesC
! ----------------------------------------------------------------------------
! ****************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine zeroVector( vec1, tracksize  )
! ----------------------------------------------------------------------------------------

integer, intent(in) :: tracksize
real(kind=p_double), dimension(tracksize), intent(inout) :: vec1
! local variables
integer :: k

do k = 1, tracksize 
  vec1(k) = 0.0d0
enddo
  
end subroutine zeroVector
! ----------------------------------------------------------------------------
! ****************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine calcCommonpart1RePart1( vecOut, dt, d, tracksize  )
! ----------------------------------------------------------------------------------------

integer, intent(in) :: tracksize
real(kind=p_double), dimension(tracksize), intent(inout) :: vecOut
real(kind=p_double), dimension(tracksize), intent(in) :: dt, d
! local variables
integer :: k

do k = 1, tracksize 
  vecOut(k) = dt(k)/(d(k)*d(k))
enddo
  
end subroutine calcCommonpart1RePart1
! ----------------------------------------------------------------------------
! ****************************************************************************


end module spec_calculations
