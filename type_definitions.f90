! ****************************************************************************
! Copyright © 2014 Instituto Superior Técnico. All rights reserved. This 
! software is the copyrighted work of Instituto Superior Técnico. 
! Reproduction, in whole or in part, on the Internet, on CD-ROM or any 
! other medium, without the prior written consent of Instituto Superior 
! Técnico is prohibited.
! ****************************************************************************

module type_definitions

use hdf5
use hdf5_utilities

implicit none

private
  
public :: t_track, t_detector_spec, t_detector_ene, t_auxvecs, t_auxvecs_pow
public :: t_input
  
!-----------------------------------------------------------------------------  
type t_track

integer, dimension(:), allocatable :: n  
real(p_double), dimension(:), allocatable :: t
real(p_double), dimension(:), allocatable :: x1
real(p_double), dimension(:), allocatable :: x2
real(p_double), dimension(:), allocatable :: x3
real(p_double), dimension(:), allocatable :: p1
real(p_double), dimension(:), allocatable :: p2
real(p_double), dimension(:), allocatable :: p3
real(p_double), dimension(:), allocatable :: dt
real(p_double), dimension(:), allocatable :: ene
real(p_double), dimension(:), allocatable :: chargetemp
real(p_double) :: charge
integer :: tracksize, ndimtrack

real(p_double) :: dtinitial
  
real(p_double), dimension(:), allocatable :: beta1
real(p_double), dimension(:), allocatable :: beta2
real(p_double), dimension(:), allocatable :: beta3
  
real(p_double), dimension(:), allocatable :: betaDot1
real(p_double), dimension(:), allocatable :: betaDot2
real(p_double), dimension(:), allocatable :: betaDot3

real(p_double), dimension(:), allocatable :: g
  
end type t_track  
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
type t_waxis

character(len=10) :: waxistype

integer, dimension(:), pointer :: firstwvec
real(p_double), dimension(:), allocatable :: waxis, waxislocal

integer :: mynw, nw
real(p_double) :: dw, wmin, wmax

integer(hsize_t), dimension(1) :: dimsfwaxis
integer(hsize_t), dimension(1) :: dims_chunkwaxis
  
end type t_waxis  
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
type t_xaxis

integer, dimension(:), pointer :: firstvec
real(p_double), dimension(:), allocatable :: xaxis, xaxislocal

integer :: mynx, nx
real(p_double) :: dx, xmin, xmax

integer(hsize_t), dimension(1) :: dimsf_xaxis
integer(hsize_t), dimension(1) :: dims_chunk_xaxis
  
end type t_xaxis  
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
type t_detector_spec

! detectors for incoherent calculation  
real(p_double), dimension(:), pointer :: spec1d
real(p_double), dimension(:,:), pointer :: spec2d
real(p_double), dimension(:,:,:), pointer :: spec3d

! detectors for coherent calculation  
real(p_double), dimension(:), pointer :: spec_1d_c1_re, spec_1d_c1_im
real(p_double), dimension(:), pointer :: spec_1d_c2_re, spec_1d_c2_im
real(p_double), dimension(:), pointer :: spec_1d_c3_re, spec_1d_c3_im

real(p_double), dimension(:,:), pointer :: spec_2d_c1_re, spec_2d_c1_im
real(p_double), dimension(:,:), pointer :: spec_2d_c2_re, spec_2d_c2_im
real(p_double), dimension(:,:), pointer :: spec_2d_c3_re, spec_2d_c3_im

real(p_double), dimension(:,:,:), pointer :: spec_3d_c1_re, spec_3d_c1_im
real(p_double), dimension(:,:,:), pointer :: spec_3d_c2_re, spec_3d_c2_im
real(p_double), dimension(:,:,:), pointer :: spec_3d_c3_re, spec_3d_c3_im

! detectors for coherent calculation in 1D detector
real(p_double), dimension(:), pointer :: e1_1d_re, e1_1d_im
real(p_double), dimension(:), pointer :: e2_1d_re, e2_1d_im
real(p_double), dimension(:), pointer :: e3_1d_re, e3_1d_im

! detectors for coherent calculation in 2D detector
real(p_double), dimension(:,:), pointer :: e1_2d_re, e1_2d_im
real(p_double), dimension(:,:), pointer :: e2_2d_re, e2_2d_im
real(p_double), dimension(:,:), pointer :: e3_2d_re, e3_2d_im

! detectors for coherent calculation in 3D detector
real(p_double), dimension(:,:,:), pointer :: e1_3d_re, e1_3d_im
real(p_double), dimension(:,:,:), pointer :: e2_3d_re, e2_3d_im
real(p_double), dimension(:,:,:), pointer :: e3_3d_re, e3_3d_im

! coherence
logical :: coherent

integer(hsize_t), dimension(:), allocatable :: dimsf, dims_chunk
integer :: ndim, ncells1, ncells2
real(p_double) :: x0, y0, z0, x1detmin, x1detmax, x2detmin, x2detmax
character(len=15) :: detector_axis

type(t_waxis) :: waxis
  
end type t_detector_spec  
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
type t_detector_ene

real(p_double), dimension(:), pointer :: pow1d
real(p_double), dimension(:,:), pointer :: pow2d
integer(hsize_t), dimension(:), allocatable :: dimsf, dims_chunk
integer :: ndim, ncells1, ncells2
real(p_double) :: x0, y0, z0
real(p_double) :: x1detmin, x1detmax, x2detmin, x2detmax
character(len=15) :: detector_axis

type(t_xaxis) :: xaxis

end type t_detector_ene  
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
type t_auxvecs
  
real(kind=p_double), dimension(:), allocatable :: a, b, c, d, n1, n2, n3, ndotr, &
      n1MinusBeta1xbetaDot2, n1MinusBeta1xbetaDot3, n2MinusBeta2xbetaDot1, &
      n2MinusBeta2xbetaDot3, n3MinusBeta3xbetaDot1, n3MinusBeta3xbetaDot2, &
      ndotbeta, eta, r, faccorr
real(kind=p_double), dimension(:), allocatable :: temp1cRe, temp1cIm, &
      temp2cRe, temp2cIm, temp3cRe, temp3cIm
real(kind=p_double), dimension(:), allocatable :: exppartRe, exppartIm, &
      commonpart1Re, commonpart2Re, commonpart2Im,temp
  
end type t_auxvecs  
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
type t_auxvecs_pow
  
real(kind=p_double), dimension(:), allocatable :: a, b, c, d, n1, n2, n3, &
      r, oneOverR, n1MinusBeta1xbetaDot2, n1MinusBeta1xbetaDot3, &
      n2MinusBeta2xbetaDot1, n2MinusBeta2xbetaDot3, n3MinusBeta3xbetaDot1, &
      n3MinusBeta3xbetaDot2
      
real(kind=p_double), dimension(:), allocatable :: temp
  
end type t_auxvecs_pow  
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
type t_input

! from namelist waxis_parameters
real(kind=p_double) :: wmin, wmax
integer :: min_dec, num_dec, ppdec, wpoints
character(len=10) :: waxistype

! from namelist track_parameters
character(len=150) :: trackfile
character(len=15) :: track_select_type
integer :: nbegin, nend, npart, ndimtrack
integer, dimension(2) :: nrange
real(kind=p_double), dimension(2) :: x1range, trange
real(kind=p_double) :: x1min, x1max, tmin, tmax, enemin 
logical :: m_weight

! from namelist detector_parameters
integer :: ndim
character(len=20) :: diag_type
logical :: endpoints
character(len=15) :: detector_axis
integer :: ncells1, ncells2
real(kind=p_double) :: x0, y0, z0
real(kind=p_double) :: x1detmin, x1detmax, x2detmin, x2detmax
character(len=8) :: emissivity
logical :: coherent

! from namelist save_parameters
character(len=150) :: filename
logical :: parallelIO

end type t_input  
!-----------------------------------------------------------------------------



end module type_definitions