! ****************************************************************************
! Copyright © 2014 Instituto Superior Técnico. All rights reserved. This 
! software is the copyrighted work of Instituto Superior Técnico. 
! Reproduction, in whole or in part, on the Internet, on CD-ROM or any 
! other medium, without the prior written consent of Instituto Superior 
! Técnico is prohibited.
! ****************************************************************************

module utilities

  use hdf5_utilities
  use type_definitions
  
implicit none
 
private
    
! #####################################

  public :: make_waxis
  public :: make_xparallelaxis
  public :: read_input_file, check_input
  public :: setup_detector, cleanup_detector
  public :: setup_detector_ene, cleanup_detector_ene
  public :: save_spec, save_ene
     
! #####################################

contains



! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine setup_detector( detector, input )
! ----------------------------------------------------------------------------------------

type(t_detector_spec), intent(inout) :: detector
type(t_input), intent(in) :: input

! local variables
integer :: nw, mynw, ndim, ncells1, ncells2

ndim = input%ndim
nw = detector%waxis%nw
mynw = detector%waxis%mynw
ncells1 = input%ncells1
ncells2 = input%ncells2

! allocation of result arrays  
allocate(detector%dimsf(ndim),detector%dims_chunk(ndim))
  
select case (ndim)   
  case (1)
    allocate(detector%spec1d(mynw))
    detector%spec1d = 0.0d0
    ! hdf5 arrays dimensions
    detector%dimsf = (/nw/)
    detector%dims_chunk = (/mynw/)
    if (input%coherent) then
      allocate(detector%spec_1d_c1_re(mynw))
      allocate(detector%spec_1d_c1_im(mynw))
      allocate(detector%spec_1d_c2_re(mynw))
      allocate(detector%spec_1d_c2_im(mynw))
      allocate(detector%spec_1d_c3_re(mynw))
      allocate(detector%spec_1d_c3_im(mynw))
      detector%spec_1d_c1_re = 0.0d0
      detector%spec_1d_c1_im = 0.0d0
      detector%spec_1d_c2_re = 0.0d0
      detector%spec_1d_c2_im = 0.0d0
      detector%spec_1d_c3_re = 0.0d0
      detector%spec_1d_c3_im = 0.0d0
    endif
  case (2)
    allocate(detector%spec2d(mynw,ncells1))
    detector%spec2d = 0.0d0
    ! hdf5 arrays dimensions
    detector%dimsf = (/nw,ncells1/)
    detector%dims_chunk = (/mynw,ncells1/)
    if (input%coherent) then
      allocate(detector%spec_2d_c1_re(mynw,ncells1))
      allocate(detector%spec_2d_c1_im(mynw,ncells1))
      allocate(detector%spec_2d_c2_re(mynw,ncells1))
      allocate(detector%spec_2d_c2_im(mynw,ncells1))
      allocate(detector%spec_2d_c3_re(mynw,ncells1))
      allocate(detector%spec_2d_c3_im(mynw,ncells1))
      detector%spec_2d_c1_re = 0.0d0
      detector%spec_2d_c1_im = 0.0d0
      detector%spec_2d_c2_re = 0.0d0
      detector%spec_2d_c2_im = 0.0d0
      detector%spec_2d_c3_re = 0.0d0
      detector%spec_2d_c3_im = 0.0d0
    endif
  case (3)
    allocate(detector%spec3d(mynw,ncells1,ncells2))
    detector%spec3d = 0.0d0
    ! hdf5 arrays dimensions
    detector%dimsf = (/nw,ncells1,ncells2/)
    detector%dims_chunk = (/mynw,ncells1,ncells2/)
    if (input%coherent) then
      allocate(detector%spec_3d_c1_re(mynw,ncells1,ncells2))
      allocate(detector%spec_3d_c1_im(mynw,ncells1,ncells2))
      allocate(detector%spec_3d_c2_re(mynw,ncells1,ncells2))
      allocate(detector%spec_3d_c2_im(mynw,ncells1,ncells2))
      allocate(detector%spec_3d_c3_re(mynw,ncells1,ncells2))
      allocate(detector%spec_3d_c3_im(mynw,ncells1,ncells2))
      detector%spec_3d_c1_re = 0.0d0
      detector%spec_3d_c1_im = 0.0d0
      detector%spec_3d_c2_re = 0.0d0
      detector%spec_3d_c2_im = 0.0d0
      detector%spec_3d_c3_re = 0.0d0
      detector%spec_3d_c3_im = 0.0d0
    endif
end select
  
detector%ndim = input%ndim
detector%ncells1 = input%ncells1
detector%ncells2 = input%ncells2
detector%x0 = input%x0
detector%y0 = input%y0
detector%z0 = input%z0
detector%x1detmin = input%x1detmin
detector%x1detmax = input%x1detmax
detector%x2detmin = input%x2detmin
detector%x2detmax = input%x2detmax
detector%detector_axis = input%detector_axis
detector%coherent = input%coherent
  
end subroutine setup_detector
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine cleanup_detector( detector )
! ----------------------------------------------------------------------------------------

type(t_detector_spec), intent(inout) :: detector

! deallocations
  
deallocate(detector%dimsf,detector%dims_chunk)
  
if (associated(detector%spec1d)) then
  deallocate(detector%spec1d)
elseif (associated(detector%spec2d)) then
  deallocate(detector%spec2d)
elseif (associated(detector%spec3d)) then
  deallocate(detector%spec3d)
endif

if (detector%coherent) then
  
  select case (detector%ndim)
    
    case (1)
    
      if (associated(detector%e1_1d_re)) deallocate(detector%e1_1d_re)
      if (associated(detector%e1_1d_im)) deallocate(detector%e1_1d_im)
      if (associated(detector%e2_1d_re)) deallocate(detector%e2_1d_re)
      if (associated(detector%e2_1d_im)) deallocate(detector%e2_1d_im)
      if (associated(detector%e3_1d_re)) deallocate(detector%e3_1d_re)
      if (associated(detector%e3_1d_im)) deallocate(detector%e3_1d_im)
      
      if (associated(detector%spec_1d_c1_re)) deallocate(detector%spec_1d_c1_re)
      if (associated(detector%spec_1d_c1_im)) deallocate(detector%spec_1d_c1_im)
      if (associated(detector%spec_1d_c2_re)) deallocate(detector%spec_1d_c2_re)
      if (associated(detector%spec_1d_c2_im)) deallocate(detector%spec_1d_c2_im)
      if (associated(detector%spec_1d_c3_re)) deallocate(detector%spec_1d_c3_re)
      if (associated(detector%spec_1d_c3_im)) deallocate(detector%spec_1d_c3_im)
      
    case (2)
    
      if (associated(detector%e1_2d_re)) deallocate(detector%e1_2d_re)
      if (associated(detector%e1_2d_im)) deallocate(detector%e1_2d_im)
      if (associated(detector%e2_2d_re)) deallocate(detector%e2_2d_re)
      if (associated(detector%e2_2d_im)) deallocate(detector%e2_2d_im)
      if (associated(detector%e3_2d_re)) deallocate(detector%e3_2d_re)
      if (associated(detector%e3_2d_im)) deallocate(detector%e3_2d_im)
      
      if (associated(detector%spec_2d_c1_re)) deallocate(detector%spec_2d_c1_re)
      if (associated(detector%spec_2d_c1_im)) deallocate(detector%spec_2d_c1_im)
      if (associated(detector%spec_2d_c2_re)) deallocate(detector%spec_2d_c2_re)
      if (associated(detector%spec_2d_c2_im)) deallocate(detector%spec_2d_c2_im)
      if (associated(detector%spec_2d_c3_re)) deallocate(detector%spec_2d_c3_re)
      if (associated(detector%spec_2d_c3_im)) deallocate(detector%spec_2d_c3_im)
      
    case (3)
    
      if (associated(detector%e1_3d_re)) deallocate(detector%e1_3d_re)
      if (associated(detector%e1_3d_im)) deallocate(detector%e1_3d_im)
      if (associated(detector%e2_3d_re)) deallocate(detector%e2_3d_re)
      if (associated(detector%e2_3d_im)) deallocate(detector%e2_3d_im)
      if (associated(detector%e3_3d_re)) deallocate(detector%e3_3d_re)
      if (associated(detector%e3_3d_im)) deallocate(detector%e3_3d_im)
      
      if (associated(detector%spec_3d_c1_re)) deallocate(detector%spec_3d_c1_re)
      if (associated(detector%spec_3d_c1_im)) deallocate(detector%spec_3d_c1_im)
      if (associated(detector%spec_3d_c2_re)) deallocate(detector%spec_3d_c2_re)
      if (associated(detector%spec_3d_c2_im)) deallocate(detector%spec_3d_c2_im)
      if (associated(detector%spec_3d_c3_re)) deallocate(detector%spec_3d_c3_re)
      if (associated(detector%spec_3d_c3_im)) deallocate(detector%spec_3d_c3_im)
    
  end select
  
endif

if (associated(detector%waxis%firstwvec)) deallocate(detector%waxis%firstwvec)
if (allocated(detector%waxis%waxis)) deallocate(detector%waxis%waxis)
if (allocated(detector%waxis%waxislocal)) deallocate(detector%waxis%waxislocal)
  
end subroutine cleanup_detector
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine setup_detector_ene( detector, input )
! ----------------------------------------------------------------------------------------

type(t_detector_ene), intent(inout) :: detector
type(t_input), intent(in) :: input

! local variables
integer :: mynx, ncells1, ncells2, ndim

ndim = input%ndim
mynx = detector%xaxis%mynx
ncells1 = input%ncells1
ncells2 = input%ncells2


! allocation of result arrays  
allocate(detector%dimsf(ndim),detector%dims_chunk(ndim))
  
select case (ndim)   
  case (2)
    allocate(detector%pow2d(mynx,ncells2))
    detector%pow2d = 0.0d0
    ! hdf5 arrays dimensions
    detector%dimsf = (/ncells1,ncells2/)
    detector%dims_chunk = (/mynx,ncells2/)
  case default
    print *, "Error; only 2D energy diagnostics are implemented"
end select
  
detector%ndim = input%ndim
detector%ncells1 = input%ncells1
detector%ncells2 = input%ncells2
detector%x0 = input%x0
detector%y0 = input%y0
detector%z0 = input%z0
detector%x1detmin = input%x1detmin
detector%x1detmax = input%x1detmax
detector%x2detmin = input%x2detmin
detector%x2detmax = input%x2detmax
detector%detector_axis = input%detector_axis
  
end subroutine setup_detector_ene
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine cleanup_detector_ene( detector )
! ----------------------------------------------------------------------------------------

type(t_detector_ene), intent(inout) :: detector

! deallocations
  
deallocate(detector%dimsf,detector%dims_chunk)
  
if (associated(detector%pow1d)) then
  deallocate(detector%pow1d)
elseif (associated(detector%pow2d)) then
  deallocate(detector%pow2d)
endif

end subroutine cleanup_detector_ene
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine read_input_file( input, comm )
! ----------------------------------------------------------------------------------------

type(t_input), intent(inout) :: input
integer, intent(in) :: comm

! local variables
integer :: ndim, ndimtrack, nbegin, nend, npart, ncells1,&
                          ncells2, min_dec, num_dec, ppdec, wpoints
character(len=20) :: diag_type
character(len=10) :: waxistype
character(len=150) :: trackfile, filename
character(len=15) :: track_select_type, detector_axis
real(kind=p_double) :: x1min, x1max, tmin, tmax, enemin       
real(kind=p_double), dimension(2) :: x1range, trange
integer, dimension(2) :: nrange
logical :: parallelIO, endpoints, coherent, m_weight
character(len=8) :: emissivity
real(kind=p_double) :: x0, y0, z0, x1detmin, x1detmax, &
                       x2detmin, x2detmax, wmin, wmax                  
integer :: mpierr, ios, myid, numprocs
logical :: stop_program

!
! wmin, wmax, wpoints = min,max, and # of points in omega-axis
! waxistype = "linear"  "log10", "linbydec"
! min_dec,num_dec, ppdec = axis for "linbydec" case
! m_weight = are the particles weighted? (TRUE/FALSE)
!
namelist /waxis_parameters/ wmin, wmax, wpoints, waxistype, min_dec, &
                            num_dec, ppdec
! trackfile = name of the trackfile                     
! ndimtrack = number of dimensions in the simulation
! track_select_type = "nrange", "nmin", "nmax", "x1min", "x1max", "x1range","tmin","tmax","trange","enemin", "npart"
! 
namelist /track_parameters/ trackfile, ndimtrack, track_select_type, &
                            nbegin, nend, nrange, x1min, x1max, x1range, &
                            tmin, tmax, trange, enemin, npart, m_weight    

                    
!
! diag_type = "energy", "standard", "farfield", "farfieldEndPoint"
! detector_axis = "x1","x2","x3"
! x0,y0,z0 = location of the detector in (x,y,z)  x=x1, y=x2, z=x3
! ncells1, ncells2 = # of cells in the detector
! coherent = flag for coherent/incoherent spectrum calculation (DEFAULT is FALSE)
! endpoints = logical flag to suggest whether an endpoint (x0,y0,z0) is specified
! emissivity = logical flag to indicate whether or not the radiation is normalized by the radius R
!
namelist /detector_parameters/ ndim, diag_type, endpoints, &
                               detector_axis, ncells1, ncells2, x0, y0, z0, &  
                               x1detmin, x1detmax, x2detmin, x2detmax, emissivity, &
                               coherent         

namelist /save_parameters/ filename, parallelIO

! Set default values

! waxis_parameters
wmin = -1.0d0
wmax = -1.0d0
wpoints=-1
waxistype = "-"
min_dec=0
num_dec=0
ppdec=-1

! track_parameters
trackfile = "-"
ndimtrack=-1
track_select_type = "none"
nbegin=-1
nend = -1
nrange(1:2) = 0
x1min = 0.0d0
x1max = 0.0d0
x1range(1:2) = 0.0d0
tmin = 0.0d0
tmax = 0.0d0
trange = (/0.0d0,0.0d0/)
enemin = 0.0d0
npart=-1
m_weight = .false.

! detector_parameters
ndim = -1
diag_type="-"
endpoints=.true.
detector_axis = ""
ncells1 = -1
ncells2=-1
x0 = -1.0d0
y0 = -1.0d0
z0 = -1.0d0
x1detmin = 0.0d0
x1detmax = 0.0d0
x2detmin = 0.0d0
x2detmax = 0.0d0
emissivity = "d2W/dwdO"
coherent = .false.

! save_parameters
filename = "out.h5"
parallelIO=.false.

call MPI_COMM_RANK( comm, myid, mpierr )
call MPI_COMM_SIZE( comm, numprocs, mpierr)
if (mpierr .ne. 0) print *, "Error getting myid"
  
!  --------------------- Read and check input data -------------------------------------
if (myid .eq. 0) then
  ! Open input file
  open (unit=2, iostat=ios, file="input_file", status="old", &
        access="sequential", action="read")
  if (ios .ne. 0) print *, "Error opening input data file"
  ! Read namelist track_parameters
  read (unit=2, nml=track_parameters)
  ! Read namelist detector_parameters
  read (unit=2, nml=detector_parameters)
  
  ! Read namelist waxis_parameters
  select case (diag_type)
    case ("energy")
    
    case default
    
      print *, "diag_type picked not available"
      
    case ("standard","farfield","farfieldEndPoints")  
    
      print *, "Reading paramters for spectrum"
    
      read (unit=2, nml=waxis_parameters)
      
      print *, "waxistype read = ", waxistype
      
  end select    
  ! Read namelist save_parameters
  read (unit=2, nml=save_parameters)
  ! Close file
  close(unit=2, iostat=ios)
  if (ios .ne. 0) print *, "Error closing input data file"
  !  --------------------- Add input data to input variable ------------------------------
  input%ndim = ndim
  input%diag_type = diag_type
  input%coherent = coherent
  input%emissivity = emissivity
  input%endpoints = endpoints
  input%trackfile = trackfile
  input%ndimtrack = ndimtrack
  input%track_select_type = track_select_type
  input%nbegin = nbegin
  input%nend = nend
  input%nrange = nrange
  input%x1min = x1min
  input%x1max = x1max
  input%x1range = x1range
  input%tmin = tmin
  input%tmax = tmax
  input%trange = trange
  input%enemin = enemin
  input%npart = npart
  input%m_weight = m_weight
  input%detector_axis = detector_axis
  input%ncells1 = ncells1
  input%ncells2 = ncells2
  input%x0 = x0
  input%y0 = y0
  input%z0 = z0
  input%x1detmin = x1detmin
  input%x1detmax = x1detmax
  input%x2detmin = x2detmin
  input%x2detmax = x2detmax
  input%wmin = wmin
  input%wmax = wmax
  input%wpoints = wpoints
  input%waxistype = waxistype
  input%min_dec = min_dec
  input%num_dec = num_dec
  input%ppdec = ppdec
  input%filename = filename
  input%parallelIO = parallelIO
		
  !  -------------------------------------------------------------------------------------
 		
  ! Check if parameters read are ok
  print *, " "
  print *, "..........................................................."
  print *, " "
  print *, "Checking input file parameters..."
  print *, " "
  call check_input( input, stop_program, numprocs)                  
  print *, "Finished checking input file parameters."
  print *, " "      
  print *, "..........................................................."
  print *, " "                    
    
  print *, " " 
  print *, "------------- Input values ---------------"
  print *, " "
  print *, "....... general parameters ........."
  print *, "ndim = ", ndim
  print *, "diag_type = ", diag_type
  select case (diag_type)
    
    case ("standard")
      
    case ("farfield","farfieldEndPoints")
    
      print *, "endpoints = ", endpoints
    
    case ("energy")
    
    case default
    
  end select
  
  print *, " "
  print *, "........ track parameters .........."
  print *, "trackfile = ", trim(trackfile)
  print *, "ndimtrack = ", ndimtrack
  print *, "track_select_type = ", track_select_type
  select case (trim(track_select_type))
    case ("nrange")
      print *, "nrange = ", nrange
    case ("nmin")  
      print *, "nmin = ", nbegin
    case ("nmax")  
      print *, "nend = ", nend
    case ("x1min")
      print *, "x1min = ", x1min
    case ("x1max")
      print *, "x1max = ", x1max
    case ("x1range")
      print *, "x1range = ", x1range
    case ("tmin")
      print *, "tmin = ", tmin
    case ("tmax")
      print *, "tmax = ", tmax
    case ("trange")
      print *, "trange = ", trange
    case ("enemin")
      print *, "enemin = ", enemin
  end select
  print *, "npart = ", npart
  if (m_weight) then
    print *, "Macroparticle weight on - only works for incoherent calculations"
  endif  
  print *, " "
  print *, "....... detector parameters ........"
  print *, "emissivity = ", emissivity
  print *, "diag_type = ", diag_type
  print *, "ndim = ", ndim
  if (coherent) print *, "coherent = ", coherent
  if (ndim .gt. 1) then 
    print *, "detector_axis = ", trim(detector_axis)
  else if (ndim .eq. 1) then
    print *, "spectrum calculated in a point (0D in space)"
  endif  
  select case (trim(detector_axis))
    case ("x1")
      print *, "y0 = ", y0
      print *, "z0 = ", z0
    case ("x2")
      print *, "x0 = ", x0
      print *, "z0 = ", z0
    case ("x3")
      print *, "x0 = ", x0
      print *, "y0 = ", y0
    case ("x1x2")
      print *, "z0 = ", z0
    case ("x1x3")
      print *, "y0 = ", y0
    case ("x2x3")
      print *, "x0 = ", x0
    case default
      print *, "x0 = ", x0
      print *, "y0 = ", y0
      print *, "z0 = ", z0  
  end select  
  if ( ((ndim .gt. 1) .and. (diag_type .ne. 'energy')) .or. &
       ((ndim .eq. 2) .and. (diag_type .eq. 'energy')) ) then
    print *, "ncells1 = ", ncells1
    print *, "x1detmin = ", x1detmin
    print *, "x1detmax = ", x1detmax
  endif   
  if ( ((ndim .gt. 2) .and. (diag_type .ne. 'energy')) .or. &
       ((ndim .eq. 2) .and. (diag_type .eq. 'energy')) ) then
    print *, "ncells2 = ", ncells2
    print *, "x2detmin = ", x2detmin
    print *, "x2detmax = ", x2detmax
  endif  
  
  if (diag_type .eq. 'energy') then
  else 
    print *, "......... waxis parameters ........."
    print *, "waxistype = ", input%waxistype
    if ((trim(input%waxistype) .eq. 'linear') .or. &
        (trim(input%waxistype) .eq. 'log10')) then
      print *, "wmin = ", input%wmin
      print *, "wmax = ", input%wmax
      print *, "wpoints = ", input%wpoints
    else
      print *, "min_dec = ", input%min_dec
      print *, "num_dec = ", input%num_dec
      print *, "ppdec = ", input%ppdec
    endif  
  endif
  print *, " "
  print *, ".......... save parameters ........."
  print *, "filename = ", trim(filename)
  print *, "parallelIO = ", parallelIO
  print *, " "
  print *, "------------------------------------------"
  print *, " "

endif
!  -------------------------------------------------------------------------------------

! Checking if program should be aborted
call MPI_BCAST(stop_program,1,MPI_LOGICAL,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting stop_program"
if  (stop_program) then
  if (myid .eq. 0) print *, "Error!! Stopping program"
  if (myid .eq. 0) print *, " "
  ! Close MPI interface
  call mpi_finalize(mpierr)
  stop
endif  

! Broadcast data (to remainder processes)
! If I'm process 0, broadcast track data to remaining processes
!........ track parameters ..........
call MPI_BCAST(input%trackfile,150,MPI_CHARACTER,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting trackfile"
call MPI_BCAST(input%ndimtrack,1,MPI_INTEGER,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting ndimtrack"
call MPI_BCAST(input%track_select_type,15,MPI_CHARACTER,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting track_select_type"
call MPI_BCAST(input%nbegin,1,MPI_INTEGER,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting nbegin"
call MPI_BCAST(input%nend,1,MPI_INTEGER,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting nend"
call MPI_BCAST(input%nrange,2,MPI_INTEGER,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting nrange"
call MPI_BCAST(input%x1min,1,MPI_DOUBLE_PRECISION,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting x1min"
call MPI_BCAST(input%x1max,1,MPI_DOUBLE_PRECISION,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting x1max"
call MPI_BCAST(input%x1range,2,MPI_DOUBLE_PRECISION,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting x1range"
call MPI_BCAST(input%tmin,1,MPI_DOUBLE_PRECISION,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting tmin"
call MPI_BCAST(input%tmax,1,MPI_DOUBLE_PRECISION,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting tmax"
call MPI_BCAST(input%trange,2,MPI_DOUBLE_PRECISION,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting trange"
call MPI_BCAST(input%enemin,1,MPI_DOUBLE_PRECISION,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting enemin"
call MPI_BCAST(input%npart,1,MPI_INTEGER,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting npart"
call MPI_BCAST(input%m_weight,1,MPI_LOGICAL,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting m_weight"
!....... detector parameters ........
call MPI_BCAST(input%emissivity,8,MPI_CHARACTER,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting emissivity"
call MPI_BCAST(input%ndim,1,MPI_INTEGER,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting ndim"
call MPI_BCAST(input%diag_type,20,MPI_CHARACTER,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting diag_type"
call MPI_BCAST(input%coherent,1,MPI_LOGICAL,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting coherent"
call MPI_BCAST(input%endpoints,1,MPI_LOGICAL,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting endpoints"
call MPI_BCAST(input%detector_axis,15,MPI_CHARACTER,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting detector_axis"
call MPI_BCAST(input%ncells1,1,MPI_INTEGER,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting ncells1"
call MPI_BCAST(input%ncells2,1,MPI_INTEGER,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting ncells2"
call MPI_BCAST(input%x0,1,MPI_DOUBLE_PRECISION,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting x0"
call MPI_BCAST(input%y0,1,MPI_DOUBLE_PRECISION,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting y0"
call MPI_BCAST(input%z0,1,MPI_DOUBLE_PRECISION,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting z0"
call MPI_BCAST(input%x1detmin,1,MPI_DOUBLE_PRECISION,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting x1detmin"
call MPI_BCAST(input%x1detmax,1,MPI_DOUBLE_PRECISION,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting x1detmax"
call MPI_BCAST(input%x2detmin,1,MPI_DOUBLE_PRECISION,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting x2detmin"
call MPI_BCAST(input%x2detmax,1,MPI_DOUBLE_PRECISION,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting x2detmax"
!......... waxis parameters .........
call MPI_BCAST(input%wmin,1,MPI_DOUBLE_PRECISION,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting wmin"
call MPI_BCAST(input%wmax,1,MPI_DOUBLE_PRECISION,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting wmax"
call MPI_BCAST(input%wpoints,1,MPI_INTEGER,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting wpoints"
call MPI_BCAST(input%waxistype,10,MPI_CHARACTER,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting waxistype"
call MPI_BCAST(input%min_dec,1,MPI_INTEGER,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting min_dec"
call MPI_BCAST(input%num_dec,1,MPI_INTEGER,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting num_dec"
call MPI_BCAST(input%ppdec,1,MPI_INTEGER,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting ppdec"
!.......... save parameters .........
call MPI_BCAST(input%filename,150,MPI_CHARACTER,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting filename"
call MPI_BCAST(input%parallelIO,1,MPI_LOGICAL,0,comm,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting parallelIO"

end subroutine read_input_file
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine check_input( input, bad_value, numprocs )               
! ----------------------------------------------------------------------------------------

type(t_input), intent(in) :: input                                   
logical, intent(inout) :: bad_value
integer, intent(in) :: numprocs
  
bad_value = .false.

! Checks to input-deck (input_file)

! Check tracks parameters
call check_track_parameters( input%trackfile, input%ndimtrack, &
                             input%npart, input%track_select_type, input%nbegin, &
                             input%nend, input%nrange, input%tmin, input%tmax, &
                             input%trange, input%x1min, input%x1max, input%x1range, &
                             input%enemin, bad_value )

select case (input%diag_type)

  case ("standard","farfield","farfieldEndPoints")

    ! Check spectrum detector parameters            
    call check_spec_detector_parameters( input%ndim, input%ncells1, &
            input%ncells2, input%min_dec, input%num_dec, input%ppdec, input%wpoints, &
            numprocs, input%diag_type, input%emissivity, input%waxistype, input%detector_axis, &
            input%x0, input%y0, input%z0, input%x1detmin, input%x1detmax, input%x2detmin,&
            input%x2detmax, input%wmin, input%wmax, bad_value )
    
  case ("energy")
  
    call check_ene_detector_parameters( input%ndim, input%ncells1, &
            input%ncells2, numprocs, input%detector_axis, &
            input%x0, input%y0, input%z0, input%x1detmin, input%x1detmax, &
            input%x2detmin, input%x2detmax, bad_value )
            
  case default
  
    print *, "Error: diag_type invalid"
    print *, "Error: no input-deck checks for this diag_type"
    print *, " "
    
end select

      
end subroutine check_input
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine check_track_parameters( trackfile, ndimtrack, npart, &
                                   track_select_type, nbegin, nend, nrange, tmin, tmax, &
                                   trange, x1min, x1max, x1range, enemin, bad_value )
! ----------------------------------------------------------------------------------------
integer, intent(in) :: ndimtrack, nbegin, nend, npart
character(len=*), intent(in) :: trackfile, track_select_type
real(kind=p_double), intent(in) :: x1min, x1max, tmin, tmax, enemin                       
real(kind=p_double), dimension(2), intent(in) :: x1range, trange
integer, dimension(2), intent(in) :: nrange               
logical, intent(inout) :: bad_value

! Check if parameters for tracks are ok

if (trackfile=="-") then
  print *, "error: missing track filename"
  bad_value = .true.
elseif ((ndimtrack .lt. 1) .or. (ndimtrack .gt. 3)) then  
  print *, "error: ndimtrack must be > 0 and < 4"
  bad_value = .true.
elseif (npart .lt. 1) then
  print *, "error: npart must be > 0"
  bad_value = .true.
endif
  
select case (trim(track_select_type))

  case ("nmin")
  
    if (nbegin .lt. 1) then
      print *, "error: nbegin must be >= 1"
      bad_value = .true.
    endif
    
  case ("nmax")
  
    if (nend .lt. 3) then
      print *, "error: nend must be >= 3"
      bad_value = .true.
    endif  
    
  case ("nrange")
  
    if ((nrange(1) .lt. 1) .or. (nrange(2) .lt. 1)) then
      print *, "error: nbegin and nend must be > 0"
      bad_value = .true.
    elseif ( nrange(2) .lt. nrange(1) ) then
      print *, "error: nbegin and nend must be > nbegin"
      bad_value = .true.
    endif  
  
  case("tmin")
    
    if (tmin .lt. 0.0d0) then
      print *, "error: tmin must be >= 0"
      bad_value = .true.
    endif
  
  case ("tmax")

    if (tmax .le. 0.0d0) then
      print *, "error: tmax must be > 0"
      bad_value = .true.
    endif  
  
  case ("trange")
  
    if (trange(2) .le. trange(1)) then
      print *, "error: trange(2) must be > trange(1)"
      bad_value = .true.
    elseif ( (trange(1) .lt. 0.0) .or. (trange(2) .le. 0.0) ) then
      print *, "error: trange(1) and trange(2) must be > 0"
      bad_value = .true. 
    endif  
        
  case("x1range")
    
    if (x1range(2) .le. x1range(1))  then
      print *, "error: x1range(2) must be > x1range(1)"
      bad_value = .true.  
    endif  
  
  case ("enemin")
    
    if (enemin .lt. 0) then
      print *, "error: enemin must be > 0"
      bad_value = .true.
    endif  
  
  case("x1min","x1max","none")
  
  case default
    print *, "error: track_select_type ",trim(track_select_type)," not available"
    bad_value = .true.  
    
end select
    
! Check for parameters present that are incompatible with track selection type

if ((tmin .ne. 0) .and. (trim(track_select_type) .ne. "tmin")) then
  print *, "error: tmin present but track_select_type is incompatible"
  bad_value = .true.
elseif ((tmax .ne. 0) .and. (trim(track_select_type) .ne. "tmax")) then
  print *, "error: tmax present but track_select_type is incompatible"
  bad_value = .true.
elseif ( ((trange(1) .ne. 0) .or. (trange(2) .ne. 0) ) .and. &
         (trim(track_select_type) .ne. "trange") ) then
  print *, "error: trange present but track_select_type is incompatible"
  bad_value = .true.      
elseif ((x1min .ne. 0) .and. (trim(track_select_type) .ne. "x1min")) then
  print *, "error: x1min present but track_select_type is incompatible"
  bad_value = .true.  
elseif ((x1max .ne. 0) .and. (trim(track_select_type) .ne. "x1max")) then
  print *, "error: x1max present but track_select_type is incompatible"
  bad_value = .true.  
elseif ( ((x1range(1) .ne. 0) .or. (x1range(2) .ne. 0) ) .and. &
         (trim(track_select_type) .ne. "x1range") ) then
  print *, "error: x1range present but track_select_type is incompatible"
  bad_value = .true.  
elseif ((enemin .ne. 0) .and. (trim(track_select_type) .ne. "enemin")) then
  print *, "error: enemin present but track_select_type is incompatible"
  bad_value = .true.
elseif ( (nbegin .ne. -1) .and. (trim(track_select_type) .ne. "iter") ) then
  print *, "error: nbegin present but track_select_type is incompatible"
  bad_value = .true.    
elseif ( (nend .ne. -1) .and. (trim(track_select_type) .ne. "iter") ) then
  print *, "error: nend present but track_select_type is incompatible"
  bad_value = .true.    
elseif ( ((nrange(1) .ne. 0) .or. (nrange(2) .ne. 0)) .and. &
         (trim(track_select_type) .ne. "nrange") ) then
  print *, "error: nbegin/nend present but track_select_type is incompatible"
  bad_value = .true.        
endif
  
end subroutine check_track_parameters
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine check_spec_detector_parameters( ndim, ncells1, ncells2, min_dec, num_dec, &
             ppdec, wpoints, numprocs, diag_type, emissivity, waxistype, &
             detector_axis, x0, y0, z0, x1detmin, x1detmax, x2detmin, x2detmax, &
             wmin, wmax, bad_value )
! ----------------------------------------------------------------------------------------
integer, intent(in) :: ndim, ncells1, ncells2
integer, intent(in) :: min_dec, num_dec, ppdec, wpoints, numprocs
character(len=*), intent(in) :: diag_type, waxistype, detector_axis, emissivity
real(kind=p_double), intent(in) :: x0, y0, z0, x1detmin, x1detmax, x2detmin, x2detmax
real(kind=p_double), intent(in) :: wmin, wmax
logical, intent(inout) :: bad_value

! Check if spectrum detector parameters are ok

select case (diag_type)
  case ("standard","farfield","farfieldEndPoints","energy")
  case default
    print *, "error: diag_type ", diag_type, " not available"
    bad_value = .true.
end select

select case (emissivity)
  case ("d2W/dwdO","d2W/dwdS")
  case default
    print *, "error: emissivity can only be d2W/dwdO or d2W/dwdS"
    bad_value = .true.
end select

if ((ndim .lt. 1) .or. (ndim .gt. 3)) then
  print *, "error: ndim must be between 1 and 3"
  bad_value = .true.
endif

if ((ndim .gt. 1) .and. (ncells1 .lt. 2)) then
  print *, "error: ndim > 1 and ncells1 < 2"
  bad_value = .true.
elseif ((ndim .gt. 2) .and. ((ncells1 .lt. 2) .or. &
        (ncells2 .lt. 2))) then 
  print *, "error: ndim > 2 and ncells1 or ncells2 < 2"
  bad_value = .true.
elseif ((x1detmax .le. x1detmin) .and. (ndim .gt. 1)) then
  print *, "error: x1detmax must be > x1detmin"
  bad_value = .true.
elseif (x2detmax .le. x2detmin .and. (ndim .gt. 2)) then
  print *, "error: x2detmax must be > x2detmin"
  bad_value = .true.
elseif (mod(wpoints,numprocs) .ne. 0) then
  print *, "error: wpoints must be a multiple of the number of processors"
  print *, "(different sized w axis partitions not implemented)"
  bad_value = .true.
elseif ((trim(waxistype) .eq. 'linear') .or. &
        (trim(waxistype) .eq. 'log10')) then
  if (wpoints .lt. 2) then 
    print *, "error: wpoints must be > 2"
    bad_value = .true.
  else if (numprocs .gt. wpoints) then
    print *, "error: wpoints must be > number of processes"
    bad_value = .true.
  else if (modulo(wpoints,numprocs) .ne. 0) then
    print *, "error: wpoints must be a multiple of the number of processes"
    print *, "       uneven w axis partitions not implemented yet"
    bad_value = .true.  
  endif
  if (wmin .lt. 0.0d0) then
    print *, "error: wmin must be > 0"
    bad_value = .true.
  endif
  if ((wmax .lt. 0.0d0) .or. (wmax .le. wmin)) then
    print *, "error: wmax must be > 0 and > wmin"
    bad_value = .true.
  endif
elseif ((trim(waxistype) .eq. 'linbydec')) then
  if ((min_dec .lt. 0) .or. (num_dec .le. 1)) then
    print *, "error: min_dec must be > 0 and num_dec must be > 0"
    bad_value = .true.
  endif  
  if (ppdec .lt. 1) then
    print *, "error: ppdec must be > 0"
    bad_value = .true.
  endif
else
  print *, "error: waxis_type ", trim(waxistype), "not available"
  bad_value = .true.
endif

select case (trim(detector_axis))
  
  case ("x1")
    if (x0 .ne. -1.0d0) then
      print *, "error: x0 should not exist for x1 detector"
      bad_value = .true.
    elseif (x1detmin .eq. x1detmax) then
      print *, "error: x1 axis limits (x1detmin and x1detmax) not properly defined"
      bad_value = .true. 
    elseif (ndim .ne. 2) then
      print *, "error: Dimensions must be 1 for an x2 detector"
      bad_value = .true.   
    endif
                       
  case ("x2")
    if (y0 .ne. -1.0d0) then
      print *, "error: y0 should not exist for x2 detector"
      bad_value = .true.
    elseif (x1detmin .eq. x1detmax) then
      print *, "error: x2 axis limits (x1detmin and x1detmax) not properly defined"
      bad_value = .true. 
    elseif (ndim .ne. 2) then
      print *, "error: Dimensions must be 2 for an x2 detector"
      bad_value = .true.   
    endif
    
  case ("x3")
    if (z0 .ne. -1.0d0) then
      print *, "error: z0 should not exist for x3 detector"
      bad_value = .true.
    elseif (x1detmin .eq. x1detmax) then
      print *, "error: x3 axis limits (x1detmin and x1detmax) not properly defined"
      bad_value = .true. 
    elseif (ndim .ne. 2) then
      print *, "error: Dimensions must be 2 for an x3 detector"
      bad_value = .true.   
    endif
    
  case ("x1x2")
    if ((x0 .ne. -1.0d0) .or. (y0 .ne. -1.0d0)) then
      print *, "error: x0 and/or y0 should not exist for x1x2 detector"
      bad_value = .true.
    elseif ((x1detmin .eq. x1detmax) .or. (x2detmin .eq. x2detmax))  then
      print *, "error: x1 and/or x2 axis limits (x1detmin,x1detmax) or (x2detmin,x2detmax)"
      print *, "       not properly defined"
      bad_value = .true. 
    elseif (ndim .ne. 3) then
      print *, "error: Dimensions must be 3 for an x1x2 detector"
      bad_value = .true.   
    endif
    
  case ("x1x3")
    if ((x0 .ne. -1.0d0) .or. (z0 .ne. -1.0d0)) then
      print *, "error: x0 and/or z0 should not exist for x1x3 detector"
      bad_value = .true.
    elseif ((x1detmin .eq. x1detmax) .or. (x2detmin .eq. x2detmax))  then
      print *, "error: x1 and/or x3 axis limits (x1detmin,x1detmax) or (x2detmin,x2detmax)"
      print *, "       not properly defined"
      bad_value = .true. 
    elseif (ndim .ne. 3) then
      print *, "error: Dimensions must be 3 for an x1x3 detector"
      bad_value = .true.   
    endif
    
  case ("x2x3")
    if ((y0 .ne. -1.0d0) .or. (z0 .ne. -1.0d0)) then
      print *, "error: y0 and/or z0 should not exist for x2x3 detector"
      bad_value = .true.
    elseif ((x1detmin .eq. x1detmax) .or. (x2detmin .eq. x2detmax))  then
      print *, "error: x2 and/or x3 axis limits (x1detmin,x1detmax) or (x2detmin,x2detmax)"
      print *, "       not properly defined"
      bad_value = .true. 
    elseif (ndim .ne. 3) then
      print *, "error: Dimensions must be 3 for an x2x3 detector"
      bad_value = .true.   
    endif
  
  case default
  
    if (ndim .gt. 1) then
      print *, "error: detector_axis not selected or not available"
      bad_value = .true.
    endif  
    
end select

end subroutine check_spec_detector_parameters
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine check_ene_detector_parameters( ndim, ncells1, ncells2, numprocs, detector_axis, &
            x0, y0, z0, x1detmin, x1detmax, x2detmin, x2detmax, bad_value )
! ----------------------------------------------------------------------------------------
integer, intent(in) :: ndim, ncells1, ncells2
integer, intent(in) :: numprocs
character(len=*), intent(in) :: detector_axis
real(kind=p_double), intent(in) :: x0, y0, z0, x1detmin, x1detmax, x2detmin, x2detmax
logical, intent(inout) :: bad_value

! Check if spectrum detector parameters are ok
if (ndim .ne. 2) then
  print *, "error: ndim must be between 2 (energy diagnostic only implemented in 2D)"
  bad_value = .true.
endif

if ( (ncells1 .lt. 2) .or. (ncells2 .lt. 2) ) then
  print *, "error: ncells1 and/or ncells2 < 2"
  bad_value = .true.
endif

if ((x1detmax .le. x1detmin) .or. (x2detmax .le. x2detmin) ) then
  print *, "error: x1detmax(x2detmax) must be > x1detmin(x2detmin)"
  bad_value = .true.
endif


if (mod(ncells1,numprocs) .ne. 0) then
  print *, "error: ncells1 must be a multiple of the number of processors"
  print *, "(different sized x axis partitions not implemented)"
  bad_value = .true.
endif

select case (trim(detector_axis))
  
    
  case ("x1x2")
    if ((x0 .ne. -1.0d0) .or. (y0 .ne. -1.0d0)) then
      print *, "error: x0 and/or y0 should not exist for x1x2 detector"
      bad_value = .true.
    elseif ((x1detmin .eq. x1detmax) .or. (x2detmin .eq. x2detmax))  then
      print *, "error: x1 and/or x2 axis limits (x1detmin,x1detmax) or (x2detmin,x2detmax)"
      print *, "       not properly defined"
      bad_value = .true. 
    endif
    
  case ("x1x3")
    if ((x0 .ne. -1.0d0) .or. (z0 .ne. -1.0d0)) then
      print *, "error: x0 and/or z0 should not exist for x1x3 detector"
      bad_value = .true.
    elseif ((x1detmin .eq. x1detmax) .or. (x2detmin .eq. x2detmax))  then
      print *, "error: x1 and/or x3 axis limits (x1detmin,x1detmax) or (x2detmin,x2detmax)"
      print *, "       not properly defined"
      bad_value = .true. 
    endif
    
  case ("x2x3")
    if ((y0 .ne. -1.0d0) .or. (z0 .ne. -1.0d0)) then
      print *, "error: y0 and/or z0 should not exist for x2x3 detector"
      bad_value = .true.
    elseif ((x1detmin .eq. x1detmax) .or. (x2detmin .eq. x2detmax))  then
      print *, "error: x2 and/or x3 axis limits (x1detmin,x1detmax) or (x2detmin,x2detmax)"
      print *, "       not properly defined"
      bad_value = .true. 
    endif
  
  case default
  
    print *, "error: detector_axis value invalid"
    bad_value = .true.
      
    
end select



end subroutine check_ene_detector_parameters
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine make_waxis( input, detector, numprocs, comm ) 
! ----------------------------------------------------------------------------------------

type(t_input), intent(in) :: input
type(t_detector_spec), intent(inout) :: detector
integer, intent(in) :: numprocs, comm              
  
! local variables  
integer :: error, div, rest, myfirstw, mylastw, myid
integer :: k, m, wcounter, wpoints
real (kind=p_double) :: dec_max, dec_min
! variables from input (temporary)
character(len=10) :: waxistype
integer, dimension(:), allocatable :: firstwvec
real(p_double), dimension(:), allocatable :: waxis, waxislocal
integer :: mynw, nw, min_dec, num_dec, ppdec
real(p_double) :: dw, wmin, wmax

! Get values from input structure
waxistype = input%waxistype
wpoints = input%wpoints
min_dec = input%min_dec
num_dec = input%num_dec
ppdec = input%ppdec
wmin = input%wmin
wmax = input%wmax
 
! execution
call MPI_COMM_RANK( comm, myid, error )
if (error .ne. 0) print *, "Error getting myid"
  
select case (trim(waxistype))
  case ('linear')
    nw = wpoints
    dw = (wmax-wmin)/nw
  case ('linbydec')
    nw = num_dec*ppdec
    print *, "nw = ", nw
  case ('log10')
    nw = wpoints
    dw = (wmax-wmin)/(nw-1)
end select
  
allocate(firstwvec(numprocs))
allocate(detector%waxis%firstwvec(numprocs))
firstwvec = 0
allocate(waxis(nw))  
allocate(detector%waxis%waxis(nw))

div = nw/numprocs
rest = nw - div*numprocs

do k = 0, numprocs-1
  if ((k .lt. (rest-1)) .or. (k .eq. (rest-1))) then
    firstwvec(k+1) = (k*(div+1) + 1)-1
  else
    firstwvec(k+1) = (rest*(div+1) + (k-rest)*div+1)-1
  endif
enddo

if ((myid .lt. (rest-1)) .or. (myid .eq. (rest-1))) then
  mynw = div+1
  myfirstw = myid*(div+1) + 1
  mylastw = myid*(div+1) + mynw
else
  mynw = div
  myfirstw = rest*(div+1) + (myid-rest)*div+1
  mylastw = rest*(div+1) + (myid-rest)*div + mynw
endif
 
!------------------- waxis calculation -------------------------
select case (trim(waxistype))
  case ('linear')
    dw = (wmax-wmin)/nw
    do k=1,nw
      waxis(k) = wmin+(k-0.5d0)*dw  
    end do
  case ('linbydec')
    wcounter = 1
    do k = min_dec, min_dec+num_dec-1
      dec_min = 1.0d1**k
      dec_max = (1.0d1)**(k+1)
      dw = (dec_max - dec_min)/ppdec
      if (myid .eq. 0) print *, "k: dec_min / dec_max = ", &
                                k, dec_min, dec_max
      do m = 1, ppdec
        waxis(wcounter) = dec_min + (0.5d0+(m-1))*dw
        wcounter = wcounter + 1
      enddo  
    enddo
  case ('log10')
    do k=1,nw
      waxis(k) = log10(wmin+(k-1)*dw)  
    enddo  
end select

!------------------- local waxis calculation -------------------
allocate(waxislocal(mynw))
allocate(detector%waxis%waxislocal(mynw))
do k=1,mynw
  waxislocal(k) = waxis(firstwvec(myid+1)+k)
end do

! Get results into detector%waxis
detector%waxis%waxistype = waxistype
detector%waxis%firstwvec = firstwvec
detector%waxis%waxis = waxis
detector%waxis%waxislocal = waxislocal
detector%waxis%mynw = mynw
detector%waxis%nw = nw
detector%waxis%dw = dw
detector%waxis%wmin = wmin
detector%waxis%wmax = wmax
detector%waxis%dimsfwaxis = nw
detector%waxis%dims_chunkwaxis = mynw

deallocate(firstwvec,waxis,waxislocal)

end subroutine make_waxis
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine make_xparallelaxis( input, detector, numprocs, comm )                                
                               
! ----------------------------------------------------------------------------------------

type(t_input), intent(in) :: input
type(t_detector_ene), intent(inout) :: detector
integer, intent(in) :: numprocs, comm        
  
! local variables  
real(kind=p_double) :: x1detmin, x1detmax
integer :: nx
real(kind=p_double) :: dx
integer, dimension(:), allocatable :: firstvec
real(kind=p_double),dimension(:),allocatable :: xaxis
real(kind=p_double),dimension(:),allocatable :: xaxislocal
integer :: error, div, rest, myfirstx, mylastx, myid, mynx
integer :: k

x1detmin = input%x1detmin
x1detmax = input%x1detmax
nx = input%ncells1
  
! execution
call MPI_COMM_RANK( comm, myid, error )
if (error .ne. 0) print *, "Error getting myid"

dx = (x1detmax-x1detmin)/nx
  
allocate(firstvec(numprocs))
allocate(detector%xaxis%firstvec(numprocs))
firstvec = 0
allocate(xaxis(nx))  
allocate(detector%xaxis%xaxis(nx))  

do k=1,nx
  xaxis(k) = x1detmin+(k-0.5d0)*dx  
end do

div = nx/numprocs
rest = nx - div*numprocs

do k = 0, numprocs-1
  if ((k .lt. (rest-1)) .or. (k .eq. (rest-1))) then
    firstvec(k+1) = (k*(div+1) + 1)-1
  else
    firstvec(k+1) = (rest*(div+1) + (k-rest)*div+1)-1
  endif
enddo

if ((myid .lt. (rest-1)) .or. (myid .eq. (rest-1))) then
  mynx = div+1
  myfirstx = myid*(div+1) + 1
  mylastx = myid*(div+1) + mynx
else
  mynx = div
  myfirstx = rest*(div+1) + (myid-rest)*div+1
  mylastx  = rest*(div+1) + (myid-rest)*div + mynx
endif

!------------------- local waxis calculation -------------------
allocate(xaxislocal(mynx))
allocate(detector%xaxis%xaxislocal(mynx))
do k=1,mynx
  xaxislocal(k) = xaxis(firstvec(myid+1)+k)
end do


detector%xaxis%firstvec = firstvec
detector%xaxis%xaxis = xaxis
detector%xaxis%xaxislocal = xaxislocal
detector%xaxis%mynx = mynx
detector%xaxis%dx = dx
detector%xaxis%xmin = x1detmin
detector%xaxis%xmax = x1detmax
detector%xaxis%dimsf_xaxis = nx
detector%xaxis%dims_chunk_xaxis = mynx

deallocate(firstvec,xaxis,xaxislocal)

end subroutine make_xparallelaxis
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine save_spec( input, detector, diag_type, dtinitial, filename, parallelIO, comm ) 
! ----------------------------------------------------------------------------------------

type(t_input), intent(in) :: input
type(t_detector_spec), intent(in) :: detector
character(len=*), intent(in) :: diag_type
real(p_double), intent(in) :: dtinitial
character(len=*), intent(in) :: filename
logical, intent(in) :: parallelIO
integer, intent(in) :: comm
! local variables  
integer(hid_t) :: filesave_id, wgr_id
real(p_double) :: x1detmin, x2detmin, x1detmax, x2detmax, dx1det, dx2det
real(p_double) :: wmin, wmax, dw
integer :: nw
integer :: error, mpierr, myid
real(p_double) :: time1, time2
character(len=4) :: axisname1, axislname1, axisname2, axislname2, temp
character(len=28) :: data_name, axis_label, axis_units

x1detmin = detector%x1detmin
x2detmin = detector%x2detmin
x1detmax = detector%x1detmax
x2detmax = detector%x2detmax

dx1det = (x1detmax-x1detmin)/(detector%ncells1)
dx2det = (x2detmax-x2detmin)/(detector%ncells2)

nw = detector%waxis%nw
dw = detector%waxis%dw
wmin = detector%waxis%wmin
wmax = detector%waxis%wmax

data_name = "spectrum"

if (input%emissivity .eq. "d2W/dwdS") then
  axis_label = "d^2 W / d\omega dS"
  axis_units = "e^2 \omega_p^2 / (\pi^2 c^3)"
else
  axis_label = "d^2 W / d\omega d\Omega"
  axis_units = "e^2/ (\pi^2 c)" 
endif

call MPI_COMM_RANK( comm, myid, mpierr )

!################### write data to hdf5 file ###################
        
if (myid .eq. 0) print *, "Before file saving"
if (myid .eq. 0) print *, " " 

 
time1 = mpi_wtime()

!########################### file setup ########################

select case (detector%ndim)   
  case (1)
    ! --- case 1D ---
    if (myid .eq. 0) print *, "--- 1D file setup ----"
    call setup_spec_file1d(filename, filesave_id, myid, parallelIO, &
      real(wmin+0.5d0*dw), real(wmax-0.5d0*dw), '\omega_p', '\omega', '\omega', &
      data_name, 0.0d0, 0, dtinitial, (/wmin/), (/wmax/), diag_type)
    ! ---------------
  case (2)
    ! --- case 2D ---
    if (myid .eq. 0) print *, "--- 2D file setup ----"
    
    axisname1 = trim(detector%detector_axis)
    temp = trim(detector%detector_axis)
    axislname1 = temp(1:1) // "_" // temp(2:2)
    
    call setup_spec_file2d(filename, filesave_id, myid, parallelIO, &
      real(wmin), real(wmax), '\omega_p', '\omega', '\omega', &
      real(x1detmin), real(x1detmax), &
      'c/\omega_p',axisname1,axislname1, data_name, 0.0d0, 0, dtinitial, &
      (/wmin, x1detmin/), (/wmax, x1detmax/),diag_type)
    ! ---------------
    
  case (3)
    ! --- case 3D ---
    if (myid .eq. 0) print *, "--- 3D file setup ----"
    
    temp = trim(detector%detector_axis)
    axisname1 = temp(1:2)
    temp = trim(detector%detector_axis)
    axislname1 = temp(1:1) // "_" // temp(2:2)
    
    axisname2 = temp(3:4)
    axislname2 = temp(3:3) // "_" // temp(4:4)
    
    call setup_spec_file3d(filename, filesave_id, myid, parallelIO, &
      real(wmin), real(wmax), '\omega_p', '\omega', '\omega', &
      real(x1detmin), real(x1detmax), &
      'c/\omega_p',axisname1,axislname1, &
      real(x2detmin), real(x2detmax), &
      'c/\omega_p',axisname2,axislname2, &
      data_name, 0.0d0, 0, dtinitial, &
      (/wmin,x1detmin,x2detmin/), &
      (/wmax,x1detmax,x2detmax/),diag_type)
    ! ---------------
end select

!###############################################################


!########################### adding dataset ########################

select case (detector%ndim)   
  case (1)
    if (myid .eq. 0) print *, "--- 1D ----"
    ! --- case 1D ---
    if (parallelIO .neqv. .true.) then
      if (myid .eq. 0) print *, "---- serial I/O ---- "
      call adddataset1d( filesave_id, data_name, detector%spec1d, &
                         detector%dimsf, detector%dims_chunk, &
                         detector%waxis%firstwvec, comm, &
                         axis_units, axis_label)
      call adddataset1d(filesave_id, "waxis", detector%waxis%waxislocal, &
                     detector%waxis%dimsfwaxis, detector%waxis%dims_chunkwaxis, &
                     detector%waxis%firstwvec, comm, '\omega_p', '\omega', &
                     detector%waxis%waxistype, '\omega')
    else
      if (myid .eq. 0) print *, "---- parallel I/O ----"
      
      call adddatasetparallel1d(filesave_id, comm, data_name, &
                                detector%dimsf, detector%dims_chunk, &
                                detector%waxis%firstwvec, detector%spec1d, &
                                axis_units, axis_label)
      
      call h5gopen_f( filesave_id, 'WAXIS', wgr_id, error) 
      if (error .ne. 0) print *, "myid, error opening waxis group = ", myid,error
      
      call adddatasetparallel1d(wgr_id, comm, "waxis", detector%waxis%dimsfwaxis, &
                        detector%waxis%dims_chunkwaxis, detector%waxis%firstwvec, &
                        detector%waxis%waxislocal, '\omega_p', '\omega', &
                        detector%waxis%waxistype, '\omega')
      
      call h5gclose_f( wgr_id, error)
      if (error .ne. 0) print *, "myid, error closing waxis group = ", myid,error
 
    endif
    ! ---------------
  case (2)
    if (myid .eq. 0) print *, "--- 2D ----"
    ! --- case 2D ---
    if (parallelIO .neqv. .true.) then
      if (myid .eq. 0) print *, "---- serial I/O ---- "
      call adddataset2d( filesave_id, detector%spec2d, detector%dimsf, &
                         detector%dims_chunk, detector%waxis%firstwvec, &
                         data_name, axis_label, axis_units, comm )
      call adddataset1d(filesave_id, "waxis", detector%waxis%waxislocal, &
                     detector%waxis%dimsfwaxis, detector%waxis%dims_chunkwaxis, &
                     detector%waxis%firstwvec, comm, '\omega_p', '\omega', &
                     detector%waxis%waxistype, '\omega')
    else
      if (myid .eq. 0) print *, "---- parallel I/O ----"
      call adddatasetparallel2d(filesave_id, comm, data_name, &
                        axis_label, axis_units, &
                        detector%dimsf, detector%dims_chunk, &
                        detector%waxis%firstwvec, detector%spec2d)
      
      call h5gopen_f( filesave_id, 'WAXIS', wgr_id, error) 
      if (error .ne. 0) print *, "myid, error opening waxis group = ", myid,error
      
      call adddatasetparallel1d(wgr_id, comm, "waxis", detector%waxis%dimsfwaxis, &
                     detector%waxis%dims_chunkwaxis, detector%waxis%firstwvec,&
                     detector%waxis%waxislocal, '\omega_p', '\omega', &
                     detector%waxis%waxistype, '\omega')
      
      call h5gclose_f( wgr_id, error)
      if (error .ne. 0) print *, "myid, error closing waxis group = ", myid,error
 
    endif
    ! ---------------
  case (3)
    if (myid .eq. 0) print *, "--- 3D ----"
    ! --- case 3D ---
    if (parallelIO .neqv. .true.) then
      if (myid .eq. 0) print *, "---- serial I/O ---- "
      call adddataset3d( filesave_id, detector%spec3d, detector%dimsf, &
                         detector%dims_chunk, detector%waxis%firstwvec, data_name, &
                         axis_label, axis_units, comm )
      
      call adddataset1d(filesave_id, "waxis", detector%waxis%waxislocal, &
                     detector%waxis%dimsfwaxis, detector%waxis%dims_chunkwaxis, &
                     detector%waxis%firstwvec, comm, '\omega_p', '\omega', &
                     detector%waxis%waxistype, '\omega')
    else
      if (myid .eq. 0) print *, "---- parallel I/O ----"
      call adddatasetparallel3d(filesave_id, comm, data_name, &
                        axis_label, axis_units, &
                        detector%dimsf, detector%dims_chunk, &
                        detector%waxis%firstwvec, detector%spec3d)
      
      call h5gopen_f( filesave_id, 'WAXIS', wgr_id, error) 
      
      call adddatasetparallel1d(wgr_id, comm, "waxis", detector%waxis%dimsfwaxis, &
                     detector%waxis%dims_chunkwaxis, detector%waxis%firstwvec, &
                     detector%waxis%waxislocal, '\omega_p', '\omega', &
                     detector%waxis%waxistype, '\omega')
      
      call h5gclose_f( wgr_id, error)
    endif
    ! ---------------
end select

!###############################################################
 
  
call closefile(filesave_id, myid, parallelIO)

time2 = mpi_wtime()

if (myid .eq. 0) then
  print *, " "
  print *, "time saving file = ", time2-time1, " s"
endif  

end subroutine save_spec
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine save_ene( detector, dtinitial, filename, parallelIO, comm ) 
! ----------------------------------------------------------------------------------------

type(t_detector_ene), intent(in) :: detector
real(p_double), intent(in) :: dtinitial
character(len=*), intent(in) :: filename
logical, intent(in) :: parallelIO
integer, intent(in) :: comm
! local variables  
integer(hid_t) :: filesave_id
real(p_double) :: x1detmin, x2detmin, x1detmax, x2detmax, dx1det, dx2det
integer :: mpierr, myid
real(p_double) :: time1, time2
character(len=4) :: detector_axis
character(len=3) :: axis_name1, axis_lname1, axis_name2, axis_lname2

x1detmin = detector%x1detmin
x2detmin = detector%x2detmin
x1detmax = detector%x1detmax
x2detmax = detector%x2detmax

detector_axis = trim(detector%detector_axis)

dx1det = (x1detmax-x1detmin)/(detector%ncells1)
dx2det = (x2detmax-x2detmin)/(detector%ncells2)

select case (trim(detector_axis))
  case ("x1x2")
    axis_name1 = "x_1"
    axis_lname1 = "x_1"
    axis_name2 = "x_2"
    axis_lname2 = "x_2"
  case ("x1x3")
    axis_name1 = "x_1"
    axis_lname1 = "x_1"
    axis_name2 = "x_3"
    axis_lname2 = "x_3"
  case ("x2x3")
    axis_name1 = "x_2"
    axis_lname1 = "x_2"
    axis_name2 = "x_3"
    axis_lname2 = "x_3"
end select

call MPI_COMM_RANK( comm, myid, mpierr )

!################### write data to hdf5 file ###################
        
 
time1 = mpi_wtime()

!########################### file setup ########################

select case (detector%ndim)   
  
  case (2)
    ! --- case 2D ---
    if (myid .eq. 0) print *, "--- 2D file setup ----"
    call setup_ene_file( filename, filesave_id, myid, parallelIO, &
        real(x1detmin), real(x1detmax), 'c/\omega_p', &
        axis_name1, axis_lname1, real(x2detmin), real(x2detmax), &
        'c/\omega_p', axis_name2, axis_lname2, "Energy", 0.0d0, 0, dtinitial)  
    ! ---------------
  case default
    print *, "Only 2D energy detectors are implemented"
    ! ---------------
end select

!###############################################################


!########################### adding dataset ########################

select case (detector%ndim)   
  
  case (2)
    if (myid .eq. 0) print *, "--- 2D ----"
    ! --- case 2D ---
    if (parallelIO .neqv. .true.) then
      if (myid .eq. 0) print *, "---- serial I/O ---- "
      call adddataset2d( filesave_id, detector%pow2d, detector%dimsf, &
                         detector%dims_chunk, detector%xaxis%firstvec, &
                         "Energy", "Energy", "e^2 \omega_p / (4\pi c)", comm )
    else
      if (myid .eq. 0) print *, "---- parallel I/O ----"
      call adddatasetparallel2d(filesave_id, comm, "Energy", &
                        "Energy", "e^2 \omega_p / (4\pi c)", &
                        detector%dimsf, detector%dims_chunk, &
                        detector%xaxis%firstvec, detector%pow2d)
    endif
    ! ---------------
  case default
  
  
end select

!###############################################################
 
  
call closefile(filesave_id, myid, parallelIO)

time2 = mpi_wtime()

end subroutine save_ene
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


end module utilities
