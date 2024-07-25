! ****************************************************************************
! Copyright © 2014 Instituto Superior Técnico. All rights reserved. This 
! software is the copyrighted work of Instituto Superior Técnico. 
! Reproduction, in whole or in part, on the Internet, on CD-ROM or any 
! other medium, without the prior written consent of Instituto Superior 
! Técnico is prohibited.
! ****************************************************************************

module track_utilities

  use hdf5_utilities
  use type_definitions
  
implicit none
 
private
    
! ##########################################################################

  public :: read_track, track_selection, cleanup_track, calc_beta, &
            calc_betaDot, setup_track
     
! ##########################################################################

contains


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine setup_track( filename,track,ndimtrack,invalidTrack,comm)
! ----------------------------------------------------------------------------------------

character(len=*), intent(in) :: filename
type(t_track), intent(inout) :: track
integer, intent(in) :: comm
integer, intent(in) :: ndimtrack
logical, intent(inout) :: invalidTrack

! local variables
integer :: mpierr, myid, error 
integer(hid_t) :: file_id, track_gid
integer :: obj_type, tracksize
character(len=20) :: track_gname ! name of track object
logical :: link_exists
  
! get myid
call MPI_COMM_RANK( comm, myid, mpierr )
if (mpierr .ne. 0) print *, "Error getting myid"

! This routine will check dimension of track (1D, 2D, 3D) and if this
! is ndimtrack <= than the track dimension

track%ndimtrack = ndimtrack

invalidTrack = .false.

  
! read file if I'm process 0
if (myid .eq. 0) then
  
  ! Open the hdf track file
  call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
    
  !----------- getting quantities from tracks -------------------
      
  ! Get track group name
  call h5gget_obj_info_idx_f(file_id, "/", 0, track_gname, obj_type, &
                             error)
    
  ! Open track group
  call h5gopen_f(file_id, track_gname, track_gid, error)

  ! Check if x2 exists
  call h5lexists_f(track_gid, "x2", link_exists, error)
  if (link_exists) then
    select case (ndimtrack) 
      case (1)
        print *, "* warning * Using 2D track as a 1D track"
        print *, "(ignoring x2)"
      case default
        print *, "[track has x2]"  
    end select 
  else if (ndimtrack .gt. 1) then
    print *, "* Error * Trying to access 2D or 3D track but found no x2"
    invalidTrack = .true.
  endif
  ! Check if x3 exists
  call h5lexists_f(track_gid, "x3", link_exists, error)
  if (link_exists) then
    select case (ndimtrack) 
      case (1)
        print *, "* warning * Using 3D track as a 1D track"
        print *, "(ignoring x2, x3)"
      case (2)
        print *, "* warning * Using 3D track as a 2D track"
        print *, "(ignoring x3)"
      case default
        print *, "[track has x3]"  
    end select
  else if (ndimtrack .gt. 2) then
    print *, "* Error * Trying to access 3D track but found no x3"  
    invalidTrack = .true.
  endif
   

  call h5gclose_f(track_gid, error)  ! closing track group
    
endif ! end of file read if I'm process 0
  
! If I'm process 0, broadcast tracksize to remaining processes
! Broadcasting already synchronizes processes
call MPI_BCAST(tracksize,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting tracksize"
      
end subroutine setup_track
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine read_track( filename,track,p,comm)
! ----------------------------------------------------------------------------------------

character(len=*), intent(in) :: filename
type(t_track), intent(inout) :: track
integer, intent(in) :: p, comm
! local variables
integer :: tracksize
integer :: mpierr, myid, error 
integer(hid_t) :: file_id, dset_id, dtype_id, space_id, track_gid
integer(hssize_t) :: tracksizetemp
integer :: obj_type
character(len=20) :: track_gname ! name of track object
  
! get myid
call MPI_COMM_RANK( comm, myid, mpierr )
if (mpierr .ne. 0) print *, "Error getting myid"

  
! read file if I'm process 0
if (myid .eq. 0) then
  
  ! Open the hdf track file
  call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
  if (error .ne. 0) print *, "Error opening trackfile, myid = ", myid
    
    
  !----------- getting quantities from tracks -------------------
      
  ! Get track group name
  call h5gget_obj_info_idx_f(file_id, "/", p-1, track_gname, obj_type, &
                             error)
  if (error .ne. 0) print *, "Error getting track name, myid = ", myid                           
    
  ! Open track group
  call h5gopen_f(file_id, track_gname, track_gid, error)
  if (error .ne. 0) print *, "Error opening track group, myid = ", myid
   
  ! time dataset
  call h5dopen_f(track_gid, "t", dset_id, error)
  call h5dget_type_f(dset_id, dtype_id, error)
  call h5dget_space_f(dset_id, space_id, error)  
  call h5sget_select_npoints_f(space_id, tracksizetemp, error)
     
  ! Allocate space for t, x1, x2, x3, p1, p2, p3 vectors 
  tracksize = int(tracksizetemp)   
  
  allocate(track%n(tracksize))    
  allocate(track%t(tracksize))
  allocate(track%x1(tracksize))
  if (track%ndimtrack .gt. 1) allocate(track%x2(tracksize))
  if (track%ndimtrack .gt. 2) allocate(track%x3(tracksize))
  allocate(track%p1(tracksize))
  allocate(track%p2(tracksize))
  allocate(track%p3(tracksize))
  allocate(track%ene(tracksize))
  allocate(track%chargetemp(tracksize))

  ! time dataset
  call h5dread_f(dset_id, H5T_NATIVE_double, track%t, &
                 (/tracksizetemp/), error)
  call h5dclose_f(dset_id, error)
  ! n dataset
  call h5dopen_f(track_gid, "n", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, track%n, &
                 (/tracksizetemp/), error)
  call h5dclose_f(dset_id, error)
  ! x1 dataset 
  call h5dopen_f(track_gid, "x1", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_double, track%x1, &
                 (/tracksizetemp/), error)
  call h5dclose_f(dset_id, error)
  if (track%ndimtrack .gt. 1) then
    ! x2 dataset
    call h5dopen_f(track_gid, "x2", dset_id, error) 
    call h5dread_f(dset_id, H5T_NATIVE_double, track%x2, &
                   (/tracksizetemp/), error)
    call h5dclose_f(dset_id, error)
  endif
  if (track%ndimtrack .gt. 2) then   
    ! x3 dataset
    call h5dopen_f(track_gid, "x3", dset_id, error) 
    call h5dread_f(dset_id, H5T_NATIVE_double, track%x3, &
                   (/tracksizetemp/), error)
    call h5dclose_f(dset_id, error)
  endif  
  ! p1 dataset
  call h5dopen_f(track_gid, "p1", dset_id, error) 
  call h5dread_f(dset_id, H5T_NATIVE_double, track%p1, &
                 (/tracksizetemp/), error)
  call h5dclose_f(dset_id, error)
  ! p2 dataset
  call h5dopen_f(track_gid, "p2", dset_id, error) 
  call h5dread_f(dset_id, H5T_NATIVE_double, track%p2, &
                 (/tracksizetemp/), error)
  call h5dclose_f(dset_id, error)
  ! p3 dataset
  call h5dopen_f(track_gid, "p3", dset_id, error) 
  call h5dread_f(dset_id, H5T_NATIVE_double, track%p3, &
                 (/tracksizetemp/), error)
  call h5dclose_f(dset_id, error)
  ! ene dataset
  call h5dopen_f(track_gid, "ene", dset_id, error) 
  call h5dread_f(dset_id, H5T_NATIVE_double, track%ene, &
                 (/tracksizetemp/), error)
  call h5dclose_f(dset_id, error)
  ! charge dataset
  call h5dopen_f(track_gid, "q", dset_id, error) 
  call h5dread_f(dset_id, H5T_NATIVE_double, track%chargetemp, &
                 (/tracksizetemp/), error)
  call h5dclose_f(dset_id, error)
  track%charge = track%chargetemp(1)
  deallocate(track%chargetemp)
  call h5gclose_f(track_gid, error)      ! closing track group
    
endif ! end of file read if I'm process 0
  
! If I'm process 0, broadcast tracksize to remaining processes
call MPI_BCAST(tracksize,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting tracksize"

  
! Remaining processes should allocate memory
if (myid .ne. 0) then
  allocate(track%n(tracksize))
  allocate(track%t(tracksize))
  allocate(track%x1(tracksize))
  if (track%ndimtrack .gt. 1) then 
    allocate(track%x2(tracksize))
    track%x2 = 0.0d0
  endif  
  if (track%ndimtrack .gt. 2) then
    allocate(track%x3(tracksize))
    track%x3 = 0.0d0
  endif  
  allocate(track%p1(tracksize),track%p2(tracksize),track%p3(tracksize))
  allocate(track%ene(tracksize))
  
  track%n = 0 
  track%t = 0.0d0 
  track%x1 = 0.0d0
  track%p1 = 0.0d0
  track%p2 = 0.0d0
  track%p3 = 0.0d0
  track%ene = 0.0d0
endif


! If I'm process 0, broadcast track data to remaining processes
call MPI_BCAST(track%n,tracksize,MPI_INTEGER,0, &
               MPI_COMM_WORLD,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting n"
call MPI_BCAST(track%t,tracksize,MPI_DOUBLE_PRECISION,0, &
               MPI_COMM_WORLD,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting t"
call MPI_BCAST(track%x1,tracksize,MPI_DOUBLE_PRECISION,0, &
               MPI_COMM_WORLD,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting x1"
if (track%ndimtrack .gt. 1) then 
  call MPI_BCAST(track%x2,tracksize,MPI_DOUBLE_PRECISION,0, &
                 MPI_COMM_WORLD,mpierr)
  if (mpierr .ne. 0) print *, "Error broadcasting x2"
endif
if (track%ndimtrack .gt. 2) then   
  call MPI_BCAST(track%x3,tracksize,MPI_DOUBLE_PRECISION,0, &
                 MPI_COMM_WORLD,mpierr)
  if (mpierr .ne. 0) print *, "Error broadcasting x3"
endif  
call MPI_BCAST(track%p1,tracksize,MPI_DOUBLE_PRECISION,0, &
               MPI_COMM_WORLD,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting p1"
call MPI_BCAST(track%p2,tracksize,MPI_DOUBLE_PRECISION,0, &
               MPI_COMM_WORLD,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting p2"
call MPI_BCAST(track%p3,tracksize,MPI_DOUBLE_PRECISION,0, &
               MPI_COMM_WORLD,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting p3"
call MPI_BCAST(track%ene,tracksize,MPI_DOUBLE_PRECISION,0, &
               MPI_COMM_WORLD,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting ene"
call MPI_BCAST(track%charge,1,MPI_DOUBLE_PRECISION,0, &
               MPI_COMM_WORLD,mpierr)
if (mpierr .ne. 0) print *, "Error broadcasting charge"


track%tracksize = tracksize
track%dtinitial = track%t(2)-track%t(1)

allocate(track%dt(tracksize))
track%dt(1:tracksize-1) = track%t(2:tracksize)-track%t(1:tracksize-1)
track%dt(tracksize) = track%t(tracksize)-track%t(tracksize-1)

!---------------------------------------------------------------
  
end subroutine read_track
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine track_selection( track, track_select_type, nbegin, nend, nrange, &
                            x1min, x1max, x1range, tmin, tmax, trange, enemin, comm )
! ----------------------------------------------------------------------------------------

type(t_track), intent(inout) :: track
character(len=15), intent(in) :: track_select_type
integer, intent(in) :: nbegin, nend
real (kind=p_double), intent(in) :: x1min, x1max, tmin, tmax, enemin
real (kind=p_double), dimension(2), intent(in) :: x1range, trange
integer, dimension(2), intent(in) :: nrange
integer, intent(in) :: comm

! local variables
real (kind=p_double), dimension(:), allocatable :: temp 
integer :: sizevec, beginindex, endindex, newtracksize
integer :: ind, track_1st_n
integer :: myid, mpierr


call MPI_COMM_RANK(comm, myid, mpierr)

! determine indexes from data read
track_1st_n = track%n(1)

select case (trim(track_select_type))

  case ('nmin')
    beginindex = maxval((/nbegin - track_1st_n + 1,1/))
    endindex = track%tracksize
    if (myid .eq. 0) then
      print *, " "
      print *, "processing iterations: ", nbegin, " to ", track%tracksize
      print *, " "
    endif  
             
  case ('nmax')
    beginindex = track_1st_n
    endindex = minval((/nend - track_1st_n + 1,track%tracksize/))
    if (myid .eq. 0) then
      print *, " "
      print *, "processing iterations: ", track_1st_n, " to ", nend
      print *, " "
    endif  
             
  case ('nrange')
    beginindex = maxval((/nrange(1) - track_1st_n + 1,1/))
    endindex = minval((/nrange(2) - track_1st_n + 1,track%tracksize/))
    if (myid .eq. 0) then
      print *, " "
      print *, "processing iterations: ", nrange(1), " to ", nrange(2)
      print *, " "
    endif  
    
  case ('x1min')
    endindex = track%tracksize
    beginindex = 1
    do ind = 1, track%tracksize
      if (track%x1(ind) .ge. x1min) then
        beginindex = ind
        exit
      endif
    enddo
    if (myid .eq. 0) then
      print *, " "
      print *, "processing track from x1min = ", x1min
      print *, " "
    endif  
  
  case ('x1max')
    endindex = track%tracksize
    beginindex = 1
    do ind = 1, track%tracksize
      if (track%x1(ind) .gt. x1max) then
        endindex = ind-1
        exit
      endif
    enddo
    if (myid .eq. 0) then
      print *, " "
      print *, "processing track in for x1 <", x1max
      print *, " "
    endif  
      
  case ('x1range')
    endindex = track%tracksize
    beginindex = 1
    do ind = 1, track%tracksize
      if (track%x1(ind) .ge. x1range(1)) then
        beginindex = ind
        exit
      endif
    enddo
    do ind = 1, track%tracksize
      if (track%x1(ind) .gt. x1range(2)) then
        endindex = ind-1
        exit
      endif
    enddo
    if (myid .eq. 0) then
      print *, " "
      print *, "processing track in the range x1 ", &
             x1range(1), " - ", x1range(2)
      print *, " "
    endif         
    
  case ('tmin')
    endindex = track%tracksize
    beginindex = 1
    do ind = 1, track%tracksize
      if (track%t(ind) .ge. tmin) then
        beginindex = ind
        exit
      endif
    enddo
    if (myid .eq. 0) then
      print *, " "
      print *, "processing track from tmin = ", tmin
      print *, " "
    endif  
  
  case ('tmax')
    endindex = track%tracksize
    beginindex = 1
    do ind = 1, track%tracksize
      if (track%t(ind) .gt. tmax) then
        endindex = ind-1
        exit
      endif
    enddo
    if (myid .eq. 0) then
      print *, " "
      print *, "processing track for t <", trange(2)
      print *, " "
    endif  
               
  case ('trange')
    endindex = track%tracksize
    beginindex = 1
    do ind = 1, track%tracksize
      if (track%t(ind) .ge. trange(1)) then
        beginindex = ind
        exit
      endif
    enddo
    do ind = 1, track%tracksize
      if (track%t(ind) .gt. trange(2)) then
        endindex = ind-1
        exit
      endif
    enddo
    if (myid .eq. 0) then
      print *, " "
      print *, "processing track in the range t ", &
             trange(1), " - ", trange(2)
      print *, " "
    endif         
    
  case ('enemin')
    endindex = track%tracksize
    beginindex = 1
    do ind = 1, track%tracksize
      if (track%ene(ind) .ge. enemin) then
        beginindex = ind
        exit
      endif
    enddo
    if (myid .eq. 0) then
      print *, " "
      print *, "processing track from enemin = ", enemin
      print *, " "
    endif  
  
  case ('none')
    beginindex = 1
    endindex = track%tracksize
  
  case default
    print *, "Error: track_select_type ", trim(track_select_type), &
             "does not exist"
    print *, "Exiting program..."
    stop
      
end select

newtracksize = endindex-beginindex+1

! Make sure selection makes sense
if (newtracksize .le. 2) then
  print *, "Error: newtracksize < 2 !!"
  print *, "Exiting program..."
  stop
endif
    
! allocate vector to hold data temporarily
sizevec = track%tracksize
allocate(temp(sizevec))

! ---------------------------------------------
temp = track%t
deallocate(track%t)
allocate(track%t(newtracksize))
track%t = temp(beginindex:endindex)

temp = track%dt
deallocate(track%dt)
allocate(track%dt(newtracksize))
track%dt = temp(beginindex:endindex)
      
temp = track%x1
deallocate(track%x1)
allocate(track%x1(newtracksize))
track%x1 = temp(beginindex:endindex)  

if (track%ndimtrack .gt. 1) then
  temp = track%x2
  deallocate(track%x2)
  allocate(track%x2(newtracksize))
  track%x2 = temp(beginindex:endindex)  
endif  

if (track%ndimtrack .gt. 2) then   
  temp = track%x3
  deallocate(track%x3)
  allocate(track%x3(newtracksize))
  track%x3 = temp(beginindex:endindex)
endif  

temp = track%p1
deallocate(track%p1)
allocate(track%p1(newtracksize))
track%p1 = temp(beginindex:endindex)
	
temp = track%p2
deallocate(track%p2)
allocate(track%p2(newtracksize))
track%p2 = temp(beginindex:endindex)
	
temp = track%p3
deallocate(track%p3)
allocate(track%p3(newtracksize))
track%p3 = temp(beginindex:endindex)
! ---------------------------------------------

track%tracksize = newtracksize

! deallocate g and sqrt1MinusBetaSquared vectors
deallocate(temp)
    
end subroutine track_selection
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine calc_beta( track )
! ----------------------------------------------------------------------------------------

type(t_track), intent(inout) :: track
! local variables
real (kind=p_double), dimension(:), allocatable :: sqrt1MinusBetaSquared 
integer :: sizevec
  
sizevec = track%tracksize
  
! allocate g and sqrt1MinusBetaSquared vectors
allocate(track%g(sizevec),sqrt1MinusBetaSquared(sizevec))
track%g = sqrt(1.0d0 + (track%p1**2 + track%p2**2 + track%p3**2))
sqrt1MinusBetaSquared = 1.0d0/track%g

! Calculate beta vector components 
allocate(track%beta1(sizevec))
allocate(track%beta2(sizevec))
allocate(track%beta3(sizevec))
track%beta1 = sqrt1MinusBetaSquared*track%p1
track%beta2 = sqrt1MinusBetaSquared*track%p2
track%beta3 = sqrt1MinusBetaSquared*track%p3
  
! deallocate and sqrt1MinusBetaSquared vectors
deallocate(sqrt1MinusBetaSquared)
    
end subroutine calc_beta
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine calc_betaDot( track )
! ----------------------------------------------------------------------------------------

type(t_track), intent(inout) :: track
! local variables
real(kind=p_double), dimension(:), allocatable :: g, sqrt1MinusBetaSquared
real(kind=p_double), dimension(:), allocatable :: oneMinusBeta2Pow3halfs 
real(kind=p_double), dimension(:), allocatable :: temp1, temp2, &
                                                  p1dot,p2dot,p3dot
integer :: tracksize
  
tracksize = track%tracksize
  
! allocate g and sqrt1MinusBetaSquared vectors
allocate(g(tracksize),sqrt1MinusBetaSquared(tracksize),&
         oneMinusBeta2Pow3halfs(tracksize))
allocate(temp1(tracksize), temp2(tracksize))
allocate(p1dot(tracksize),p2dot(tracksize),p3dot(tracksize))
  
! Calculate betaDot vector components 
allocate(track%betaDot1(tracksize),track%betaDot2(tracksize),&
         track%betaDot3(tracksize))
  
g = sqrt(1.0d0 + (track%p1**2 + track%p2**2 + track%p3**2))
sqrt1MinusBetaSquared = 1.0d0/g
oneMinusBeta2Pow3halfs = 1.0d0/(g**3)

! Calculate pdot
temp1 = cshift(track%p1,1)
temp2 = cshift(track%p1,-1)
p1dot(2:tracksize-1) = (temp1(2:tracksize-1)-temp2(2:tracksize-1) )/ &
                       (2.0*track%dt(2:tracksize-1))
p1dot(1) = (track%p1(2)-track%p1(1))/track%dt(1)
p1dot(tracksize) = (track%p1(tracksize)-track%p1(tracksize-1))/ &
                    track%dt(tracksize)
		
temp1 = cshift(track%p2,1)
temp2 = cshift(track%p2,-1)
p2dot(2:tracksize-1) = (temp1(2:tracksize-1)-temp2(2:tracksize-1) )/ &
                       (2.0*track%dt(2:tracksize-1))
p2dot(1) = (track%p2(2)-track%p2(1))/track%dt(1)
p2dot(tracksize) = (track%p2(tracksize)-track%p2(tracksize-1))/ &
                    track%dt(tracksize)
		
temp1 = cshift(track%p3,1)
temp2 = cshift(track%p3,-1)
p3dot(2:tracksize-1) = (temp1(2:tracksize-1)-temp2(2:tracksize-1) )/ &
                       (2.0*track%dt(2:tracksize-1))
p3dot(1) = (track%p3(2)-track%p3(1))/track%dt(1)
p3dot(tracksize) = (track%p3(tracksize)-track%p3(tracksize-1))/ &
                    track%dt(tracksize) 

! Calculate beta time derivative vector components
track%betaDot1 = sqrt1MinusBetaSquared*p1dot - oneMinusBeta2Pow3halfs * &
        ( track%p1*p1dot + track%p2*p2dot + track%p3*p3dot )*track%p1
track%betaDot2 = sqrt1MinusBetaSquared*p2dot - oneMinusBeta2Pow3halfs * &
        ( track%p1*p1dot + track%p2*p2dot + track%p3*p3dot )*track%p2
track%betaDot3 = sqrt1MinusBetaSquared*p3dot - oneMinusBeta2Pow3halfs * &
        ( track%p1*p1dot + track%p2*p2dot + track%p3*p3dot )*track%p3
  
! deallocate g and sqrt1MinusBetaSquared vectors
deallocate(g,sqrt1MinusBetaSquared,oneMinusBeta2Pow3halfs,temp1,temp2)
deallocate(p1dot,p2dot,p3dot)
    
end subroutine calc_betaDot
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine cleanup_track( track )
! ----------------------------------------------------------------------------------------

type(t_track), intent(inout) :: track
  
! deallocate all pointers allocated
if (allocated(track%n)) deallocate(track%n)
if (allocated(track%t)) deallocate(track%t)
if (allocated(track%x1)) deallocate(track%x1)
if (allocated(track%x2)) deallocate(track%x2)
if (allocated(track%x3)) deallocate(track%x3)
if (allocated(track%p1)) deallocate(track%p1)
if (allocated(track%p2)) deallocate(track%p2)
if (allocated(track%p3)) deallocate(track%p3)
if (allocated(track%ene)) deallocate(track%ene)
if (allocated(track%dt)) deallocate(track%dt)
if (allocated(track%chargetemp)) deallocate(track%chargetemp) 
  
if (allocated(track%beta1)) deallocate(track%beta1) 
if (allocated(track%beta2)) deallocate(track%beta2) 
if (allocated(track%beta3)) deallocate(track%beta3) 
    
if (allocated(track%g)) deallocate(track%g)    
    
if (allocated(track%betaDot1)) deallocate(track%betaDot1) 
if (allocated(track%betaDot2)) deallocate(track%betaDot2) 
if (allocated(track%betaDot3)) deallocate(track%betaDot3)  
    
end subroutine cleanup_track
! ----------------------------------------------------------------------------------------
! ****************************************************************************************




end module track_utilities
