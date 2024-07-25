! ****************************************************************************
! Copyright © 2014 Instituto Superior Técnico. All rights reserved. This 
! software is the copyrighted work of Instituto Superior Técnico. 
! Reproduction, in whole or in part, on the Internet, on CD-ROM or any 
! other medium, without the prior written consent of Instituto Superior 
! Técnico is prohibited.
! ****************************************************************************

program jrad

! Program that reads a track file and determines the power and/or
! emissivity of the radiation it emits
  
! #################################################################
  
use hdf5
use utilities
use hdf5_utilities
use track_utilities
use type_definitions
use rad_calculations
use spec_calculations
use ene_calculations
    
implicit none
  
! variables
integer :: error
integer(hid_t) :: file_id       ! File identifier 
integer :: num_tracks           ! Number of objects in file

! variables to count time taken by the program
real (kind=p_double) :: starttime, endtime, &
                        timea1, timea2, time1st, time2nd, time3rd
                        

! variables for parallel calculations
integer :: mpierr, myid, numprocs

! Track variables
type(t_track) :: track

! Detector variables
type(t_detector_spec) :: detector_spec
type(t_detector_ene) :: detector_ene

! Input variables (for namelists)
type(t_input) :: input


!--------------------------------------------------------------------

! Initialize MPI interface
call MPI_INIT(mpierr)
if (mpierr .ne. 0) print *, "mpierr = ", mpierr
call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, mpierr)
if (mpierr .ne. 0) print *, "mpierr = ", mpierr
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, mpierr)

time1st = 0
time2nd = 0
time3rd = 0

starttime = mpi_wtime()

!--------------------------------------------------------------------

call read_input_file( input, MPI_COMM_WORLD )

!-------------- results arrays allocation ---------------------

timea1 = mpi_wtime()

select case (input%diag_type)

  case ("standard","farfield","farfieldEndPoints")
  
    !------------------------- make w axis------------------------

    if (myid .eq. 0) print *, "Setting up waxis and detector"

    call make_waxis( input , detector_spec, numprocs, MPI_COMM_WORLD )
  
    call setup_detector( detector_spec, input )
  
  case ("energy")
  
    !------------------------- setup parallel axis ------------------------

    if (myid .eq. 0) print *, "Setting up xaxis and detector"

    call make_xparallelaxis( input, detector_ene, numprocs, MPI_COMM_WORLD )

    call setup_detector_ene( detector_ene, input )


  case default
  
    print *, ""               

end select

timea2 = mpi_wtime()

timea2 = timea2-timea1                   

!--------------------------------------------------------------------

! Initialize FORTRAN interface.
call h5open_f(error) 

if (myid .eq. 0) then
    
  ! Open the hdf track file
  call h5fopen_f(input%trackfile, H5F_ACC_RDONLY_F, file_id, error)
  if (error .ne. 0) print *, "error opening trackfile = ", error
  call h5gn_members_f(file_id, "/", num_tracks, error)
  if (error .ne. 0) print *, "error reading trackfile group members = ", error
  print *, " "
  print *, "# tracks = ", num_tracks
  call h5fclose_f(file_id, error)
  
  if (input%npart .gt. num_tracks) then
    print *, "error: npart must be <= number of tracks in file"
    stop
  endif

endif  

! Put barrier to ensure synchronization between processors
call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
if (mpierr .ne. 0) print *, "Error in barrier 0"
  
!####################################################################  
!--------------------------------------------------------------------
! Calculate radiation: energy, spectra

call calc_radiation( input, track, detector_ene, detector_spec, myid, MPI_COMM_WORLD ) 

!--------------------------------------------------------------------
!####################################################################


select case (input%diag_type)

  case ("standard","farfield","farfieldEndPoints")
  
    ! Save spectrum
    call save_spec( input, detector_spec, input%diag_type, track%dtinitial, &
                    input%filename, input%parallelIO, MPI_COMM_WORLD )                 

    if (myid .eq. 0) then
      print *, " "
      print *, "After saving spectrum file"
    endif  
    
    call cleanup_detector( detector_spec )
    
    if (myid .eq. 0) then
      print *, " "
      print *, "After cleanup detector" 
    endif  
    
  case ("energy")
  
    ! Save energy
    call save_ene( detector_ene, track%dtinitial, &
                   input%filename, input%parallelIO, MPI_COMM_WORLD ) 
    if (myid .eq. 0) then
      print *, " "
      print *, "After saving energy file"
    endif
      
    call cleanup_detector_ene( detector_ene )
    
    if (myid .eq. 0) then
      print *, " "
      print *, "After cleanup detector" 
      print *, " "
    endif  
               
  case default
  
    if (myid .eq. 0) print *, "diag_type not available: nothing to save"                 

end select

!###############################################################


! Close FORTRAN interface.
call h5close_f(error)

endtime = mpi_wtime()

if (myid .eq. 0) then 
  print *, " "
  print *, "total time = ", endtime-starttime, " s"
endif  

! Close MPI interface
call mpi_finalize(mpierr)

!###############################################################
   
end program jrad
