! ****************************************************************************
! Copyright © 2014 Instituto Superior Técnico. All rights reserved. This 
! software is the copyrighted work of Instituto Superior Técnico. 
! Reproduction, in whole or in part, on the Internet, on CD-ROM or any 
! other medium, without the prior written consent of Instituto Superior 
! Técnico is prohibited.
! ****************************************************************************

module hdf5_utilities

use hdf5

implicit none

include 'mpif.h'
  
! ############ inputs ############
  
! data from track file
integer, parameter, public :: p_double = kind(1.0d0)
integer, parameter, public :: p_single = kind(1.0e0)
integer, parameter, public :: p_k_is    = 4
integer, public, parameter :: p_int64 = selected_int_kind(8)

! parameters  
real (kind=p_double), parameter :: twopi = 6.28318530717959d0
real (kind=p_double), parameter :: pi    = 3.14159265358979d0
    
! #####################################

interface freemem
    module procedure freemem_1d_db
    module procedure freemem_2d_db
    module procedure freemem_3d_db
end interface

interface setup_spec_file
    module procedure setup_spec_file1d
    module procedure setup_spec_file2d
    module procedure setup_spec_file3d
end interface
  
interface setup_ene_file
    module procedure setup_ene_file2d
end interface
     
! #####################################
  
public :: addstringattribute, addstringarrayattribute, &
          addfloatattribute, addfloatarrayattribute, &
          addintegerattribute, addintegerarrayattribute, &
          adddoubleattribute, adddoublearrayattribute, &
          createfile, closefile, addgroup, &
          adddataset1d, adddataset2d, &
          adddatasetparallel1d, adddatasetparallel2d, freemem     
public :: setup_spec_file, setup_ene_file
  
contains

!-----------------------------------------------------------------------------

  
!--------------- subroutine to create attribute--------------------------
subroutine addstringattribute(object_id, att_name, att_data)
	  
! subroutine arguments
integer(hid_t), intent(in) :: object_id
character(len=*), intent(in) :: att_data 
character(len=*), intent(in) :: att_name
	
! local variables
! attribute, att. space and att. type identifier
integer(hid_t) :: att_id
integer(hid_t) :: att_space_id, att_type_id
! attribute rank, dims, length, data and data dims
integer :: att_rank
integer(hsize_t), dimension(1) :: att_dims
integer(size_t) :: att_len
! Debug variables
integer :: error
	 
att_rank = 1 
att_dims(1) = 1 
att_len = len(att_data)
	  
! Create scalar data space for the attribute. 
call h5screate_simple_f(att_rank, att_dims, att_space_id, error)

! Create datatype for the attribute.
call h5tcopy_f(H5T_NATIVE_CHARACTER, att_type_id, error)
if (error .ne. 0) print *, "error in h5tcopy_f", error
call h5tset_size_f(att_type_id, att_len, error)
if (error .ne. 0) print *, "error in h5tset_size_f", error
     
! Create dataset attribute.
call h5acreate_f(object_id, att_name, att_type_id, att_space_id, &
                 att_id, error)
if (error .ne. 0) print *, "error in h5acreate_f", error                 

! Write the attribute data.
call h5awrite_f(att_id, att_type_id, att_data, att_dims, error)
if (error .ne. 0) print *, "error in h5awrite_f", error

! Close the attribute. 
call h5aclose_f(att_id, error)
    
! Close datatype and dataspace
call h5tclose_f( att_type_id, error )
call h5sclose_f( att_space_id, error )
  	
end subroutine addstringattribute
!--------------------------- end of subroutine --------------------------
  
!--------------- subroutine to create attribute--------------------------
subroutine addstringarrayattribute(object_id, att_name, att_data)
	  
! subroutine arguments
integer(hid_t), intent(in) :: object_id
character(len=*), dimension(:), intent(in) :: att_data 
character(len=*), intent(in) :: att_name
	
! local variables
! attribute, att. space and att. type identifier
integer(hid_t) :: att_id
integer(hid_t) :: att_space_id, att_type_id
! attribute rank, dims, length, data and data dims
integer :: att_rank
integer(hsize_t), dimension(1) :: att_dims
integer(size_t) :: att_len
integer :: j
! Debug variables
integer :: error
	 
att_rank = 1 
att_dims(1) = size(att_data) 
att_len = 0
	
do j = 1, att_dims(1)
  if ( len(att_data(j)) .gt. att_len) att_len = len(att_data(j))
enddo
	  
! Create scalar data space for the attribute. 
call h5screate_simple_f(att_rank, att_dims, att_space_id, error)

! Create datatype for the attribute.
call h5tcopy_f(H5T_NATIVE_CHARACTER, att_type_id, error)
call h5tset_size_f(att_type_id, att_len, error)
     
! Create dataset attribute.
call h5acreate_f(object_id, att_name, att_type_id, att_space_id, &
                 att_id, error)

! Write the attribute data.
call h5awrite_f(att_id, att_type_id, att_data, att_dims, error)

! Close the attribute. 
call h5aclose_f(att_id, error)
    
! Close datatype and dataspace
call h5tclose_f( att_type_id, error )
call h5sclose_f( att_space_id, error )
  	
end subroutine addstringarrayattribute
!--------------------------- end of subroutine --------------------------
    
!--------------- subroutine to create attribute--------------------------
subroutine addfloatattribute(object_id, att_name, att_data)
	  
! subroutine arguments
integer(hid_t), intent(in) :: object_id
real (kind=p_single), intent(in) ::  att_data 
character(len=*), intent(in) :: att_name        ! Attribute name
	
! attribute, att. space and att. type identifier
integer(hid_t) :: att_id
integer(hid_t) :: att_space_id
! attribute rank, dims, length, data and data dims
integer :: att_rank
integer(hsize_t), dimension(1) :: att_dims
! Debug variables
integer :: error
	
att_rank = 1
att_dims(1) = 1
error = 0
	
! Create scalar data space for the attribute. 
call h5screate_simple_f(att_rank, att_dims, att_space_id, error)
if (error .ne. 0) print *, "error creating dataspace in addfloatattribute"
    
! Create dataset attribute.
call h5acreate_f(object_id, att_name, H5T_NATIVE_REAL, att_space_id, &
                 att_id, error)
if (error .ne. 0) print *, "error creating attribute in addfloatattribute"                 

! Write the attribute data.
call h5awrite_f(att_id, H5T_NATIVE_REAL, att_data, att_dims, error)
if (error .ne. 0) print *, "error writing attribute in addfloatattribute"
    
! Close the attribute. 
call h5aclose_f(att_id, error)
if (error .ne. 0) print *, "error closing attribute in addfloatattribute"
  	
call h5sclose_f( att_space_id, error )
    	
end subroutine addfloatattribute
!--------------------------- end of subroutine --------------------------
  
!--------------- subroutine to create attribute--------------------------
subroutine addfloatarrayattribute(object_id, att_name, att_data)
	  
! subroutine arguments
integer(hid_t), intent(in) :: object_id
real (kind=p_single), dimension(:), intent(in) ::  att_data 
character(len=*), intent(in) :: att_name     ! Attribute name
	
! attribute, att. space and att. type identifier
integer(hid_t) :: att_id
integer(hid_t) :: att_space_id
! attribute rank, dims, length, data and data dims
integer :: att_rank
integer(hsize_t), dimension(1) :: att_dims
! Debug variables
integer :: error
	
att_rank = 1
att_dims(1) = size(att_data)
	
! Create scalar data space for the attribute. 
call h5screate_simple_f(att_rank, att_dims, att_space_id, error)

! Create dataset attribute.
call h5acreate_f(object_id, att_name, H5T_NATIVE_REAL, att_space_id, &
                 att_id, error)

! Write the attribute data.
call h5awrite_f(att_id, H5T_NATIVE_REAL, att_data, att_dims, error)
     
! Close the attribute. 
call h5aclose_f(att_id, error)

call h5sclose_f( att_space_id, error )
  	
end subroutine addfloatarrayattribute
!--------------------------- end of subroutine --------------------------
  
  
!--------------- subroutine to create attribute--------------------------
subroutine adddoubleattribute(object_id, att_name, att_data)
	  
! subroutine arguments
integer(hid_t), intent(in) :: object_id
real (kind=p_double), intent(in) ::  att_data 
character(len=*), intent(in) :: att_name      ! Attribute name
	
! attribute, att. space and att. type identifier
integer(hid_t) :: att_id
integer(hid_t) :: att_space_id
! attribute rank, dims, length, data and data dims
integer :: att_rank
integer(hsize_t), dimension(1) :: att_dims
! Debug variables
integer :: error
	
att_rank = 1
att_dims(1) = 1
error = 0
	
! Create scalar data space for the attribute. 
call h5screate_simple_f(att_rank, att_dims, att_space_id, error)
if (error .ne. 0) print *, "error creating dataspace in addfloatattribute"
    
! Create dataset attribute.
call h5acreate_f(object_id, att_name, H5T_NATIVE_double, att_space_id, &
                 att_id, error)
if (error .ne. 0) print *, "error creating attribute in addfloatattribute"                 

! Write the attribute data.
call h5awrite_f(att_id, H5T_NATIVE_double, att_data, att_dims, error)
if (error .ne. 0) print *, "error writing attribute in addfloatattribute"
    
! Close the attribute. 
call h5aclose_f(att_id, error)

call h5sclose_f( att_space_id, error )
  	
end subroutine adddoubleattribute
!--------------------------- end of subroutine --------------------------
  
  
!--------------- subroutine to create attribute--------------------------
subroutine adddoublearrayattribute(object_id, att_name, att_data)
	  
! subroutine arguments
integer(hid_t), intent(in) :: object_id
real (kind=p_double), dimension(:), intent(in) ::  att_data 
character(len=*), intent(in) :: att_name    ! Attribute name
	
! attribute, att. space and att. type identifier
integer(hid_t) :: att_id
integer(hid_t) :: att_space_id
! attribute rank, dims, length, data and data dims
integer :: att_rank
integer(hsize_t), dimension(1) :: att_dims
! Debug variables
integer :: error
	
att_rank = 1
att_dims(1) = size(att_data)
	
! Create scalar data space for the attribute. 
call h5screate_simple_f(att_rank, att_dims, att_space_id, error)

! Create dataset attribute.
call h5acreate_f(object_id, att_name, H5T_NATIVE_double, att_space_id, &
                 att_id, error)

! Write the attribute data.
call h5awrite_f(att_id, H5T_NATIVE_double, att_data, att_dims, error)
     
! Close the attribute. 
call h5aclose_f(att_id, error)

call h5sclose_f( att_space_id, error )
  	
end subroutine adddoublearrayattribute
!--------------------------- end of subroutine --------------------------
  
!--------------- subroutine to create attribute--------------------------
subroutine addintegerattribute(object_id, att_name, att_data)
	  
! subroutine arguments
integer(hid_t), intent(in) :: object_id
character(len=*), intent(in) :: att_name     ! Attribute name
integer, intent(in) ::  att_data  
! attribute, att. space and att. type identifier
integer(hid_t) :: att_id
integer(hid_t) :: att_space_id
! attribute rank, dims, length, data and data dims
integer :: att_rank
integer(hsize_t), dimension(1) :: att_dims
! Debug variables
integer :: error
	
att_rank = 1
att_dims(1) = 1
! Create scalar data space for the attribute. 
call h5screate_simple_f(att_rank, att_dims, att_space_id, error)

! Create dataset attribute.
call h5acreate_f(object_id, att_name, H5T_NATIVE_INTEGER, att_space_id, &
                 att_id, error)

! Write the attribute data.
call h5awrite_f(att_id, H5T_NATIVE_INTEGER, att_data, att_dims, error)
     
! Close the attribute. 
call h5aclose_f(att_id, error)

call h5sclose_f( att_space_id, error )
  	
end subroutine addintegerattribute
!--------------------------- end of subroutine --------------------------

  
!--------------- subroutine to create attribute--------------------------
subroutine addintegerarrayattribute(object_id, att_name, att_data)
	  
! subroutine arguments
integer(hid_t), intent(in) :: object_id
character(len=*), intent(in) :: att_name    ! Attribute name
integer, dimension(:), intent(in) ::  att_data  
! attribute, att. space and att. type identifier
integer(hid_t) :: att_id
integer(hid_t) :: att_space_id
! attribute rank, dims, length, data and data dims
integer :: att_rank
integer(hsize_t), dimension(1) :: att_dims
! Debug variables
integer :: error
	
att_rank = 1
att_dims(1) = size(att_data)
    
! Create scalar data space for the attribute. 
call h5screate_simple_f(att_rank, att_dims, att_space_id, error)

! Create dataset attribute.
call h5acreate_f(object_id, att_name, H5T_NATIVE_INTEGER, att_space_id, &
                 att_id, error)

! Write the attribute data.
call h5awrite_f(att_id, H5T_NATIVE_INTEGER, att_data, att_dims, error)
     
! Close the attribute. 
call h5aclose_f(att_id, error)

call h5sclose_f( att_space_id, error )
  	
end subroutine addintegerarrayattribute
!--------------------------- end of subroutine --------------------------
    
!--------------- subroutine to create group -----------------------------
subroutine addgroup(object_id, groupname, group_id)
	  
! subroutine arguments
! id of object to which the group will belong to
integer(hid_t), intent(in) :: object_id
integer(hid_t), intent(out) :: group_id
character(len=18), intent(in) :: groupname ! Group name

integer :: error
	  
! Create a group in the file.
call h5gcreate_f(object_id, groupname, group_id, error)

! Close the group.
call h5gclose_f(group_id, error)

end subroutine addgroup
!--------------------------- end of subroutine --------------------------

   
!------------------ subroutine to shift vector right --------------------
subroutine createfile(filename, file_id, myid, parallel)

character(len=*), intent(in) :: filename
integer(hid_t), intent(inout) :: file_id
integer, intent(in) :: myid
logical, intent(in) :: parallel
! local variables
character(len=150) :: filesavename
integer(hid_t) :: plist_id
integer :: fileInfo
! use independent or collective transfers
integer :: transferMode 
integer :: error
logical :: gpfs
  
transferMode = H5FD_MPIO_COLLECTIVE_F
fileInfo = MPI_INFO_NULL
gpfs = .false.

filesavename = trim(filename) 
   
    
if (parallel .neqv. .true.) then
    
   
  if (myid .eq. 0) then
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    if (error .ne. 0) print *, "error creating plist, myid = ", myid
        
    ! Create the file to save result in slabs.
    call h5fcreate_f(filesavename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
    if (error .ne. 0) print *, "error creating file, myid = ", myid
        
    ! Close the property list
    call h5pclose_f(plist_id, error)
    if (error .ne. 0) print *, "error closing file property list, myid = ", myid
  endif
      
else
    
    
  ! Setup file property list for parallel I/O (with MPI)
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  if (error .ne. 0) print *, "error creating plist, myid = ", myid
  call MPI_INFO_CREATE( fileInfo, error )
  if (error .ne. 0) print *, "error in MPI_INFO_CREATE, myid = ", myid

  if ( gpfs ) then
    call MPI_INFO_SET( fileInfo, "IBM_largeblock_io", "true", error)
    if (error .ne. 0) print *, "error setting gpfs option, myid = ", myid
  endif
    
  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
  if (error .ne. 0) print *, "error setting plist, myid = ", myid
      
  ! Create the file to save result collectively.
  call h5fcreate_f(filesavename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
  if (error .ne. 0) print *, "error creating file, myid = ", myid
      
  call h5pclose_f(plist_id, error)
  if (error .ne. 0) print *, "error closing plist, myid = ", myid
      
  ! free fileInfo if necessary
  if ( fileInfo /= MPI_INFO_NULL ) then
    call MPI_INFO_FREE( fileInfo, error )
    if (error .ne. 0) print *, "error freeing fileInfo, myid = ", myid
  endif
   
endif
   
    
end subroutine createfile
!--------------------------- end of subroutine --------------------------

  
!------------------ subroutine to shift vector right --------------------
subroutine closefile(file_id, myid, parallel)

integer(hid_t), intent(inout) :: file_id
integer, intent(in) :: myid
logical, intent(in) :: parallel
integer :: error
    
if (parallel .neqv. .true.) then
  if (myid .eq. 0) then
    ! Create the file to save result in slabs.
    call h5fclose_f(file_id, error)
    if (error .ne. 0) print *, "error closing file, myid = ", myid
  endif
else
  call h5fclose_f(file_id, error)
  if (error .ne. 0) print *, "error closing file, myid = ", myid
endif
    
end subroutine closefile
!--------------------------- end of subroutine --------------------------
  
  
!------------------ subroutine to write dataset in root --------------------
subroutine adddataset3d( file_id, dataset, dimsf, dims_chunk, firstwvec, &
                         name, lname, units, comm )
    
integer(hid_t), intent(in) :: file_id
real (kind=p_double), dimension(:,:,:), intent(in), target :: dataset
integer(hsize_t), dimension(3), intent(in) :: dimsf, dims_chunk
integer, dimension(:), intent(in) :: firstwvec
character(len=*), intent(in) :: name, lname, units
integer, intent(in) :: comm
! local variables
integer :: dset_rank
real (kind=p_double), dimension(:,:,:), pointer :: write_data
real (kind=p_double), dimension(:,:,:), pointer :: comm_buffer_1 => null(), comm_buffer_2 => null()
integer :: source_id, numprocs, error, myid
integer :: ping, ping_handle, comm_handle, comm_size
integer(hid_t) :: filespace_id, dset_id, memspace_id, dcplID
integer, dimension(MPI_STATUS_SIZE) :: stat
integer :: comm_tag
integer(hsize_t), dimension(3) :: offset, count, block, stride
logical :: dataspaceExists
    
! We will assume all datasets have the same size, given by dims_chunk
    
comm_tag = 0
    
dset_rank = 3
    
dataspaceExists = .false.
   
if ( comm == MPI_COMM_NULL ) then
  
  myid = 0
  numprocs = 1

else

  call MPI_COMM_RANK( comm, myid, error )
  call MPI_COMM_SIZE( comm, numprocs, error )

endif  
    
comm_size = dims_chunk(1)*dims_chunk(2)*dims_chunk(3)
    
if (myid == 0) then

  ! crete property list for datset creation
  call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, error)
  
  ! set chunk size
  call h5pset_chunk_f(dcplID, dset_rank, dims_chunk, error)
      
  ! ----------------------- Dataset creation --------------------------
  ! Create the data space for the  dataset. 
  call h5screate_simple_f(dset_rank, dimsf, filespace_id, error)
  if (error .ne. 0) print *, "error creating filespace, myid = ", myid
  ! Create dataset.
  call h5dcreate_f(file_id, name, H5T_NATIVE_double, filespace_id, &
                   dset_id, error, dcplID)
  ! SPECTRUM dataset attributes
  call addstringattribute(dset_id, 'UNITS', units)
  call addstringattribute(dset_id, 'LONG_NAME', lname)     
                            
  if (error .ne. 0) print *, "error creating chunked dataset, myid = ", myid                 
  call h5sclose_f(filespace_id, error)
  if (error .ne. 0) print *, "error closing filespace, myid = ", myid
  dataspaceExists = .true.
  ! -------------------------------------------------------------------

  write_data => dataset

  ! --------------------- Dataset chunk writing -----------------------
  do source_id = 0, numprocs-1
    
    if ( source_id < numprocs - 1 ) then
	  ! notify node that we are ready
      call mpi_isend( ping, 1, MPI_INTEGER, source_id + 1, comm_tag, comm, ping_handle, error )

	  ! post receive for data from node
      if ( mod( source_id, 2 ) == 0 ) then
        comm_buffer_1 => null()
        allocate( comm_buffer_1(dims_chunk(1),dims_chunk(2),dims_chunk(3)), stat = error )
        if (error .ne. 0) print *, "my id, comm_buffer1, error = ", myid, error
	         
        call mpi_irecv( comm_buffer_1, comm_size, MPI_DOUBLE_PRECISION, source_id + 1, &
                        comm_tag, comm, comm_handle, error )
    else
        comm_buffer_2 => null()
        allocate( comm_buffer_2(dims_chunk(1),dims_chunk(2),dims_chunk(3)), stat = error )
        if (error .ne. 0) print *, "my id, comm_buffer2, error = ", myid, error
        call mpi_irecv( comm_buffer_2, comm_size, MPI_DOUBLE_PRECISION, source_id + 1, &
                        comm_tag, comm, comm_handle, error )
      endif
    endif 
       
    ! write available data
    stride(1) = 1
    stride(2) = 1
    stride(3) = 1 
    count(1) =  1
    count(2) =  1
    count(3) =  1 
    block(1) = dims_chunk(1)
    block(2) = dims_chunk(2)
    block(3) = dims_chunk(3)
    if (source_id .eq. 0) then
      offset(1) = 0
    else  
      offset(1) = firstwvec(source_id+1)
    endif
    offset(2) = 0
    offset(3) = 0
        
    ! create memory dataspace
    call h5screate_simple_f(dset_rank, block, memspace_id, error)
	   
	! select hyperslab in the file	   
    call h5dget_space_f(dset_id, filespace_id, error)
	   
    call h5sselect_hyperslab_f( filespace_id, H5S_SELECT_SET_F, offset, count, error, &
                                stride, block)
	
	! write data
    call h5dwrite_f(dset_id, H5T_NATIVE_double, write_data, dimsf, error, &
                    file_space_id = filespace_id, mem_space_id = memspace_id)                    

	! close resources
    call h5sclose_f(filespace_id, error)
    call h5sclose_f(memspace_id, error)
	    
    if ( source_id < numprocs-1 ) then
      ! wait for messages to complete
      call mpi_wait( ping_handle, stat, error )
      call mpi_wait( comm_handle, stat, error )
    endif
	   
	! free available data
    if ( mod( source_id, 2 ) == 0 ) then
      call freemem( comm_buffer_2 )
      comm_buffer_2 => null()
      write_data => comm_buffer_1
    else
      call freemem( comm_buffer_1 )
      comm_buffer_1 => null()
      write_data => comm_buffer_2
    endif

  enddo
	  
  ! close the dataset
  call h5dclose_f( dset_id, error )
  !call h5pclose_f( dcplID, error )
	  
  ! debug
  if ( associated( comm_buffer_1) .or. associated( comm_buffer_2 ) ) then
    write(0,*) '(*error*) comm buffers are still allocated '
    stop
  endif
	  
else
    
      
  ! wait for root node to be ready
  call mpi_irecv( ping, 1, MPI_INTEGER, 0, comm_tag, comm, ping_handle, error )
  call mpi_wait( ping_handle, stat, error )


  ! send data to root node
  call mpi_isend( dataset, size(dataset), MPI_DOUBLE_PRECISION, 0, comm_tag, comm, &
                  comm_handle, error )
  call mpi_wait( comm_handle, stat, error )   
      
      
endif
    
end subroutine adddataset3d
!--------------------------- end of subroutine --------------------------
  
  
!------------------ subroutine to write dataset in root --------------------
subroutine adddataset2d( file_id, dataset, dimsf, dims_chunk, firstwvec, &
                         name, lname, units, comm )
    
integer(hid_t), intent(in) :: file_id
real (kind=p_double), dimension(:,:), intent(in), target :: dataset
integer(hsize_t), dimension(2), intent(in) :: dimsf, dims_chunk
integer, dimension(:), intent(in) :: firstwvec
character(len=*), intent(in) :: name, lname, units
integer, intent(in) :: comm
! local variables
integer :: dset_rank
real (kind=p_double), dimension(:,:), pointer :: write_data
real (kind=p_double), dimension(:,:), pointer :: comm_buffer_1 => null(), comm_buffer_2 => null()
integer :: source_id, numprocs, error, myid
integer :: ping, ping_handle, comm_handle, comm_size
integer(hid_t) :: filespace_id, dset_id, memspace_id, dcplID
integer, dimension(MPI_STATUS_SIZE) :: stat
integer :: comm_tag
integer(hsize_t), dimension(2) :: offset, count, block, stride
logical :: dataspaceExists
    
! We will assume all datasets have the same size, given by dims_chunk
    
comm_tag = 0
    
dset_rank = 2
    
dataspaceExists = .false.
    
if ( comm == MPI_COMM_NULL ) then
  myid = 0
  numprocs = 1
else
  call MPI_COMM_RANK( comm, myid, error )
  call MPI_COMM_SIZE( comm, numprocs, error )
endif  
    
comm_size = dims_chunk(1)*dims_chunk(2)
    
if (myid == 0) then

  ! crete property list for datset creation
  call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, error)
  
  ! set chunk size
  call h5pset_chunk_f(dcplID, dset_rank, dims_chunk, error)
    
  ! ----------------------- Dataset creation --------------------------
  ! Create the data space for the  dataset. 
  call h5screate_simple_f(dset_rank, dimsf, filespace_id, error)
  if (error .ne. 0) print *, "error creating filespace, myid = ", myid
  ! Create dataset.
  call h5dcreate_f(file_id, name, H5T_NATIVE_double, filespace_id, &
                   dset_id, error, dcplID)
  ! SPECTRUM dataset attributes
  call addstringattribute(dset_id, 'UNITS', units)
  call addstringattribute(dset_id, 'LONG_NAME', lname)     
                            
  if (error .ne. 0) print *, "error creating chunked dataset, myid = ", myid                 
  call h5sclose_f(filespace_id, error)
  if (error .ne. 0) print *, "error closing filespace, myid = ", myid
  dataspaceExists = .true.
  ! -------------------------------------------------------------------

  write_data => dataset

  ! --------------------- Dataset chunk writing -----------------------
  do source_id = 0, numprocs-1
    
    if ( source_id < numprocs - 1 ) then
	  ! notify node that we are ready
      call mpi_isend( ping, 1, MPI_INTEGER, source_id + 1, comm_tag, comm, ping_handle, error )

	  ! post receive for data from node
      if ( mod( source_id, 2 ) == 0 ) then
        comm_buffer_1 => null()
        allocate( comm_buffer_1(dims_chunk(1),dims_chunk(2)), stat = error )
        if (error .ne. 0) print *, "my id, error = ", myid, error
	         
        call mpi_irecv( comm_buffer_1, comm_size, MPI_DOUBLE_PRECISION, source_id + 1, &
                        comm_tag, comm, comm_handle, error )
      else
        comm_buffer_2 => null()
        allocate( comm_buffer_2(dims_chunk(1),dims_chunk(2)), stat = error )
        if (error .ne. 0) print *, "my id, error = ", myid, error
        call mpi_irecv( comm_buffer_2, comm_size, MPI_DOUBLE_PRECISION, source_id + 1, &
                        comm_tag, comm, comm_handle, error )
      endif
    endif 
       
    ! write available data
    stride(1) = 1
    stride(2) = 1 
    count(1) =  1
    count(2) =  1 
    block(1) = dims_chunk(1)
    block(2) = dims_chunk(2)
    if (source_id .eq. 0) then
      offset(1) = 0
    else  
      offset(1) = firstwvec(source_id+1)
    endif
    offset(2) = 0

	! create memory dataspace
    call h5screate_simple_f(dset_rank, block, memspace_id, error)
	   
	! select hyperslab in the file	   
    call h5dget_space_f(dset_id, filespace_id, error)
	   
    call h5sselect_hyperslab_f( filespace_id, H5S_SELECT_SET_F, offset, count, error, &
                                stride, block)
	
	! write data
    call h5dwrite_f( dset_id, H5T_NATIVE_double, write_data, dimsf, error, &
                     file_space_id = filespace_id, mem_space_id = memspace_id)                    

    ! close resources
    call h5sclose_f(filespace_id, error)
    call h5sclose_f(memspace_id, error)
	    
    if ( source_id < numprocs-1 ) then
      ! wait for messages to complete
      call mpi_wait( ping_handle, stat, error )
      call mpi_wait( comm_handle, stat, error )
    endif
	   
    ! free available data
    if ( mod( source_id, 2 ) == 0 ) then
      call freemem( comm_buffer_2 )
      comm_buffer_2 => null()
      write_data => comm_buffer_1
    else
      call freemem( comm_buffer_1 )
      comm_buffer_1 => null()
      write_data => comm_buffer_2
    endif

  enddo
	  
  ! close the dataset
  call h5dclose_f( dset_id, error )
  !call h5pclose_f( dcplID, error )
	  
  ! debug
  if ( associated( comm_buffer_1) .or. associated( comm_buffer_2 ) ) then
    write(0,*) '(*error*) comm buffers are still allocated '
    stop
  endif
	  
else
    
      
  ! wait for root node to be ready
  call mpi_irecv( ping, 1, MPI_INTEGER, 0, comm_tag, comm, ping_handle, error )
  call mpi_wait( ping_handle, stat, error )


  ! send data to root node
  call mpi_isend( dataset, size(dataset), MPI_DOUBLE_PRECISION, 0, comm_tag, comm, &
                  comm_handle, error )
  call mpi_wait( comm_handle, stat, error )   
      
      
endif
    
end subroutine adddataset2d
!--------------------------- end of subroutine --------------------------
  
  
!------------------ subroutine to shift vector right --------------------
subroutine adddataset1d( file_id, name, dataset, dimsf, dims_chunk, firstwvec, comm, &
                         units_att, lname_att, type_att, name_att)
    
integer(hid_t), intent(in) :: file_id
character(len=*), intent(in) :: name
real (kind=p_double), dimension(:), intent(in), target :: dataset
integer(hsize_t), dimension(1), intent(in) :: dimsf, dims_chunk
integer, dimension(:), intent(in) :: firstwvec
integer, intent(in) :: comm
character(len=*), intent(in) :: units_att, lname_att
character(len=*), optional, intent(in) :: type_att, name_att
    
! local variables
integer :: dset_rank
real (kind=p_double), dimension(:), pointer :: write_data
real (kind=p_double), dimension(:), pointer :: comm_buffer_1 => null(), &
                                                   comm_buffer_2 => null()
integer :: source_id, numprocs, error, myid
integer :: ping, ping_handle, comm_handle, comm_size
integer(hid_t) :: filespace_id, dset_id, memspace_id, wgr_id, dcplID
integer, dimension(MPI_STATUS_SIZE) :: stat
integer :: comm_tag
integer(hsize_t), dimension(1) :: offset, count, block, stride
    
! We will assume all datasets have the same size, given by dims_chunk
    
comm_tag = 0
dset_rank = 1
    
call MPI_COMM_RANK( comm, myid, error )
call MPI_COMM_SIZE( comm, numprocs, error )
    
comm_size = dims_chunk(1)
    
if (myid .eq. 0) then

  ! crete property list for datset creation
  call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, error)
  
  ! set chunk size
  call h5pset_chunk_f(dcplID, dset_rank, dims_chunk, error)
    
  ! ----------------------- Dataset creation --------------------------
  ! Create the data space for the  dataset. 
  call h5screate_simple_f(dset_rank, dimsf, filespace_id, error)
  if (error .ne. 0) print *, "error creating filespace, myid = ", myid
      
  ! Create dataset.
  if ((present(name_att)) .and. (name_att .eq. "\omega")) then
    ! ---------------------------
    ! add waxis information
    call h5gcreate_f( file_id, 'WAXIS', wgr_id, error) 
    if (error .ne. 0) print *, "error creating WAXIS group, myid = ", myid
    call h5dcreate_f(wgr_id, name, H5T_NATIVE_double, filespace_id, &
                     dset_id, error, dcplID)
    if (error .ne. 0) print *, "error creating WAXIS dataset, myid = ", myid                  
    ! ---------------------------
  else
    call h5dcreate_f(file_id, name, H5T_NATIVE_double, filespace_id, &
                     dset_id, error)
  endif             
            
  ! SPECTRUM dataset attributes
  if (present(type_att)) then
    call addstringattribute( dset_id, 'TYPE', type_att ) 
  endif  
  call addstringattribute( dset_id, 'UNITS', units_att ) 
  if (present(name_att)) then
    call addstringattribute( dset_id, 'NAME', name_att )
  endif   
  call addstringattribute( dset_id, 'LONG_NAME', lname_att )
                                  
  if (error .ne. 0) print *, "error creating chunked dataset, myid = ", myid                 
  call h5sclose_f(filespace_id, error)
  if (error .ne. 0) print *, "error closing filespace, myid = ", myid
  ! -------------------------------------------------------------------

  write_data => dataset

  ! --------------------- Dataset chunk writing -----------------------
  do source_id = 0, numprocs-1
    
    if ( source_id < numprocs - 1 ) then
	  ! notify node that we are ready
      call mpi_isend( ping, 1, MPI_INTEGER, source_id + 1, comm_tag, comm, &
                      ping_handle, error )

	  ! post receive for data from node
      if ( mod( source_id, 2 ) == 0 ) then
        comm_buffer_1 => null()
        allocate( comm_buffer_1(dims_chunk(1)), stat = error )
        if (error .ne. 0) print *, "my id, error = ", myid, error
	         
        call mpi_irecv( comm_buffer_1, comm_size, MPI_DOUBLE_PRECISION, &
                        source_id + 1, comm_tag, comm, comm_handle, error )
      else
        comm_buffer_2 => null()
        allocate( comm_buffer_2(dims_chunk(1)), stat = error )
        if (error .ne. 0) print *, "my id, error = ", myid, error
        call mpi_irecv( comm_buffer_2, comm_size, MPI_DOUBLE_PRECISION, &
                        source_id + 1, comm_tag, comm, comm_handle, error )
      endif
    endif 
       
    ! write available data
    stride(1) = 1
    count(1) =  1
    block(1) = dims_chunk(1)
    offset(1) = firstwvec(source_id+1)

	! create memory dataspace
    call h5screate_simple_f(dset_rank, block, memspace_id, error)
	   
    ! select hyperslab in the file	   
    call h5dget_space_f(dset_id, filespace_id, error)
  
    call h5sselect_hyperslab_f( filespace_id, H5S_SELECT_SET_F, offset, count, &
                                error, stride, block)
	
	! write data
    call h5dwrite_f(dset_id, H5T_NATIVE_double, write_data, dimsf, error, &
                    file_space_id = filespace_id, mem_space_id = memspace_id)                    

	! close resources
    call h5sclose_f(filespace_id, error)
    call h5sclose_f(memspace_id, error)
	    
    if ( source_id < numprocs-1 ) then
      ! wait for messages to complete
      call mpi_wait( ping_handle, stat, error )
      call mpi_wait( comm_handle, stat, error )
    endif
	   
	! free available data
    if ( mod( source_id, 2 ) == 0 ) then
      call freemem( comm_buffer_2 )
      comm_buffer_2 => null()
      write_data => comm_buffer_1
    else
      call freemem( comm_buffer_1 )
      comm_buffer_1 => null()
      write_data => comm_buffer_2
    endif

  enddo
	  
  ! close the dataset
  call h5dclose_f( dset_id, error )
	  
else
  
  ! wait for root node to be ready
  call mpi_irecv( ping, 1, MPI_INTEGER, 0, comm_tag, comm, ping_handle, error )
  call mpi_wait( ping_handle, stat, error )

  ! send data to root node
  call mpi_isend( dataset, size(dataset), MPI_DOUBLE_PRECISION, 0, comm_tag, comm, &
                  comm_handle, error )
  call mpi_wait( comm_handle, stat, error )   
  
endif
    
if ((present(name_att)) .and. (name_att .eq. "\omega") .and. (myid .eq. 0)) then
  call h5gclose_f( wgr_id, error)                 
endif    
    
end subroutine adddataset1d
!--------------------------- end of subroutine --------------------------
  
  
!------------------ subroutine to write 3d dataset chunk --------------------
subroutine adddatasetparallel3d(file_id, comm, name, lname, units, dimsf, dims_chunk, &
                                firstwvec, datasetchunk)
  
integer(hid_t), intent(in) :: file_id
integer, intent(in) :: comm
character(len=*), intent(in) :: name, lname, units
integer(hsize_t), dimension(3), intent(in) :: dimsf, dims_chunk
integer, dimension(:), intent(in) :: firstwvec
real (kind=p_double), dimension(:,:,:), intent(in) :: datasetchunk
! local variables
integer :: dset_rank
integer(hid_t) :: filespace_id, dset_id, memspace_id
integer :: error, myid
integer(hid_t) :: plist_id      ! Property list identifier 
integer(hsize_t), dimension(3) :: count, offset, stride, block
    
dset_rank = 3
    
call MPI_COMM_RANK( comm, myid, error )
    
    
! Create the data space for the  dataset. 
call h5screate_simple_f(dset_rank, dimsf, filespace_id, error)
if (error .ne. 0) print *, "error creating filespace, myid = ", myid
call h5screate_simple_f(dset_rank, dims_chunk, memspace_id, error)
if (error .ne. 0) print *, "error creating memspace, myid = ", myid
  
! Create chunked dataset.
call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
if (error .ne. 0) print *, "error creating chunked dataset plist, myid = ", myid
call h5pset_chunk_f(plist_id, dset_rank, dims_chunk, error)
if (error .ne. 0) print *, "error setting chunk dims, myid = ", myid
call h5dcreate_f(file_id, name, H5T_NATIVE_double, filespace_id, &
                 dset_id, error, plist_id)
if (error .ne. 0) print *, "error creating chunked dataset, myid = ", myid                 
call h5sclose_f(filespace_id, error)
if (error .ne. 0) print *, "error closing filespace, myid = ", myid
  
call addstringattribute( dset_id, 'UNITS', units ) 
call addstringattribute( dset_id, 'LONG_NAME', lname ) 
  
! Each process defines dataset in memory and writes it to the hyperslab in the file. 
stride(1) = 1 
stride(2) = 1 
stride(3) = 1 
count(1) =  1 
count(2) =  1 
count(3) =  1 
block(1) = dims_chunk(1)
block(2) = dims_chunk(2)
block(3) = dims_chunk(3)
if (myid .eq. 0) then
  offset(1) = 0
else  
  offset(1) = firstwvec(myid+1)
endif  
offset(2) = 0
offset(3) = 0
  
! Select hyperslab in the file.
call h5dget_space_f(dset_id, filespace_id, error)
if (error .ne. 0) print *, "error getting filespace_id, myid = ", myid
call h5sselect_hyperslab_f (filespace_id, H5S_SELECT_SET_F, offset, count, &
                            error, stride, block)
if (error .ne. 0) print *, "error selecting hyperslab in the file, myid = ", myid
                              
! Create property list for collective dataset write
call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
if (error .ne. 0) print *, "error creating transfer plist, myid = ", myid
call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
if (error .ne. 0) print *, "error setting transfer plist, myid = ", myid
  
    
! Write the dataset collectively. 
call h5dwrite_f(dset_id, H5T_NATIVE_double, datasetchunk, dimsf, error, &
                file_space_id = filespace_id, mem_space_id = memspace_id, &
                xfer_prp = plist_id)
if (error .ne. 0) print *, "error writing dataset 2D collectively, myid = ", myid                
  
! Close dataspaces.
call h5sclose_f(filespace_id, error)
if (error .ne. 0) print *, "error closing filespace, myid = ", myid 
call h5sclose_f(memspace_id, error)
if (error .ne. 0) print *, "error closing memspace, myid = ", myid 
    
! Close the dataset.
call h5dclose_f(dset_id, error)
if (error .ne. 0) print *, "error closing dataset, myid = ", myid
! Close the property list.
call h5pclose_f(plist_id, error)
if (error .ne. 0) print *, "error closing plist, myid = ", myid
    
end subroutine adddatasetparallel3d
!--------------------------- end of subroutine --------------------------
  
  
!------------------ subroutine to shift vector right --------------------
subroutine adddatasetparallel2d(file_id, comm, name, lname, units, dimsf,&
                                dims_chunk, firstwvec, datasetchunk)
  
integer(hid_t), intent(in) :: file_id
integer, intent(in) :: comm
character(len=*), intent(in) :: name, lname, units
integer(hsize_t), dimension(2), intent(in) :: dimsf, dims_chunk
integer, dimension(:), intent(in) :: firstwvec
real (kind=p_double), dimension(:,:), intent(in) :: datasetchunk
! local variables
integer :: dset_rank
integer(hid_t) :: filespace_id, dset_id, memspace_id
integer :: error, myid
integer(hid_t) :: plist_id      ! Property list identifier 
integer(hsize_t), dimension(2) :: count, offset, stride, block
    
dset_rank = 2
    
call MPI_COMM_RANK( comm, myid, error )
    
    
! Create the data space for the  dataset. 
call h5screate_simple_f(dset_rank, dimsf, filespace_id, error)
if (error .ne. 0) print *, "error creating filespace, myid = ", myid
call h5screate_simple_f(dset_rank, dims_chunk, memspace_id, error)
if (error .ne. 0) print *, "error creating memspace, myid = ", myid
  
! Create chunked dataset.
call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
if (error .ne. 0) print *, "error creating chunked dataset plist, myid = ", myid
call h5pset_chunk_f(plist_id, dset_rank, dims_chunk, error)
if (error .ne. 0) print *, "error setting chunk dims, myid = ", myid
call h5dcreate_f(file_id, name, H5T_NATIVE_double, filespace_id, &
                 dset_id, error, plist_id)
if (error .ne. 0) print *, "error creating chunked dataset, myid = ", myid                 
call h5sclose_f(filespace_id, error)
if (error .ne. 0) print *, "error closing filespace, myid = ", myid
  
call addstringattribute( dset_id, 'UNITS', units ) 
call addstringattribute( dset_id, 'LONG_NAME', lname ) 
  
! Each process defines dataset in memory and writes it to the hyperslab in the file. 
stride(1) = 1 
stride(2) = 1 
count(1) =  1 
count(2) =  1 
block(1) = dims_chunk(1)
block(2) = dims_chunk(2)
if (myid .eq. 0) then
  offset(1) = 0
else  
  offset(1) = firstwvec(myid+1)
endif  
offset(2) = 0
    
! Select hyperslab in the file.
call h5dget_space_f(dset_id, filespace_id, error)
if (error .ne. 0) print *, "error getting filespace_id, myid = ", myid
call h5sselect_hyperslab_f (filespace_id, H5S_SELECT_SET_F, offset, count, error, &
                            stride, block)
if (error .ne. 0) print *, "error selecting hyperslab in the file, myid = ", myid
                              
! Create property list for collective dataset write
call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
if (error .ne. 0) print *, "error creating transfer plist, myid = ", myid
call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
if (error .ne. 0) print *, "error setting transfer plist, myid = ", myid
    
   
! Write the dataset collectively. 
call h5dwrite_f(dset_id, H5T_NATIVE_double, datasetchunk, dimsf, error, &
                file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
if (error .ne. 0) print *, "error writing dataset 2D collectively, myid = ", myid                
  
! Close dataspaces.
call h5sclose_f(filespace_id, error)
if (error .ne. 0) print *, "error closing filespace, myid = ", myid 
call h5sclose_f(memspace_id, error)
if (error .ne. 0) print *, "error closing memspace, myid = ", myid 
    
! Close the dataset.
call h5dclose_f(dset_id, error)
if (error .ne. 0) print *, "error closing dataset, myid = ", myid
! Close the property list.
call h5pclose_f(plist_id, error)
if (error .ne. 0) print *, "error closing plist, myid = ", myid
    
end subroutine adddatasetparallel2d
!--------------------------- end of subroutine --------------------------
  
  
! ****************************************************************************************  
!------------------ subroutine to shift vector right --------------------
subroutine adddatasetparallel1d(file_id, comm, name, dimsf, dims_chunk, &
                                firstwvec, datasetchunk, &
                                units_att, lname_att, type_att, name_att)
  
integer(hid_t), intent(in) :: file_id
integer, intent(in) :: comm
character(len=*), intent(in) :: name
integer(hsize_t), dimension(1), intent(in) :: dimsf, dims_chunk
integer, dimension(:), intent(in) :: firstwvec
real(kind=p_double), dimension(:), intent(in) :: datasetchunk
character(len=*), intent(in) :: units_att, lname_att
character(len=*), optional, intent(in) :: type_att, name_att
    
! local variables
integer :: dset_rank
integer(hid_t) :: filespace_id, dset_id, memspace_id
integer :: error, myid
integer(hid_t) :: plist_id      ! Property list identifier 
integer(hsize_t), dimension(1) :: count, offset, stride, block
    
dset_rank = 1
    
call MPI_COMM_RANK( comm, myid, error )
    
! Create the data space for the  dataset. 
call h5screate_simple_f(dset_rank, dimsf, filespace_id, error)
if (error .ne. 0) print *, "error creating filespace, myid = ", myid
call h5screate_simple_f(dset_rank, dims_chunk, memspace_id, error)
if (error .ne. 0) print *, "error creating memspace, myid = ", myid
  
! Create chunked dataset.
call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
if (error .ne. 0) print *, "error creating chunked dataset plist, myid = ", myid
call h5pset_chunk_f(plist_id, dset_rank, dims_chunk, error)
if (error .ne. 0) print *, "error setting chunk dims, myid = ", myid
call h5dcreate_f(file_id, name, H5T_NATIVE_double, filespace_id, &
                 dset_id, error, plist_id)
if (error .ne. 0) print *, "error creating chunked dataset, myid = ", myid                 
call h5sclose_f(filespace_id, error)
if (error .ne. 0) print *, "error closing filespace, myid = ", myid
  
! Each process defines dataset in memory and writes it to the hyperslab in the file. 
stride = 1 
count =  1 
block = dims_chunk
offset = firstwvec(myid+1)
   
    
! Select hyperslab in the file.
call h5dget_space_f(dset_id, filespace_id, error)
if (error .ne. 0) print *, "error getting filespace_id, myid = ", myid
call h5sselect_hyperslab_f (filespace_id, H5S_SELECT_SET_F, offset, count, error, &
                            stride, block)
if (error .ne. 0) print *, "error selecting hyperslab in the file, myid = ", myid
                              
! Create property list for collective dataset write
call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
if (error .ne. 0) print *, "error creating transfer plist, myid = ", myid
call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
if (error .ne. 0) print *, "error setting transfer plist, myid = ", myid
     
! Write the dataset collectively. 
call h5dwrite_f(dset_id, H5T_NATIVE_double, datasetchunk, dimsf, error, &
                file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
if (error .ne. 0) print *, "error writing dataset 1D collectively, myid = ", myid                
  
! Close dataspaces.
call h5sclose_f(filespace_id, error)
if (error .ne. 0) print *, "error closing filespace, myid = ", myid 
call h5sclose_f(memspace_id, error)
if (error .ne. 0) print *, "error closing memspace, myid = ", myid 
    
if (present(type_att)) then
  call addstringattribute( dset_id, 'TYPE', type_att ) 
endif  
call addstringattribute( dset_id, 'UNITS', units_att ) 
if (present(name_att)) then
  call addstringattribute( dset_id, 'NAME', name_att )
endif   
call addstringattribute( dset_id, 'LONG_NAME', lname_att )
    
! Close the dataset.
call h5dclose_f(dset_id, error)
if (error .ne. 0) print *, "error closing dataset, myid = ", myid
! Close the property list.
call h5pclose_f(plist_id, error)
if (error .ne. 0) print *, "error closing plist, myid = ", myid
    
end subroutine adddatasetparallel1d
!--------------------------- end of subroutine --------------------------
! ****************************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine setup_spec_file1d( filesavename, filesave_id, myid, &
        parallelIO, axis_min1, axis_max1, axis_units1, axis_name1, &
        axis_lname1, name, time, iter, dt, xmin, xmax, diag_type)
! ----------------------------------------------------------------------------------------
  
character(len=*), intent(in) :: filesavename
integer(hid_t), intent(inout) :: filesave_id
integer, intent(in) :: myid
real(kind=p_single), intent(in) :: axis_min1, axis_max1
real(kind=p_double), intent(in) :: time, dt
integer, intent(in) :: iter
real(kind=p_double), dimension(:), intent(in) :: xmin, xmax
character(len=*), intent(in) :: axis_units1, axis_name1, axis_lname1
character(len=*), intent(in) :: name
character(len=*), intent(in) :: diag_type
    
real(p_double), dimension(2) :: axis_range
integer(hid_t) :: rootID, axisGroupID, dataspaceID, datasetID, plistID
integer(hsize_t), dimension(1) :: dims
integer :: error
  
integer, parameter :: izero = ichar('0')
logical :: parallelIO

! For serial I/O only node 0 creates the file and sets atributes;
! however for parallel I/O this must be done by all nodes.
! ---------------------- PARALLEL IO -------------------------------------------

if (parallelIO) then
   
  ! create the file
  call createfile(filesavename, filesave_id, myid, parallelIO)

  ! this is required for older versions of hdf5 that don't allow 
  ! setting atributes for the file, only the root group
  call h5gopen_f( filesave_id, '/', rootID, error )
 
  ! add name property
  call addstringattribute(rootID, 'NAME', name)

  ! add file attributes
  ! create a collective write property list for axis values
  call h5pcreate_f(H5P_DATASET_XFER_F, plistID, error) 
  call h5pset_dxpl_mpio_f(plistID, H5FD_MPIO_COLLECTIVE_F, error)

  call addstringattribute( rootID, 'TYPE', 'grid' ) 
	  
  ! add axis information
  call h5gcreate_f( rootID, 'AXIS', axisGroupID, error) 
	  
  dims(1) = 2
  ! ---------- AXIS 1 ---------
  call h5screate_simple_f(1, dims, dataspaceID, error ) 
  call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+1), &
                    H5T_NATIVE_DOUBLE, dataspaceID, &
                    datasetID, error )
  call addstringattribute( datasetID, 'TYPE', 'linear' ) 
  call addstringattribute( datasetID, 'UNITS', axis_units1 ) 
  call addstringattribute( datasetID, 'NAME', axis_name1 ) 
  call addstringattribute( datasetID, 'LONG_NAME', axis_lname1 ) 

  axis_range(1) = axis_min1
  axis_range(2) = axis_max1
		
  call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, axis_range, dims, &
                   error, xfer_prp = plistID )
		
  call h5dclose_f( datasetID, error )
  
  ! ---------------------------
  call h5gclose_f(axisGroupID, error) 

  ! ---------------------------
  !! add waxis information
  call h5gcreate_f( rootID, 'WAXIS', axisGroupID, error) 
  call h5gclose_f(axisGroupID, error) 
  ! ---------------------------

  call h5pclose_f(plistID, error)
  
  ! time information 
  call adddoubleattribute( rootID, 'TIME', time ) 
  call addintegerattribute( rootID, 'ITER', iter ) 
  
  ! simulation information
  ! possibly move to a separate group
  call adddoubleattribute( rootID, 'DT', dt )
  call addstringattribute( rootID, 'TIME UNITS', '1/\omega_0' ) 
  
  call adddoublearrayattribute( rootID, 'XMIN', xmin ) 
  call adddoublearrayattribute( rootID, 'XMAX', xmax ) 
  call addintegerarrayattribute(filesave_id, 'PERIODIC', (/0/))
  call addintegerarrayattribute(filesave_id, 'MOVE C', (/0/))
  
  call addstringattribute( rootID, 'SPEC_METHOD', diag_type )
    
  call h5gclose_f( rootID, error )

! ---------------------- SERIAL IO -------------------------------------------
else

  ! create the file
  call createfile(filesavename, filesave_id, myid, parallelIO)

  if ( myid == 0 ) then
  
    ! this is required for older versions of hdf5 that don't allow 
    ! setting atributes for the file, only the root group
    call h5gopen_f( filesave_id, '/', rootID, error )
 
    ! add name property
    call addstringattribute(rootID, 'NAME', name)

    ! add file attributes

    ! set a default write property list for axis values
    plistID = H5P_DEFAULT_F

    call addstringattribute( rootID, 'TYPE', 'grid' ) 
	  
    ! add axis information
    call h5gcreate_f( rootID, 'AXIS', axisGroupID, error) 
	  
    dims(1) = 2
    ! ---------- AXIS 1 ---------
    call h5screate_simple_f(1, dims, dataspaceID, error ) 
    call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+1), &
                      H5T_NATIVE_DOUBLE, dataspaceID, &
                      datasetID, error )
    call addstringattribute( datasetID, 'TYPE', 'linear' ) 
    call addstringattribute( datasetID, 'UNITS', axis_units1 ) 
    call addstringattribute( datasetID, 'NAME', axis_name1 ) 
    call addstringattribute( datasetID, 'LONG_NAME', axis_lname1 ) 

    axis_range(1) = axis_min1
    axis_range(2) = axis_max1
		
    call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, axis_range, dims, &
                     error, xfer_prp = plistID )
		
    call h5dclose_f( datasetID, error )
  
    ! ---------------------------
	
    call h5gclose_f( axisGroupID, error ) 
	  
    ! time information 
    call adddoubleattribute( rootID, 'TIME', time ) 
    call addintegerattribute( rootID, 'ITER', iter ) 
  
    ! simulation information
    ! possibly move to a separate group
    call adddoubleattribute( rootID, 'DT', dt )
    call addstringattribute( rootID, 'TIME UNITS', '1/\omega_0' ) 
  
    call adddoublearrayattribute( rootID, 'XMIN', xmin ) 
    call adddoublearrayattribute( rootID, 'XMAX', xmax ) 
    call addintegerarrayattribute(filesave_id, 'PERIODIC', (/0/))
    call addintegerarrayattribute(filesave_id, 'MOVE C', (/0/))
    
    call addstringattribute( rootID, 'SPEC_METHOD', diag_type ) 
    
    call h5gclose_f( rootID, error )

  endif
endif
  
end subroutine setup_spec_file1d
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************  
! ----------------------------------------------------------------------------------------
subroutine setup_spec_file2d( filesavename, &
        filesave_id, myid, parallelIO, &
        axis_min1, axis_max1, axis_units1, axis_name1, axis_lname1, &
        axis_min2, axis_max2, axis_units2, axis_name2, axis_lname2, &
        name, time, iter, dt, xmin, xmax, diag_type)
! ----------------------------------------------------------------------------------------
  
character(len=*), intent(in) :: filesavename
integer(hid_t), intent(inout) :: filesave_id
integer, intent(in) :: myid
real(kind=p_single), intent(in) :: axis_min1, axis_max1
real(kind=p_single), intent(in) :: axis_min2, axis_max2
real(kind=p_double), intent(in) :: time, dt
integer, intent(in) :: iter
real(kind=p_double), dimension(:), intent(in) :: xmin, xmax
character(len=*), intent(in) :: axis_units1, axis_name1, axis_lname1, &
                                axis_units2, axis_name2, axis_lname2, &
                                name
character(len=*), intent(in) :: diag_type
    
real(p_double), dimension(2) :: axis_range
integer(hid_t) :: rootID, axisGroupID, dataspaceID, datasetID, plistID
integer(hsize_t), dimension(1) :: dims
integer :: error
  
integer, parameter :: izero = ichar('0')
logical :: parallelIO


! For serial I/O only node 0 creates the file and sets atributes;
! however for parallel I/O this must be done by all nodes.
! ---------------------- PARALLEL IO -------------------------------------------

if (parallelIO) then
   
   
  ! create the file
  call createfile(filesavename, filesave_id, myid, parallelIO)

  ! this is required for older versions of hdf5 that don't allow setting 
  ! atributes for the file, only the root group
  call h5gopen_f( filesave_id, '/', rootID, error )
 
  ! add name property
  call addstringattribute(rootID, 'NAME', name)

  ! add file attributes
  ! create a collective write property list for axis values
  call h5pcreate_f(H5P_DATASET_XFER_F, plistID, error) 
  call h5pset_dxpl_mpio_f(plistID, H5FD_MPIO_COLLECTIVE_F, error)

  call addstringattribute( rootID, 'TYPE', 'grid' ) 
	  
  ! add axis information
  call h5gcreate_f( rootID, 'AXIS', axisGroupID, error) 
	  
  dims(1) = 2
  ! ---------- AXIS 1 ---------
  call h5screate_simple_f(1, dims, dataspaceID, error ) 
  call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+1), &
                    H5T_NATIVE_DOUBLE, dataspaceID, &
                    datasetID, error )
  call addstringattribute( datasetID, 'TYPE', 'linear' ) 
  call addstringattribute( datasetID, 'UNITS', axis_units1 ) 
  call addstringattribute( datasetID, 'NAME', axis_name1 ) 
  call addstringattribute( datasetID, 'LONG_NAME', axis_lname1 ) 

  axis_range(1) = axis_min1
  axis_range(2) = axis_max1
		
  call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, axis_range, dims, &
                   error, xfer_prp = plistID )
		
  call h5dclose_f( datasetID, error )
  
  ! ---------- AXIS 2 ---------
  call h5screate_simple_f(1, dims, dataspaceID, error ) 
  call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+2), &
                    H5T_NATIVE_DOUBLE, dataspaceID, &
                    datasetID, error )
  call addstringattribute( datasetID, 'TYPE', 'linear' ) 
  call addstringattribute( datasetID, 'UNITS', axis_units2 ) 
  call addstringattribute( datasetID, 'NAME', axis_name2 ) 
  call addstringattribute( datasetID, 'LONG_NAME', axis_lname2 ) 

  axis_range(1) = axis_min2
  axis_range(2) = axis_max2
		
  call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, axis_range, dims, &
                   error, xfer_prp = plistID )
	
  call h5sclose_f( dataspaceID, error )
  call h5dclose_f( datasetID, error )
  ! ---------------------------
  call h5gclose_f(axisGroupID, error) 

  ! ---------------------------
  call h5gcreate_f( rootID, 'WAXIS', axisGroupID, error) 
  call h5gclose_f(axisGroupID, error) 
  ! ---------------------------

  call h5pclose_f(plistID, error)
  
  ! time information 
  call adddoubleattribute( rootID, 'TIME', time ) 
  call addintegerattribute( rootID, 'ITER', iter ) 
  
  ! simulation information
  ! possibly move to a separate group
  call adddoubleattribute( rootID, 'DT', dt )
  call addstringattribute( rootID, 'TIME UNITS', '1/\omega_0' ) 
  
  call adddoublearrayattribute( rootID, 'XMIN', xmin ) 
  call adddoublearrayattribute( rootID, 'XMAX', xmax ) 
  call addintegerarrayattribute(filesave_id, 'PERIODIC', (/0,0/))
  call addintegerarrayattribute(filesave_id, 'MOVE C', (/0,0/))
  
  call addstringattribute( rootID, 'SPEC_METHOD', diag_type )
    
  call h5gclose_f( rootID, error )

! ---------------------- SERIAL IO -------------------------------------------
else


  ! create the file
  call createfile(filesavename, filesave_id, myid, parallelIO)

  if ( myid == 0 ) then
  
    ! this is required for older versions of hdf5 that don't allow 
    ! setting atributes for the file, only the root group
    call h5gopen_f( filesave_id, '/', rootID, error )
 
    ! add name property
    call addstringattribute(rootID, 'NAME', name)

    ! add file attributes

    ! set a default write property list for axis values
    plistID = H5P_DEFAULT_F

    call addstringattribute( rootID, 'TYPE', 'grid' ) 
	  
    ! add axis information
    call h5gcreate_f( rootID, 'AXIS', axisGroupID, error) 
	  
    dims(1) = 2
    ! ---------- AXIS 1 ---------
    call h5screate_simple_f(1, dims, dataspaceID, error ) 
    call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+1), &
                      H5T_NATIVE_DOUBLE, dataspaceID, &
                      datasetID, error )
    call addstringattribute( datasetID, 'TYPE', 'linear' ) 
    call addstringattribute( datasetID, 'UNITS', axis_units1 ) 
    call addstringattribute( datasetID, 'NAME', axis_name1 ) 
    call addstringattribute( datasetID, 'LONG_NAME', axis_lname1 ) 

    axis_range(1) = xmin( 1 )
    axis_range(2) = xmax( 1 )
		
    call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, axis_range, dims, &
                     error, xfer_prp = plistID )
		
    call h5dclose_f( datasetID, error )
  
    ! ---------- AXIS 2 ---------
    call h5screate_simple_f(1, dims, dataspaceID, error ) 
    call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+2), &
                      H5T_NATIVE_DOUBLE, dataspaceID, &
                      datasetID, error )
    call addstringattribute( datasetID, 'TYPE', 'linear' ) 
    call addstringattribute( datasetID, 'UNITS', axis_units2 ) 
    call addstringattribute( datasetID, 'NAME', axis_name2 ) 
    call addstringattribute( datasetID, 'LONG_NAME', axis_lname2 ) 

    axis_range(1) = xmin( 2 )
    axis_range(2) = xmax( 2 )
		
    call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, axis_range, dims, &
                     error, xfer_prp = plistID )
		
    call h5dclose_f( datasetID, error )
    ! ---------------------------
	
    call h5gclose_f( axisGroupID, error ) 
	
    ! time information 
    call adddoubleattribute( rootID, 'TIME', time ) 
    call addintegerattribute( rootID, 'ITER', iter ) 
  
    ! simulation information
    ! possibly move to a separate group
    call adddoubleattribute( rootID, 'DT', dt )
    call addstringattribute( rootID, 'TIME UNITS', '1/\omega_0' ) 
  
    call adddoublearrayattribute( rootID, 'XMIN', xmin ) 
    call adddoublearrayattribute( rootID, 'XMAX', xmax ) 
    call addintegerarrayattribute(filesave_id, 'PERIODIC', (/0,0/))
    call addintegerarrayattribute(filesave_id, 'MOVE C', (/0,0/))
    
    call addstringattribute( rootID, 'SPEC_METHOD', diag_type ) 
    
    call h5gclose_f( rootID, error )

  endif
endif
  
end subroutine setup_spec_file2d
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine setup_spec_file3d( filesavename, &
        filesave_id, myid, parallelIO, &
        axis_min1, axis_max1, axis_units1, axis_name1, axis_lname1, &
        axis_min2, axis_max2, axis_units2, axis_name2, axis_lname2, &
        axis_min3, axis_max3, axis_units3, axis_name3, axis_lname3, &
        name, time, iter, dt, xmin, xmax, diag_type )
! ----------------------------------------------------------------------------------------
  
character(len=*), intent(in) :: filesavename
integer(hid_t), intent(inout) :: filesave_id
integer, intent(in) :: myid
real(kind=p_single), intent(in) :: axis_min1, axis_max1, &
                                   axis_min2, axis_max2, &
                                   axis_min3, axis_max3
real(kind=p_double), intent(in) :: time, dt
integer, intent(in) :: iter
real(kind=p_double), dimension(:), intent(in) :: xmin, xmax
character(len=*), intent(in) :: axis_units1, axis_name1, axis_lname1, &
                                axis_units2, axis_name2, axis_lname2, &
                                axis_units3, axis_name3, axis_lname3
character(len=*), intent(in) :: name
character(len=*), intent(in) :: diag_type
    
real(p_double), dimension(2) :: axis_range
integer(hid_t) :: rootID, axisGroupID, dataspaceID, datasetID, plistID
integer(hsize_t), dimension(1) :: dims
integer :: error
  
integer, parameter :: izero = ichar('0')
logical :: parallelIO

! For serial I/O only node 0 creates the file and sets atributes;
! however for parallel I/O this must be done by all nodes.
! ---------------------- PARALLEL IO -------------------------------------------

if (parallelIO) then
   
   
  ! create the file
  call createfile(filesavename, filesave_id, myid, parallelIO)


  ! this is required for older versions of hdf5 that don't allow setting atributes
  ! for the file, only the root group
  call h5gopen_f( filesave_id, '/', rootID, error )
 
 
  ! add name property
  call addstringattribute(rootID, 'NAME', name)

  ! add file attributes
  ! create a collective write property list for axis values
  call h5pcreate_f(H5P_DATASET_XFER_F, plistID, error) 
  call h5pset_dxpl_mpio_f(plistID, H5FD_MPIO_COLLECTIVE_F, error)

  call addstringattribute( rootID, 'TYPE', 'grid' ) 
	  
  ! add axis information
  call h5gcreate_f( rootID, 'AXIS', axisGroupID, error) 
	  
  dims(1) = 2
  ! ---------- AXIS 1 ---------
  call h5screate_simple_f(1, dims, dataspaceID, error ) 
  call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+1), &
                    H5T_NATIVE_DOUBLE, dataspaceID, &
                    datasetID, error )
  call addstringattribute( datasetID, 'TYPE', 'linear' ) 
  call addstringattribute( datasetID, 'UNITS', axis_units1 ) 
  call addstringattribute( datasetID, 'NAME', axis_name1 ) 
  call addstringattribute( datasetID, 'LONG_NAME', axis_lname1 ) 

  axis_range(1) = axis_min1
  axis_range(2) = axis_max1
		
  call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, axis_range, dims, &
                   error, xfer_prp = plistID )
		
  call h5dclose_f( datasetID, error )
  
  ! ---------- AXIS 2 ---------
  call h5screate_simple_f(1, dims, dataspaceID, error ) 
  call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+2), &
                    H5T_NATIVE_DOUBLE, dataspaceID, &
                    datasetID, error )
  call addstringattribute( datasetID, 'TYPE', 'linear' ) 
  call addstringattribute( datasetID, 'UNITS', axis_units2 ) 
  call addstringattribute( datasetID, 'NAME', axis_name2 ) 
  call addstringattribute( datasetID, 'LONG_NAME', axis_lname2 ) 

  axis_range(1) = axis_min2
  axis_range(2) = axis_max2
		
  call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, axis_range, dims, &
                   error, xfer_prp = plistID )
		
  call h5dclose_f( datasetID, error )

  ! ---------- AXIS 3 ---------
  call h5screate_simple_f(1, dims, dataspaceID, error ) 
  call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+3), &
                    H5T_NATIVE_DOUBLE, dataspaceID, &
                    datasetID, error )
  call addstringattribute( datasetID, 'TYPE', 'linear' ) 
  call addstringattribute( datasetID, 'UNITS', axis_units3 ) 
  call addstringattribute( datasetID, 'NAME', axis_name3 ) 
  call addstringattribute( datasetID, 'LONG_NAME', axis_lname3 ) 

  axis_range(1) = axis_min3
  axis_range(2) = axis_max3
		
  call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, axis_range, dims, &
                   error, xfer_prp = plistID )
	
  call h5sclose_f( dataspaceID, error )
  call h5dclose_f( datasetID, error )
  ! ---------------------------
  call h5gclose_f( axisGroupID, error ) 

  ! ---------------------------
  if (myid .eq. 0) print *, "before WAXIS group"
  ! add waxis information
  call h5gcreate_f( rootID, 'WAXIS', axisGroupID, error) 
  call h5gclose_f(axisGroupID, error) 
  ! ---------------------------

  call h5pclose_f(plistID, error)
	 
  
  ! time information 
  call adddoubleattribute( rootID, 'TIME', time ) 
  call addintegerattribute( rootID, 'ITER', iter ) 
  
  ! simulation information
  ! possibly move to a separate group
  call adddoubleattribute( rootID, 'DT', dt )
  call addstringattribute( rootID, 'TIME UNITS', '1/\omega_0' ) 
  
  call adddoublearrayattribute( rootID, 'XMIN', xmin ) 
  call adddoublearrayattribute( rootID, 'XMAX', xmax ) 
  call addintegerarrayattribute(filesave_id, 'PERIODIC', (/0,0,0/))
  call addintegerarrayattribute(filesave_id, 'MOVE C', (/0,0,0/))
  
  call addstringattribute( rootID, 'SPEC_METHOD', diag_type ) 
    
  call h5gclose_f( rootID, error )

! ---------------------- SERIAL IO -------------------------------------------
else


  if ( myid == 0 ) then
   
    ! create the file
    call createfile(filesavename, filesave_id, myid, parallelIO)

    ! this is required for older versions of hdf5 that don't allow 
    ! setting atributes for the file, only the root group
    call h5gopen_f( filesave_id, '/', rootID, error )
 
    ! add name property
    call addstringattribute(rootID, 'NAME', name)

    ! add file attributes

    ! set a default write property list for axis values
    plistID = H5P_DEFAULT_F

    call addstringattribute( rootID, 'TYPE', 'grid' ) 
	  
    ! add axis information
    call h5gcreate_f( rootID, 'AXIS', axisGroupID, error) 
	  
    dims(1) = 2
    ! ---------- AXIS 1 ---------
    call h5screate_simple_f(1, dims, dataspaceID, error ) 
    call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+1), &
                      H5T_NATIVE_DOUBLE, dataspaceID, &
                      datasetID, error )
    call addstringattribute( datasetID, 'TYPE', 'linear' ) 
    call addstringattribute( datasetID, 'UNITS', axis_units1 ) 
    call addstringattribute( datasetID, 'NAME', axis_name1 ) 
    call addstringattribute( datasetID, 'LONG_NAME', axis_lname1 ) 

    axis_range(1) = xmin( 1 )
    axis_range(2) = xmax( 1 )
		
    call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, axis_range, dims, &
                     error, xfer_prp = plistID )
		
    call h5dclose_f( datasetID, error )
  
    ! ---------- AXIS 2 ---------
    call h5screate_simple_f(1, dims, dataspaceID, error ) 
    call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+2), &
                      H5T_NATIVE_DOUBLE, dataspaceID, &
                      datasetID, error )
    call addstringattribute( datasetID, 'TYPE', 'linear' ) 
    call addstringattribute( datasetID, 'UNITS', axis_units2 ) 
    call addstringattribute( datasetID, 'NAME', axis_name2 ) 
    call addstringattribute( datasetID, 'LONG_NAME', axis_lname2 ) 

    axis_range(1) = xmin( 2 )
    axis_range(2) = xmax( 2 )
		
    call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, axis_range, dims, &
                     error, xfer_prp = plistID )
		
    call h5dclose_f( datasetID, error )
    
    ! ---------- AXIS 3 ---------
    call h5screate_simple_f(1, dims, dataspaceID, error ) 
    call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+3), &
                      H5T_NATIVE_DOUBLE, dataspaceID, &
                      datasetID, error )
    call addstringattribute( datasetID, 'TYPE', 'linear' ) 
    call addstringattribute( datasetID, 'UNITS', axis_units3 ) 
    call addstringattribute( datasetID, 'NAME', axis_name3 ) 
    call addstringattribute( datasetID, 'LONG_NAME', axis_lname3 ) 

    axis_range(1) = xmin( 3 )
    axis_range(2) = xmax( 3 )
		
    call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, axis_range, dims, &
                     error, xfer_prp = plistID )
		
    call h5dclose_f( datasetID, error )
    
    ! ---------------------------
	
    call h5gclose_f( axisGroupID, error ) 
	 
    ! time information 
    call adddoubleattribute( rootID, 'TIME', time ) 
    call addintegerattribute( rootID, 'ITER', iter ) 
  
    ! simulation information
    ! possibly move to a separate group
    call adddoubleattribute( rootID, 'DT', dt )
    call addstringattribute( rootID, 'TIME UNITS', '1/\omega_0' ) 
  
    call adddoublearrayattribute( rootID, 'XMIN', xmin ) 
    call adddoublearrayattribute( rootID, 'XMAX', xmax ) 
    call addintegerarrayattribute(filesave_id, 'PERIODIC', (/0,0,0/))
    call addintegerarrayattribute(filesave_id, 'MOVE C', (/0,0,0/))
    
    call addstringattribute( rootID, 'SPEC_METHOD', diag_type ) 
    
    call h5gclose_f( rootID, error )

  endif
endif

  
end subroutine setup_spec_file3d
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************  
! ----------------------------------------------------------------------------------------
subroutine setup_ene_file2d( filesavename, &
        filesave_id, myid, parallelIO, &
        axis_min1, axis_max1, axis_units1, axis_name1, axis_lname1, &
        axis_min2, axis_max2, axis_units2, axis_name2, axis_lname2, &
        name, time, iter, dt)
! ----------------------------------------------------------------------------------------
  
character(len=*), intent(in) :: filesavename
integer(hid_t), intent(inout) :: filesave_id
integer, intent(in) :: myid
real(kind=p_single), intent(in) :: axis_min1, axis_max1
real(kind=p_single), intent(in) :: axis_min2, axis_max2
real(kind=p_double), intent(in) :: time, dt
integer, intent(in) :: iter
character(len=*), intent(in) :: axis_units1, axis_name1, axis_lname1, &
                                axis_units2, axis_name2, axis_lname2, &
                                name
!character(len=*), intent(in) :: diag_type
    
real(p_double), dimension(2) :: axis_range
integer(hid_t) :: rootID, axisGroupID, dataspaceID, datasetID, plistID
integer(hsize_t), dimension(1) :: dims
integer :: error
  
integer, parameter :: izero = ichar('0')
logical :: parallelIO


! For serial I/O only node 0 creates the file and sets atributes;
! however for parallel I/O this must be done by all nodes.
! ---------------------- PARALLEL IO -------------------------------------------

if (parallelIO) then
   
  ! create the file
  call createfile(filesavename, filesave_id, myid, parallelIO)

  ! this is required for older versions of hdf5 that don't allow setting 
  ! atributes for the file, only the root group
  call h5gopen_f( filesave_id, '/', rootID, error )
  if (error .ne. 0) print *, "error opening group"
 
  ! add name property
  call addstringattribute(rootID, 'NAME', name)

  ! add file attributes
  ! create a collective write property list for axis values
  call h5pcreate_f(H5P_DATASET_XFER_F, plistID, error) 
  call h5pset_dxpl_mpio_f(plistID, H5FD_MPIO_COLLECTIVE_F, error)

  call addstringattribute( rootID, 'TYPE', 'grid' ) 
	  
  ! add axis information
  call h5gcreate_f( rootID, 'AXIS', axisGroupID, error) 
	  
  dims(1) = 2
  ! ---------- AXIS 1 ---------
  call h5screate_simple_f(1, dims, dataspaceID, error ) 
  call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+1), &
                    H5T_NATIVE_DOUBLE, dataspaceID, &
                    datasetID, error )
  call addstringattribute( datasetID, 'TYPE', 'linear' ) 
  call addstringattribute( datasetID, 'UNITS', axis_units1 ) 
  call addstringattribute( datasetID, 'NAME', axis_name1 ) 
  call addstringattribute( datasetID, 'LONG_NAME', axis_lname1 ) 

  axis_range(1) = axis_min1
  axis_range(2) = axis_max1
		
  call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, axis_range, dims, &
                   error, xfer_prp = plistID )
		
  call h5dclose_f( datasetID, error )
  
  ! ---------- AXIS 2 ---------
  call h5screate_simple_f(1, dims, dataspaceID, error ) 
  call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+2), &
                    H5T_NATIVE_DOUBLE, dataspaceID, &
                    datasetID, error )
  call addstringattribute( datasetID, 'TYPE', 'linear' ) 
  call addstringattribute( datasetID, 'UNITS', axis_units2 ) 
  call addstringattribute( datasetID, 'NAME', axis_name2 ) 
  call addstringattribute( datasetID, 'LONG_NAME', axis_lname2 ) 

  axis_range(1) = axis_min2
  axis_range(2) = axis_max2
		
  call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, axis_range, dims, &
                   error, xfer_prp = plistID )
	
  call h5sclose_f( dataspaceID, error )
  call h5dclose_f( datasetID, error )
  ! ---------------------------
  call h5gclose_f(axisGroupID, error) 

  call h5pclose_f(plistID, error)
  
  ! time information 
  call adddoubleattribute( rootID, 'TIME', time ) 
  call addintegerattribute( rootID, 'ITER', iter ) 
  
  ! simulation information
  ! possibly move to a separate group
  call adddoubleattribute( rootID, 'DT', dt )
  call addstringattribute( rootID, 'TIME UNITS', '1/\omega_0' ) 
  
  call addfloatarrayattribute( rootID, 'XMIN', (/axis_min1,axis_min2/) ) 
  call addfloatarrayattribute( rootID, 'XMAX', (/axis_max1,axis_max2/) ) 
  call addintegerarrayattribute(filesave_id, 'PERIODIC', (/0,0/))
  call addintegerarrayattribute(filesave_id, 'MOVE C', (/0,0/))
    
  call h5gclose_f( rootID, error )

! ---------------------- SERIAL IO -------------------------------------------
else


  ! create the file
  call createfile(filesavename, filesave_id, myid, parallelIO)

  if ( myid == 0 ) then
  
    ! this is required for older versions of hdf5 that don't allow 
    ! setting atributes for the file, only the root group
    call h5gopen_f( filesave_id, '/', rootID, error )
 
    ! add name property
    call addstringattribute(rootID, 'NAME', name)

    ! add file attributes

    ! set a default write property list for axis values
    plistID = H5P_DEFAULT_F

    call addstringattribute( rootID, 'TYPE', 'grid' ) 
  
    ! add axis information
    call h5gcreate_f( rootID, 'AXIS', axisGroupID, error) 
	  
    dims(1) = 2
    ! ---------- AXIS 1 ---------
    call h5screate_simple_f(1, dims, dataspaceID, error ) 
    call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+1), &
                      H5T_NATIVE_DOUBLE, dataspaceID, &
                      datasetID, error )
    call addstringattribute( datasetID, 'TYPE', 'linear' ) 
    call addstringattribute( datasetID, 'UNITS', axis_units1 ) 
    call addstringattribute( datasetID, 'NAME', axis_name1 ) 
    call addstringattribute( datasetID, 'LONG_NAME', axis_lname1 ) 

    axis_range(1) = axis_min1
    axis_range(2) = axis_max1
		
    call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, axis_range, dims, &
                     error, xfer_prp = plistID )
		
    call h5dclose_f( datasetID, error )
  
    ! ---------- AXIS 2 ---------
    call h5screate_simple_f(1, dims, dataspaceID, error ) 
    call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+2), &
                      H5T_NATIVE_DOUBLE, dataspaceID, &
                      datasetID, error )
    call addstringattribute( datasetID, 'TYPE', 'linear' ) 
    call addstringattribute( datasetID, 'UNITS', axis_units2 ) 
    call addstringattribute( datasetID, 'NAME', axis_name2 ) 
    call addstringattribute( datasetID, 'LONG_NAME', axis_lname2 ) 

    axis_range(1) = axis_min2
    axis_range(2) = axis_max2
		
    call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, axis_range, dims, &
                     error, xfer_prp = plistID )
		
    call h5dclose_f( datasetID, error )
    ! ---------------------------
	
    call h5gclose_f( axisGroupID, error ) 
	
    ! time information 
    call adddoubleattribute( rootID, 'TIME', time ) 
    call addintegerattribute( rootID, 'ITER', iter ) 
  
    ! simulation information
    ! possibly move to a separate group
    call adddoubleattribute( rootID, 'DT', dt )
    call addstringattribute( rootID, 'TIME UNITS', '1/\omega_0' ) 
  
    call addfloatarrayattribute( rootID, 'XMIN', (/axis_min1,axis_min2/) ) 
    call addfloatarrayattribute( rootID, 'XMAX', (/axis_max1,axis_max2/) ) 
    call addintegerarrayattribute(filesave_id, 'PERIODIC', (/0,0/))
    call addintegerarrayattribute(filesave_id, 'MOVE C', (/0,0/))
    
    call h5gclose_f( rootID, error )

  endif
endif
  
end subroutine setup_ene_file2d
! ----------------------------------------------------------------------------------------
! ****************************************************************************************


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine freemem_1d_db( ptr2free ) 
! ----------------------------------------------------------------------------------------

real(kind=p_double),dimension(:),pointer :: ptr2free

! local variables  
integer :: error

if ( associated(ptr2free) ) then

  deallocate(ptr2free, stat=error)
  if (error .ne. 0) print *, "Error deallocating 1d pointer"
  
  ptr2free => null()

endif

end subroutine freemem_1d_db
! ----------------------------------------------------------------------------------------
! ****************************************************************************************   


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine freemem_2d_db( ptr2free ) 
! ----------------------------------------------------------------------------------------

real(kind=p_double),dimension(:,:),pointer :: ptr2free

! local variables  
integer :: error

if ( associated(ptr2free) ) then

  deallocate(ptr2free, stat=error)
  if (error .ne. 0) print *, "Error deallocating 2d pointer"
  
  ptr2free => null()

endif

end subroutine freemem_2d_db
! ----------------------------------------------------------------------------------------
! ****************************************************************************************      


! ****************************************************************************************
! ----------------------------------------------------------------------------------------
subroutine freemem_3d_db( ptr2free ) 
! ----------------------------------------------------------------------------------------

real(kind=p_double),dimension(:,:,:),pointer :: ptr2free

! local variables  
integer :: error

if ( associated(ptr2free) ) then

  deallocate(ptr2free, stat=error)
  if (error .ne. 0) print *, "Error deallocating 3d pointer"
  
  ptr2free => null()

endif

end subroutine freemem_3d_db
! ----------------------------------------------------------------------------------------
! ****************************************************************************************      

end module hdf5_utilities
