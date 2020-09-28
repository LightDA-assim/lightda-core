module advect1d_model_interface_tests

  use advect1d_assimilate_interfaces, ONLY: advect1d_interface
  use random_integer, ONLY: randint
  use hdf5
  use exceptions, ONLY: error_status, throw, new_exception
  use hdf5_exceptions, ONLY: new_hdf5_exception
  use util, only: str

  implicit none

  integer, parameter::istep = 1

contains

  subroutine run_all_advect1d(iface)
    class(advect1d_interface)::iface

  end subroutine run_all_advect1d

  function generate_ensemble(istep, state_size, n_ensemble) result(ensemble)

    integer, allocatable::obs_positions(:)
    integer, intent(in)::istep, state_size, n_ensemble
    real(kind=8)::ensemble(n_ensemble, state_size)

    call random_number(ensemble)

  end function generate_ensemble

  subroutine write_ensemble(istep, ensemble, status)

    use iso_fortran_env, ONLY: output_unit

    !! Write the ensemble to disk

    ! Arguments
    integer, intent(in)::istep
        !! Iteration number
    real(kind=8)::ensemble(:, :)
        !! ensemble state
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    integer::state_size, n_ensemble
    character(len=80)::preassim_filename
    integer(HID_T)::h5file_h, dset_h, dataspace, memspace
    integer(HSIZE_T)::offset(2), count(2), stride(2), data_block(2)
    logical::exists
    integer::ierr

    n_ensemble = size(ensemble, 1)
    state_size = size(ensemble, 2)

    stride = (/1, 1/)
    offset = (/0, 0/)
    count = (/1, 1/)
    data_block = (/n_ensemble, state_size/)

    call h5open_f(ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception(ierr, &
                                            'HDF5 initialization error', &
                                            procedure='write_ensemble'))
      return
    end if

    ! Set the HDF5 filename
    ! Set the HDF5 filename
    write (preassim_filename, "(A,I0,A)") &
      'ensembles/', istep, '/preassim.h5'

    call system('mkdir -p ensembles/'//str(istep))

    ! Create the file
    call h5fcreate_f(preassim_filename, H5F_ACC_TRUNC_F, h5file_h, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception(ierr, &
                                            'HDF5 error creating output file', &
                                            filename=preassim_filename, &
                                            procedure='write_ensemble'))
      return
    end if

    ! Define a dataspace
    call h5screate_simple_f(2, data_block, dataspace, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception(ierr, &
                                            'HDF5 error creating dataspace', &
                                            filename=preassim_filename, &
                                            procedure='write_ensemble'))
      return
    end if

    ! Create the dataset
    call h5dcreate_f(h5file_h, 'ensemble_state', H5T_IEEE_F64LE, &
                     dataspace, dset_h, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception(ierr, &
                                            'HDF5 error creating dataset', &
                                            filename=preassim_filename, &
                                            procedure='write_ensemble'))
      return
    end if

    ! Write the data
    call h5dwrite_f( &
      dset_h, H5T_IEEE_F64LE, ensemble, data_block, ierr, &
      file_space_id=dataspace)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception(ierr, &
                                            'HDF5 error writing data', &
                                            filename=preassim_filename, &
                                            procedure='write_ensemble'))
      return
    end if

    ! Close the dataspace
    call h5sclose_f(dataspace, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception(ierr, &
                                            'HDF5 error closing dataspace', &
                                            filename=preassim_filename, &
                                            procedure='write_ensemble'))
      return
    end if

    ! Close the dataset
    call h5dclose_f(dset_h, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception(ierr, &
                                            'HDF5 error closing dataset', &
                                            filename=preassim_filename, &
                                            procedure='write_ensemble'))
      return
    end if

    ! Close the file
    call h5fclose_f(h5file_h, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception(ierr, &
                                            'HDF5 error closing file', &
                                            filename=preassim_filename, &
                                            procedure='write_ensemble'))
      return
    end if

  end subroutine write_ensemble

end module advect1d_model_interface_tests

program test_advect1d_model_interface
  use model_interface_tests, ONLY: run_all
  use advect1d_assimilate_interfaces, ONLY: advect1d_interface, &
                                            new_advect1d_interface
  use advect1d_model_interface_tests, ONLY: &
    write_ensemble, generate_ensemble, run_all_advect1d
  use system_mpi
  implicit none
  type(advect1d_interface)::model_interface
  integer, parameter::istep = 1
  integer, parameter::n_observations = 5
  integer, parameter::state_size = 100
  integer, parameter::n_ensemble = 15
  real(kind=8)::ensemble(n_ensemble, state_size)
  integer::rank, ierr

  call mpi_init(ierr)

  call mpi_comm_rank(mpi_comm_world, rank, ierr)

  ensemble = generate_ensemble(istep, state_size, n_ensemble)

  call MPI_Bcast(ensemble, state_size*n_ensemble, MPI_DOUBLE_PRECISION, &
                 0, mpi_comm_world, ierr)

  if (rank == 0) call write_ensemble(istep, ensemble)

  call mpi_barrier(mpi_comm_world, ierr)

  model_interface = new_advect1d_interface( &
                    istep, n_ensemble, state_size, mpi_comm_world)

  call run_all(model_interface)
  call run_all_advect1d(model_interface)

  call mpi_finalize(ierr)
end program test_advect1d_model_interface
