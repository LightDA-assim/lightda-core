module advect1d_model_interface_tests

  use advect1d_assimilate_interfaces, ONLY: advect1d_interface
  use random_integer, ONLY: randint
  use hdf5

  implicit none

  integer, parameter::istep = 1

contains

  subroutine write_observations(observations, obs_errors, obs_positions)
    real(kind=8), intent(in)::observations(:), obs_errors(:)
    integer, intent(in)::obs_positions(:)
    character(len=50)::obs_filename
    integer(HID_T)::h5file_h, dset_h, dataspace
    integer(HSIZE_T)::data_block(2)
    integer::ierr
    character(len=10)::istepstr

    data_block = (/1, size(observations)/)

    ! Convert istep to string
    write (istepstr, '(I10)') istep

    ! Create parent directory
    call execute_command_line('mkdir -p ensembles/'//adjustl(trim(istepstr)))

    ! Set the HDF5 filename
    write (obs_filename, "(A,I0,A)") &
      'ensembles/', istep, '/observations.h5'

    ! Delete any preexisting file
    call execute_command_line('rm -f '//adjustl(trim(obs_filename)))

    ! Create the file
    call h5fcreate_f(obs_filename, H5F_ACC_TRUNC_F, h5file_h, ierr)

    if (ierr /= 0) error stop

    ! Define a dataspace for the observations
    call h5screate_simple_f(1, (/int(size(observations), 8)/), dataspace, ierr)

    if (ierr /= 0) error stop

    ! Create the observations dataset
    call h5dcreate_f(h5file_h, 'observations', H5T_IEEE_F64LE, &
                     dataspace, dset_h, ierr)

    if (ierr /= 0) error stop

    ! Write the data
    call h5dwrite_f(dset_h, H5T_NATIVE_DOUBLE, observations, &
                    data_block, ierr, file_space_id=dataspace)

    if (ierr /= 0) error stop

    ! Close the dataset
    call h5dclose_f(dset_h, ierr)

    if (ierr /= 0) error stop

    ! Create the obs_errors dataset
    call h5dcreate_f(h5file_h, 'obs_errors', H5T_IEEE_F64LE, &
                     dataspace, dset_h, ierr)

    if (ierr /= 0) error stop

    ! Write the data
    call h5dwrite_f( &
      dset_h, H5T_NATIVE_DOUBLE, obs_errors, data_block, &
      ierr, file_space_id=dataspace)

    if (ierr /= 0) error stop

    ! Close the dataset
    call h5dclose_f(dset_h, ierr)

    ! Create the obs_positions dataset
    call h5dcreate_f(h5file_h, 'obs_positions', H5T_NATIVE_INTEGER, &
                     dataspace, dset_h, ierr)

    if (ierr /= 0) error stop

    ! Write the data
    call h5dwrite_f(dset_h, H5T_NATIVE_INTEGER, obs_positions, &
                    data_block, ierr, file_space_id=dataspace)

    if (ierr /= 0) error stop

    ! Close the dataset
    call h5dclose_f(dset_h, ierr)

    ! Close the dataspace
    call h5sclose_f(dataspace, ierr)

    ! Close the file
    call h5fclose_f(h5file_h, ierr)

  end subroutine write_observations

  subroutine generate_inputs(iface)

    class(advect1d_interface)::iface
    real(kind=8), allocatable::observations(:), obs_errors(:)
    integer, allocatable::obs_positions(:)
    integer::state_size, n_obs, i

    state_size = iface%get_state_size(istep)
    n_obs = iface%get_subset_obs_count(istep, 0, state_size)

    allocate (observations(n_obs), obs_errors(n_obs), obs_positions(n_obs))

    call random_number(observations)
    call random_number(obs_errors)

    do i = 1, n_obs
      obs_positions(i) = randint(state_size/2)
    end do

    call write_observations(observations, obs_errors, obs_positions)

  end subroutine generate_inputs

  subroutine test_localization(iface)
    class(advect1d_interface)::iface
    integer::istep, subset_offset, subset_size, imodel, &
              iobs1, iobs2, state_size, n_obs
    integer, allocatable::obs_positions(:)
    integer::i, model_pos, obs_pos1, obs_pos2, domain_size
    real::weight

    state_size = iface%get_state_size(istep)
    n_obs = iface%get_subset_obs_count(istep, 0, state_size)

    domain_size = state_size/2

    allocate (obs_positions(n_obs))

    obs_positions = iface%get_obs_positions()

    do i = 1, 100
      imodel = randint(state_size)
      iobs1 = randint(n_obs)
      iobs2 = randint(n_obs)

      model_pos = mod(imodel - 1, domain_size)
      obs_pos1 = obs_positions(iobs1)
      obs_pos2 = obs_positions(iobs2)

      weight = iface%get_weight_obs_obs(istep, iobs1, iobs2)

      if (obs_pos1 /= obs_pos2 .and. weight >= 1) then
        print *, 'Weight should be less than 1 for differently &
             &placed observations'
        error stop
      end if

      weight = iface%get_weight_model_obs(istep, imodel, iobs1)
      if (model_pos /= obs_pos1 .and. &
          abs(model_pos - obs_pos1) /= domain_size .and. weight == 1) then
        print *, 'Weight should be less than 1 when model position &
             &is different from observation position'
        error stop
      end if

      weight = iface%get_weight_model_obs(istep, obs_pos1 + 1, iobs1)

      if (weight /= 1) then
        print *, 'Weight should equal one when model position &
             &is the same as observing position'
        error stop
      end if

    end do

  end subroutine test_localization

  subroutine run_all_advect1d(iface)
    class(advect1d_interface)::iface

    call test_localization(iface)

  end subroutine run_all_advect1d

end module advect1d_model_interface_tests

program test_advect1d_model_interface
  use model_interface_tests, ONLY: run_all
  use advect1d_assimilate_interfaces, ONLY: advect1d_interface, &
                                            new_advect1d_interface
  use advect1d_model_interface_tests, ONLY: generate_inputs, run_all_advect1d
  use system_mpi
  implicit none
  type(advect1d_interface)::model_interface
  integer, parameter::n_observations = 5
  integer, parameter::state_size = 100
  integer, parameter::n_ensemble = 15
  integer::rank, ierr

  call mpi_init(ierr)

  model_interface = new_advect1d_interface( &
                    n_ensemble, n_observations, state_size, mpi_comm_world)

  call mpi_comm_rank(mpi_comm_world, rank, ierr)

  if (rank == 0) call generate_inputs(model_interface)

  call mpi_barrier(mpi_comm_world, ierr)

  call model_interface%read_observations(1)

  call run_all(model_interface)
  call run_all_advect1d(model_interface)

  call mpi_finalize(ierr)
end program test_advect1d_model_interface
