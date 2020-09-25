#include "mpi_types.h"

module assimilation_manager_tests

  USE dummy_model_interfaces, ONLY: dummy_model_interface, new_dummy_model
  use mod_assimilation_manager, ONLY: assimilation_manager, &
                                      new_assimilation_manager
  use random_observations, ONLY: &
    random_observation_set, new_random_observation_set
  use localization, ONLY: base_localizer
  use mod_assimilation_filter, ONLY: assimilation_filter
  use mod_dummy_model_forward_operator, ONLY: &
    dummy_model_forward_operator, new_dummy_model_forward_operator
  use system_mpi

  implicit none

contains

  subroutine test_assimilate_dummy()

    type(dummy_model_interface), target::model_interface
    type(assimilation_manager)::assim_mgr
    type(assimilation_filter)::filter
    type(dummy_model_forward_operator)::forward_operator
    type(random_observation_set)::observation_sets(1)
    type(base_localizer)::localizer
    integer::comm_size, rank, ierr
    integer, parameter::state_size = 100
    integer, parameter::batch_size = 15
    integer, parameter::istep = 1
    integer, parameter::n_ensemble = 15
    integer, parameter::n_observations = 20

    call mpi_comm_size(mpi_comm_world, comm_size, ierr)

    call mpi_comm_rank(mpi_comm_world, rank, ierr)

    ! Load observations
    observation_sets(1) = new_random_observation_set( &
                          n_observations, mpi_comm_world)

    ! Initialize i/o interface for accessing data for assimilation
    model_interface = new_dummy_model( &
                      n_ensemble, n_observations, state_size, mpi_comm_world)

    ! Initialize forward operator
    forward_operator = new_dummy_model_forward_operator(model_interface)

    assim_mgr = new_assimilation_manager( &
                model_interface, istep, n_ensemble, forward_operator, &
                observation_sets, batch_size, localizer, filter, mpi_comm_world)

    ! Run the assimilation
    call assim_mgr%assimilate(istep)

  end subroutine test_assimilate_dummy

end module assimilation_manager_tests

program test_assimilatiion_manager

  use assimilation_manager_tests
  use system_mpi

  call mpi_init(ierr)

  call test_assimilate_dummy()

  call mpi_finalize(ierr)

end program test_assimilatiion_manager
