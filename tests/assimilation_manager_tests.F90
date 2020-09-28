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
  use distributed_array, ONLY: darray, darray_segment
  use exceptions, ONLY: error_status, throw, new_exception
  use util, ONLY: str

  implicit none

contains

  subroutine test_assimilate_dummy(status)

    class(error_status), intent(out), allocatable, optional::status
        !! Error status

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
    type(darray), target::state_before_filtering(n_ensemble)
    type(darray), target::state_after_filtering(n_ensemble)
    type(darray_segment), pointer:: &
      segment_before_filtering, segment_after_filtering
        !! Pointers to darray segments

    integer::imember, isegment ! Loop counters

    call mpi_comm_size(mpi_comm_world, comm_size, ierr)

    call mpi_comm_rank(mpi_comm_world, rank, ierr)

    ! Load observations
    observation_sets(1) = new_random_observation_set( &
                          n_observations, mpi_comm_world)

    ! Initialize i/o interface for accessing data for assimilation
    model_interface = new_dummy_model( &
                      n_ensemble, n_observations, state_size, mpi_comm_world)

    do imember = 1, n_ensemble
      state_before_filtering(imember) = model_interface%get_state_darray( &
                                        istep, imember)
    end do

    ! Initialize forward operator
    forward_operator = new_dummy_model_forward_operator(model_interface)

    assim_mgr = new_assimilation_manager( &
                model_interface, istep, n_ensemble, forward_operator, &
                observation_sets, batch_size, localizer, filter, mpi_comm_world)

    ! Run the assimilation
    call assim_mgr%assimilate(istep)

    do imember = 1, n_ensemble
      state_after_filtering(imember) = model_interface%get_state_darray( &
                                       istep, imember)
    end do

    do imember = 1, n_ensemble

      if (size(state_before_filtering(imember)%segments) /= &
          size(state_after_filtering(imember)%segments)) then
        call throw( &
          status, new_exception( &
          'Number of state segments changed after filtering on member '// &
          str(imember)//'. Was '// &
          str(size(state_before_filtering(imember)%segments))// &
          ', now '//str(size(state_after_filtering(imember)%segments)), &
          'test_assimilate_dummy'))
        return
      end if

      do isegment = 1, size(state_before_filtering(imember)%segments)
        segment_before_filtering => &
          state_before_filtering(imember)%segments(isegment)
        segment_after_filtering => &
          state_after_filtering(imember)%segments(isegment)

        if (segment_before_filtering%length /= &
            segment_after_filtering%length) then
          call throw( &
            status, new_exception( &
            'Size changed for segment ' &
            //str(isegment)//' on member '//str(imember)//'. Was '// &
            str(segment_before_filtering%length)// &
            ', now '//str(segment_after_filtering%length), &
            'test_assimilate_dummy'))
          return
        end if

        if (segment_before_filtering%rank /= segment_after_filtering%rank) then
          call throw( &
            status, &
            new_exception( &
            'Rank changed for segment ' &
            //str(isegment)//' on member '//str(imember)//'. Was '// &
            str(segment_before_filtering%rank)// &
            ', now '//str(segment_after_filtering%rank), &
            'test_assimilate_dummy'))
          return
        end if

        if (segment_before_filtering%rank /= segment_after_filtering%rank) then
          call throw( &
            status, &
            new_exception('Rank changed for segment ' &
                          //str(isegment)//' on member '//str(imember)// &
                          '. Was '//str(segment_before_filtering%rank)// &
                          ', now '//str(segment_after_filtering%rank), &
                          'test_assimilate_dummy'))
          return
        end if

        if (segment_before_filtering%rank == rank) then

          if (size(segment_before_filtering%data) /= &
              segment_before_filtering%length) then
            call throw( &
                 status, new_exception( &
                 'Wrong buffer size for segment ' &
                 //str(isegment)//' on member '//str(imember)//'. &
                 &Expected '//str(segment_before_filtering%length)// &
                 ', got '//str(size(segment_before_filtering%data)), &
                 'test_assimilate_dummy'))
            return
          end if

          if (size(segment_after_filtering%data) /= &
              segment_after_filtering%length) then
            call throw( &
                 status, new_exception( &
                 'Wrong buffer size for segment ' &
                 //str(isegment)//' on member '//str(imember)//'. &
                 &Expected '//str(segment_after_filtering%length)// &
                 ', got '//str(size(segment_after_filtering%data)), &
                 'test_assimilate_dummy'))
            return
          end if

          if (any(segment_before_filtering%data /= &
                  segment_after_filtering%data)) &
            then
            call throw(status, new_exception('Data changed for segment ' &
                                             //str(isegment)//' on member '// &
                                             str(imember)//'.', &
                                             'test_assimilate_dummy'))
            return
          end if
        end if

      end do

    end do

  end subroutine test_assimilate_dummy

end module assimilation_manager_tests

program test_assimilatiion_manager

  use assimilation_manager_tests
  use system_mpi

  call mpi_init(ierr)

  call test_assimilate_dummy()

  call mpi_finalize(ierr)

end program test_assimilatiion_manager
