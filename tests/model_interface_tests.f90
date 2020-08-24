module model_interface_tests
  use assimilation_model_interface, ONLY: base_model_interface
  use system_mpi
  use random_integer, ONLY: randint
  implicit none
contains

  subroutine test_localization(iface)

    class(base_model_interface)::iface
    integer::istep, subset_offset, subset_size, imodel, &
              iobs1, iobs2, state_size, n_obs, i
    real::weight

    istep = 1

    state_size = iface%get_state_size(istep)
    n_obs = iface%get_subset_obs_count(istep, 0, state_size)

    do i = 1, 100
      imodel = randint(state_size)
      iobs1 = randint(n_obs)
      iobs2 = randint(n_obs)

      weight = iface%get_weight_obs_obs(istep, iobs1, iobs2)

      if (weight > 1 .or. weight < 0) then
        print *, 'Weight', weight, 'out of range.'
        error stop
      end if

      weight = iface%get_weight_obs_obs(istep, iobs1, iobs1)

      if (weight /= 1) then
        print *, 'Weight should equal one for identical observations'
        error stop
      end if

      weight = iface%get_weight_model_obs(istep, imodel, iobs1)

      if (weight > 1 .or. weight < 0) then
        print *, 'Weight', weight, 'out of range.'
        error stop
      end if

    end do

  end subroutine test_localization

  subroutine test_io_counts(iface)
    class(base_model_interface)::iface
    real(kind=8), allocatable::buf(:)
    integer, allocatable::ranks(:), counts(:), offsets(:)
    integer::i, istep

    istep = 1

    call iface%get_io_ranks(istep, 1, ranks, counts, offsets)

    if (sum(counts) /= iface%get_state_size(istep)) then
      print *, 'Expected io counts to add up to', &
        iface%get_state_size(istep), 'got sum of', sum(counts)
      error stop
    end if

  end subroutine test_io_counts

  subroutine test_buffer_readwrite(iface)
    class(base_model_interface)::iface
    real(kind=8), allocatable::buf(:)
    integer::length, ierr, i, irank
    integer, parameter::shift = 1
    integer, allocatable::ranks(:), counts(:), offsets(:)
    integer::rank, istep

    istep = 1

    call mpi_comm_rank(mpi_comm_world, rank, ierr)

    call iface%get_io_ranks(istep, 1, ranks, counts, offsets)

    do irank = 1, size(ranks)

      if (ranks(irank) /= rank) cycle

      ! Set requested buffer length
      length = min(iface%get_state_size(istep), counts(irank), 5)

      allocate (buf(length))

      if (length > 0) then

        ! Get model state subset
        buf = iface%get_state_subset(istep, 1, offsets(irank), length)

        ! Write to buffer
        do i = 1, size(buf)
          buf(i) = i
        end do

        ! Write model state subset
        call iface%set_state_subset(istep, 1, offsets(irank), length, buf)

        ! Get new buffer, shifted relative to original
        buf = iface%get_state_subset( &
              istep, 1, offsets(irank) + shift, length - 1)

        ! Check values in new buffer
        do i = 1, size(buf) - 1
          if (buf(i) /= i + shift) then
            print *, 'Wrote value of', i + shift, &
              'read back value of', buf(i + shift)
            error stop
          end if
        end do

      end if

    end do

  end subroutine test_buffer_readwrite

  subroutine run_all(iface)
    class(base_model_interface)::iface

    call test_io_counts(iface)
    call test_buffer_readwrite(iface)

  end subroutine run_all

end module model_interface_tests
