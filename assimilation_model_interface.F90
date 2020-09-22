#include "mpi_types.h"

module assimilation_model_interface

  use exceptions, ONLY: throw, new_exception, error_status
  use system_mpi

  implicit none

  type, abstract::base_model_interface

    !! Base class for model interfaces

    integer::n_ensemble
    MPI_COMM_TYPE::comm
  contains
    procedure(I_get_state_size), deferred::get_state_size
    procedure(I_set_state_subset), deferred::set_state_subset
    procedure(I_get_state_darray), deferred::get_state_darray
    procedure::read_state
    procedure::write_state
  end type base_model_interface

  abstract interface

    function I_get_state_darray(this, istep, imember, status) &
      result(state_darray)

       !! Get the requested ensemble member state as a darray

      use exceptions, ONLY: error_status
      use distributed_array, ONLY: darray
      import base_model_interface

      ! Arguments
      class(base_model_interface)::this
           !! Model interface
      integer, intent(in)::istep
           !! Iteration number
      integer, intent(in)::imember
           !! Ensemble member index
      class(error_status), intent(out), allocatable, optional::status
           !! Error status

      ! Returns
      type(darray)::state_darray
           !! State array represented as a darray object

    end function I_get_state_darray

    function I_get_state_size(this, istep, status) result(size)

       !! Returns the number elements in the model state array

      use exceptions, ONLY: error_status
      import base_model_interface

      implicit none

      ! Arguments
      class(base_model_interface)::this
           !! Model interface
      integer, intent(in)::istep
           !! Iteration number
      class(error_status), intent(out), allocatable, optional::status
           !! Error status

      ! Returns
      integer::size
           !! State size

    end function I_get_state_size

    subroutine I_set_state_subset( &
      this, istep, imember, subset_offset, subset_size, subset_state, status)

       !! Returns model state values for a given subset of the model state
       !! Note that ability to write data may depend on the rank of the calling
       !! processor, so this subroutine should only be called after a call to
       !! get_io_ranks, and the values of subset_offset and subset_size should
       !! be chosen to fall within a continguous segment of data stored on the
       !! calling processor rank as determined by the output from get_io_ranks.
       !! Failure to satisfy this requirement may result in an error.

      use exceptions, ONLY: error_status
      import base_model_interface

      implicit none

      ! Arguments
      class(base_model_interface)::this
           !! Model interface
      integer, intent(in)::istep
           !! Iteration number
      integer, intent(in)::imember
           !! Ensemble member index
      integer, intent(in)::subset_offset
           !! Offset of subset from start of state array
      integer, intent(in)::subset_size
           !! Size of subset
      real(kind=8), intent(in)::subset_state(subset_size)
           !! Values of the model state in the requested subset
      class(error_status), intent(out), allocatable, optional::status
           !! Error status
    end subroutine I_set_state_subset

  end interface

contains

  subroutine read_state(this, istep, status)

    !! Load the model state from disk for a given iteration `istep`. This
    !! subroutine is to be called before any calls to get_state_subset for the
    !! specified iteration `istep`. The default implementation is a no-op,
    !! enabling a model interface to load the model state incrementally as
    !! ensemble state is requested.

    ! Arguments
    class(base_model_interface)::this
        !! Model interface
    integer, intent(in)::istep
        !! Iteration number
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

  end subroutine read_state

  subroutine write_state(this, istep, status)

    !! Record the model state after assimilation of `istep`. This
    !! subroutine is to be called after all calls to set_state_subset have
    !! completed for the specified iteration `istep`. The default
    !! implementation is a no-op, enabling a model interface to load the model
    !! state incrementally as ensemble state is received from the assimilation
    !! workers.

    ! Arguments
    class(base_model_interface)::this
        !! Model interface
    integer, intent(in)::istep
        !! Iteration number
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

  end subroutine write_state

end module assimilation_model_interface
