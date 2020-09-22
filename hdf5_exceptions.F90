module hdf5_exceptions

  use hdf5
  use exceptions, ONLY: exception, error_status, no_error

  implicit none

  type, extends(exception)::hdf5_exception
    integer::code
    character(:), allocatable::filename
  contains
    procedure::print => hdf5_print_error

    ! Some versions of gfortran fail to compile polymorphic code with
    ! overriden finalizers, so we use a macro to control whether the finalizer
    ! is used
#ifdef OVERRIDABLE_FINALIZERS
    final::hdf5_exception_finalize
#endif
  end type hdf5_exception

contains

  function new_hdf5_exception(code, message, filename, procedure) &
    result(hdf_result)

    !! Create a new HDF5 exception, or a no_error object if the error code is 0.

    ! Arguments
    integer, intent(in)::code
        !! HDF5 error code
    character(*), intent(in), optional::message
        !! Custom error message
    character(*), intent(in), optional::procedure
        !! Procedure where error occured
    character(*), intent(in), optional::filename

    ! Result
    class(error_status), allocatable::hdf_result

    type(no_error)::success    ! No-error object in case the code is 0

    type(hdf5_exception)::exc  ! HDF5 exception object

    if (code == 0) then
      ! No error
      hdf_result = success
      return
    end if

    exc%code = code
    if (present(procedure)) exc%procedure = procedure
    if (present(filename)) exc%filename = filename

    if (present(message)) then
      exc%message = message
    else
      exc%message = 'HDF5 exception'
    end if

#ifndef OVERRIDABLE_FINALIZERS
    if (present(filename)) then
      exc%message = exc%message//' on file '//filename
    end if
#endif

    hdf_result = exc

  end function new_hdf5_exception

  subroutine hdf5_print_error(this)

    !! Print the HDF5 exception information

    use iso_fortran_env, ONLY: error_unit

    ! Arguments
    class(hdf5_exception), intent(in)::this
        !! Exception

    integer::ierr
        !! HDF5 error code

    write (error_unit, *) this%message

    if (allocated(this%filename)) write (error_unit, *) 'on file '//this%filename

    call h5eprint_f(ierr)

  end subroutine hdf5_print_error

  subroutine hdf5_exception_finalize(this)

    !! Finalize an exception. If this%handled is false, print error message
    !! and terminate. Otherwise do nothing.

    use iso_fortran_env, ONLY: error_unit

    ! Arguments
    type(hdf5_exception), intent(inout)::this
        !! Exception object

    call this%default_handler()

  end subroutine hdf5_exception_finalize

end module hdf5_exceptions
