module observations

  implicit none

  type, abstract::observation_set
   contains
     procedure(get_size), deferred::get_size
     procedure(get_values), deferred::get_values
  end type observation_set

  abstract interface

     function get_size(this,status) result(size)

       use exceptions, ONLY: error_status

       import observation_set

       ! Arguments
       class(observation_set)::this
           !! Observation set
       class(error_status), intent(out), allocatable, optional::status
           !! Error status

       ! Result
       integer::size
           !! Number of observation values in set

     end function get_size

     function get_values(this,status) result(values)

       use exceptions, ONLY: error_status

       import observation_set

       ! Arguments
       class(observation_set)::this
           !! Observation set
       class(error_status), intent(out), allocatable, optional::status
           !! Error status

       ! Result
       real(kind=8), allocatable::values(:)
           !! Observation values

     end function get_values

  end interface

end module observations
