module advect1d_model_interface_tests

  use advect1d_assimilate_interfaces,ONLY:advect1d_interface
  use random_integer, ONLY: randint
  implicit none

  integer, parameter::istep=1

contains

  subroutine test_localization()
    type(advect1d_interface)::iface
    integer::istep,subset_offset,subset_size,imodel,iobs1,iobs2,state_size,n_obs
    integer,allocatable::obs_positions(:)
    integer::i,model_pos,obs_pos1,obs_pos2
    real::weight

    state_size=iface%get_state_size()
    n_obs=iface%get_subset_obs_count(istep,0,state_size)

    allocate(obs_positions(n_obs))

    obs_positions=iface%get_obs_positions()

    do i=1,100
       imodel=randint(state_size)
       iobs1=randint(n_obs)
       iobs2=randint(n_obs)

       model_pos=mod(imodel,state_size/2)
       obs_pos1=obs_positions(iobs1)
       obs_pos2=obs_positions(iobs2)

       weight=iface%get_weight_obs_obs(istep,iobs1,iobs2)

       if(obs_pos1/=obs_pos2 .and. weight>=1) then
          print *,'Weight should be less than 1 for differently placed observations'
          error stop
       end if

       weight=iface%get_weight_model_obs(istep,imodel,iobs1)
       if(weight==1) then
          print *,'Weight be less than 1 when model position is different from observation position'
          error stop
       end if

       weight=iface%get_weight_model_obs(istep,obs_pos1,iobs1)

       if(weight/=1) then
          print *,'Weight should equal one when model position is the same as observing position'
          error stop
       end if

    end do

  end subroutine test_localization

end module advect1d_model_interface_tests

program test_advect1d_model_interface
  use model_interface_tests, ONLY: run_all
  use advect1d_assimilate_interfaces,ONLY:advect1d_interface,new_advect1d_interface
  use system_mpi
  implicit none
  type(advect1d_interface)::model_interface
  integer,parameter::n_observations=5
  integer,parameter::state_size=100
  integer,parameter::n_ensemble=15
  integer::ierr

  call mpi_init(ierr)

  model_interface=new_advect1d_interface(n_ensemble,n_observations,state_size,mpi_comm_world)

  call run_all(model_interface)
  call mpi_finalize(ierr)
end program test_advect1d_model_interface
