module dummy_assimilator

  implicit none

contains

  ! Copyright (c) 2004-2018 Lars Nerger
  !
  ! This file is part of PDAF.
  !
  ! PDAF is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU Lesser General Public License
  ! as published by the Free Software Foundation, either version
  ! 3 of the License, or (at your option) any later version.
  !
  ! PDAF is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ! GNU Lesser General Public License for more details.
  !
  ! You should have received a copy of the GNU Lesser General Public
  ! License along with PDAF.  If not, see <http://www.gnu.org/licenses/>.
  !
  !$Id: PDAF-D_lenkf_analysis_rsm.F90 29 2018-03-09 20:06:46Z lnerger $
  !BOP

  ! !ROUTINE: PDAF_lenkf_analysis_rsm --- Perform LEnKF analysis step
  !
  ! !INTERFACE:
  SUBROUTINE dummy_assimilator_assimilate( &
    step, ind_p, dim_p, dim_obs_p, dim_obs, dim_ens, rank_ana, &
    state_p, ens_p, predictions, innovations, U_add_obs_err, &
    U_localize, forget, flag, info)

    ! !DESCRIPTION:
    ! Analysis step of ensemble Kalman filter with
    ! representer-type formulation.  In this version
    ! HP is explicitly computed.  This variant is
    ! optimal if the number of observations is
    ! smaller than or equal to half of the ensemble
    ! size.
    ! The final ensemble update uses a block
    ! formulation to reduce memory requirements.
    !
    ! Variant for domain decomposition.
    !
    ! !  This is a core routine of PDAF and
    !    should not be changed by the user   !
    !
    ! !REVISION HISTORY:
    ! 2003-10 - Lars Nerger - Initial code
    ! Later revisions - see svn log
    !
    ! !USES:
    ! Include definitions for real type of different precision
    ! (Defines BLAS/LAPACK routines and MPI_REALTYPE)

    USE iso_c_binding

    IMPLICIT NONE

    abstract INTERFACE
      SUBROUTINE add_obs_err(step, ind_p, dim_obs, HPH, info)
        ! Add observation error covariance matrix
        USE iso_c_binding
        INTEGER(c_int32_t), INTENT(in), value :: step, ind_p, dim_obs
        REAL(c_double), INTENT(inout) :: HPH(dim_obs, dim_obs)
        class(*), intent(in)::info
      END SUBROUTINE add_obs_err
      SUBROUTINE localize(step, ind_p, dim_p, dim_obs, HP_p, HPH, info)
        ! Apply localization to HP and HPH^T
        USE iso_c_binding
        INTEGER(c_int32_t), INTENT(in), value :: step, ind_p, dim_p, dim_obs
        REAL(c_double), INTENT(inout) :: HP_p(dim_obs, dim_p)
        REAL(c_double), INTENT(inout) :: HPH(dim_obs, dim_obs)
        class(*), intent(in)::info
      END SUBROUTINE localize

    END INTERFACE

    ! !ARGUMENTS:
    INTEGER(c_int32_t), INTENT(in), value :: step      ! Iteration number
    INTEGER(c_int32_t), INTENT(in), value :: ind_p      ! Integer index of PE
    INTEGER(c_int32_t), INTENT(in), value  :: dim_p     ! PE-local dimension of model state
    INTEGER(c_int32_t), INTENT(in), value :: dim_obs_p  ! PE-local dimension of observation vector
    INTEGER(c_int32_t), INTENT(in), value :: dim_obs    ! Global dimension of observation vector
    INTEGER(c_int32_t), INTENT(in), value :: dim_ens   ! Size of state ensemble
    INTEGER(c_int32_t), INTENT(in), value :: rank_ana  ! Rank to be considered for inversion of HPH
    REAL(c_double), INTENT(inout)  :: state_p(dim_p)        ! PE-local ensemble mean state
    REAL(c_double), INTENT(inout)  :: ens_p(dim_p, dim_ens) ! PE-local state ensemble
    REAL(c_double), INTENT(inout)  :: predictions(dim_obs, dim_ens) ! PE-local state ensemble
    REAL(c_double), INTENT(in)     :: innovations(dim_obs, dim_ens) ! Global array of innovations
    REAL(c_double), INTENT(in), value     :: forget    ! Forgetting factor
    INTEGER(c_int32_t), INTENT(out) :: flag    ! Status flag
    class(*), intent(inout)::info

    procedure(add_obs_err) :: U_add_obs_err
    procedure(localize) :: U_localize

  END SUBROUTINE dummy_assimilator_assimilate

end module dummy_assimilator
