from cffi import FFI
ffi=FFI()
ffi.cdef("""
void lenkf_rsm(double*,double*,double*,double*,double*,double*,double*,int32_t,int32_t,int32_t);
""")

import libsuffix
lib=ffi.dlopen('./libenkf'+libsuffix.suffix)

def lenkf_rsm(ensemble_state,forward_operator,observations,obs_errors,localization_obs_obs=None,localization_obs_model=None):
    import numpy as np
    
    assert ensemble_state.flags['C_CONTIGUOUS'], \
        "ensemble_state is not contiguous in memory (C order)"
    assert forward_operator.flags['C_CONTIGUOUS'], \
        "forward_operator is not contiguous in memory (C order)"
    assert observations.flags['C_CONTIGUOUS'], \
        "observations is not contiguous in memory (C order)"
    assert obs_errors.flags['C_CONTIGUOUS'], \
        "obs_errors is not contiguous in memory (C order)"

    assert ensemble_state.dtype==np.float64
    assert forward_operator.dtype==np.float64
    assert observations.dtype==np.float64
    assert obs_errors.dtype==np.float64

    model_size,n_ensemble=ensemble_state.shape
    n_observations=len(observations)

    assert len(obs_errors)==n_observations
    assert forward_operator.shape==(n_observations,model_size)

    if localization_obs_obs is None:
        localization_obs_obs=np.ones([n_observations,n_observations])
    assert localization_obs_obs.dtype==np.float64
    assert localization_obs_obs.shape==(n_observations,n_observations)
    assert localization_obs_obs.flags['C_CONTIGUOUS']

    if localization_obs_model is None:
        localization_obs_model=np.ones([n_observations,model_size])
    assert localization_obs_model.dtype==np.float64
    assert localization_obs_model.shape==(n_observations,model_size)
    assert localization_obs_model.flags['C_CONTIGUOUS']

    new_ensemble_state=np.empty((model_size,n_ensemble))
    new_ensemble_state_ptr=ffi.cast('double*',new_ensemble_state.ctypes.data)
    ensemble_state_ptr=ffi.cast('double*',ensemble_state.ctypes.data)
    forward_operator_ptr=ffi.cast('double*',forward_operator.ctypes.data)
    observations_ptr=ffi.cast('double*',observations.ctypes.data)
    obs_errors_ptr=ffi.cast('double*',obs_errors.ctypes.data)
    localization_obs_obs_ptr=ffi.cast('double*',localization_obs_obs.ctypes.data)
    localization_obs_model_ptr=ffi.cast('double*',localization_obs_model.ctypes.data)
    lib.enkf_analysis_from_obs_wrapper(
        new_ensemble_state_ptr,
        ensemble_state_ptr,
        forward_operator_ptr,
        observations_ptr,
        obs_errors_ptr,
        localization_obs_obs_ptr,
        localization_obs_model_ptr,
        model_size,
        n_ensemble,
        n_observations
    )

    return new_ensemble_state
