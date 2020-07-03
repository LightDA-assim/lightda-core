from cffi import FFI
ffi=FFI()
ffi.cdef("""
typedef void (*U_localize)(int,int,int,int,double*,double*);
typedef void (*U_add_obs_err)(int,int,int,double*);

void lenkf_analysis_rsm_c(int32_t,int32_t,int32_t,int32_t,int32_t,int32_t,int32_t,double*,double*,double*,double*,U_add_obs_err,U_localize, double,int32_t*);

void enkf_analysis_from_obs_wrapper(double*,double*,double*,double*,double*,double*,double*,int32_t,int32_t,int32_t);
void enkf_analysis_from_innovations_wrapper(double*,double*,double*,double*,double*,double*,int32_t,int32_t,int32_t);
""")
import numpy as np

import libsuffix
lib=ffi.dlopen('./libassim'+libsuffix.suffix)

ctype2dtype = {}

# Integer types
for prefix in ('int', 'uint'):
    for log_bytes in range(4):
        ctype = '%s%d_t' % (prefix, 8 * (2**log_bytes))
        dtype = '%s%d' % (prefix[0], 2**log_bytes)
        # print( ctype )
        # print( dtype )
        ctype2dtype[ctype] = np.dtype(dtype)

# Floating point types
ctype2dtype['float'] = np.dtype('f4')
ctype2dtype['double'] = np.dtype('f8')

def ptr_to_array(ptr,shape,**kwargs):
    import numpy as np
    length = np.prod(shape)
    # Get the canonical C type of the elements of ptr as a string.
    T = ffi.getctype(ffi.typeof(ptr).item)
    # print( T )
    # print( ffi.sizeof( T ) )

    if T not in ctype2dtype:
        raise RuntimeError("Cannot create an array for element type: %s" % T)

    a = np.frombuffer(ffi.buffer(ptr, length * ffi.sizeof(T)), ctype2dtype[T])\
          .reshape(shape, **kwargs)
    return a

def lenkf_rsm(step,ind_p,ensemble_state,predictions,innovations,add_obs_err,localization,forget):

    import numpy as np
    
    assert ensemble_state.flags['F_CONTIGUOUS'], \
        "ensemble_state is not contiguous in memory (F order)"
    assert predictions.flags['F_CONTIGUOUS'], \
        "resid is not contiguous in memory (F order)"
    assert innovations.flags['F_CONTIGUOUS'], \
        "resid is not contiguous in memory (F order)"

    assert ensemble_state.dtype==np.float64
    assert predictions.dtype==np.float64
    assert innovations.dtype==np.float64

    model_size,n_ensemble=ensemble_state.shape
    n_observations,n_ensemble_resid=predictions.shape

    assert(n_ensemble==n_ensemble_resid)

    assert(innovations.shape==predictions.shape)
    
    mean_state=np.empty([model_size])
    
    ensemble_state_ptr=ffi.cast('double*',ensemble_state.ctypes.data)
    predictions_ptr=ffi.cast('double*',predictions.ctypes.data)
    innovations_ptr=ffi.cast('double*',innovations.ctypes.data)
    mean_state_ptr=ffi.cast('double*',mean_state.ctypes.data)

    flag_ptr=ffi.new('int32_t*')

    @ffi.callback("void(*)(int,int,int,int,double*,double*)")
    def localization_cb(step,ind_p,dim_p,dim_obs,HP_p_ptr,HPH_ptr):
        HP_p=ptr_to_array(HP_p_ptr,(dim_obs,dim_p),order='F')
        HPH=ptr_to_array(HPH_ptr,(dim_obs,dim_obs))
        localization(step,ind_p,HP_p,HPH)

    @ffi.callback("void(*)(int,int,int,double*)")
    def add_obs_err_cb(step,ind_p,dim_obs,HPH_ptr):
        HPH=ptr_to_array(HPH_ptr,(dim_obs,dim_obs),order='F')
        add_obs_err(step,ind_p,HPH)

    lib.lenkf_analysis_rsm_c(
        step,
        ind_p,
        model_size,
        n_observations,
        n_observations,
        n_ensemble,
        0,
        mean_state_ptr,
        ensemble_state_ptr,
        predictions_ptr,
        innovations_ptr,
        add_obs_err_cb,
        localization_cb,
        forget,
        flag_ptr
    )

    return ensemble_state
