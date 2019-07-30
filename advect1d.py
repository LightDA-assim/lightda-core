from cffi import FFI
ffi=FFI()
ffi.cdef("""
void advect1d_step(double*,double*,double,double,int32_t,int32_t);
""")

import libsuffix
lib=ffi.dlopen('./libadvect1d'+libsuffix.suffix)

def advance(u,a,dx,dt,limiter):
    import numpy as np
    
    assert u.flags['F_CONTIGUOUS'], \
        "u is not contiguous in memory (Fortran order)"
    assert a.flags['F_CONTIGUOUS'], \
        "a is not contiguous in memory (Fortran order)"

    assert u.dtype==np.float64
    assert a.dtype==np.float64

    u_ptr=ffi.cast('double*',u.ctypes.data)
    a_ptr=ffi.cast('double*',a.ctypes.data)
    lib.advect1d_step(u_ptr,a_ptr,dx,dt,limiter,len(u))
