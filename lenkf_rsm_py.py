import numpy as np

from cffi import FFI
ffi=FFI()
ffi.cdef("""
void compute_residual(int32_t,int32_t,double*,double*);
""")
import numpy as np

import libsuffix
lib=ffi.dlopen('./libassim'+libsuffix.suffix)

def compute_residual(a):

    assert a.flags['F_CONTIGUOUS'], \
        "ensemble_state is not contiguous in memory (F order)"

    n,m=a.shape
    resid=np.asfortranarray(np.empty(a.shape))
    a_ptr=ffi.cast('double*',a.ctypes.data)
    resid_ptr=ffi.cast('double*',resid.ctypes.data)
    lib.compute_residual(n,m,a_ptr,resid_ptr)

    return resid

def call_dgemm(transa,transb,alpha,beta,a,b):

    from dgemm import dgemm

    a=np.asfortranarray(a)
    b=np.asfortranarray(b)

    if transa in ('n','N'):
        m,k=a.shape
        lda=m
    else:
        k,m=a.shape
        lda=k

    if transb in ('n','N'):
        k,n=b.shape
        ldb=k
    else:
        n,k=b.shape
        ldb=n

    c=np.asfortranarray(np.empty([m,n]))
    ldc=m

    dgemm('n','t',m,n,k,alpha,a,b,beta,c,lda,ldb,ldc)

    return c

def analysis_with_pseudoinverse(HP,HPH,ens,innovations,rank):

    from scipy.linalg import eigh
    eigvals,eigvecs=eigh(HPH,lower=False)
    n_obs=innovations.shape[0]
    iupper=n_obs
    ilower=n_obs-rank
    eigvals=eigvals[ilower:iupper]
    eigvecs=eigvecs[:,-rank:]
    eigvecs_scaled=eigvecs/eigvals
    HPH_inv=np.dot(eigvecs,eigvecs_scaled.T)
    repres=np.dot(HPH_inv,innovations)
    ens_new=ens+np.dot(HP.T,repres)
    return ens_new

def analysis_direct_solve(HP,HPH,ens,innovations):

    repres=np.linalg.solve(HPH,innovations)
    ens_new=ens+np.dot(HP.T,repres)
    return ens_new

def lenkf_rsm_py(ens,predictions,innovations,localize=None,add_obs_errors=None):
    resid_pred=compute_residual(predictions)
    resid_ens=compute_residual(ens)
    n_ens=ens.shape[1]
    alpha=1.0/float(n_ens-1)
    HP=np.dot(resid_pred,resid_ens.T)/float(n_ens-1)
    HPH=np.dot(resid_pred,resid_pred.T)/float(n_ens-1)
    if localize:
        localize(0,0,HP,HPH)
        add_obs_errors(0,0,HPH)
    return analysis_direct_solve(HP,HPH,ens,innovations)
    
