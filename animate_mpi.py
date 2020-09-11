from matplotlib import pyplot as plt
import matplotlib.animation
import numpy as np
from advect1d import advance
from scipy.signal import sawtooth, square, hann
from scipy.ndimage.filters import convolve1d
from enkf import lenkf_rsm

def build_ensemble(n,n_ensemble,nghost=1):
    """
    Generate an initial ensemble state.

    Arguments:
    - n: Number of points in model grid
    - n_ensemble: Number of ensemble members
    - nghost: Number of ghost cells
    """

    # u is a normally distributed random variable
    u=np.random.normal(0,2,[n,n_ensemble])

    # a consists of a normally distributed random constant value for each
    # ensemble member, plus a normally distributed random component within
    # each ensemble member
    a=np.random.normal(np.random.normal(0,2,[1,n_ensemble]),2,[n,n_ensemble])

    # Update ghost cells
    updateboundary(u,a,nghost)

    # Combine u and a into a single array and return
    return np.vstack([u,a])

def updateboundary(u,a,nghost=1):
    """
    Update ghost cells for u and a

    Arguments:
    - u: Passively advected variable
    - a: Flow velocity
    - nghost: Number of ghost cells
    """
    
    # Update ghost cells
    for i in range(nghost):
        u[i]=u[-nghost*2+i]
        a[i]=a[-nghost*2+i]
        u[i-nghost]=u[i+nghost]
        a[i-nghost]=a[i+nghost]

def advance_to_time(u,a,dx,dt,limiter):
    """
    Advance the model state to the time increment specified by dt

    Arguments:
    - u: Passively advected variable
    - a: Flow velocity
    - dx: Grid cell size
    - dt: Time increment
    - limiter: Limiter type (integer from 0 to 5)
    """
    
    # Set initial time
    t=0

    # Set CFL
    cfl=0.5
    
    while t<dt:
        # Find maximum stable time step
        dt_max=cfl*np.abs(dx/a).min()

        # Time step to move toward dt but not cross it
        this_dt=min(dt_max,dt-t)

        # Advance u and a
        advance(u,a,dx,this_dt,limiter)
        a_old=a.copy()
        advance(a,a_old,dx,this_dt,limiter)

        # Update boundary conditions
        updateboundary(u,a)

        t+=this_dt

def advance_ensemble_to_time(ensemble,dx,dt,limiter):
    """
    Advance all ensemble members to the time specified by dt

    Arguments:
    - ensemble: Ensemble state array
    - dx: Grid cell size
    - dt: Time increment
    - limiter: Limiter type (integer from 0 to 5)
    """

    model_size,n_ensemble=ensemble.shape
    n=int(model_size/2)
    
    for imember in range(n_ensemble):
        u=ensemble[:n,imember].flatten()
        a=ensemble[n:,imember].flatten()
        advance_to_time(u,a,dx,dt,limiter)
        ensemble[:n,imember]=u
        ensemble[n:,imember]=a

def get_distances(x1,x2):
    """
    Get the distances between the points in two arrays x1 and x2 as a 2-D matrix
    """

    # Displacements between all points
    displacements=x2[:]-x1[:,np.newaxis]
    
    if len(x1.shape)==1:
        # 1-D points, so distances are just the absolute value of the
        # displacements
        distances=np.abs(displacements)
    else:
        # Untested, but in general np.linalg.norm should return the
        # Cartesian distance between points of arbitrary dimension
        distances=np.linalg.norm(displacements,axis=range(len(x1.shape)-1))

    return distances

def get_distances_periodic(x1,x2,periodic_dimensions):
    """
    Get distances between points in a periodic domain
    """
    
    displacements=x2[:]-x1[:,np.newaxis]
    if len(x1.shape)==1:
        distances=np.minimum(np.abs(displacements),
                             np.abs(periodic_dimensions-np.abs(displacements)))
    else:
        raise NotImplementedError('get_displacements_periodic only works with 1-dimensional values')
    
    return distances

def gaspari_cohn_mid(z,c):
    """
    Gaspari-Cohn correlation function for middle distances (between c and 2*c)
    
    Arguments:
    - z: Points to be evaluated
    - c: Cutoff value
    """
    return 1./12*(z/c)**5 - 0.5*(z/c)**4 + 5./8*(z/c)**3 \
        + 5./3*(z/c)**2 - 5*z/c - 2./3*c/z + 4

def gaspari_cohn_close(z,c):
    """
    Gaspari-Cohn correlation function for middle distances (z<=c)
    
    Arguments:
    - z: Points to be evaluated
    - c: Cutoff value
    """
    return -0.25*(z/c)**5 + 0.5*(z/c)**4 + 5./8*(z/c)**3 - 5./3*(z/c)**2 + 1

def localize_gaspari_cohn(dist,c):
    """
    Gaspari-Cohn correlation function
    
    Arguments:
    - z: Points to be evaluated
    - c: Cutoff value
    """

    # Initialize localization array
    localization=np.zeros(dist.shape)

    # Mask for mid-distance points
    mid_mask=(dist<=2*c)

    # c for mid-distance points (array indexed by mid_mask, or
    # fallback to scalar)
    try:
        c_mid=c[mid_mask]
    except TypeError:
        c_mid=c

    # Compute weights for mid-distance points
    close_weights=gaspari_cohn_mid(dist[mid_mask],c_mid)
    localization[mid_mask]=close_weights
    
    # Mask for close points
    close_mask=(dist<=c)&(dist>0)
    
    # c for mid-distance points (array indexed by close_mask, or
    # fallback to scalar)
    try:
        c_close=c[close_mask]
    except TypeError:
        c_close=c
    localization[close_mask]=gaspari_cohn_close(dist[close_mask],c_close)
    localization[dist==0]=1
    return localization

def inflate_ensemble(ensemble,inflation_factor):
    """
    Apply covariance inflation to the ensemble state

    Arguments:
    - ensemble: Ensemble state array
    - inflation_factor: Inflation factor
    """

    ensemble_mean=np.mean(ensemble,axis=1).reshape([ensemble.shape[0],1])
    
    inflated_state=(
        (ensemble-ensemble_mean)*inflation_factor[:,np.newaxis]
        +ensemble_mean)
    
    return inflated_state

def update_inflation_factor_localized(ensemble,forward_operator,observations,obs_errors,prior_inflation_factor,inflation_variance,weight):
    """
    Adaptive inflation as described in Anderson (2009)
    """

    predictions=np.dot(forward_operator,ensemble)
    ensemble_var=np.var(predictions,axis=1)

    actual_distance=np.abs(np.mean(predictions,axis=1)-observations)
    expected_lambda=(1+weight*(np.sqrt(prior_inflation_factor)-1))**2
    n_observations=len(observations)
    model_size=ensemble.shape[0]

    thetabar=np.sqrt(
        expected_lambda*inflation_variance
        +obs_errors.reshape([n_observations,1])**2)
    thetabar=np.maximum(thetabar,1e-10)
    lbar=np.sqrt(2*np.pi*thetabar)*np.exp(
        -0.5*actual_distance.reshape(
            [n_observations,1])**2*thetabar**-2)
    dtheta_dlambda=0.5*inflation_variance*weight*(
        1-weight+weight*np.sqrt(prior_inflation_factor))/ \
        (thetabar/np.sqrt(prior_inflation_factor))
    lprime=np.maximum(
        (lbar*(actual_distance.reshape([n_observations,1])**2*thetabar**-2-1)
         /thetabar)*dtheta_dlambda,1e-10)
    
    coeffs=[
        np.ones([n_observations,model_size]),
        lbar/lprime-2*prior_inflation_factor,
        (prior_inflation_factor**2-inflation_variance
         -lbar*prior_inflation_factor/lprime)
    ]

    coeffs=np.array(coeffs)
    posterior_inflation_factor=np.empty([n_observations,model_size])

    # Array of positive and negative for use in quadratic equation
    pm=np.array([1,-1])
    
    for i in range(n_observations):
        for j in range(model_size):
            a,b,c=coeffs[:,i,j]
            # Find roots using quadratic equation
            roots=(-b+pm*np.sqrt(b**2-4*a*c))/(2*a)

            # Pick the root closest to the prior inflation factor
            closest_root_ind=np.argmin(
                np.abs(roots-prior_inflation_factor[i,j]))
            posterior_inflation_factor[i,j]=roots[closest_root_ind]

    return posterior_inflation_factor

class ensemble_animator(object):

    def add_obs_err(self,step,ind_p,HPH):
        for i in range(self.n_obs):
            HPH[i,i]+=self.obs_errors[i]**2

    def localize(self,step,ind_p,HP_p,HPH):
        HP_p[:,:]=HP_p*self.localization_obs_model[:,self.batch_boundaries[ind_p]:self.batch_boundaries[ind_p+1]]
        HPH[:,:]=HPH*self.localization_obs_obs

    def __init__(self,n=100,n_ensemble=15,n_obs=80,cutoff=0.35):
    
        self.n=n
        nghost=1
        domain_width=2*np.pi
        self.x=np.linspace(-float(nghost)/n*domain_width,domain_width*(n+nghost)/n,n)
        self.u_true=np.sin(8*self.x)
        self.a_true=np.ones(self.x.shape)+np.sin(self.x)*0.5
        self.cfl=0.8
        self.dx=np.diff(self.x)[0]
        self.dt=0.01
        self.assimilate_every=5
        
        self.limiter=5

        self.n_ensemble=n_ensemble
        self.n_obs=n_obs

        self.ensemble=np.asfortranarray(build_ensemble(self.n,self.n_ensemble))

        self.inflation_factor=np.ones([n_obs,n*2])*1.13
        self.inflation_variance=0.1
        self.inflation_max=2
        self.inflation_min=1.0

        self.obs_positions=np.random.randint(0,self.n-1,self.n_obs)
        self.obs_locations=self.x[self.obs_positions]
        self.forward_operator=np.zeros([self.n_obs,self.n*2])
        self.forward_operator[np.arange(self.n_obs),self.obs_positions]=1

        assim_batchsize=50

        batch_boundaries=np.arange(0,self.ensemble.shape[0],assim_batchsize)
        self.batch_boundaries=np.concatenate([batch_boundaries,[self.ensemble.shape[0]]])
        print(self.batch_boundaries)

        cutoff_u_a=0.4

        obs_model_distances=get_distances_periodic(
            self.obs_locations,np.tile(self.x,[2]),np.max(self.x))
        self.localization_obs_model=localize_gaspari_cohn(
            obs_model_distances,
            np.hstack([np.ones([n_obs,n])*cutoff,np.ones([n_obs,n])*cutoff_u_a]))
        
        obs_obs_distances=get_distances_periodic(
            self.obs_locations,self.obs_locations,np.max(self.x))
        self.localization_obs_obs=localize_gaspari_cohn(
            obs_obs_distances,
            cutoff)

        self.obs_errors=np.random.lognormal(-3,1,self.n_obs)

        self.create_figure()

    def create_figure(self):
        self.fig,axes=plt.subplots(2,1,sharex='all',figsize=[8,4.5])
        for ax in axes:
            ax.set_xlim(0,2*np.pi)

        axes[0].set_ylim(-4,4)
        axes[1].set_ylim(-3,5)

        self.artists={}

        self.artists['u']={}
        self.artists['a']={}

        self.artists['u']['ensemble']=axes[0].plot(
            self.x,self.ensemble[:self.n,:])
        self.artists['a']['ensemble']=axes[1].plot(
            self.x,self.ensemble[self.n:,:])
        self.artists['a']['true'],=axes[1].plot(
            self.x,self.a_true,linewidth=6,color='fuchsia',alpha=0.5)
        self.artists['u']['true'],=axes[0].plot(
            self.x,self.u_true,linewidth=6,color='fuchsia',alpha=0.5)
        self.artists['u']['mean'],=axes[0].plot(
            self.x,
            np.mean(self.ensemble[self.n:],axis=1),
            linewidth=2,color='k')
        self.artists['a']['mean'],=axes[1].plot(
            self.x,
            np.mean(self.ensemble[:self.n],axis=1),
            linewidth=2,color='k')
        
        obs=np.dot(self.forward_operator,np.hstack([self.u_true,self.a_true]))
        self.artists['u']['obs'],=axes[0].plot(self.obs_locations,obs,linestyle='',marker='o')

        predictions=np.dot(self.forward_operator,self.ensemble)
        self.artists['u']['predictions']=axes[0].plot(self.obs_locations,predictions,linestyle='',marker='.')


    def apply_inflation(self,predictions,obs_w_errors):
        # Compute residuals
        resid_pred=predictions-np.mean(predictions,axis=1)[:,np.newaxis]
        resid_ens=self.ensemble-np.mean(self.ensemble,axis=1)[:,np.newaxis]

        # Weight factors for inflation
        hp=resid_pred.dot(resid_ens.T)
        inflation_weight=self.localization_obs_model*hp

        # Update inflation factor
        self.inflation_factor=update_inflation_factor_localized(self.ensemble,self.forward_operator,obs_w_errors,self.obs_errors,self.inflation_factor,self.inflation_variance,inflation_weight)

        # Apply inflation limits
        self.inflation_factor=np.minimum(
            self.inflation_factor,self.inflation_max)
        self.inflation_factor=np.maximum(
            self.inflation_factor,self.inflation_min)

        # Get new inflation variance
        self.inflation_variance=np.var(self.inflation_factor)

        # Inflate ensemble
        self.ensemble=inflate_ensemble(self.ensemble,np.mean(self.inflation_factor,axis=0))

    def render_frame(self,i):

        if not np.all(np.isfinite(self.ensemble)):
            raise ValueError('nan or inf found in ensemble state array')

        advance_to_time(self.u_true,self.a_true,self.dx,self.dt,self.limiter)
        advance_ensemble_to_time(self.ensemble,self.dx,self.dt,self.limiter)

        if not np.all(np.isfinite(self.ensemble)):
            raise ValueError('nan or inf found in ensemble state array')

        obs_true=np.dot(self.forward_operator,np.hstack([self.u_true,self.a_true]))
        obs_w_errors=np.random.normal(obs_true,self.obs_errors)

        predictions=np.dot(self.forward_operator,self.ensemble)

        do_assimilation=(i%self.assimilate_every==self.assimilate_every-1)
        
        if do_assimilation:

            # Covariance inflation
            #self.apply_inflation(predictions,obs_w_errors)

            # Get observation perturbations
            obs_perturbations=np.random.normal(
                0,self.obs_errors,[self.n_ensemble,self.n_obs]).T

            # Compute innovations
            innovations=np.asfortranarray(obs_w_errors[:,np.newaxis]+obs_perturbations-predictions)

            import tables
            import os

            path='ensembles/{}'.format(i)

            if not os.path.exists(path):
                os.makedirs(path)

            with tables.open_file(os.path.join(path,'observations.h5'),'w') as h5file:
                h5file.create_array(h5file.root,'observations',obs_w_errors)
                h5file.create_array(h5file.root,'obs_positions',self.obs_positions)
                h5file.create_array(h5file.root,'obs_errors',self.obs_errors)

            with tables.open_file(os.path.join(path,'preassim.h5'),'w') as h5file:
                h5file.create_array(h5file.root,'ensemble_state',self.ensemble)
                h5file.create_array(h5file.root,'predictions',predictions)
                h5file.create_array(h5file.root,'innovations',innovations)

            import subprocess

            print(os.getcwd())

            cmd=['mpirun','-n','4','./advect1d_assimilate',
                                '--istep',str(i),
                                '--state_size',str(self.ensemble.shape[0]),
                                '--batch_size',str(int(
                                    self.ensemble.shape[0]/10)),
                                '--n_ensemble',str(self.ensemble.shape[1]),
                                '--n_observations',str(len(self.obs_positions))]

            print(' '.join(cmd))

            p=subprocess.Popen(cmd)
            retcode=p.wait()

            if retcode!=0:
                import sys
                sys.exit()

            for i in range(self.ensemble.shape[1]):
                with tables.open_file(
                        os.path.join(path,'postassim_{}.h5'.format(i))
                ) as h5file:
                    self.ensemble[:,i]=h5file.root.ensemble_state

            # Update predictions
            predictions=np.dot(self.forward_operator,self.ensemble)

        # Array of artists that were updated in this frame
        modified_artists=[]

        # Update plots of ensemble members
        for i,variable in enumerate(['u','a']):
            for j,member in enumerate(self.artists[variable]['ensemble']):
                member.set_data(
                    self.x,
                    self.ensemble[i*self.n:(i+1)*self.n,j].T)
                modified_artists.append(member)

        # Update plots of ensemble predictions
        for j,member in enumerate(self.artists['u']['predictions']):
            member.set_data(self.obs_locations,predictions[:,j])
            modified_artists.append(member)

        # Update plots of observations
        if do_assimilation:
            self.artists['u']['obs'].set_data(self.x[self.obs_positions],obs_w_errors)

        # Update plots of true state
        self.artists['u']['true'].set_data(self.x,self.u_true)
        self.artists['a']['true'].set_data(self.x,self.a_true)

        # Update plots of ensemble mean
        self.artists['u']['mean'].set_data(
            self.x,np.mean(self.ensemble[:self.n],axis=1))
        self.artists['a']['mean'].set_data(
            self.x,np.mean(self.ensemble[self.n:],axis=1))

        modified_artists.extend([
            self.artists['u']['true'],self.artists['a']['true'],
            self.artists['u']['obs'],
            self.artists['u']['mean'],self.artists['a']['mean']
        ])

        return modified_artists

def plot_localization(localization=None):
    """
    Plot a localization array
    """
    fig=plt.figure()
    ax=fig.add_subplot(1,1,1)

    if localization is None:
        x=np.linspace(0,1,100)
        y=np.sort(np.random.rand(40))

        c=0.1
    
        localization=localize_gaspari_cohn(get_distances_periodic(x,y,1),c)

    img=ax.imshow(localization)
    plt.colorbar(img)

def plot_correlation():
    fig=plt.figure()
    ax=fig.add_subplot(1,1,1)

    x=np.linspace(0,1,200)
    y=localize_gaspari_cohn(x,0.5)
    ax.plot(x,y)
    
if __name__=='__main__':
    animator=ensemble_animator()
    sort_inds=np.argsort(animator.obs_locations)
    #plot_localization(animator.localization_obs_model[sort_inds,:])
    #plot_correlation()
    ani = matplotlib.animation.FuncAnimation(animator.fig, animator.render_frame)

    
    plt.show()
