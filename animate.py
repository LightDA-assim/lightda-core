from matplotlib import pyplot as plt
import matplotlib.animation
import numpy as np
from advect1d import advance
from scipy.signal import sawtooth, square, hann
from scipy.ndimage.filters import convolve1d
from enkf import lenkf_rsm

def smooth(y,window_size):
    #win = np.ones(window_size)/window_size
    win=hann(window_size)
    y_smooth = convolve1d(y, win, mode='wrap')
    
    return y_smooth

def build_ensemble(n,n_ensemble,nghost=1):
    u=np.random.normal(0,1,[n,n_ensemble])
    a=np.random.normal(0,1,[n,n_ensemble])
    #for i in range(n_ensemble):
    #    a[nghost:-nghost,i]=smooth(a[nghost:-nghost,i],64)
    #for i in range(n_ensemble):
    #    u[nghost:-nghost,i]=smooth(u[nghost:-nghost,i],64)

    updateboundary(u,a,nghost)
    
    return np.vstack([u,a])

def updateboundary(u,a,nghost=1):
        
    # Update ghost cells
    for i in range(nghost):
        u[i]=u[-nghost*2+i]
        a[i]=a[-nghost*2+i]
        u[i-nghost]=u[i+nghost]
        a[i-nghost]=a[i+nghost]

def advance_to_time(u,a,dx,dt,limiter):
    t=0
    cfl=0.8
    while t<dt:
        dt_max=cfl*np.abs(dx/a).min()
        this_dt=min(dt_max,dt-t)
        advance(u,a,dx,this_dt,limiter)
        a_old=a.copy()
        advance(a,a_old,dx,this_dt,limiter)
        #a[:]=np.maximum(np.minimum(a,100),-100)[:]

        updateboundary(u,a)

        t+=this_dt

def advance_ensemble_to_time(ensemble,dx,dt,limiter):

    model_size,n_ensemble=ensemble.shape
    n=model_size/2
    
    for imember in range(n_ensemble):
        u=ensemble[:n,imember].flatten()
        a=ensemble[n:,imember].flatten()
        advance_to_time(u,a,dx,dt,limiter)
        ensemble[:n,imember]=u
        ensemble[n:,imember]=a

def get_distances(x1,x2):
    displacements=x2[:]-x1[:,np.newaxis]
    if len(x1.shape)==1:
        distances=np.abs(displacements)
    else:
        distances=np.linalg.norm(displacements,axis=range(len(x1.shape)-1))
    return distances

def get_distances_periodic(x1,x2,periodic_dimensions):
    displacements=x2[:]-x1[:,np.newaxis]
    if len(x1.shape)==1:
        distances=np.minimum(np.abs(displacements),
                             periodic_dimensions-np.abs(displacements))
    else:
        raise NotImplementedError('get_displacements_periodic only works with 1-dimensional values')
    return distances

def gaspari_cohn_mid(z,c):
    return 1./12*(z/c)**5 - 0.5*(z/c)**4 + 5./8*(z/c)**3 \
        + 5./3*(z/c)**2 - 5*z/c - 2./3*c/z + 4

def gaspari_cohn_close(z,c):
    return -0.25*(z/c)**5 + 0.5*(z/c)**4 + 5./8*(z/c)**3 - 5./3*(z/c)**2 + 1

def localize_gaspari_cohn(dist,c):
    localization=np.zeros(dist.shape)
    localization[dist<=2*c]=gaspari_cohn_mid(dist[dist<=2*c],c)
    localization[(dist<=c)&(dist>0)]=gaspari_cohn_close(dist[(dist<=c)&(dist>0)],c)
    localization[dist==0]=1
    return localization

class ensemble_animator(object):

    def __init__(self,n=100,n_ensemble=15,n_obs=80,cutoff=1.0):
    
        self.n=n
        self.x=np.linspace(0,2*np.pi*(n+1)/n,n)
        self.u_true=np.sin(8*self.x)
        self.a_true=np.ones(self.x.shape)+np.sin(self.x)*0.2
        self.cfl=0.8
        self.dx=np.diff(self.x)[0]
        self.dt=0.01
        self.assimilate_every=5
        
        self.limiter=5

        self.n_ensemble=n_ensemble
        self.n_obs=n_obs

        self.ensemble=build_ensemble(self.n,self.n_ensemble)

        self.inflation_factor=np.ones([n_obs,n*2])*1.12
        self.inflation_variance=0.03
        self.inflation_max=1.4
        self.inflation_min=1.0

        self.fig,axes=plt.subplots(2,1,sharex='all')
        self.fig.set_size_inches([16,9])
        for ax in axes:
            ax.set_xlim(0,2*np.pi)

        axes[0].set_ylim(-2,2)
        axes[1].set_ylim(-2,4)

        self.artists={}

        self.artists['u']={}
        self.artists['a']={}

        self.artists['u']['ensemble']=axes[0].plot(
            self.x,self.ensemble[:self.n,:])
        self.artists['a']['ensemble']=axes[1].plot(
            self.x,self.ensemble[self.n:,:])
        self.artists['a']['true'],=axes[1].plot(self.x,self.a_true,linewidth=6,color='fuchsia',alpha=0.5)
        self.artists['u']['true'],=axes[0].plot(self.x,self.u_true,linewidth=6,color='fuchsia',alpha=0.5)

        self.obs_positions=np.random.randint(0,self.n-1,self.n_obs)
        self.obs_locations=self.x[self.obs_positions]
        self.forward_operator=np.zeros([self.n_obs,self.n*2])
        self.forward_operator[np.arange(self.n_obs),self.obs_positions]=1
        obs=np.dot(self.forward_operator,np.hstack([self.u_true,self.a_true]))
        self.artists['u']['obs'],=axes[0].plot(self.obs_locations,obs,linestyle='',marker='o')

        self.localization_obs_model=localize_gaspari_cohn(
            get_distances_periodic(self.obs_locations,np.tile(self.x,[2]),np.max(self.x)),
            cutoff)
        self.localization_obs_obs=localize_gaspari_cohn(
            get_distances_periodic(self.obs_locations,self.obs_locations,np.max(self.x)),
            cutoff)

        self.obs_errors=np.random.lognormal(-3,1,self.n_obs)
        #self.obs_errors=np.ones([self.n_obs])*0.01

    def render_frame(self,i):

        if not np.all(np.isfinite(self.ensemble)):
            raise(ValueError,'nan or inf found in ensemble state array')

        print 'u:',np.max(np.abs(self.ensemble[:self.n,:]))
        print 'du:',np.max(np.abs(self.ensemble[1:self.n+1,:]-self.ensemble[:self.n,:]))
        print 'a:',np.max(np.abs(self.ensemble[self.n:,:]))
        print 'da:',np.max(np.abs(self.ensemble[self.n+1:,:]-self.ensemble[self.n:-1,:]))

        print 'advancing'
        advance_to_time(self.u_true,self.a_true,self.dx,self.dt,self.limiter)
        advance_ensemble_to_time(self.ensemble,self.dx,self.dt,self.limiter)
        
        print 'u:',np.max(np.abs(self.ensemble[:self.n,:]))
        print 'du:',np.max(np.abs(self.ensemble[1:self.n+1,:]-self.ensemble[:self.n,:]))
        print 'a:',np.max(np.abs(self.ensemble[self.n:,:]))
        print 'da:',np.max(np.abs(self.ensemble[self.n+1:,:]-self.ensemble[self.n:-1,:]))
        
        if not np.all(np.isfinite(self.ensemble)):
            raise(ValueError,'nan or inf found in ensemble state array')

        obs_true=np.dot(self.forward_operator,np.hstack([self.u_true,self.a_true]))
        obs_w_errors=np.random.normal(obs_true,self.obs_errors)

        if i%self.assimilate_every==self.assimilate_every-1:
            pass

        modified_artists=[]

        for i,variable in enumerate(['u','a']):
            for j,member in enumerate(self.artists[variable]['ensemble']):
                member.set_data(self.x,self.ensemble[i*self.n:(i+1)*self.n,j].T)
                modified_artists.append(member)
        self.artists['u']['obs'].set_data(self.x[self.obs_positions],obs_w_errors)

        self.artists['u']['true'].set_data(self.x,self.u_true)
        self.artists['a']['true'].set_data(self.x,self.a_true)

        modified_artists.extend([
            self.artists['u']['true'],self.artists['a']['true'],self.artists['u']['obs']])

        return modified_artists

def plot_localization(localization=None):
    fig=plt.figure()
    ax=fig.add_subplot(1,1,1)

    if localization is None:
        x=np.linspace(0,1,100)
        y=np.sort(np.random.rand(40))

        c=0.1
    
        localization=localize_gaspari_cohn(get_distances_periodic(x,y,1),c)

    img=ax.imshow(localization)
    plt.colorbar(img)


    
if __name__=='__main__':
    animator=ensemble_animator()
    #animator.render_frame(0)
    sort_inds=np.argsort(animator.obs_locations)
    plot_localization(animator.localization_obs_model[sort_inds,:])
    ani = matplotlib.animation.FuncAnimation(animator.fig, animator.render_frame)

    
    plt.show()
