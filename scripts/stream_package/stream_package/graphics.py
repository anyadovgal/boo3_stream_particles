from galpy.orbit import Orbit
from galpy import potential
import numpy as np
from astropy import units
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation
from IPython.display import HTML
from galpy.potential import MWPotential2014,ChandrasekharDynamicalFrictionForce,HernquistPotential,MovingObjectPotential
from galpy.util import conversion
from streamtools.df import streamspraydf

#Galpy internal units scaling factors
ro = 8. #distance to Galactic centre in kpc, scales distances
vo = 220. #circular velocity at solar circle
to=conversion.time_in_Gyr(ro=ro,vo=vo)
mo=conversion.mass_in_msol(ro=ro,vo=vo)

class StellarStream3D(object):
    def __init__(self, gorbit, starorbits, dt, tdisrupt, pot=MWPotential2014, ro=8., vo=220.):
        # Initializing the stream object that can be animated
        self.o = gorbit
        self.oe = starorbits
        self.tesc = -1.*dt
        self.nstar = len(self.oe)/2
        self.tdisrupt = tdisrupt
        self.pot = pot
        
        # since we're dealing w/ galpy units:
        self.ro = ro
        self.vo = vo
        self.to = conversion.time_in_Gyr(ro=ro, vo=vo)
        
    def _init_fig(self, xlim=(-50,50), ylim=(-50,50), zlim=(-50,50)):
        self.fig = plt.figure(figsize = (10,10))
        #self.ax = self.fig.add_subplot(projection='3d', xlim=xlim, ylim=ylim, zlim=zlim)
        self.ax = p3.Axes3D(self.fig, auto_add_to_figure=False)
        self.fig.add_axes(self.ax)
        
        self.ax.set_xlim3d(xlim)
        self.ax.set_ylim3d(ylim)
        self.ax.set_zlim3d(zlim)
        
        self.ax.set_xlabel('X (kpc)')
        self.ax.set_ylabel('Y (kpc)')
        self.ax.set_zlabel('Z (kpc)')
        self.txt_title=self.ax.set_title('')
        self.pt, = self.ax.plot([],[],[],'.')
        self.line, = self.ax.plot([], [], [], lw=2, c='k')
        
        self.ax.plot(0,0,0, 'g*')
        
        self.ax.view_init(20,80)
        
    def _set_data(self, gdata, sdata):
        self.gdata = gdata
        self.sdata = sdata
        
    def _ani_init(self):
        self.line.set_data([], [])
        self.line.set_3d_properties([])
        
        self.pt.set_data([], [])
        self.pt.set_3d_properties([])
        
        return self.line, self.pt
    
    def _ani_update(self, i):
        
#         if i < 5:
#             x = self.gdata[0:i+1,0]
#             y = self.gdata[0:i+1,1]
#             z = self.gdata[0:i+1,2]
#         else:
#             x = self.gdata[i-5:i+1,0]
#             y = self.gdata[i-5:i+1,1]
#             z = self.gdata[i-5:i+1,2]

        x = self.gdata[0:i+1,0]
        y = self.gdata[0:i+1,1]
        z = self.gdata[0:i+1,2]
        
        self.line.set_data(x,y)
        self.line.set_3d_properties(z)
        
        #escape index
        escindx = self.tesc/self.to <= self.ts[i]
        
        if np.sum(escindx) > 0:
            self.pt.set_data(self.sdata[i][0][escindx],self.sdata[i][1][escindx])
            self.pt.set_3d_properties(self.sdata[i][2][escindx])
        else:
            self.pt.set_data([],[])
            self.pt.set_3d_properties([])
            
        self.txt_title.set_text('%s' % str (self.ts[i]*self.to))
        
        return self.line, self.pt
    
    def animate(self, frames=300, interval=50, xlim=(-50,50), ylim=(-50,50), zlim=(-50,50), solarmotion=[-11.1, 24.0, 7.25]):
        
        self._init_fig(xlim, ylim, zlim)
        
        self.ts=np.linspace(-1.*self.tdisrupt/self.to,0.,frames)
        
        tsint=np.linspace(0,-1.*self.tdisrupt/self.to,1000)
        
        self.o.integrate(tsint,self.pot)

        gdata=np.zeros(shape=(frames,3))
        
        for i in range(0, frames):
            gdata[i] = [self.o.x(self.ts[i])*self.ro, self.o.y(self.ts[i])*self.ro, self.o.z(self.ts[i])*self.ro]
            
        sdata = np.zeros(shape=(frames, 3, int(2*self.nstar)))
        self.oe.integrate(tsint, self.pot)
        
        for i in range(0, frames):
            sdata[i] = [self.oe.x(self.ts[i]), self.oe.y(self.ts[i]), self.oe.z(self.ts[i])]
            
        self._set_data(gdata, sdata)
        
        self.anim = animation.FuncAnimation(self.fig, self._ani_update, init_func=self._ani_init, frames=frames, interval=interval, blit=False)

        
def streamorbits(mass, o, tdisrupt, pot=MWPotential2014, nstar=100):
    #Distribution function for the leading tail
    #Note input values are scaled to galpy units
    spdf= streamspraydf(mass/mo,
                   progenitor=o,
                   pot=pot,
                   tdisrupt=tdisrupt/to)
    #Distribution function for the trailing tail
    spdft= streamspraydf(mass/mo,
                   progenitor=o,
                   pot=pot,
                   tdisrupt=tdisrupt/to,
                   leading=False)
    
    #Sample the distribution functions to create nstar
    #returndt lets you know the time the star was ejected 
    #integrate=True returns the stars integrated 
    #position and velocity at time = tdisrupt
    #All values in galpy units
    RvR,dt= spdf.sample(n=nstar,returndt=True,integrate=True)
    RvRt,dtt= spdft.sample(n=nstar,returndt=True,integrate=True)
    
    vxvv=np.column_stack([RvR[0],RvR[1],RvR[2],RvR[3],RvR[4],RvR[5]])
    vxvvt=np.column_stack([RvRt[0],RvRt[1],RvRt[2],RvRt[3],RvRt[4],RvRt[5]])

    vxvva=np.column_stack([np.append(RvR[0],RvRt[0]),
                           np.append(RvR[1],RvRt[1]),
                           np.append(RvR[2],RvRt[2]),
                           np.append(RvR[3],RvRt[3]),
                           np.append(RvR[4],RvRt[4]),
                           np.append(RvR[5],RvRt[5])])

    oleading=Orbit(vxvv,ro=ro,vo=vo,solarmotion=[-11.1, 24.0, 7.25])
    otrailing=Orbit(vxvvt,ro=ro,vo=vo,solarmotion=[-11.1, 24.0, 7.25])

    oall=Orbit(vxvva,ro=ro,vo=vo,solarmotion=[-11.1, 24.0, 7.25])
    dtall=np.append(dt,dtt)
    
    return oall, oleading, otrailing, dtall, dt, dtt
    
def streamorbitslmc(mass, o, tdisrupt, pot=MWPotential2014, nstar=100):
    mass_lmc=1.0e11 #solar masses
    rscale_lmc=10.2 #kpc

    #Initialize and integrate the orbit of the LMC
    #Note orbit has to be integrated back 5 Gyr
    #Note we assume the LMC experienced dynamical friction due to MW
    o_lmc = Orbit.from_name('LMC', ro=ro, vo=vo, solarmotion=[-11.1, 24.0, 7.25])
    ts= np.linspace(0.,-tdisrupt/to,1001)
    cdf= ChandrasekharDynamicalFrictionForce(GMs=mass_lmc/mo, rhm=rscale_lmc/ro, dens=pot, ro=ro,vo=vo)
    o_lmc.integrate(ts,pot+cdf)

    #Setup a moving Hernquist potential to represent the LMC
    pot_lmc = HernquistPotential(mass_lmc/mo,rscale_lmc/ro,ro=ro,vo=vo)
    moving_pot_lmc = MovingObjectPotential(o_lmc, pot_lmc,ro=ro,vo=vo)

    #Add the moving Hernquest potential to the MW
    total_pot = [pot]
    total_pot += [moving_pot_lmc]
    
    #Now when you intitialize the distribution functions
    #you have to specify rtpot, which is the potential
    #used to calculate the cluster's tidal radius. The
    #calculation can't be done when the MovingPotential
    #is included in total_pot

    #Distribution function for the leading tail
    #Note input values are scaled to galpy units
    spdf_lmc= streamspraydf(mass/mo,
                       progenitor=o,
                       pot=total_pot,
                       tdisrupt=tdisrupt/to,
                       rtpot=pot)
    #Distribution function for the trailing tail
    spdft_lmc= streamspraydf(mass/mo,
                       progenitor=o,
                       pot=total_pot,
                       tdisrupt=tdisrupt/to,
                       rtpot=pot,
                       leading=False)
    
    #Sample the distribution functions to ,create nstar
    #returndt lets you know the time the star was ejected 
    #integrate=True returns the stars integrated 
    #position and velocity at time = tdisrupt
    #All values in galpy units

    RvR_lmc,dt_lmc= spdf_lmc.sample(n=nstar,returndt=True,integrate=True)
    RvRt_lmc,dtt_lmc= spdft_lmc.sample(n=nstar,returndt=True,integrate=True)
    
    vxvva_lmc=np.column_stack([np.append(RvR_lmc[0],RvRt_lmc[0]),
                           np.append(RvR_lmc[1],RvRt_lmc[1]),
                           np.append(RvR_lmc[2],RvRt_lmc[2]),
                           np.append(RvR_lmc[3],RvRt_lmc[3]),
                           np.append(RvR_lmc[4],RvRt_lmc[4]),
                           np.append(RvR_lmc[5],RvRt_lmc[5])])

    oall_lmc=Orbit(vxvva_lmc,ro=ro,vo=vo,solarmotion=[-11.1, 24.0, 7.25])
    dtall_lmc=np.append(dt_lmc,dtt_lmc)
    
    return o_lmc, oall_lmc, dtall_lmc



def lmc_potential_wacc(tdisrupt, pot=MWPotential2014):
    olmc= Orbit.from_name('LMC')
    cdf= ChandrasekharDynamicalFrictionForce(GMs=10.**11.*units.Msun,rhm=5.*units.kpc,
                                                 dens=pot)
    ts= np.linspace(0.,-tdisrupt,1001)*units.Gyr
    olmc.integrate(ts,pot+cdf)
    
    lmcpot= HernquistPotential(amp=2*10.**11.*units.Msun,
                               a=5.*units.kpc/(1.+np.sqrt(2.))) #rhm = (1+sqrt(2)) a
    moving_lmcpot= MovingObjectPotential(olmc,pot=lmcpot)
    
    from galpy.potential import (evaluateRforces, evaluatephitorques,
                                 evaluatezforces)
    loc_origin= 1e-4 # Small offset in R to avoid numerical issues
    ax= lambda t: evaluateRforces(moving_lmcpot,loc_origin,0.,phi=0.,t=t,
                                      use_physical=False)
    ay= lambda t: evaluatephitorques(moving_lmcpot,loc_origin,0.,phi=0.,t=t,
                                        use_physical=False)/loc_origin
    az= lambda t: evaluatezforces(moving_lmcpot,loc_origin,0.,phi=0.,t=t,
                                      use_physical=False)
    
    t_intunits= olmc.time(use_physical=False)[::-1] # need to reverse the order for interp
    ax4int= np.array([ax(t) for t in t_intunits])
    ax_int= lambda t: np.interp(t,t_intunits,ax4int)
    ay4int= np.array([ay(t) for t in t_intunits])
    ay_int= lambda t: np.interp(t,t_intunits,ay4int)
    az4int= np.array([az(t) for t in t_intunits])
    az_int= lambda t: np.interp(t,t_intunits,az4int)
    
    from galpy.potential import NonInertialFrameForce
    nip= NonInertialFrameForce(a0=[ax_int,ay_int,az_int])
    
    return moving_lmcpot + pot + nip
    
def streamorbitslmc_acc_lt(mass, o, tdisrupt, pot=MWPotential2014, nstar=100):
    o_lmc = Orbit.from_name('LMC', ro=ro, vo=vo, solarmotion=[-11.1, 24.0, 7.25])
    total_pot = lmc_potential_wacc(tdisrupt, pot)
    
    #Now when you intitialize the distribution functions
    #you have to specify rtpot, which is the potential
    #used to calculate the cluster's tidal radius. The
    #calculation can't be done when the MovingPotential
    #is included in total_pot

    #Distribution function for the leading tail
    #Note input values are scaled to galpy units
    spdf_lmc= streamspraydf(mass/mo,
                       progenitor=o,
                       pot=total_pot,
                       tdisrupt=tdisrupt/to,
                       rtpot=pot)
    #Distribution function for the trailing tail
    spdft_lmc= streamspraydf(mass/mo,
                       progenitor=o,
                       pot=total_pot,
                       tdisrupt=tdisrupt/to,
                       rtpot=pot,
                       leading=False)
    
    #Sample the distribution functions to ,create nstar
    #returndt lets you know the time the star was ejected 
    #integrate=True returns the stars integrated 
    #position and velocity at time = tdisrupt
    #All values in galpy units

    RvR_lmc,dt_lmc= spdf_lmc.sample(n=nstar,returndt=True,integrate=True)
    RvRt_lmc,dtt_lmc= spdft_lmc.sample(n=nstar,returndt=True,integrate=True)
    
    vxvva_lmc=np.column_stack([np.append(RvR_lmc[0],RvRt_lmc[0]),
                           np.append(RvR_lmc[1],RvRt_lmc[1]),
                           np.append(RvR_lmc[2],RvRt_lmc[2]),
                           np.append(RvR_lmc[3],RvRt_lmc[3]),
                           np.append(RvR_lmc[4],RvRt_lmc[4]),
                           np.append(RvR_lmc[5],RvRt_lmc[5])])
    
    vxvv=np.column_stack([RvR_lmc[0],RvR_lmc[1],RvR_lmc[2],RvR_lmc[3],RvR_lmc[4],RvR_lmc[5]])
    vxvvt=np.column_stack([RvRt_lmc[0],RvRt_lmc[1],RvRt_lmc[2],RvRt_lmc[3],RvRt_lmc[4],RvRt_lmc[5]])

    oleading_lmc=Orbit(vxvv,ro=ro,vo=vo,solarmotion=[-11.1, 24.0, 7.25])
    otrailing_lmc=Orbit(vxvvt,ro=ro,vo=vo,solarmotion=[-11.1, 24.0, 7.25])

    oall_lmc=Orbit(vxvva_lmc,ro=ro,vo=vo,solarmotion=[-11.1, 24.0, 7.25])
    dtall_lmc=np.append(dt_lmc,dtt_lmc)
    
    return o_lmc, oall_lmc, oleading_lmc, otrailing_lmc, dtall_lmc
    
def streamorbitslmc_acc(mass, o, tdisrupt, pot=MWPotential2014, nstar=100):
    o_lmc = Orbit.from_name('LMC', ro=ro, vo=vo, solarmotion=[-11.1, 24.0, 7.25])
    total_pot = lmc_potential_wacc(tdisrupt, pot)
    
    #Now when you intitialize the distribution functions
    #you have to specify rtpot, which is the potential
    #used to calculate the cluster's tidal radius. The
    #calculation can't be done when the MovingPotential
    #is included in total_pot

    #Distribution function for the leading tail
    #Note input values are scaled to galpy units
    spdf_lmc= streamspraydf(mass/mo,
                       progenitor=o,
                       pot=total_pot,
                       tdisrupt=tdisrupt/to,
                       rtpot=pot)
    #Distribution function for the trailing tail
    spdft_lmc= streamspraydf(mass/mo,
                       progenitor=o,
                       pot=total_pot,
                       tdisrupt=tdisrupt/to,
                       rtpot=pot,
                       leading=False)
    
    #Sample the distribution functions to ,create nstar
    #returndt lets you know the time the star was ejected 
    #integrate=True returns the stars integrated 
    #position and velocity at time = tdisrupt
    #All values in galpy units

    RvR_lmc,dt_lmc= spdf_lmc.sample(n=nstar,returndt=True,integrate=True)
    RvRt_lmc,dtt_lmc= spdft_lmc.sample(n=nstar,returndt=True,integrate=True)
    
    vxvva_lmc=np.column_stack([np.append(RvR_lmc[0],RvRt_lmc[0]),
                           np.append(RvR_lmc[1],RvRt_lmc[1]),
                           np.append(RvR_lmc[2],RvRt_lmc[2]),
                           np.append(RvR_lmc[3],RvRt_lmc[3]),
                           np.append(RvR_lmc[4],RvRt_lmc[4]),
                           np.append(RvR_lmc[5],RvRt_lmc[5])])

    oall_lmc=Orbit(vxvva_lmc,ro=ro,vo=vo,solarmotion=[-11.1, 24.0, 7.25])
    dtall_lmc=np.append(dt_lmc,dtt_lmc)
    
    return o_lmc, oall_lmc, dtall_lmc
    
    
def orbit_plots(o, pot, tint=5):
    delt = np.linspace(0,-tint/to,1000)
    deltfwd = np.linspace(0, tint/to, 1000)
    
    oint_p1 = o()
    oint_p1.integrate(delt, pot=pot)
    ointfwd_p1 = o()
    ointfwd_p1.integrate(deltfwd, pot=pot)

    fig, axs = plt.subplots(2, 2, figsize=(10,10))

    # Orbit in X-Y
    ax1 = axs[0,0]
    ax1.plot(oint_p1.x(delt)*ro, oint_p1.y(delt)*ro, "b--")
    ax1.plot(ointfwd_p1.x(deltfwd)*ro, ointfwd_p1.y(deltfwd)*ro, "b-")
    ax1.set_xlabel('X (kpc)')
    ax1.set_ylabel('Y (kpc)')

    # Orbit in X-Z
    ax2 = axs[0,1]
    ax2.plot(oint_p1.x(delt)*ro, oint_p1.z(delt)*ro, "b--")
    ax2.plot(ointfwd_p1.x(deltfwd)*ro, ointfwd_p1.z(deltfwd)*ro, "b-")
    ax2.set_xlabel('X (kpc)')
    ax2.set_ylabel('Z (kpc)')

    # Orbit in Y-Z
    ax3 = axs[1,0]
    ax3.plot(oint_p1.y(delt)*ro, oint_p1.z(delt)*ro, "b--")
    ax3.plot(ointfwd_p1.y(deltfwd)*ro, ointfwd_p1.z(deltfwd)*ro, "b-")
    ax3.set_xlabel('Y (kpc)')
    ax3.set_ylabel('Z (kpc)')

    # Orbit radius vs. time
    ax4 = axs[1,1]
    ax4.plot(delt*to, oint_p1.r(delt)*ro, "k--")
    ax4.plot(deltfwd*to, ointfwd_p1.r(deltfwd)*ro, "k-")
    ax4.set_xlabel('time (Gyr)')
    ax4.set_ylabel('Radius (kpc)')

    print('Pericentre is at', min(oint_p1.r(delt)*ro), 'kpc')
    print('Apocentre is at', max(oint_p1.r(delt)*ro), 'kpc')

    return fig