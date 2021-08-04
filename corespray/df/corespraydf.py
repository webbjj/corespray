""" The corespraydf class
  
"""

__author__ = "Steffani Grondin & Jeremy J Webb"

__all__ = [
    "corespraydf",
]

from galpy.orbit import Orbit
from galpy.potential import MWPotential2014,PlummerPotential,KingPotential,MovingObjectPotential
from galpy.util import bovy_conversion,conversion,bovy_plot
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import animation

from clustertools import cart_to_cyl


class corespraydf(object):

	""" A class for initializing a distribution function for stars that are ejected from the core of a globular cluster

	Parameters
	----------

	gcorbit : string or galpy orbit instance
		Name of Galactic Globular Cluster from which to simulate core ejection or a Galpy orbit instance
	pot : galpy potential
		Potentional to be used for orbit integration (default: MWPotential2014)
	mu0 : float
		average 1D velocity in the core (default: 0 km/s)
	sig0 : float
		avarege 1D velocity disperions in the core (default 10.0 km/s)
	vesc0 : float
		escape velocity from the core (default: 10.0 km/s)
	rho0 : float
		core density (default: 1 Msun/pc^3)
	mgc : float
		globular cluster mass - needed if cluster's potential is to be included in orbit integration of escapers (default: None)
	rgc : float
		scale radius of globular cluster (assuming Plummer potential) or tidal radius of globular cluster (assuming King potential) (default: None)
	W0 : float
		King central potential parameter (default: None, which results in cluster potential taken to be a Plummer)
	mmin : float
		minimum stellar mass in core (default (0.1))
	mmax : float
		maximum stellar mass in the core (default: 1.4)
	alpha : float
		slope of the stellar mass function in the core (default: -1.35)
	emin : float
		minimum binary energy (default: 10.0**36 J)
	emax : float
		maximum binary energy (default: 10.0**40 J)
	ro : float
		galpy length scaling parameter (default: 8.)
	vo : float
		galpy veocicity scaling parameter (default: 220.)
	q : float
		exponenet for calculating probability of stellar escape from three-body system (#Equation 7.23)


	History
	-------
	2021 - Written - Grandin (UofT)

	"""

	def __init__(self,gcorbit,pot=MWPotential2014,mu0=0.,sig0=10.0,vesc0=10.0,rho0=1.,mgc=None,rgc=None,W0=None,mmin=0.1,mmax=1.4,alpha=-1.35,emin=10.0**36,emax=10.0**40.0,ro=8.,vo=220.,q=-3):

		if isinstance(gcorbit,str):
			self.gcname=gcorbit
			self.o = Orbit.from_name(self.gcname, ro=ro, vo=vo, solarmotion=[-11.1, 24.0, 7.25])
		else:
			self.gcname='unknown'
			self.o=gcorbit

		self.mu0,self.sig0,self.vesc0,self.rho0=mu0,sig0,vesc0,rho0

		self.mmin,self.mmax,self.alpha=mmin,mmax,alpha

		#Mean separation of star's in the core equal to twice the radius of a sphere that contains one star
		#Assume all stars in the core have mass equal to the mean mass
		masses=self._power_law_distribution_function(1000, self.alpha, self.mmin, self.mmax)
		self.mbar=np.mean(masses)
		self.rsep=((self.mbar/self.rho0)/(4.*np.pi/3.))**(1./3.)

		#Limits of binary energy distribution
		self.emin=emin
		self.emax=emax

		self.ro,self.vo=ro,vo
		self.to=conversion.time_in_Gyr(ro=self.ro,vo=self.vo)*1000.
		self.mo=bovy_conversion.mass_in_msol(ro=self.ro,vo=self.vo)


		self.q=q

		self.mwpot=pot

		if mgc is None:
			self.gcpot=None
		else:
			if W0 is None:
				ra=rgc/1.3
				self.gcpot=PlummerPotential(mgc/self.mo,ra/self.ro,ro=self.ro,vo=self.vo)
			else:
				self.gcpot=KingPotential(W0,mgc/self.mo,rgc/self.ro,ro=self.ro,vo=self.vo)

	def sample(self,tdisrupt=1000.,rate=1.,nstar=None,npeak=5.,binaries=False,verbose=False):
		""" A function for sampling the core ejection distribution function

		Parameters
		----------

		tdisrupt : float
			time over which sampling begins (Myr)
		rate : float
			ejection rate (default 1 per Myr)
		nstar : float
			if set, nstar stars will be ejected randomly from tdisrupt to 0 Myr. Rate is recalculated. (default : None)
		npeak : float
			when sampling kick velocity distribution function, sampling range will be from 0 to npeak*vpeak, where vpeak is the peak in the distribution function (default: 5)
		binaries : bool
			keep track of binaries that receive recoil kicks greater than the cluster's escape velocity (default : False)
		verbose : bool
			print additional information to screen (default: False)

		Returns
		----------



		History
		-------
		2021 - Written - Grandin/Webb (UofT)

		"""

		grav=4.302e-3 #pc/Msun (km/s)^2

		self.tdisrupt=tdisrupt


		#Select escape times
		#If nstar is not None, randomly select escapers between tstart and tend
		if nstar is not None:
			self.nstar=nstar
			self.rate=nstar/self.tdisrupt
		else:
			self.rate=rate
			self.nstar=self.tdisrupt*rate
			
		self.tesc=-1.*self.tdisrupt*np.random.rand(self.nstar)

		ts=np.linspace(0.,-1.*self.tdisrupt/self.to,1000)
		self.o.integrate(ts,self.mwpot)

		if self.gcpot is None:
			self.pot=self.mwpot
		else:
			moving_pot=MovingObjectPotential(self.o,self.gcpot,ro=self.ro,vo=self.vo)
			self.pot=[self.mwpot,moving_pot]
	    
	    #Generate kick velocities for escaped stars and binaries
		vxkick=np.zeros(self.nstar)
		vykick=np.zeros(self.nstar)
		vzkick=np.zeros(self.nstar)
	    
		vxkickb=np.zeros(self.nstar)
		vykickb=np.zeros(self.nstar)
		vzkickb=np.zeros(self.nstar)

		self.vesc=np.array([])

		nescape=0

		self.mstar=np.zeros(self.nstar)
		self.mb1=np.zeros(self.nstar)
		self.mb2=np.zeros(self.nstar)
		self.eb=np.zeros(self.nstar)


		while nescape < self.nstar:
			ms,m_a,m_b=self._power_law_distribution_function(3, self.alpha, self.mmin, self.mmax)   
			mb=m_a+m_b
			M=ms+mb

			prob=self._prob_three_body_escape(ms,m_a,m_b,self.q)

			if np.random.rand() < prob:
			    
				vxs,vys,vzs=np.random.normal(self.mu0,self.sig0,3)
				vstar=np.sqrt(vxs**2.+vys**2.+vzs**2.)
				vxb,vyb,vzb=np.random.normal(self.mu0,self.sig0,3)
				vbin=np.sqrt(vxb**2.+vyb**2.+vzb**2.)

				rdot=np.sqrt((vxs-vxb)**2.+(vys-vyb)**2.+(vzs-vzb)**2.)

				ebin,semi=self._sample_binding_energy(m_a,m_b,-1,self.emin,self.emax)

				e0=0.5*(mb*ms/M)*(rdot**2.)-grav*ms*mb/self.rsep + ebin

				vs=self._sample_escape_velocity(e0,ms,mb,npeak)

				if vs >self.vesc0:  

				    self.vesc=np.append(self.vesc,vs)

				    vxkick[nescape]=vs*(vxs/vstar)
				    vykick[nescape]=vs*(vys/vstar)
				    vzkick[nescape]=vs*(vzs/vstar)

				    #Check to see if recoil binary will also escape
				    #Binary kick velocity is calculated assuming total linear momentum of system sums to zero
				    vsb=vs*ms/mb

				    if vsb > self.vesc0:
					    vxkickb[nescape]=-vxkick[nescape]*ms/mb
					    vykickb[nescape]=-vykick[nescape]*ms/mb
					    vzkickb[nescape]=-vzkick[nescape]*ms/mb
				    else:
					    vxkickb[nescape]=-9999.
					    vykickb[nescape]=-9999.
					    vzkickb[nescape]=-9999.

				    self.mstar[nescape]=ms
				    self.mb1[nescape]=m_a
				    self.mb2[nescape]=m_b
				    self.eb[nescape]=ebin

				    nescape+=1

				if verbose: print('DEBUG: ',nescape,prob,vs,self.vesc0)
		
		Re0, phie0, ze0, vRe0, vTe0, vze0=np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])

		#Initial and final positions and velocities
		vxvv_i=[]
		vxvv_f=[]

		for i in range(0,self.nstar):
			xi,yi,zi=self.o.x(self.tesc[i]/self.to),self.o.y(self.tesc[i]/self.to),self.o.z(self.tesc[i]/self.to)
			vxi=vxkick[i]+self.o.vx(self.tesc[i]/self.to)
			vyi=vykick[i]+self.o.vy(self.tesc[i]/self.to)
			vzi=vzkick[i]+self.o.vz(self.tesc[i]/self.to)

			if verbose: print(i,self.tesc[i],xi,yi,zi,vxi,vyi,vzi)

			#Save initial positions and velocities
			Ri, phii, zi, vRi, vTi, vzi=cart_to_cyl(xi,yi,zi,vxi,vyi,vzi)
			vxvv_i.append([Ri/self.ro, vRi/self.vo, vTi/self.vo, zi/self.ro, vzi/self.vo, phii])

			#Integrate orbit from tesc to 0.
			os=Orbit(vxvv_i[-1],ro=self.ro,vo=self.vo,solarmotion=[-11.1, 24.0, 7.25])
			ts=np.linspace(self.tesc[i]/self.to,0.,1000)
			os.integrate(ts,self.pot)

			#Save final positions and velocities
			vxvv_f.append([os.R(0.)/self.ro,os.vR(0.)/self.vo,os.vT(0.)/self.vo,os.z(0.)/self.ro,os.vz(0.)/self.vo,os.phi(0.)])

		#Save initial and final positions and velocities of kicked stars at t=0 in orbit objects
		self.oi=Orbit(np.column_stack(vxvv_i),ro=self.ro,vo=self.vo,solarmotion=[-11.1, 24.0, 7.25])
		self.of=Orbit(np.column_stack(vxvv_f),ro=self.ro,vo=self.vo,solarmotion=[-11.1, 24.0, 7.25])

		if binaries:
			#Integrate orbits of binary star with kick velocities greater than the cluster's escape speed:
			#Initial and final positions and velocities
			vxvvb_i=[]
			vxvvb_f=[]

			for i in range(0,self.nstar):
				xi,yi,zi=self.o.x(self.tesc[i]/self.to),self.o.y(self.tesc[i]/self.to),self.o.z(self.tesc[i]/self.to)
				vxi=vxkickb[i]+self.o.vx(self.tesc[i]/self.to)
				vyi=vykickb[i]+self.o.vy(self.tesc[i]/self.to)
				vzi=vzkickb[i]+self.o.vz(self.tesc[i]/self.to)

				if verbose: print(i,self.tesc[i],xi,yi,zi,vxi,vyi,vzi)

				#Save initial positions and velocities
				Ri, phii, zi, vRi, vTi, vzi=cart_to_cyl(xi,yi,zi,vxi,vyi,vzi)
				vxvvb_i.append([Ri/self.ro, vRi/self.vo, vTi/self.vo, zi/self.ro, vzi/self.vo, phii])

				#Integrate orbit from tesc to 0. if kick velocity is higher than cluster's escape velocity

				if vxkickb[i]>-9999.:
					os=Orbit(vxvv_i[-1],ro=self.ro,vo=self.vo,solarmotion=[-11.1, 24.0, 7.25])
					ts=np.linspace(self.tesc[i]/self.to,0.,1000)
					os.integrate(ts,self.pot)
					vxvvb_f.append([os.R(0.)/self.ro,os.vR(0.)/self.vo,os.vT(0.)/self.vo,os.z(0.)/self.ro,os.vz(0.)/self.vo,os.phi(0.)])

				else:
					vxvvb_f.append(vxvvb_i[-1])

			self.obi=Orbit(np.column_stack(vxvvb_i),ro=self.ro,vo=self.vo,solarmotion=[-11.1, 24.0, 7.25])
			self.obf=Orbit(np.column_stack(vxvvb_f),ro=self.ro,vo=self.vo,solarmotion=[-11.1, 24.0, 7.25])

		return self.of

	def _prob_three_body_escape(self,ms,m_a,m_b,q):

		#Equation 7.23
		prob=(ms**q)/(ms**q+m_a**q+m_b**q)
		return prob

	def _sample_binding_energy(self,mb1,mb2,alpha,emin,emax):
	    #Opik's Law
		#Default binding energy distribution is:
		# power law of slope -1 
		# Between 10.0**36 and 10.0**40 J

		grav=4.302e-3 #pc/Msun (km/s)^2


		if isinstance(mb1,float):
			n=1
		else:
			n=len(mb1)

		ebin_si=self._power_law_distribution_function(n,alpha,emin,emax) #Joules = kg (m/s)^2
		ebin=ebin_si/1.9891e30 # Msun (m/s)^2
		ebin/=(1000.0*1000.0) #Msun (km/s)^2


		#Log normal a:
		semi_pc=(0.5*grav*mb1*mb2)/ebin
		semi_au=semi_pc/4.84814e-6    

		semi=semi_pc

		return ebin,semi

	def _sample_escape_velocity(self,e0,ms,mb,npeak=5):
		#randomly sample between npeak*vs_peak

		vs_peak=self._escape_velocity_distribution_peak(e0,ms,mb)
		match=False

		while not match:
			vstemp=np.random.rand()*npeak*vs_peak
			amptemp=np.random.rand()*vs_peak

			if amptemp < self._escape_velocity_distribution(vstemp,e0,ms,mb):
				vs=vstemp
				match=True

		return vs


	def _escape_velocity_distribution(self,vs,e0,ms,mb):
		#Equation 7.19
		M=ms+mb
		fv=(3.5*(np.fabs(e0)**(7./2.))*ms*M/mb)*vs/((np.fabs(e0)+0.5*(ms*M/mb)*(vs**2.))**(9./2.))
		return fv

	def _escape_velocity_distribution_peak(self,e0,ms,mb):
		M=ms+mb
		vs_peak=0.5*np.sqrt((M-ms)/(ms*M))*np.sqrt(np.fabs(e0))

		return vs_peak

	def _power_law_distribution_function(self,n,alpha,xmin,xmax):

	    eta=alpha+1.
	    
	    if xmin==xmax:
	        x=xmin
	    elif alpha==0:
	        x=xmin+np.random.random(n)*(xmax-xmin)
	    elif alpha>0:
	        x=xmin+np.random.power(eta,n)*(xmax-xmin)
	    elif alpha<0 and alpha!=-1.:
	        x=(xmin**eta + (xmax**eta - xmin**eta)*np.random.rand(n))**(1./eta)
	    elif alpha==-1:
	        x=np.log10(xmin)+np.random.random(n)*(np.log10(xmax)-np.log10(xmin))
	        x=10.0**x
	        
	    if n==1:
	        return x
	    else:      
	        return np.array(x)

	def _init_fig(self,xlim=(-20,20),ylim=(-20,20)):
	    self.fig = plt.figure()
	    self.ax = plt.axes(xlim=xlim, ylim=ylim)
	    self.ax.set_xlabel('X (kpc)')
	    self.ax.set_ylabel('Y (kpc)')
	    self.txt_title=self.ax.set_title('')
	    self.line, = self.ax.plot([], [], lw=2)
	    self.pt, = self.ax.plot([],[],'.')

	def _set_data(self,gcdata,sdata):
	    self.gcdata = gcdata
	    self.sdata=sdata

	def _ani_init(self):
	    self.line.set_data([], [])
	    self.pt.set_data([],[])
	    return self.line,self.pt

	def _ani_update(self, i):

		if i < 5:
		    x = self.gcdata[0:i+1,0]
		    y = self.gcdata[0:i+1,1]
		else:
			x = self.gcdata[i-5:i+1,0]
			y = self.gcdata[i-5:i+1,1]	
		self.line.set_data(x, y)

		escindx=self.tesc/self.to <= self.ts[i]


		if np.sum(escindx)>0:
			self.pt.set_data(self.sdata[i][0][escindx],self.sdata[i][1][escindx])
		else:
			self.pt.set_data([],[])

		self.txt_title.set_text('%s' % str (self.ts[i]*self.to))

		return self.line,self.pt


	def animate(self,frames=100,interval=50,xlim=(-20,20),ylim=(-20,20)):

		"""Animate the ejection of stars from the cluster's core
		
		Parameters
    	----------

    	frames : int
    		number of frames to use for animation (default:100)
		interval : float
			time intercal between frames (default: 50 Myr)
		xlim : tuple
			xlimits for figure
		ylim : tuple
			ylimts for figure

	    History
   		-------
	    2021 - Written - Webb (UofT)

	    """

		self._init_fig(xlim,ylim)

		self.ts=np.linspace(-1.*self.tdisrupt/self.to,0.,frames)
		tsint=np.linspace(0.,-1.*self.tdisrupt/self.to,1000)
		self.oe.integrate(tsint,self.pot)

		gcdata=np.zeros(shape=(frames,2))

		for i in range(0,frames):
		    gcdata[i]=[self.o.x(self.ts[i]),self.o.y(self.ts[i])]

		sdata=np.zeros(shape=(frames,2,self.nstar))

		for i in range(0,frames):
			sdata[i]=[self.oe.x(self.ts[i]),self.oe.y(self.ts[i])]

		self._set_data(gcdata,sdata)

		self.anim = animation.FuncAnimation(self.fig, self._ani_update, init_func=self._ani_init, frames=frames, interval=interval, blit=False)

	def snapout(self,filename='corespray.dat'):
		"""Output present day positions, velocities, escape times, and escape velocities of stars
		
		Parameters
    	----------

    	filename: str
    		file name to write data to (default: corespray.dat)

	    History
   		-------
	    2021 - Written - Webb (UofT)

	    """
		R=np.append(self.o.R(0.),self.oe.R(0.))
		vR=np.append(self.o.vR(0.),self.oe.vR(0.))
		vT=np.append(self.o.vT(0.),self.oe.vT(0.))
		z=np.append(self.o.z(0.),self.oe.z(0.))
		vz=np.append(self.o.vz(0.),self.oe.vz(0.))
		phi=np.append(self.o.phi(0.),self.oe.phi(0.))

		vesc=np.append(0.,self.vesc)
		tesc=np.append(0.,self.tesc)

		np.savetxt(filename,np.column_stack([R,vR,vT,z,vz,phi,vesc,tesc]))
