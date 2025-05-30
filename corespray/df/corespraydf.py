""" The corespraydf class
  
"""

__author__ = "Steffani M Grondin & Jeremy J Webb"

__all__ = [
    "corespraydf",
]

from galpy.orbit import Orbit
from galpy.potential import MWPotential2014,PlummerPotential,KingPotential,MovingObjectPotential
from galpy.util import conversion,coords
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import animation

import time
from tqdm import tqdm


class corespraydf(object):

    """ A class for initializing a distribution function for stars that are ejected from the core of a globular cluster
    -- If emin and emax are None, assume limits are between twice the hard-soft boundary and twice the contact boundary between two solar mass stars
    Parameters
    ----------

    gcorbit : string or galpy orbit instance
        Name of Galactic Globular Cluster from which to simulate core ejection or a Galpy orbit instance
    pot : galpy potential
        Potentional to be used for orbit integration (default: MWPotential2014)
    mgc : float
        globular cluster mass - needed if cluster's potential is to be included in orbit integration of escapers (default: None)
    rgc : float
        half-mass radius of globular cluster (assuming Plummer potential) or tidal radius of globular cluster (assuming King potential) (default: None)
    W0 : float
        King central potential parameter (default: None, which results in cluster potential taken to be a Plummer)
    ro : float
        galpy length scaling parameter (default: 8.)
    vo : float
        galpy veocicity scaling parameter (default: 220.)
    verbose : bool
        print additional information to screen (default: False)
    timing : bool
        print timing information to screen (default: False)
    show_progress : bool
        show progress bar when integrating orbits (default: False)


    History
    -------
    2021 - Written - Grondin (UofT)
    2024 - Updated - Evans (UofT)

    """

    def __init__(self,gcorbit,pot=MWPotential2014,mgc=None,rgc=None,W0=None,ro=8.,vo=220.,verbose=False, timing=False, show_progress=False):

        if isinstance(gcorbit,str):
            self.gcname=gcorbit
            self.o = Orbit.from_name(self.gcname, ro=ro, vo=vo, solarmotion=[-11.1, 24.0, 7.25])
        else:
            self.gcname='unknown'
            self.o=gcorbit

        self.ro,self.vo=ro,vo
        self.to=conversion.time_in_Gyr(ro=self.ro,vo=self.vo)*1000.
        self.mo=conversion.mass_in_msol(ro=self.ro,vo=self.vo)

        self.mwpot=pot

        if mgc is None:
            self.gcpot=None
        else:
            if W0 is None:
                ra=rgc/1.3
                self.gcpot=PlummerPotential(mgc/self.mo,ra/self.ro,ro=self.ro,vo=self.vo)
            else:
                self.gcpot=KingPotential(W0,mgc/self.mo,rgc/self.ro,ro=self.ro,vo=self.vo)

        self.binaries=False

        self.timing=timing
        self.show_progress=show_progress

    def _get_K(self, mass):
        """ Calculate the Hurley type K for a given mass """

        if mass > 0.7:
            return 1
        else:
            return 0
        
    def sample_three_body(self,tdisrupt=1000.,rate=1.,nstar=None,mu0=0.,sig0=10.0,vesc0=10.0,rho0=1.,mmin=0.1,mmax=1.4,alpha=-1.35,seps=None, masses=None,m1=None,m2a=None,m2b=None,emin=None,emax=None,balpha=-1,q=-3, npeak=5.,binaries=False,verbose=False, show_progress=None, **kwargs):
        """ A function for sampling the three-body interaction core ejection distribution function

        Parameters
        ----------

        tdisrupt : float
            time over which sampling begins (Myr)
        rate : float
            ejection rate (default 1 per Myr)
        nstar : float
            if set, nstar stars will be ejected randomly from tdisrupt to 0 Myr. Rate is recalculated. (default : None)
        mu0 : float
            average 1D velocity in the core (default: 0 km/s)
        sig0 : float
            avarege 1D velocity disperions in the core (default 10.0 km/s)
        vesc0 : float
            escape velocity from the core (default: 10.0 km/s)
        rho0 : float
            core density (default: 1 Msun/pc^3)
        mgc : float
            globular cluster mass in solar masses - needed if cluster's potential is to be included in orbit integration of escapers (default: None)
        rgc : float
            half-mass radius of globular cluster (assuming Plummer potential) or tidal radius of globular cluster (assuming King potential) in kpc (default: None)
        W0 : float
            King central potential parameter (default: None, which results in cluster potential taken to be a Plummer)
        mmin : float
            minimum stellar mass in core (default (0.1 Msun))
        mmax : float
            maximum stellar mass in the core (default: 1.4 Msun)
        alpha : float
            slope of the stellar mass function in the core (default: -1.35)
        masses : None or array of floats (Msun)
            array of masses to be used instead of drawing for a power-law mass function (default: None)
            Note : mmin, mmax, and alpha will be overridden
        m1 : None or float or array of floats (Msun)
            fixed mass or array of masses for drawing single kickedstar.
            Overrides masses (for the star picking) if provided (default: None)
            Note : (mmin, mmax, alpha) or (masses) must still be provided to determine the mean mass in the core
        m2a, m2b : None or floats or arrays of floats
            Binary star A mass and star B mass (Msun)  (default: None)
                - If None, will be sampled from masses or power-law distribution
                - If floats, will be the fixed masses for the binary stars
                - If arrays, binary masses will be randomly sampled from the arrays  (ararys must be same length)
            Note : (mmin, mmax, alpha) or (masses) must still be provided to determine the mean mass in the core
        seps : None or float or array of floats
            semi-major axis of the binary at time of ejection (default: None)
                - If None, will be sampled from the hard-soft boundary to the contact boundary
                - If float, will be the fixed semi-major axis for the binary stars
                - If array, will be randomly sampled from the array (array must be same length as m2a and m2b) 
            Note : if seps is not None, then emin, emax, a_min, and a_max do not matter
        emin : float
            minimum binary energy in Joules (default: None)
        emax : float
            maximum binary energy in Joules (default: None)
        balpha : float
            power-law slope of initial binary binding energy distribution (default: -1)
        q : float
            exponenet for calculating probability of stellar escape from three-body system (#Equation 7.23) (default: -3)
        npeak : float
            when sampling kick velocity distribution function, sampling range will be from 0 to npeak*vpeak, where vpeak is the peak in the distribution function (default: 5)
        binaries : bool
            keep track of binaries that receive recoil kicks greater than the cluster's escape velocity (default : False)
        verbose : bool
            print additional information to screen (default: False)

        Key Word Arguments
        ----------
        nrandom : int
            Nunber of random numbers to sample in a given batch
        dt : float
            Fixed time step for orbit integration (default: 0.1 Myr)
        rsample : bool
            Sample separation between single star and binary within core (default: False)
        nrsep : float
            Number of mean separations to sample out to when sampling separation between single and binary stars (default : 2)
        initialize : bool
            initialize orbits only, do not integrate (default:False)
            Note if initialize == True then the initial, un-integrated orbits, will be returned
        escape_only : bool
            Keep track of kicked single stars with vkick > vesc only (default:True)
        method : str
            Method to use for orbit integration (default: None, which uses galpy's default method)
        K1, K2a, K2b : arrays of integers
            Hurley integer stellar types (https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract) associated with the m1, m2a, and m2b arrays, respectively. Helps track the kinds of stars involved in the interaction.
            If provided, must be the same length as the corresponding mass array.
            If not provided, will be calculated based on the mass of the star assuming it's main sequence.

        Returns
        ----------
        of : orbit
            galpy orbit instance for kicked stars

        if binaries:
            obf : orbit
                galpy orbit instance for recoil binary stars

        Also sets the following attributes of the corespraydf instance, one 
        for each escaped star (or kicked star if escape_only is False):
        vesc : array
            escape velocities of stars from the cluster core (km/s)
        dr : array
            pericentre separation between single star and binary at time of ejection (pc)
        mstar, mb1, mb2 : arrays
            masses of the escaped star and each binary star (Msun)
        Kstar, Kb1, Kb2 : arrays
            Hurley integer stellar types (https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract) of the escaped star and each binary star
        semi, semif : arrays
            semi-major axes of the binary at time of ejection and after recoil kick (pc)
        eb, ebf : arrays
            binding energies of the binary at time of ejection and after recoil kick (Msun (km/s)^2)
        e0 : array
            initial energy of the whole three-body system at time of ejection (Msun (km/s)^2)
        sindx : array
            boolean array indicating whether the star escaped with a kick velocity greater than the escape velocity of the core
        bindx : array
            boolean array indicating whether the binary escaped with a kick velocity greater than the escape velocity of the core (only if binaries=True)
        vescb : array
            escape velocities of binaries from the cluster core (km/s) (only if binaries=True)
        binx : array
            boolean array indicating whether the binary received a recoil kick greater than the escape velocity of the core (only if binaries=True)

        History
        -------
        2021 - Written - Grondin/Webb (UofT)
        2024 - Updated - Evans (UofT)
        """

        grav=4.302e-3 #pc/Msun (km/s)^2
        msolar=1.9891e30

        nrandom=kwargs.get('nrandom',1000)
        self.timing=kwargs.get('timing',self.timing)
        self.show_progress=kwargs.get('show_progress',self.show_progress)

        if show_progress is not None:
            self.show_progress=show_progress

        self.tdisrupt=tdisrupt

        rsample=kwargs.get('rsample',False)
        initialize=kwargs.get('initialize',False)

        #Select escape times
        #If nstar is not None, randomly select escapers between tstart and tend
        if nstar is not None:
            self.nstar=nstar
            self.rate=nstar/self.tdisrupt
        else:
            self.rate=rate
            self.nstar=self.tdisrupt*rate
            
        self.tesc=-1.*self.tdisrupt*np.random.rand(self.nstar)

        ntstep=kwargs.get('ntstep',10000)
        dt = kwargs.get('dt',0.1)
        nrsep=kwargs.get('nrsep',1)
        method=kwargs.get('method',None)

        ntstep = np.ceil(((self.tdisrupt/self.to)/dt)).astype(int)

        ts=np.linspace(0.,-1.*self.tdisrupt/self.to,ntstep)

        if method is not None:
            self.o.integrate(ts,self.mwpot,method=method)
        else:
            self.o.integrate(ts,self.mwpot)

        if self.gcpot is None:
            self.pot=self.mwpot
        else:
            moving_pot=MovingObjectPotential(self.o,self.gcpot,ro=self.ro,vo=self.vo)
            self.pot=[self.mwpot,moving_pot]
        
        self.mu0,self.sig0,self.vesc0,self.rho0=mu0,sig0,vesc0,rho0

        if masses is None:
            self.mmin,self.mmax,self.alpha=mmin,mmax,alpha
            #Mean separation of star's in the core equal to twice the radius of a sphere that contains one star
            #Assume all stars in the core have mass equal to the mean mass
            self.masses=self._power_law_distribution_function(1000, self.alpha, self.mmin, self.mmax)
        else:
            self.masses=masses

        self.mbar=np.mean(self.masses)
        self.rsep=((self.mbar/self.rho0)/(4.*np.pi/3.))**(1./3.)

        #Limits of binary energy distribution
        #If emin and emax are None, assume limits are between twice the hard-soft boundary and twice the contact boundary between two solar mass stars
        #If seps is not None, then the semi-major axis of the binary will be sampled from seps and a_min/a_max/e_min/e_max do not matter

        a_hs=grav*self.mbar/(sig0**2.) # pc
        a_max=2.*a_hs
        a_min=(4.0/215.032)*4.84814e-6 #pc

        if emin is None:
            e_min=grav*(self.mbar**2.)/(2.0*a_max) #Msun (km/s)**2
            e_min*=(1000.0**2.)
            self.emin=e_min*msolar
        else:
            self.emin=emin

        if emax is None:
            e_max=grav*(self.mbar**2.)/(2.0*a_min) #Msun (km/s)**2
            e_max*=(1000.0**2.)
            self.emax=e_max*msolar
        else:
            self.emax=emax

        if verbose:
            print('Sample Binary Energies between: ',self.emin,' and ',self.emax,' J')

        self.q=q

        #Generate kick velocities for escaped stars and binaries
        vxkick=np.zeros(self.nstar)
        vykick=np.zeros(self.nstar)
        vzkick=np.zeros(self.nstar)
        
        vxkickb=np.zeros(self.nstar)
        vykickb=np.zeros(self.nstar)
        vzkickb=np.zeros(self.nstar)

        self.vesc=np.array([])

        self.dr=np.array([])

        nescape=0

        self.mstar=np.zeros(self.nstar)
        self.mb1=np.zeros(self.nstar)
        self.mb2=np.zeros(self.nstar)
        self.Kstar=np.zeros(self.nstar)
        self.Kb1=np.zeros(self.nstar)
        self.Kb2=np.zeros(self.nstar)
        self.semi=np.zeros(self.nstar)
        self.semif=np.zeros(self.nstar)
        self.eb=np.zeros(self.nstar)
        self.ebf=np.zeros(self.nstar)
        self.vperi = np.zeros(self.nstar)

        self.e0=np.zeros(self.nstar)

        self.sindx=np.ones(self.nstar,dtype=bool)

        if binaries:
            self.binaries=True
            self.bindx=np.zeros(self.nstar,dtype=bool)
            self.vescb=np.array([])
        else:
            self.binaries=False


        if self.timing: dttime=time.time()

        escape_only=kwargs.get('escape_only',True)

        while nescape < self.nstar:

                if m1 is None:
                    if masses is None:
                        ms=self._power_law_distribution_function(1, self.alpha, self.mmin, self.mmax)
                    else:
                        ms=np.random.choice(self.masses,1)
                    K1 = self._get_K(ms)
                elif isinstance(m1,float) or isinstance(m1,int):
                    ms=m1
                    K1 = self._get_K(ms)
                else:
                    ind = np.random.choice(range(0,len(m1)),1)
                    ms = m1[ind]
                    
                    if kwargs.get('K1',None) is not None:
                        if len(kwargs['K1']) != len(m1):
                            raise ValueError('K1 must be the same length as m1 if it is an array')
                        K1 = kwargs['K1'][ind]
                    else:
                        K1 = self._get_K(m1[ind])

                    

                if isinstance(m2a,np.ndarray) and isinstance(m2b,np.ndarray):
                    if len(m2a) != len(m2b):
                        raise ValueError('m2a and m2b must be the same length if they are arrays')
                        
                    randind = np.random.choice(range(0,len(m2b)),1)
                    m_a = m2a[randind]
                    m_b = m2b[randind]

                    if kwargs.get('K2a',None) is not None:
                        if len(kwargs.get('K2a')) != len(m2a):
                            raise ValueError('K2a must be the same length as m2a if it is an array')
                        K2_a = kwargs['K2a'][randind]
                    else:
                        K2_a = self._get_K(m_a)

                    if kwargs.get('K2b',None) is not None:
                        if len(kwargs.get('K2b')) != len(m2b):
                            raise ValueError('K2b must be the same length as m2b if it is an array')
                        K2_b = kwargs['K2b'][randind]

                    else:
                        K2_b = self._get_K(m_b)
                else:
                    if m2a is None:
                        if masses is None:
                            m_a=self._power_law_distribution_function(1, self.alpha, self.mmin, self.mmax)
                        else:
                            m_a=np.random.choice(self.masses,1)

                        K2_a = self._get_K(m_a)                               

                    elif isinstance(m2a,float) or isinstance(m2a,int):
                        m_a=m2a
                        K2_a = self._get_K(m_a)
                    else:
                        ind = np.random.choice(range(0,len(m2a)),1)
                        m_a = m2a[ind]
                        if kwargs.get('K2a',None) is not None:
                            if len(kwargs['K2a']) != len(m2a):
                                raise ValueError('K2a must be the same length as m2a if it is an array')
                            K2_a = kwargs['K2a'][ind]
                        else:
                            K2_a = self._get_K(m_a)

                    if m2b is None:
                        if masses is None:
                            m_b=self._power_law_distribution_function(1, self.alpha, self.mmin, self.mmax)
                        else:
                            m_b=np.random.choice(self.masses,1)
                        K2_b = self._get_K(m_b)
                    elif isinstance(m2b,float) or isinstance(m2b,int) :
                        m_b=m2b
                        K2_b = self._get_K(m_b)
                    else:
                        ind = np.random.choice(range(0,len(m2b)),1)
                        m_b = m2b[ind]
                        if kwargs.get('K2b',None) is not None:
                            if len(kwargs['K2b']) != len(m2b):
                                raise ValueError('K2b must be the same length as m2b if it is an array')
                            K2_b = kwargs['K2b'][ind]
                        else:
                            K2_b = self._get_K(m_b)

                mb=m_a+m_b
                M=ms+mb

                prob=self._prob_three_body_escape(ms,m_a,m_b,self.q)

                if np.random.rand() < prob:
                
                    vxs,vys,vzs=np.random.normal(self.mu0,self.sig0,3)
                    vstar=np.sqrt(vxs**2.+vys**2.+vzs**2.)
                    vxb,vyb,vzb=np.random.normal(self.mu0,self.sig0,3)
                    vbin=np.sqrt(vxb**2.+vyb**2.+vzb**2.)

                    rdot=np.sqrt((vxs-vxb)**2.+(vys-vyb)**2.+(vzs-vzb)**2.)

                    if seps is None:
                        ebin,semi=self._sample_binding_energy(m_a,m_b,balpha,self.emin,self.emax)
                    elif 'randind' in locals():
                        semi = seps[randind] * 4.84814e-6 #AU to pc
                        ebin = grav * m_a * m_b / semi #Msun (km/s)^2
                    elif isinstance(seps,float) or isinstance(seps,int):
    
                        semi = seps * 4.84814e-6
                        ebin = grav * m_a * m_b / semi #Msun (km/s)^2
                    else:
                        semi = np.random.choice(seps,1) * 4.84814e-6
                        ebin = grav * m_a * m_b / semi #Msun (km/s)^2

                    #ebin,semi=self._sample_binding_energy(m_a,m_b,balpha,self.emin,self.emax)

                    if rsample:
                        #Sample between semi-major axis of the binary and nrsep * mean separation in the core
                        rs=semi+np.random.rand()*(nrsep*(self.rsep-semi))
                        phis=2.0*np.pi*np.random.rand()
                        thetas=np.arccos(1.0-2.0*np.random.rand())

                        xs=rs*np.cos(phis)*np.sin(thetas)
                        ys=rs*np.sin(phis)*np.sin(thetas)
                        zs=rs*np.cos(thetas)

                        rb=a_max/2.+np.random.rand()*(nrsep*self.rsep/2.-a_max/2.)
                        phib=2.0*np.pi*np.random.rand()
                        thetab=np.arccos(1.0-2.0*np.random.rand())

                        xb=rb*np.cos(phib)*np.sin(thetab)
                        yb=rb*np.sin(phib)*np.sin(thetab)
                        zb=rb*np.cos(thetab)

                        dr=np.sqrt((xs-xb)**2.+(ys-yb)**2.+(zs-zb)**2.)

                        e0=0.5*(mb*ms/M)*(rdot**2.)-grav*ms*mb/dr - ebin
                    else:
                        e0=0.5*(mb*ms/M)*(rdot**2.)-grav*ms*mb/self.rsep - ebin
                        dr=self.rsep

                    if kwargs.get('e0bug',False):
                        print('UNFIXING EB')
                        e0+=(2.0*ebin)

                    vs=self._sample_escape_velocity(e0,ms,mb,npeak,nrandom)

                    if (vs >= self.vesc0 or not escape_only):
                        
                        if vs<self.vesc0:
                            self.sindx[nescape]=False

                        self.vesc=np.append(self.vesc,vs)
                        self.dr=np.append(self.dr,dr)
                        vxkick[nescape]=vs*(vxs/vstar)
                        vykick[nescape]=vs*(vys/vstar)
                        vzkick[nescape]=vs*(vzs/vstar)

                        pxi=ms*vxs+mb*vxb
                        pyi=ms*vys+mb*vyb
                        pzi=ms*vzs+mb*vzb

                        vxkickb[nescape]=(pxi-ms*vxkick[nescape])/mb
                        vykickb[nescape]=(pyi-ms*vykick[nescape])/mb
                        vzkickb[nescape]=(pzi-ms*vzkick[nescape])/mb

                        vsb=np.sqrt(vxkickb[nescape]**2.+ vykickb[nescape]**2.+ vzkickb[nescape]**2.)

                        if binaries:
                            #Check to see if recoil binary will also escape
                            #Binary kick velocity is calculated assuming total linear momentum of system sums to zero

                            self.vescb=np.append(self.vescb,vsb)

                            if vsb > self.vesc0:
                                self.bindx[nescape]=True

                        self.mstar[nescape]=ms
                        self.mb1[nescape]=m_a
                        self.mb2[nescape]=m_b
                        self.Kstar[nescape]=K1
                        self.Kb1[nescape]=K2_a
                        self.Kb2[nescape]=K2_b
                        self.eb[nescape]=ebin
                        self.e0[nescape]=e0
                        self.semi[nescape]=semi

                        #Final relative velocity between single and binary
                        rdotf=vs+vsb
                        #Final binding energy and semi major axis of the binary
                        ebf=e0-0.5*(mb*ms/M)*(rdotf**2.)
                        self.ebf[nescape]=ebf

                        semif = -1*(0.5*grav*m_a*m_b)/ebf

                        self.semif[nescape] = semif

                        #seps[randind] = semif / 4.84814e-6 #AU to pc
                   

                        #Make sure new semi-major axis is less than original
                        #print(e0,rdot,ebin,ebin/e0,self.semi[nescape],rdotf,ebf,ebf/e0,self.semif[nescape])
                        #assert self.semif[nescape] <= self.semi[nescape]

                        nescape+=1

                    if verbose: print('Sampling: ',nescape,prob,vs,self.vesc0)
        
        if self.timing: print(nescape,' three body encounters simulated in ', time.time()-dttime,' s')

        if initialize:
            self.oi=self._initialize_orbits(vxkick,vykick,vzkick,False,verbose,**kwargs)
            self.of=None
            if binaries:
                self.obi=self._initialize_orbits(vxkickb,vykickb,vzkickb,False,verbose,**kwargs)
                self.obf=None

            if binaries:
                return self.oi,self.obi
            else:
                return self.oi


        else:
            self.oi,self.of=self._integrate_orbits(vxkick,vykick,vzkick,False,verbose,**kwargs)

            if binaries:
                self.obi,self.obf=self._integrate_orbits(vxkickb,vykickb,vzkickb,binaries,verbose,**kwargs)

            if binaries:
                return self.of,self.obf
            else:
                return self.of

    def _integrate_orbits(self,vxkick,vykick,vzkick,binaries=False,verbose=False,**kwargs):
        
        #Set integration method (see https://docs.galpy.org/en/v1.7.2/orbit.html)
        method=kwargs.get('method',None)
        #Add a minimal offset to prevent stars from being initialized at r=0 in a the cluster.
        offset=kwargs.get('offset',1e-9)
        xoffsets=np.random.normal(0.0,offset,len(vxkick))
        yoffsets=np.random.normal(0.0,offset,len(vykick))
        zoffsets=np.random.normal(0.0,offset,len(vzkick))

        Re0, phie0, ze0, vRe0, vTe0, vze0=np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])

        #Initial and final positions and velocities
        vxvv_i=[]
        vxvv_f=[]

        if self.timing: dttottime=time.time()

        if binaries:
            tqdmlabel='Integrating binary orbits...'
        else:
            tqdmlabel='Integrating single star orbits...'

        iterator = tqdm(range(0,len(self.tesc)), desc=tqdmlabel) \
                    if self.show_progress else range(0,len(self.tesc))

        for i in iterator:
            if self.timing: dttime=time.time()

            xi,yi,zi=self.o.x(self.tesc[i]/self.to)+xoffsets[i],self.o.y(self.tesc[i]/self.to)+yoffsets[i],self.o.z(self.tesc[i]/self.to)+zoffsets[i]
            vxi=vxkick[i]+self.o.vx(self.tesc[i]/self.to)
            vyi=vykick[i]+self.o.vy(self.tesc[i]/self.to)
            vzi=vzkick[i]+self.o.vz(self.tesc[i]/self.to)

            if verbose: print(i,self.tesc[i],xi,yi,zi,vxi,vyi,vzi)

            #Save initial positions and velocities

            Ri, phii, zi = coords.rect_to_cyl(xi, yi, zi)
            vRi, vTi, vzi = coords.rect_to_cyl_vec(vxi, vyi, vzi, xi, yi, zi)

            vxvv_i.append([Ri/self.ro, vRi/self.vo, vTi/self.vo, zi/self.ro, vzi/self.vo, phii])

            if not binaries or (binaries and self.bindx[i]):
                #Integrate orbit from tesc to 0.
                os=Orbit(vxvv_i[-1],ro=self.ro,vo=self.vo,solarmotion=[-11.1, 24.0, 7.25])

                ntstep=kwargs.get('ntstep',10000)
                dt = kwargs.get('dt',0.01)

                ntstep = np.ceil(((-self.tesc[i]/self.to)/dt)).astype(int)

                if ntstep < 100:
                    ntstep = 100

                ts=np.linspace(self.tesc[i]/self.to,0.,ntstep)

                if method is None:
                    os.integrate(ts,self.pot) 
                else:
                    os.integrate(ts,self.pot,method=method,progressbar=True)

                #Save final positions and velocities
                vxvv_f.append([os.R(0.)/self.ro,os.vR(0.)/self.vo,os.vT(0.)/self.vo,os.z(0.)/self.ro,os.vz(0.)/self.vo,os.phi(0.)])

                if self.timing: print('ORBIT ',i,' INTEGRATED FROM %f WITH VK= %f KM/S IN' % (self.tesc[i],self.vesc[i]),time.time()-dttime,' s (Rp=%f)' % os.rperi())
                
            elif binaries and not self.bindx[i]:
                vxvv_f.append([self.o.R(0.)/self.ro,self.o.vR(0.)/self.vo,self.o.vT(0.)/self.vo,self.o.z(0.)/self.ro,self.o.vz(0.)/self.vo,self.o.phi(0.)])



        if self.timing: print('ALL ORBITS INTEGRATED IN',time.time()-dttottime,' s' )

        #Save initial and final positions and velocities of kicked stars at t=0 in orbit objects
        oi=Orbit(vxvv_i,ro=self.ro,vo=self.vo,solarmotion=[-11.1, 24.0, 7.25])
        of=Orbit(vxvv_f,ro=self.ro,vo=self.vo,solarmotion=[-11.1, 24.0, 7.25])

        return oi,of

    def _initialize_orbits(self,vxkick,vykick,vzkick,binaries=False,verbose=False,**kwargs):
        
        #Set integration method (see https://docs.galpy.org/en/v1.7.2/orbit.html)
        method=kwargs.get('method',None)
        #Add a minimal offset to prevent stars from being initialized at r=0 in a the cluster.
        offset=kwargs.get('offset',1e-9)
        xoffsets=np.random.normal(0.0,offset,len(vxkick))
        yoffsets=np.random.normal(0.0,offset,len(vykick))
        zoffsets=np.random.normal(0.0,offset,len(vzkick))

        Re0, phie0, ze0, vRe0, vTe0, vze0=np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])

        #Initial and final positions and velocities
        vxvv_i=[]
        vxvv_f=[]

        if self.timing: dttottime=time.time()


        for i in range(0,self.nstar):
            if self.timing: dttime=time.time()

            xi,yi,zi=self.o.x(self.tesc[i]/self.to)+xoffsets[i],self.o.y(self.tesc[i]/self.to)+yoffsets[i],self.o.z(self.tesc[i]/self.to)+zoffsets[i]
            vxi=vxkick[i]+self.o.vx(self.tesc[i]/self.to)
            vyi=vykick[i]+self.o.vy(self.tesc[i]/self.to)
            vzi=vzkick[i]+self.o.vz(self.tesc[i]/self.to)

            if verbose: print(i,self.tesc[i],xi,yi,zi,vxi,vyi,vzi)

            #Save initial positions and velocities

            Ri, phii, zi = coords.rect_to_cyl(xi, yi, zi)
            vRi, vTi, vzi = coords.rect_to_cyl_vec(vxi, vyi, vzi, xi, yi, zi)

            vxvv_i.append([Ri/self.ro, vRi/self.vo, vTi/self.vo, zi/self.ro, vzi/self.vo, phii])

        #Save initial and final positions and velocities of kicked stars at t=0 in orbit objects
        oi=Orbit(vxvv_i,ro=self.ro,vo=self.vo,solarmotion=[-11.1, 24.0, 7.25])

        return oi

    def _prob_three_body_escape(self,ms,m_a,m_b,q):

        #Equation 7.23
        prob=(ms**q)/(ms**q+m_a**q+m_b**q)
        return prob

    def _sample_binding_energy(self,mb1,mb2,balpha,emin,emax):
        #Opik's Law
        #Default binding energy distribution is:
        # power law of slope -1 
        # Between 10.0**36 and 10.0**40 J

        grav=4.302e-3 #pc/Msun (km/s)^2


        if isinstance(mb1,float) or isinstance(mb1,int):
            n=1
        else:
            n=len(mb1)

        ebin_si=self._power_law_distribution_function(n,balpha,emin,emax) #Joules = kg (m/s)^2
        ebin=ebin_si/1.9891e30 # Msun (m/s)^2
        ebin/=(1000.0*1000.0) #Msun (km/s)^2


        #Log normal a:
        semi_pc=(0.5*grav*mb1*mb2)/ebin
        semi_au=semi_pc/4.84814e-6    

        semi=semi_pc

        return ebin,semi

    def _sample_escape_velocity(self,e0,ms,mb,npeak=5,nrandom=1000):
        #randomly sample between npeak*vs_peak

        vs_peak=self._escape_velocity_distribution_peak(e0,ms,mb)

        vstemp=np.random.rand(nrandom)*npeak*vs_peak
        probs = self._escape_velocity_distribution(np.ones(nrandom)*vstemp,np.ones(nrandom)*e0,np.ones(nrandom)*ms,np.ones(nrandom)*mb)
 
        if np.sum(probs)==0:

            return -99
        else:
            v = np.random.choice(vstemp,1,p=probs/np.sum(probs))
            return v[0]

    def _escape_velocity_distribution(self,vs,e0,ms,mb):
        #Equation 7.19
        M=ms+mb

        fv=(3.5*(np.fabs(e0)**(7./2.))*ms*M/mb)*vs/((np.fabs(e0)+0.5*(ms*M/mb)*(vs**2.))**(9./2.))
        
        return fv

    def _escape_velocity_distribution_peak(self,e0,ms,mb):
        M=ms+mb
        vs_peak=0.5*np.sqrt((M-ms)/(ms*M))*np.sqrt(np.fabs(e0))

        return vs_peak

    def sample_uniform(self,tdisrupt=1000.,rate=1.,nstar=None,vmin=0.,vmax=500.,verbose=False, **kwargs):
        """ A function for sampling a uniform core ejection distribution function

        Parameters
        ----------

        tdisrupt : float
            time over which sampling begins (Myr)
        rate : float
            ejection rate (default 1 per Myr)
        nstar : float
            if set, nstar stars will be ejected randomly from tdisrupt to 0 Myr. Rate is recalculated. (default : None)
        vmin : float
            minimum kick velocity
        vmax : float
            maximum kick velocity
        verbose : bool
            print additional information to screen (default: False)

        Returns
        ----------
        of : orbit
            galpy orbit instance for kicked stars


        History
        -------
        2021 - Written - Grondin/Webb (UofT)

        """

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

        ntstep=kwargs.get/('ntstep',10000)
        ts=np.linspace(0.,-1.*self.tdisrupt/self.to,ntstep)
        self.o.integrate(ts,self.mwpot)

        if self.gcpot is None:
            self.pot=self.mwpot
        else:
            moving_pot=MovingObjectPotential(self.o,self.gcpot,ro=self.ro,vo=self.vo)
            self.pot=[self.mwpot,moving_pot]
        
        #Generate kick velocities for escaped stars and binaries
        self.vesc=vmin+(vmax-vmin)*np.random.rand(self.nstar)


        #Assume a random direction
        vxs=np.random.normal(0,1.,self.nstar)
        vys=np.random.normal(0,1.,self.nstar)
        vzs=np.random.normal(0,1.,self.nstar)

        vstar=np.sqrt(vxs**2.+vys**2.+vzs**2.)

        vxkick=self.vesc*(vxs/vstar)
        vykick=self.vesc*(vys/vstar)
        vzkick=self.vesc*(vzs/vstar)

        self.oi,self.of=self._integrate_orbits(vxkick,vykick,vzkick,**kwargs)

        return self.of

    def sample_gaussian(self,tdisrupt=1000.,rate=1.,nstar=None,vmean=100.,vsig=10.,verbose=False, **kwargs):
        """ A function for sampling a uniform core ejection distribution function

        Parameters
        ----------

        tdisrupt : float
            time over which sampling begins (Myr)
        rate : float
            ejection rate (default 1 per Myr)
        nstar : float
            if set, nstar stars will be ejected randomly from tdisrupt to 0 Myr. Rate is recalculated. (default : None)
        vmean : float
            average kick velocity
        vsig : float
            standard deviation of kick velocity distribution
        verbose : bool
            print additional information to screen (default: False)

        Returns
        ----------
        of : orbit
            galpy orbit instance for kicked stars


        History
        -------
        2022 - Written - Grondin/Webb (UofT)

        """

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

        ntstep=kwargs.get('ntstep',10000)
        ts=np.linspace(0.,-1.*self.tdisrupt/self.to,ntstep)
        self.o.integrate(ts,self.mwpot)

        if self.gcpot is None:
            self.pot=self.mwpot
        else:
            moving_pot=MovingObjectPotential(self.o,self.gcpot,ro=self.ro,vo=self.vo)
            self.pot=[self.mwpot,moving_pot]
        
        #Generate kick velocities for escaped stars and binaries
        self.vesc=np.random.normal(vmean,vsig,nstar)

        #Assume a random direction
        vxs=np.random.normal(vmean,vsig,self.nstar)
        vys=np.random.normal(vmean,vsig,self.nstar)
        vzs=np.random.normal(vmean,vsig,self.nstar)

        vstar=np.sqrt(vxs**2.+vys**2.+vzs**2.)

        vxkick=self.vesc*(vxs/vstar)
        vykick=self.vesc*(vys/vstar)
        vzkick=self.vesc*(vzs/vstar)

        self.oi,self.of=self._integrate_orbits(vxkick,vykick,vzkick,**kwargs)
    
        return self.of


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
        self.pt2, = self.ax.plot([],[],'.')

    def _set_data(self,gcdata,sdata,bdata):
        self.gcdata = gcdata
        self.sdata=sdata
        self.bdata=bdata

    def _ani_init(self):
        self.line.set_data([], [])
        self.pt.set_data([],[])
        self.pt2.set_data([],[])

        return self.line,self.pt,self.pt2

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

        if self.binaries:
            if np.sum(escindx)>0:
                self.pt2.set_data(self.bdata[i][0][escindx*self.bindx],self.bdata[i][1][escindx*self.bindx])
            else:
                self.pt2.set_data([],[])    

        self.txt_title.set_text('%s' % str (self.ts[i]*self.to))

        return self.line,self.pt,self.pt2


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
        self.of.integrate(tsint,self.pot)
        if self.binaries:
            self.obf.integrate(tsint,self.pot)

        gcdata=np.zeros(shape=(frames,2))

        for i in range(0,frames):
            gcdata[i]=[self.o.x(self.ts[i]),self.o.y(self.ts[i])]

        sdata=np.zeros(shape=(frames,2,self.nstar))

        for i in range(0,frames):
            sdata[i]=[self.of.x(self.ts[i]),self.of.y(self.ts[i])]

        if self.binaries:
            bdata=np.zeros(shape=(frames,2,self.nstar))
            for i in range(0,frames):
                bdata[i]=[self.obf.x(self.ts[i]),self.obf.y(self.ts[i])]
        else:
            bdata=None

        self._set_data(gcdata,sdata,bdata)

        self.anim = animation.FuncAnimation(self.fig, self._ani_update, init_func=self._ani_init, frames=frames, interval=interval, blit=False)

    def snapout(self,filename='corespray.dat',filenameb='coresprayb.dat'):
        """Output present day positions, velocities, escape times, and escape velocities of stars
        
        Parameters
        ----------

        filename: str
            file name to write data to (default: corespray.dat)
        filenameb: str
            file name to write binary data to (default: corespray.dat)
        History
           -------
        2021 - Written - Webb (UofT)

        """
        R=np.append(self.o.R(0.),self.of.R(0.))
        vR=np.append(self.o.vR(0.),self.of.vR(0.))
        vT=np.append(self.o.vT(0.),self.of.vT(0.))
        z=np.append(self.o.z(0.),self.of.z(0.))
        vz=np.append(self.o.vz(0.),self.of.vz(0.))
        phi=np.append(self.o.phi(0.),self.of.phi(0.))

        vesc=np.append(0.,self.vesc)
        tesc=np.append(0.,self.tesc)

        np.savetxt(filename,np.column_stack([R,vR,vT,z,vz,phi,vesc,tesc]))

        if self.binaries:
            R=self.obf.R(0.)
            vR=self.obf.vR(0.)
            vT=self.obf.vT(0.)
            z=self.obf.z(0.)
            vz=self.obf.vz(0.)
            phi=self.obf.phi(0.)

            vesc=self.vescb
            tesc=self.tesc

            np.savetxt(filenameb,np.column_stack([R,vR,vT,z,vz,phi,vesc,tesc,self.bindx]))
