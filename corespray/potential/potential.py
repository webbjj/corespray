import numpy as np

from galpy import potential
try:
    from galpy.util import conversion
except:
    import galpy.util.bovy_conversion as conversion

from galpy.potential import MWPotential2014

import copy

def load_galpy_potential(mwpot=MWPotential2014, bar=False,dw_arm=False,trans_arm=False,ro=8.,vo=220., **kwargs):
    """ A function for loading a galpy potentail to be used as the Galactic potential
    -- Note that the potential is designed to being at a negative time and end at time=0

    Parameters
    ----------

    mwpot : galpy potential
        Underlying Galactic potential (default : MWPotential2014)
    bar : bool
        Option to include bar potential (default : False)
    dw_arm : bool
        Option to include density wave spiral arm potential (default : False)
    tran_arm : bool
        Option to include transient spiral arm potential (default : False)
    ro : float
        galpy length scaling parameter (default: 8.)
    vo : float
        galpy veocicity scaling parameter (default: 220.)

    Key Word Arguments
    ----------

    if bar:
        tform_bar : float
            start of bar growth / bar period (default : -1.0e-9)
        tsteady : float
            time from tform at which bar is full grown / bar period (default : 5)
        omegab : float
            pattern speed of the bar (default : 35.75 km s−1 kpc−1)
        barphi : float
            angle between line connecting the sun with the Galactic centre and the bar's major axix in radians (default : 25 degrees )
        rb : float
            length of the bar in kpc (default : 5 kpc )
        Af : float
            strength of the bar in kpc (default : 645.3 km^2/s^2 )

    if dw_arm:
        N : int
            number of spiral arms (default : 4)
        amp : float
            amplitude to be applied to the potential (default : 1)
        phi_ref : float
            reference angle (default : 45 degrees)
        alpha : float
            pitch angle of the logarithmic spiral arms in radians (default : 12. degrees)
        omega : float
            rotational pattern speed of the spiral arms (default : 21.725 km s−1 kpc−1)
        tform_arm : float
            start of arm growth in Myr (default : tform_bar if bar else 0 Myr)

    if trans_arm:
        N : int
            number of spiral arms (default : 2)
        amp : float
            amplitude to be applied to the potential (default : 0.75)
        phi_ref : float
            reference angle in radians (default : 25 degrees)
        alpha : float
            pitch angle of the logarithmic spiral arms in radians (default : 25. degrees)
        tstart : float
            start time of transient spiral arms in Myr (default : None)
        vpo : float
            amplitude of circular-velocity curve in km/s (default : 220.0 km/s)
        beta : float
            power-kaw amplitude of the circular-velocity curve (default : -0.1)
        to : float
            reference time at which the potential reaches maximum in Myr (default : -500 Myr)
        pa : float (HERE)
            position angle (default : if bar omega*to, else 1.3*to)
        sigma : float
            standard deviation of Gaussian Wrapper potential (default : 1.3)
        transwrap : bool
            Use TransientWrapperPotential for transient spiral arm potential (default : True)
            Note - this is the preferred method if there are arms that don't start until later in the simulation
            or have fully decayed before the simulation finishes

    Notes
    -------
        if trans_arm and tstart is not given, will add a series of three transient winding spirals where
        to is peak (--> Initialize 3 per to). The lifetime is a gaussian that grows and decays.

        if trans_arm and tstart is given, will add a series of transient winding spirals dating back to tstart where
        to is peak (--> Initialize 3 per to). The lifetime is a gaussian that grows and decays.


    History
    -------
    2022 - Written - Grandin (UofT)

    """

    pot=copy.deepcopy(mwpot)

    tcon=conversion.time_in_Gyr(ro=ro,vo=vo)

    if bar:
        # Add bar potential (https://docs.galpy.org/en/stable/reference/potentialdehnenbar.html)
        tform_bar=kwargs.get('tform_bar',-1.0e9) 
        tsteady=kwargs.get('tsteady', 5.)
        omegab= kwargs.get('omegab',35.75) * (ro/vo)
        barphi=kwargs.get('barphi',np.deg2rad(25.)) 
        rb=kwargs.get('rb',5.)/ro
        Af=kwargs.get('Af',(220.*220.0/75.))/(vo**2.)

        dp= potential.DehnenBarPotential(omegab=omegab,rb=rb,Af=Af,tform=tform_bar,tsteady=tsteady,barphi=barphi)

        pot.append(dp)

    if dw_arm:
        # Add density wave spiral arm potential (http://galpy.readthedocs.io/en/latest/reference/potentialspiralarms.html)
        # Wrapper (https://docs.galpy.org/en/stable/reference/potentialdehnensmoothwrapper.html)

        N=kwargs.get('N',4)
        amp=kwargs.get('amp',1)
        phi_ref=kwargs.get('phi_ref',np.deg2rad(45.))
        alpha=kwargs.get('alpha',np.deg2rad(12.))
        omega=kwargs.get('omega',21.725)*(ro/vo)

        if bar:
            tform_arm=kwargs.get('tform_arm',None)

            if tform_arm is None:
                tform_arm=dp.tform()
            else:
                tform_arm=(tform_arm/1000.0)/tcon
        else:
            tform_arm=kwargs.get('tform_arm',0.0)
            tform_arm=(tform_arm/1000.0)/tcon

        density_wave_spiral_arm=potential.DehnenSmoothWrapperPotential(pot=potential.SpiralArmsPotential(N=N,amp=amp,phi_ref=phi_ref,alpha=alpha,omega=omega),tform=tform_arm)
        pot.append(density_wave_spiral_arm)

    elif trans_arm:
        # Add transient spiral arm potential (http://galpy.readthedocs.io/en/latest/reference/potentialspiralarms.html)
        # Wrapper (https://docs.galpy.org/en/stable/reference/potentialcorotwrapper.html)
        # Wrapper (https://docs.galpy.org/en/stable/reference/potentialgaussampwrapper.html)

        N=kwargs.get('N',2)
        amp=kwargs.get('amp',0.75)
        phi_ref=kwargs.get('phi_ref',np.deg2rad(25.))
        alpha=kwargs.get('alpha',np.deg2rad(25.))

        sp= potential.SpiralArmsPotential(N=N,amp=amp,phi_ref=phi_ref,alpha=alpha)


        tstart=kwargs.get('tstart',None)

        vpo=kwargs.get('vpo',220.0)/vo
        beta=kwargs.get('beta',-0.1)
        to=kwargs.get('to',-500.0)/(1000.0*tcon)

        if bar:
            pa=(kwargs.get('omegab',35.75) * (ro/vo))*to
        else:
            pa=kwargs.get('pa',35.75*(ro/vo)*to)

        sigma = kwargs.get('sigma',1.3)

        if tstart==None:
            
            csp= potential.GaussianAmplitudeWrapperPotential(\
                   pot=potential.CorotatingRotationWrapperPotential(pot=sp,vpo=vpo,to=to,beta=beta,pa=pa),
                   to=to,sigma=sigma)
            pot.append(csp)
            csp3= potential.GaussianAmplitudeWrapperPotential(\
                   pot=potential.CorotatingRotationWrapperPotential(pot=sp,vpo=vpo,to=to/2.,beta=beta,pa=pa),
                   to=to/2.,sigma=sigma)
            pot.append(csp3)
            csp5= potential.GaussianAmplitudeWrapperPotential(\
                   pot=potential.CorotatingRotationWrapperPotential(pot=sp,vpo=vpo,to=0.,beta=beta,pa=pa),
                   to=0.,sigma=sigma)
            pot.append(csp5)
        else:
            #Assume tstart is in Myr and setup so there are three transient winding spirals per 500 Myr
            tstart=(tstart/1000.)/tcon
            t=0.
            tpot=[]
            while t > tstart:
                csp= potential.GaussianAmplitudeWrapperPotential(\
                       pot=potential.CorotatingRotationWrapperPotential(pot=sp,vpo=1.,to=t,beta=beta,pa=pa),
                       to=t,sigma=sigma)
                tpot.append(csp)
                t+=to/2.
            
            if kwargs.get('transwrap',True):
                pot.append(TransientWrapperPotential(amp=1.,pot=tpot,_init=True))
            else:
                pot=pot+tpot


    return pot

class TransientWrapperPotential(potential.Potential):
    def __init__(self,amp=1.,pot=None,ro=None,vo=None,_init=None,**kwargs):
        """
        NAME:

           __init__

        PURPOSE:

           initialize a TransientWrapperPotential

        INPUT:

           amp - amplitude to be applied to the potential (default: 1.)

           pot - Potential instance or list thereof; the amplitude of this will be grown by this wrapper

           ro - galpy scaling parameter (default: None)

           vo - galpy scaling parameter (default: None)


        OUTPUT:

           (none)

        HISTORY:

           2022 - Webb (UofT)

        """
        if not _init: return None # Don't run __init__ at the end of setup
        potential.Potential.__init__(self,amp=amp,ro=ro,vo=vo)
        self._pot= pot
        self.isNonAxi= potential._isNonAxi(self._pot)
        self.ro=ro
        self.vo=vo
        
        self.times=np.array([])
        self.tend=np.array([])
        for i in range(0,len(self._pot)):
            self.times=np.append(self.times,self._pot[i]._to-10.*np.sqrt(self._pot[i]._sigma2))
            self.tend=np.append(self.tend,self._pot[i]._to+10.*np.sqrt(self._pot[i]._sigma2))

        self.isNonAxi = True

    def __getattr__(self,attribute):
        if attribute == '_evaluate' \
                or attribute == '_Rforce' or attribute == '_zforce' \
                or attribute == '_phiforce' \
                or attribute == '_R2deriv' or attribute == '_z2deriv' \
                or attribute == '_Rzderiv' or attribute == '_phi2deriv' \
                or attribute == '_Rphideriv' or attribute == '_dens' \
                or attribute == 'isNonAxi':
            return lambda R,Z,phi=0.,t=0.: \
                self._wrap(attribute,R,Z,phi=phi,t=t)
        else:
            return super(TransientWrapperPotential,self).__getattr__(attribute)

    def _wrap_pot_func(self,attribute):
        if attribute == '_evaluate':
            return evaluatePotentials
        elif attribute == '_dens':
            return evaluateDensities
        elif attribute == '_Rforce':
            return evaluateRforces
        elif attribute == '_zforce':
            return evaluatezforces
        elif attribute == '_phiforce':
            return evaluatephiforces
        elif attribute == '_R2deriv':
            return evaluateR2derivs
        elif attribute == '_z2deriv':
            return evaluatez2derivs
        elif attribute == '_Rzderiv':
            return evaluateRzderivs
        elif attribute == '_phi2deriv':
            return lambda p,R,Z,phi=0.,t=0.: \
                evaluatePotentials(p,R,Z,phi=phi,t=t,dphi=2)
        elif attribute == '_Rphideriv':
            return lambda p,R,Z,phi=0.,t=0.: \
                evaluatePotentials(p,R,Z,phi=phi,t=t,dR=1,dphi=1)
        else: #pragma: no cover
            raise AttributeError("Attribute %s not found in for this TransientWrapperPotential" % attribute)
    
    def _evaluate(self,R,z,phi=0.,t=0.):
        mindx=self.match_time(t)
        if len(mindx)>0:
            return potential.evaluatePotentials(self._pot[mindx[0]:mindx[-1]+1],R,z,phi=phi,t=t)
        else:
            return 0.

    def _Rforce(self,R,z,phi=0.,t=0.):
        mindx=self.match_time(t)
        if len(mindx)>0:
            return potential.evaluateRforces(self._pot[mindx[0]:mindx[-1]+1],R,z,phi=phi,t=t)
        else:
            return 0.

    def _zforce(self,R,z,phi=0.,t=0.):
        mindx=self.match_time(t)
        if len(mindx)>0:
            return potential.evaluatezforces(self._pot[mindx[0]:mindx[-1]+1],R,z,phi=phi,t=t)
        else:
            return 0.

    def _phiforce(self,R,z,phi=0.,t=0.):
        mindx=self.match_time(t)
        if len(mindx)>0:
            return potential.evaluatephiforces(self._pot[mindx[0]:mindx[-1]+1],R,z,phi=phi,t=t)
        else:
            return 0.

    def match_time(self,t):
        
        mindx=(self.times<=t) * (self.tend>=t)
        indx=np.linspace(0,len(self.times)-1,len(self.times),dtype=int)

        #print((self.times<=t))
        #print((self.tend>=t))
        #print('MINDX: ',np.ndarray.tolist(indx[mindx]))
        return np.ndarray.tolist(indx[mindx])
