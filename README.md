# corespray
A python package for sampling a distribution function for stars that have been ejected from a star cluster's core. If you use corespray in your research, please cite [Grondin et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.tmp.3150G/abstract) and link to https://github.com/webbjj/corespray/. Please also see the requirements below, as corespray relies heavily on other packages (e.g. [galpy](https://docs.galpy.org/en/)) that should be cited as well.

# Installation 
To install corespray from GitHub, clone the repository and install via setup tools:

`git clone https://github.com/webbjj/corespray.git`  
`cd corespray`  
`python setup.py install`  

Please note that if you don’t have permission to write files to the default install location (often /usr/local/lib), you will either need to run:

`sudo python setup.py install`

or

`python setup.py install --prefix='PATH'`

where ‘PATH’ is a directory that you do have permission to write in and is in your PYTHONPATH.

# Requirements  

corespray requires the following python packages: [galpy](https://docs.galpy.org/en/v1.8.1/), [matplotlib](https://matplotlib.org/), [numpy](https://numpy.org/)

# Sample

To initialize a corespray distribution function:

`from corespray import corespraydf`
`cspray = corespraydf(gcorbit, pot, mgc, rgc, W0, ro, vo)`

where the parameters are defined as follows

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

There are then options to sample a kick velocity distribution due to three-body encounters, a uniform kick velocity distribution, or a gaussian kick velocity distribution. 

To sample a kick velocity distribution due to three-body encounters, run

`os,ob=cspray.sample_three_body(tdisrupt=1000.,rate=1.,nstar=None,mu0=0.,sig0=10.0,vesc0=10.0,rho0=1.,mmin=0.1,mmax=1.4,alpha=-1.35,masses=None,m1=None,m2a=None,m2b=None,emin=None,emax=None,balpha=-1,q=-3, npeak=5.,binaries=False,verbose=False, show_progress=False, **kwargs)`

where the parameters are defined as follows

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
        show_progress : bool
            show a progress bar during the orbital integration

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

Alternatively to sample a gaussian kick velocity distribution

`os = cspray.sample_gaussian(tdisrupt=1000.,rate=1.,nstar=None,vmin=0.,vmax=500.)`

where the parameters are as follows:

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
		Returns
		----------
		of : orbit
			galpy orbit instance for kicked stars

Finally, it is possible to sample a uniform distribution of kick velocities

`os = cspray.sample_uniform(tdisrupt=1000.,rate=1.,nstar=None,vmin=0.,vmax=500.)`

which has parameters

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
		Returns
		----------
		of : orbit
			galpy orbit instance for kicked stars

