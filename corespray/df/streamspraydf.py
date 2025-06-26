from galpy.df import chen24spraydf, fardal15spraydf
from galpy.util import conversion
from galpy.orbit import Orbit
import numpy as np

def streamspraydf(mgc,gcorbit,gcpot,galpot,tdisrupt,nstar,method='Fardal15',integrate=True,**kwargs):


	ro=gcorbit._ro
	vo=gcorbit._vo
	zo=gcorbit._zo
	solarmotion=gcorbit._solarmotion



	if method=='Fardal15':
		spdf = fardal15spraydf(mgc,progenitor=gcorbit,pot=galpot,tdisrupt=tdisrupt,**kwargs)
		spdft = fardal15spraydf(mgc,progenitor=gcorbit,pot=galpot,leading=False,tdisrupt=tdisrupt,**kwargs)

	elif method=='Chen24':
		spdf = chen24spraydf(mgc,progenitor=gcorbit,pot=galpot,tdisrupt=tdisrupt, progpot=gcpot,**kwargs)
		spdft = chen24spraydf(mgc,progenitor=gcorbit,pot=galpot,leading=False,tdisrupt=tdisrupt, progpot=gcpot,**kwargs)

	oil,dtl = spdf.sample(n=nstar,returndt=True,integrate=False,**kwargs)
	oit,dtt = spdft.sample(n=nstar,returndt=True,integrate=False,**kwargs)

	vkl=np.zeros(len(oil))
	vkt=np.zeros(len(oil))


	for i in range(0,len(oil)):
		vkl[i]=np.sqrt((oil.vx()[i]-gcorbit.vx(-dtl[i]))**2.+(oil.vy()[i]-gcorbit.vy(-dtl[i]))**2.+(oil.vz()[i]-gcorbit.vz(-dtl[i]))**2.)
		vkt[i]=np.sqrt((oit.vx()[i]-gcorbit.vx(-dtt[i]))**2.+(oit.vy()[i]-gcorbit.vy(-dtt[i]))**2.+(oit.vz()[i]-gcorbit.vz(-dtt[i]))**2.)


	if integrate:
		out = np.empty((6, len(oil)))
		outt = np.empty((6, len(oit)))

		# Now integrate the orbits
		for i in range(0,len(oil)):
			o = Orbit([oil.R()[i]/ro, oil.vR()[i]/vo, oil.vT()[i]/vo, oil.z()[i]/ro, oil.vz()[i]/vo, oil.phi()[i]],ro=ro,vo=vo,zo=zo,solarmotion=solarmotion)
			o.integrate(np.linspace(-dtl[i], 0.0, 10001), galpot)
			o = o(0.0)
			out[:, i] = [o.R()/ro, o.vR()/vo, o.vT()/vo, o.z()/ro, o.vz()/vo, o.phi()]

			o = Orbit([oit.R()[i]/ro, oit.vR()[i]/vo, oit.vT()[i]/vo, oit.z()[i]/ro, oit.vz()[i]/vo, oit.phi()[i]],ro=ro,vo=vo,zo=zo,solarmotion=solarmotion)
			o.integrate(np.linspace(-dtt[i], 0.0, 10001), galpot)
			o = o(0.0)
			outt[:, i] = [o.R()/ro, o.vR()/vo, o.vT()/vo, o.z()/ro, o.vz()/vo, o.phi()]


		ofl = Orbit(
		    vxvv=out.T,
		    ro=ro,
		    vo=vo,
		    zo=zo,
		    solarmotion=solarmotion,
		)

		oft = Orbit(
		    vxvv=outt.T,
		    ro=ro,
		    vo=vo,
		    zo=zo,
		    solarmotion=solarmotion,
		)

	return oil,ofl,dtl,vkl,oit,oft,dtt,vkt






