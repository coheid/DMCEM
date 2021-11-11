from .component import *


## EnergySector
## =================================================
class EnergySector(Component):
	""" Model of an individual energy production sector """

	## __init__
	## ---------------------------------------------
	def __init__(self, m, i):
		super(EnergySector, self).__init__(m)
		self.m.setComp(self, i)
		self.i        = i    ## energy sector index
		self.obs      = []   ## list of observables
		self._alpha   = None ## \alpha_i(t) == capital elasticity [cfg]
		self._ecc     = None ## c_i         == constant per-unit extraction costs per resource [cfg]
		self._kappa   = None ## \kappa_i    == relative productivity of this sector [cfg]
		self._nu      = None ## \nu_i(t)    == energy elasticity [cfg]
		self._r       = None ## R_i         == initial resource stock [cfg]
		self._rcrit   = None ## R_i^{crit}  == critical resource stock [cfg] 
		self._evv     = None ## v_i(t)      == resource price [cfg if dothrow=False]
		self._evvprev = None ## v_i(t)      == resource price of previous step
		self._zeta    = None ## \zeta_i     == carbon content per unit resource [cfg]
		setPars(self)

	## init
	## ---------------------------------------------
	def init(self):
		""" Extract all variables from cfg """
		super(EnergySector, self).init(["alpha","ecc","evv","kappa","nu","r","rcrit","zeta"])

	## push
	## ---------------------------------------------
	def push(self):
		""" Computes the resource prices of the timestep """
		self._evvprev = self._evv

	## runEvv
	## ---------------------------------------------
	def runEvv(self):
		""" Computes the resource prices of the timestep """
		if not self.valid: return 
		self._evv = self._ecc + self.m._err*(self._evvprev - self._ecc)  ## v_{i,t}, Eqn (10) 

	## start
	## ---------------------------------------------
	def start(self, par=None, var="central"):
		""" (Re-)Initialize parameters of this object and reset cache """
		if not self.valid: return 
		super(EnergySector, self).start(par, var)
		self._evvprev = self._evv


