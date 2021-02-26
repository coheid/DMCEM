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
		self.i      = i    ## energy sector index
		self._alpha = None ## \alpha_i(t) == capital elasticity [cfg]
		self._ecc   = None ## c_i         == constant per-unit extraction costs per resource [cfg]
		self._kappa = None ## \kappa_i    == relative productivity of this sector [cfg]
		self._nu    = None ## \nu_i(t)    == energy elasticity [cfg]
		self._r     = None ## R_i         == initial resource stock [cfg]
		self._rcrit = None ## R_i^{crit}  == critical resource stock [cfg] 
		self._evv   = None ## v_i(t)      == resource price
		self._zeta  = None ## \zeta_i     == carbon content per unit resource [cfg]

	## init
	## ---------------------------------------------
	def init(self):
		""" Extract all variables from cfg """
		super(EnergySector, self).init(["alpha","ecc","kappa","nu","r","rcrit","zeta"])

	## run
	## ---------------------------------------------
	def run(self):
		""" Computations per step, solving optimization problem """
		if not self.valid: return 
		if self.t == self.m.t: return ## sync to market time!
		self._evv = self._ecc + self.m._err*(self._evv - self._ecc) ## v_{i,t}, Eqn (10) 
		self.t  = self.m.t


