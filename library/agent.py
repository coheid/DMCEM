from .component import *


## EnergyAgent
## =================================================
class EnergyAgent(Component):
	""" Model of an individual energy production sector in a region """

	## __init__
	## ---------------------------------------------
	def __init__(self, m, l, i):
		super(EnergyAgent, self).__init__(m)
		self.m.setComp(self, l, i)
		self.l      = l    ## region index
		self.i      = i    ## energy sector index
		self.obs    = []   ## list of observables
		self._e     = None ## E_i^l(t)      == energy production
		self._ekk   = None ## k_i^l(t)      == capital share
		self._enn   = None ## n_i^l(t)      == labor share
		self._epp   = None ## p_i^l(t)      == energy price
		self._eta   = None ## \eta_i^l(t)   == energy production share [cfg]
		self._k     = None ## K_i^l(t)      == capital allocation
		self._n     = None ## N_i^l(t)      == labor allocation
		self._pi    = None ## \Pi_i^l       == lifetime profit
		self._q     = None ## Q_i^l         == productivity parameter [cfg]
		self._r     = None ## R_i^l         == resource stock [cfg]
		self._x     = None ## X_i^l(t)      == resource input
		self._xbar  = None ## \bar{X}_i^l   == initial resource consumption [cfg]
		self._z     = None ## Z_i^l(t)      == carbon emissions
		setPars(self)

	## init
	## ---------------------------------------------
	def init(self):
		""" Extract all variables from cfg """
		super(EnergyAgent, self).init(["eta", "q", "r", "xbar"])

	## run
	## ---------------------------------------------
	def run(self):
		""" Computations per step, solving optimization problem """
		if not self.valid: return 
		r         = self.m.getComp("Region"      , self.l)
		es        = self.m.getComp("EnergySector", self.i)
		self._e   = self._q * self._k**es._alpha * self._n**(1-es._alpha-es._nu)  ## E_{i,t}^l, Eqs (4,7)
		self._e  *= self._x**es._nu if self._x>0 else 1 ## E_{i,t}^l, Eqs (4,7)

	## start
	## ---------------------------------------------
	def start(self, par=None, var="central"):
		""" (Re-)Initialize parameters of this object and reset cache """
		if not self.valid: return 
		super(EnergyAgent, self).start(par, var)
		self._x = self._xbar  ## starting value for X_{t,l}^i


