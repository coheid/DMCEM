import math
from .component import *


## Region
## =================================================
class Region(Component):
	""" Model of an individual region """

	## __init__
	## ---------------------------------------------
	def __init__(self, m, l):
		super(Region, self).__init__(m)
		self.m.setComp(self, l)
		self.l         = l    ## region index
		self.obs       = []   ## list of observables
		self._c        = None ## C^l(t)         == consumption
		self._d        = None ## D^l(t)         == climate damage
		self._eaa      = None ## a^l(t)         == labor productivity
		self._eaa0     = None ## a_0^l          == initial labor productivity [cfg]
		self._egga     = None ## g_a^l(t)       == specific growth rate for labor productivity
		self._egga0    = None ## g_a^l          == specific growth rate for labor productivity [cfg]
		self._eggbara0 = None ## \bar{g}_a^l    == mean growth rate for labor productivity [cfg]
		self._eggbarn0 = None ## \bar{g}_n^l    == mean growth rate for population size [cfg]
		self._eggn     = None ## g_n^l(t)       == specific growth rate for population size
		self._eggn0    = None ## g_n^l          == specific growth rate for population size [cfg]
		self._ekk0     = None ## k_0^l(t)       == capital share
		self._ekk0new  = None ## k_0^l(t+1)     == capital share (temporarily, for next iteration)
		self._enn      = None ## n^l(t)         == population size
		self._ennbar0  = None ## \bar{n}_{0}^l  == population size in 2010 [cfg]
#		self._ennbar9  = None ## \bar{n}_{9}^l  == population size in 2100 [cfg]
#		self._ennbar19 = None ## \bar{n}_{19}^l == population size in 2200 [cfg]
		self._eww      = None ## w^l(t)         == wages
		self._gamma    = None ## gamma^l        == climate damage parameter [cfg]
		self._k        = None ## K^l(t)         == capital allocation in the final sector
		self._knew     = None ## K^l(t+1)       == capital allocation in the final sector (temporarily, for next iteration)
		self._k0       = None ## K^l(t=0)       == initial capital allocation in the final sector
		self._kinit    = None ## K_0^l          == initial share of capital allocation [cfg]
		self._mu       = None ## \mu^l(t)       == consumption share
		self._n        = None ## \bar{N}^l(t)   == total labor supply
		self._n0       = None ## N_0^l(t)       == labor supply of final sector
		self._pi       = None ## \Pi^l          == lifetime profit 
		self._t        = None ## T^l(t)         == transfer to be received by consumers
		self._tprev    = None ## T^l(t)         == transfer to be received by consumers of previous step
		self._theta    = None ## theta^l        == transfer policy
		self._w        = None ## W^l(t)         == lifetime labor income
		self._wprev    = None ## W^l(t)         == lifetime labor income of previous step
		self._y        = None ## Y^l(t)         == total regional production [cfg]
		setPars(self)

	## init
	## ---------------------------------------------
	def init(self):
		""" Extract all variables from cfg """
		super(Region, self).init(["eaa0", "egga0", "eggbara0", "eggbarn0", "eggn0", "ennbar0", "gamma", "kinit", "y"])

	## run
	## ---------------------------------------------
	def run(self):
		""" Computations per step, solving optimization problems """
		if not self.valid: return 
		## produced energy
		self._e = 0
		for i in range(self.m._nSectors):
			es       = self.m.getComp("EnergySector",         i)
			ea       = self.m.getComp("EnergyAgent" , self.l, i)
			ea.run() ## make sure that ea._e = E_{i,t}^l are up-to-date
			self._e += es._kappa*(ea._e**self.m.rho)  ## E_t^l without exponent, Eqn (2)
		self._e = math.pow(self._e, 1./self.m.rho) ## E_t^l, Eqn (2)
		## climate damage
		cm        = self.m.getComp("ClimateModel")
		self._d   = 1.0 - math.exp(-self._gamma*(cm._s-cm._sbar))  ## D_t^l, Eqn (15)
		## output
		self._y   = (1-self._d) * self._k**self.m._alpha0 * self._n0**(1-self.m._alpha0-self.m._nu0) * self._e**self.m._nu0  ## Y_t^l, Eqn (1)

	## push
	## ---------------------------------------------
	def push(self):
		""" Pushes parameters of labor supply and other variables to next iteration """
		if not self.valid: return 
		## labor supply
		self._eggn  = self._eggbarn0+self._eggn0*math.exp(-self.m._deltan*(self.m.t+1)) ## g_{n,t}^l, Eqn (41)
		self._egga  = self._eggbara0+self._egga0*math.exp(-self.m._deltaa*(self.m.t+1)) ## g_{a,t}^l, Eqn (41)
		self._enn   = (1+self._eggn)*self._enn  ## n_t^l, Eqn (40)
		self._eaa   = (1+self._egga)*self._eaa  ## a_t^l, Eqn (40)
		## transfer
		self._wprev = self._w
		self._tprev = self._t

	## start
	## ---------------------------------------------
	def start(self, par=None, var="central"):
		""" (Re-)Initialize parameters of this object and reset cache """
		if not self.valid: return 
		super(Region, self).start(par, var)
		self._enn   = self._ennbar0  ## starting value for n_t^l
		self._eaa   = self._eaa0     ## starting value for a_t^l
		self._mu    = 0.0            ## starting value for \mu^l
		self._w     = 0.0            ## starting value for W_t^l
		self._wprev = 0.0            ## starting value for W_t^l
		self._t     = 0.0            ## starting value for T_t^l
		self._tprev = 0.0            ## starting value for T_t^l
		self._theta = 0.0            ## starting value for \theta^l


