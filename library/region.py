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
		self._enn      = None ## n^l(t)         == population size
		self._ennbar0  = None ## \bar{n}_{0}^l  == population size in 2010 [cfg]
#		self._ennbar9  = None ## \bar{n}_{9}^l  == population size in 2100 [cfg]
#		self._ennbar19 = None ## \bar{n}_{19}^l == population size in 2200 [cfg]
		self._eww      = None ## w^l(t)         == wages
		self._gamma    = None ## gamma^l        == climate damage parameter [cfg]
		self._k        = None ## K^l(t)         == capital allocation in the final sector
		self._k0       = None ## K^l(t=0)       == initial capital allocation in the final sector
		self._kinit    = None ## K_0^l          == initial share of capital allocation [cfg]
		self._mu       = None ## \mu^l(t)       == consumption share
		self._n        = None ## \bar{N}^l(t)   == total labor supply
		self._n0       = None ## N_0^l(t)       == labor supply of final sector
		self._pi       = None ## \Pi^l          == lifetime profit 
		self._t        = None ## T^l(t)         == transfer to be received by consumers
		self._theta    = None ## theta^l        == transfer policy
#		self._tau       = None ## tau^l(t)       == tax policy
		self._w        = None ## W^l(t)         == lifetime labor income
		self._y        = None ## Y^l(t)         == total regional production [cfg]

	## init
	## ---------------------------------------------
	def init(self):
		""" Extract all variables from cfg """
		super(Region, self).init(["eaa0", "egga0", "eggbara0", "eggbarn0", "eggn0", "ennbar0", "gamma", "kinit", "y"])

	## initEtas
	## ---------------------------------------------
	def initEtas(self):
		""" Computes the initial energy mixes """
		if not self.valid: return 
		if self._c == 0: return
		## compute total initial consumption
		self._xbar = 0 ## sum of \bar{X}_i^l over i
		for i in range(self.m._nSectors):
			ea          = self.m.getComp("EnergyAgent", self.l, i)
			self._xbar += ea._xbar
		## compute difference between C^l and \bar{X}^l and distribute it to the renewables
		diff = self._c-self._xbar 
		for i in range(self.m._nSectors):
			ea = self.m.getComp("EnergyAgent", self.l, i)
			if ea._xbar==0: ea._xbar = diff; break ## FIXME: works with only one renewable sector!
		## compute energy mixes
		for i in range(self.m._nSectors):
			ea      = self.m.getComp("EnergyAgent", self.l, i)
			ea._eta = ea._xbar/self._xbar

	## run
	## ---------------------------------------------
	def run(self):
		""" Computations per step, solving optimization problems """
		if not self.valid: return 
		if self.t == self.m.t: return ## sync to market time!
		## produced energy
		self._e = 0
		for i in range(self.m._nSectors):
			es       = self.m.getComp("EnergySector",         i)
			ea       = self.m.getComp("EnergyAgent" , self.l, i)
			ea.run() ## make sure that ea._e = E_{i,t}^l are up-to-date
			self._e += es._kappa*(ea._e)**self.m.rho   ## E_t^l without exponent, Eqn (2)
		self._e = math.pow(self._e, 1./self.m.rho)     ## E_t^l, Eqn (2)
		## climate damage and total production
		cm        = self.m.getComp("ClimateModel")
		self._d   = 1 - math.exp(-self._gamma*(cm._s-cm._sbar))  ## D_t^l, Eqn (15)
		self._y   = (1-self._d) * self._k**self.m._alpha0 * self._n**(1-self.m._alpha0-self.m._nu0) * self._e**self.m._nu0  ## Y_t^l, Eqn (1)

	## runEnn
	## ---------------------------------------------
	def runEnn(self):
		""" Computations of labor supply """
		if not self.valid: return 
##		#print(self.t,self.m.t)
##		#if self.t == self.m.t: return ## sync to market time!
		self._eggn     = self._eggbarn0+self._eggn0*math.exp(-self.m._deltan*self.m.t*self.m.tstep) ## g_{n,t}^l, Eqn (41)
		self._egga     = self._eggbara0+self._egga0*math.exp(-self.m._deltaa*self.m.t*self.m.tstep) ## g_{a,t}^l, Eqn (41)
		self._enn      = (1+self._eggn)*self._enn       ## n_t^l, Eqn (40)
		self._eaa      = (1+self._egga)*self._eaa       ## a_t^l, Eqn (40)
		self._n        = self._enn*self._eaa            ## N_t^{l,s} = \bar{N}_t^l, page 19

	## start
	## ---------------------------------------------
	def start(self, par=None, var="central"):
		""" (Re-)Initialize parameters of this object and reset cache """
		if not self.valid: return 
		super(Region, self).start(par, var)
		self._enn = self._ennbar0  ## starting value for n_t^l
		self._eaa = self._eaa0     ## starting value for a_t^l
		self._w   = 0              ## starting value for W_t^l
		self._t   = 0              ## starting value for T_t^l
		self._k   = self.m._kbar*self._k0       ## starting value for K_t^l
		self._theta = 1 # FIXME


