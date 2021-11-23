import math
from .component import *



## ClimateModel
## =================================================
class ClimateModel(Component):
	""" Model of the climate, i.e. CO2 concentrations """

	## init
	## ---------------------------------------------
	def __init__(self, m):
		""" Constructor """
		super(ClimateModel, self).__init__(m)
		self.m.setComp(self)
		self.t       = 0    ## time value for sync
		self.obs     = []   ## list of observables
		self._phil   = None ## \varphi_L == share of emissions becoming permanent [cfg]
		self._phi0   = None ## \varphi_0 == share of emissions becoming non-permanent [cfg]
		self._phi    = None ## \varphi   == decay rate of non-permanent CO2 [cfg]
		self._s      = None ## S(t)      == current total amount of CO2
		self._s1     = None ## S_1(t)    == current amount of permanent CO2 [cfg]
		self._s1prev = None ## S_1(t)    == current amount of permanent CO2 previous step
		self._s2     = None ## S_2(t)    == current amount of non-permanent CO2 [cfg]
		self._s2prev = None ## S_2(t)    == current amount of non-permanent CO2 previous step
		self._sbar   = None ## \bar{S}   == total pre-industrial concentration [cfg]
		self._t      = None ## TEMP(t)   == current global temperature [cfg]
		setPars(self)

	## init
	## ---------------------------------------------
	def init(self):
		""" Extract all variables from cfg """
		super(ClimateModel, self).init(["phil", "phi0", "phi", "s1", "s2", "sbar", "t"])

	## push
	## ---------------------------------------------
	def push(self):
		""" Push S1 and S2 to next step """
		self._s1prev = self._s1
		self._s2prev = self._s2

	## run
	## ---------------------------------------------
	def run(self):
		""" Computations per step, solving optimization problems """
		if not self.valid: return 
		self._s1 =               self._s1prev + self._phil               *self.m._z ## S_{1,t}, Eqn (14a)
		self._s2 = (1-self._phi)*self._s2prev + (1-self._phil)*self._phi0*self.m._z ## S_{2,t}, Eqn (14b)
		self._s  = self._s1 + self._s2                                              ## S_t    , page 7
		#self._s  = self._s * 1.09 ## FIXME: for emissions test
		self._t  = 3*math.log(self._s/self._sbar)/math.log(2)                       ## TEMP_t , Eqn (31)

	## start
	## ---------------------------------------------
	def start(self, par=None, var="central"):
		""" (Re-)Initialize parameters of this object and reset cache """
		if not self.valid: return 
		super(ClimateModel, self).start(par, var)
		self._s1prev = self._s1 ## starting value for S1
		self._s2prev = self._s2 ## starting value for S2



