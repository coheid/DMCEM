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
		self._s      = None ## S(t)      == current total amount of CO2
		self._s1     = None ## S_1(t)    == current amount of permanent CO2 [cfg]
		self._s2     = None ## S_2(t)    == current amount of non-permanent CO2 [cfg]
		self._phil   = None ## \varphi_L == share of emissions becoming permanent [cfg]
		self._phi0   = None ## \varphi_0 == share of emissions becoming non-permanent [cfg]
		self._phi    = None ## \varphi   == decay rate of non-permanent CO2 [cfg]
		self._sbar   = None ## \bar{S}   == total pre-industrial concentration [cfg]

	## init
	## ---------------------------------------------
	def init(self):
		""" Extract all variables from cfg """
		super(ClimateModel, self).init(["s1", "s2", "phil", "phi0", "phi", "sbar"])

	## run
	## ---------------------------------------------
	def run(self):
		""" Computations per step, solving optimization problems """
		if not self.valid: return 
		if self.t == self.m.t: return ## sync to market time!
		self._s1 =               self._s1 + self._phil               *self.m._z ## S_{1,t}, Eqn (14a)
		self._s2 = (1-self._phi)*self._s2 + (1-self._phil)*self._phi0*self.m._z ## S_{2,t}, Eqn (14b)
		self._s  = self._s1 + self._s2                                          ## S_t    , page 7
		self.t   = self.m.t

