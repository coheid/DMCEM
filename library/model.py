import numpy  as np 
import pandas as pd 
import os

from core.config    import *
from core.functions import *
from core.master    import *
from core.tui       import *

from .component     import endPars,extractPars,startPars,storePars
from .agent         import *
from .climate       import *
from .region        import *
from .sector        import *


## throw
## -------------------------------------------------
def throw(mean, variance, minimum=None, maximum=None, direction=None, incl=True):
	""" Helper function to throw initial value for a 
	variable based on certain constraints """
	while True:
		x = np.random.normal(mean, variance)
		if minimum is not None:
			if     incl and x<=minimum: continue
			if not incl and x< minimum: continue
		if maximum is not None:
			if     incl and x>=maximum: continue
			if not incl and x> maximum: continue 
		if direction=="up" and x<=mean: continue
		if direction=="dn" and x>=mean: continue
		return x
	return 0


## Model
## =================================================
class Model(Master):
	""" Core class for the simulation and market """

	## main methods
	## =============================================

	## __init__
	## ---------------------------------------------
	def __init__(self, name, path, base):
		""" Constructor """
		super(Model, self).__init__(name, [path])
		self.base      = base
		self.out       = "%s/output/%s/"%(self.base, os.path.splitext(os.path.basename(path))[0])
		mkdir(self.out)
		self.tui       = Tui(self)
		self.ui        = self.tui
		self.idx1      = None ## mimick component
		self.idx2      = None
		self.it        = 0  ## counter for successful iterations
		self.t         = 0  ## time value
		self.objs      = {} ## pointers to the components
		self.pars      = {} ## parameter values and systematic variations
		self.obs       = [] ## list of observables
		self.cache     = {} ## cache for observable values per timestep
		self._alpha0   = None ## \alpha_0      == capital elasticity [cfg]
		self._beta     = None ## \beta         == utility parameter [cfg]
		self._ctry0    = None ##               == mean     for throwing C [cfg]
		self._ctry1    = None ##               == variance for throwing C [cfg]
		self._ccrit    = None ## C^{crit}      == critical value for C [cfg]
		self._deltaa   = None ## \delta_a      == convergence for labor productivity [cfg]
		self._deltan   = None ## \delta_n      == convergence for population size [cfg]
		self._ekk      = None ## k(t)          == sum of k_{i,t}^l over i,l
		self._enn      = None ## n_0^l(t)      == labor cost share of final sector
		self._eqq      = None ## q(t)          == discount factor [cfg]
		self._err      = None ## r(t)          == capital price
		self._err0     = None ## r(0)          == capital price
		self._kbar     = None ## \bar{K}(t)    == sum of K_0^l(t) over l [cfg]
		self._nbar     = None ## \bar{N}(t)    == sum of N^l(t) over l 
		self._nRegions = None ## L             == number of regions [cfg]
		self._nSectors = None ## I             == number of sectors [cfg]
		self._nu0      = None ## \nu_0         == energy elasticity [cfg]
		self._sigma    = None ## \sigma        == utility parameter [cfg]
		self._rho      = None ## \rho          == substitution elasticity [cfg]
		self._tau      = None ## \tau^{eff}(t) == carbon tax revenue
		self._taubar   = None ## \bar{\tau}(t) == carbon tax 
		self._vtry1    = None ##               == variance for throwing v [cfg]
		self._y        = None ## Y(t)          == sum of Y^l(t) over l 
		self._z        = None ## Z(t)          == sum of Z_i^l(t) over l,i, total emissions 

	## __getattr__
	## ---------------------------------------------
	def __getattr__(self, key):
		""" Easy-access to config variables """
		if self.cfg.has(key): return self.cfg.get(key).value
		return None

	## __hasattr__
	## ---------------------------------------------
	def __hasattr__(self, key):
		""" Easy-access to config variables """
		return self.cfg.has(key)
		
	## procedure
	## ---------------------------------------------
	def procedure(self):
		""" Main function, systematics loop, everything to be done """
		## start
		self.initPars () ## initializes parameters
		self.load     () ## builds components
		self.initComps() ## initializes components
		## central loop
		self.start() ## start sequence, sets initial values
		self.run  () ## run sequence, iterate through the individual steps of the simulation
		self.end  () ## end sequence, produce and store results
		## loop over systematic variations
		varsToSys = [p for p in self.pars.keys() if len(self.pars[p].items())>1]
		if len(varsToSys)==0: return
		for p in varsToSys:
			for unc,val in self.pars[p].items():
				if unc=="central": continue
				self.start(p, unc) ## start sequence
				self.run  (      ) ## run sequence
				self.end  (p, unc) ## end sequence


	## helper methods
	## =============================================

	## dump
	## ---------------------------------------------
	def dump(self, name):
		""" Save everything of this iteration to disk """
		## set up output folder
		out = "%s%s/"%(self.out, "tmp")
		mkdir(out)
		## prepare the giant csv
		csv = {}
		csv["Model"]      = endPars(self,self)
		csv["Model"]["t"] = [x*self.tstep for x in self.tsteps]
		## climate model
		cm = self.getComp("ClimateModel")
		csv[cm.name] = cm.end()
		## energy sectors
		for i in range(self._nSectors):
			es = self.getComp("EnergySector", i)
			csv[es.name] = es.end()
		## regions and energy agents
		for l in range(self._nRegions):
			r = self.getComp("Region", l)
			csv[r.name] = r.end()
			for i in range(self._nSectors):
				ea = self.getComp("EnergyAgent", l, i)
				csv[ea.name] = ea.end()
		## export the giant csv
		df = pd.DataFrame()
		for comp,vals in csv.items():
			cols = ["%s_%s"%(comp,x) for x in list(vals.keys())]
			df1 = pd.DataFrame(columns=cols)
			for obs,data in vals.items():
				key = "%s_%s"%(comp, obs)
				df1.insert(cols.index(key), key, data, True)
			df.join(df1)
		df.to_csv("%s/%s.csv"%(out, name))

	## end
	## ---------------------------------------------
	def end(self, par=None, var="central"):
		""" End sequence, save everything to disk """
		self.dump("final")
		name = "loop%d_%s_%s"%(self.it, par, var) if par and var else "loop%d"%self.it
		tmp  = "%s%s/"%(self.out, "tmp")
		out  = "%s%s/"%(self.out, name )
		mv(tmp, out)

	## getComp
	## ---------------------------------------------
	def getComp(self, objname, idx1=None, idx2=None):
		""" Retrieve a stored component """
		k = (objname, idx1, idx2)
		if not k in self.objs.keys(): return None
		return self.objs[k]

	## getPar
	## ---------------------------------------------
	def getPar(self, objname, parname, varname="central"):
		""" Retrieve a stored parameter value """
		if not objname in self.pars                  .keys(): return None
		if not parname in self.pars[objname]         .keys(): return None
		if not varname in self.pars[objname][parname].keys(): return None
		return self.pars[objname][parname][varname]

	## initComps
	## ---------------------------------------------
	def initComps(self):
		""" Initialize components """
		## build components
		self.getComp("ClimateModel").init()
		for i in range(self._nSectors):
			self.getComp("EnergySector", i).init()
		for l in range(self._nRegions):
			self.getComp("Region", l).init()
			for i in range(self._nSectors):
				self.getComp("EnergyAgent", l, i).init()

	## initPars
	## ---------------------------------------------
	def initPars(self):
		""" Initialize parameters """
		## set model parameters
		self.tsteps = [t for t in range(1, self.nTsteps+1)]
		extractPars(self, self, ["alpha0","beta","ctry0","ctry1","ccrit","deltaa","deltan","eqq","kbar","nSectors","nRegions","nu0","rho","sigma","vtry1"])
		startPars  (self, self) ## FIXME: not good; should do it with start of others but need nRegions et al BEFORE loading components, so keep it here for now (no syst var though)

	## initVals
	## ---------------------------------------------
	def initVals(self):
		""" Compute initial values, e.g. the initial lifetime resource profit per region """
		## initial capital
		for l in range(self._nRegions):
			r     = self.getComp("Region", l)
			r._k  = r._kinit*self._kbar  ## K_{0,t}^l, page 42
			r._k0 = r._k
		## efficient tax
		self._g      = 0.01 ## FIXME
		cm           = self.getComp("ClimateModel")
		self._taubar = cm._phil/(1-self._beta*(1+self._g)**(1-self._sigma))+cm._phi0*(1-cm._phil)/(1-self._beta*(1+self._g)**(1-self._sigma)*(1-cm._phi)) ## \bar{\tau}^{eff}, Eqn (28)
		## lifetime profit of resource sectors
		yt = 0
		for l in range(self._nRegions):
			r     = self.getComp("Region", l)
			r._pi = 0
			yt   += r._gamma*r._y
			for i in range(self._nSectors):
				es     = self.getComp("EnergySector", i)
				ea     = self.getComp("EnergyAgent" , l, i)
				ea._pi = (es._evv - es._ecc)*ea._r ## Pi_i^l, Eqn (16)
				r._pi += ea._pi
		## initial capital price and tax
		r0         = self.getComp("Region", 0)
		self._err  = self._alpha0*r0._y/r0._k0 ## r_t, Eqn (3)
		self._err0 = self._err ## r_0, initial value for r_t
		self._eqq *= 1./self._err ## q_t, page 6
		self._tau  = self._taubar*yt ## \tau_t^{eff}, Eqn (28)

	## load
	## ---------------------------------------------
	def load(self):
		""" Load the model by building all objects """
		ClimateModel(self)
		for i in range(self._nSectors):
			EnergySector(self, i)
		for l in range(self._nRegions):
			Region(self, l)
			for i in range(self._nSectors):
				EnergyAgent(self, l, i)

	## run
	## ---------------------------------------------
	def run(self):
		""" Core sequence, simulation loop """
		## the simulation loop
		it           = 0
		self.rndcorr = {}
		while self.tui.statlvl!=ERROR and it<self.nAttempts:
			it += 1
			self.tui.msg("Doing iteration %02d"%it)
			## throwing initial values
			self.throw   () ## throw variables
			self.initVals() ## compute initial values
## FIXME: store initial values too..? 
			## main temporal iteration
			v = False
			for self.t in self.tsteps:
				#if any([not c.valid for c in self.objs.values()]): break ## stop if a component failed
				self.tui.msg("Computing t=%02d"%self.t)
				self.runStep()
				self.store()
				v = self.verifyStep()
				if not v: break
				self.runPush()
			## if verification failed, restart with better throwing; else leave
			if not v: 
				self.tui.msg("Per-step verification failed! Dumping and going to restart..")
				self.dump("iteration%02d"%it)
				continue
			## end-of-loop verification	
			if not self.verify(): 
				self.tui.msg("End-of-loop verification failed! Dumping and going to restart..")
				self.dump("iteration%02d"%it)
			else:
				break

	## runPush
	## ---------------------------------------------
	def runPush(self):
		""" Pushes all relevant observables from one timestep to the next """
		## new capital prices
		r0         = self.getComp("Region", 0) ## first region, l=1
		self._err  = self._alpha0*r0._y/r0._k  ## r_t, Eqn (3)
		self._eqq *= 1./self._err              ## q_t, page 6
		## new energy prices
		for l in range(self._nRegions):
			r = self.getComp("Region", l)
			for i in range(self._nSectors):
				es      = self.getComp("EnergySector",    i)
				ea      = self.getComp("EnergyAgent" , l, i)
				ea._epp = self._nu0*r._y*(es._kappa/ea._e)*(ea._e/r._e)**self._rho  ## p_{i,t}^l, Eqn (3)
		## new resource prices
		for i in range(self._nSectors):
			es = self.getComp("EnergySector", i)
			es.run()  ## compute v_{i,t}, Eqn (10)
		## new consumption
		self._c = math.pow(self._beta*self._err, 1./self._sigma)*self._c) ## \bar{C}_t, Eqn (20)

	## runStep
	## ---------------------------------------------
	def runStep(self):
		""" Computation per time step """
		## step Ia: labor allocation
		self._nbar = 0
		for l in range(self._nRegions):
			r   = self.getComp("Region", l)
			r.runEnn()
			self._nbar += r._n ## \bar{N}_t = sum of \bar{N}_t^l over l, page 35 implicit
		self._enn = 1-self._alpha0-self._nu0 ## n_{0,t}^l, page 35
		for l in range(self._nRegions):
			r      = self.getComp("Region", l)
			totenn = self._enn  ## sum of n_{i,t}^l over i, page 35 implicit
			for i in range(self._nSectors):
				es      = self.getComp("EnergySector", i   )
				ea      = self.getComp("EnergyAgent" , l, i)
				ea._enn = (1-es._alpha-es._nu)*self._nu0*ea._eta ## n_{i,t}^l, page 35
				totenn += ea._enn
			r._n0  = self._enn / totenn * self._n   ## N_{0,t}^l, page 35
			for i in range(self._nSectors):
				ea      = self.getComp("EnergyAgent" , l, i)
				ea._n   = ea._enn / totenn * self._n   ## N_{i,t}^l, page 35
		## regional consumption
		tot = 0
		for l in range(self._nRegions):
			r      = self.getComp("Region", l)
			r._eww = self._enn*r._y/r._n0 ## w_t^l, Eqn (3)
			r._w  += self._eqq*r._eww*r._n ## W_t^l, page 7
			zx     = 0
			for i in range(self._nSectors):
				es  = self.getComp("EnergySector", i)
				ea  = self.getComp("EnergyAgent" , l, i)
				zx += es._zeta*ea._x
			r._t  += self._eqq*self._tau*zx ## T^l(t) without \theta^l, Eqn (12)
## FIXME: what about \theta^l?
			tot   += self._err0*r._k0 + r._w + r._pi + r._theta*r._t 
		for l in range(self._nRegions):
			r        = self.getComp("Region", l)
			r._mu    = (self._err0*r._k0 + r._w + r._pi + r._theta*r._t) / tot ## \mu^l, Eqn (19)
			r._c     = r._mu*ctot ## C_t^l, Eqn (19)
			self._c += r._c       ## \bar{C}_t
		## step Ib: capital allocation
		totekk = 0 ## sum of k_{i,t}^l over i,l, page 36
		for l in range(self._nRegions):
			r       = self.getComp("Region", l)
			r._ekk0 = self._alpha0*r._y ## k_{0,t}^l, page 36
			totekk += r._ekk0
			for i in range(self._nSectors):
				es      = self.getComp("EnergySector", i)
				ea      = self.getComp("EnergyAgent" , l, i)
				ea._ekk = es._alpha*self._nu0*ea._eta*r._y ## k_{i,t}^l, page 36
				totekk += ea._ekk  ## k_t, page 36
		kprev      = self._kbar ## \bar{K}_t of previous step, page 36
		self._kbar = 0          ## \bar{K}_t
		for l in range(self._nRegions):
			r           = self.getComp("Region", l)
			r._k        = r._ekk0 / totekk * kprev  ## K_{0,t}^l, page 36
			self._kbar += r._k  ## \bar{K}_t
			for i in range(self._nSectors):
				ea          = self.getComp("EnergyAgent", l, i)
				ea._k       = ea._ekk / totekk * kprev  ## K_{i,t}^l, page 36
				self._kbar += ea._k  ## \bar{K}_t, page 35
		## step Ic: resource allocation
		yt = 0  ## old sum of gamma^l*Y_t^l (recomputed below)
		for l in range(self._nRegions):
			r   = self.getComp("Region", l)
			yt += r._gamma*r._y
		r0      = self.getComp("Region", 0) ## first region, l=1
		self._z = 0  ## sum over Z_{i,t}^l, total emissions
		for l in range(self._nRegions):
			r = self.getComp("Region", l)
			for i in range(self._nSectors):
				es       = self.getComp("EnergySector", i)
				if es._ecc==0: continue 
				ea       = self.getComp("EnergyAgent" , l, i)
				ea._x    = self._nu0*es._nu*r._y * ea._eta / (es._ecc + self._alpha0*r0._y/r0._k*(es._evv-es._ecc) + es._zeta*self._tau*yt)          ## X_{i,t}^l, page 36
				ea._z    = es._zeta*ea._x  ## Z_{i,t}^l, Eqn (5)
				self._z += ea._z           ## Z_t
		## step IIa: climate model
		self.getComp("ClimateModel").run() ## computing S_t
		## step IIb: energy outputs (total E_t^l determined in IIc below)
		for l in range(self._nRegions):
			for i in range(self._nSectors):
				ea = self.getComp("EnergyAgent", l, i)
				ea.run() ## computing E_{t,l}^i
		## step IIc: total production (including E_t^l, tax)
		self._y = 0
		yt      = 0
		for l in range(self._nRegions):
			r = self.getComp("Region", l)
			r.run() ## compute E_t^l and Y_t^l
			self._y += r._y
			yt      += r._gamma*r._y 
			for i in range(self._nSectors):
				es = self.getComp("EnergySector", i)
				ea = self.getComp("EnergyAgent" , l, i)
				ea._eta = es._kappa * (ea._e/r._e)**self._rho  ## \eta_{i,t}^l, Eqn (32)
		self._tau = self._taubar*yt ## \tau_t^{eff}, Eqn (28)

	## setComp
	## ---------------------------------------------
	def setComp(self, obj, idx1=None, idx2=None):
		""" Store a component """
		self.objs[(obj.__class__.__name__, idx1, idx2)] = obj

	## setPar
	## ---------------------------------------------
	def setPar(self, objname, value, parname, varname="central"):
		""" Store a parameter value """
		if not objname in self.pars         .keys(): self.pars[objname]          = {}
		if not parname in self.pars[objname].keys(): self.pars[objname][parname] = {}
		self.pars[objname][parname][varname] = value

	## start
	## ---------------------------------------------
	def start(self, parname=None, varname="central"):
		""" Start sequence, set initial values of model and component """
		## increment counter for successful attempts
		self.it+=1
		## initialize components
		self.getComp("ClimateModel").start(par=parname, var=varname)
		for i in range(self._nSectors):
			self.getComp("EnergySector", i).start(par=parname, var=varname)
		for l in range(self._nRegions):
			self.getComp("Region", l).start(par=parname, var=varname)
			for i in range(self._nSectors):
				self.getComp("EnergyAgent", l, i).start(par=parname, var=varname)		

	## store
	## ---------------------------------------------
	def store(self):
		""" Save a evolvable observable in buffer """
		storePars(self)
		self.getComp("ClimateModel").store()
		for i in range(self._nSectors):
			self.getComp("EnergySector", i).store()
		for l in range(self._nRegions):
			self.getComp("Region", l).store()
			for i in range(self._nSectors):
				self.getComp("EnergyAgent", l, i).store()

	## throw
	## ---------------------------------------------
	def throw(self):
		""" Throw initial values for variables """
		self._c = throw(self._c if self._c else self.ctry0, \
		                self.ctry1, \
		                direction=self.rndcorr["c"] if "c" in self.rndcorr.keys() else None) ## \bar{C}
		for i in range(self._nSectors):
			es      = self.getComp("EnergySector", i)
			es._evv = es._ecc
			if es._r==0: continue
			idx     = "v%d"%i
			es._evv = throw(es._ecc, self.vtry1, minimum=es._ecc, incl=False, \
			                direction=self.rndcorr[idx] if idx in self.rndcorr.keys() else None) ## v_i

	## verify
	## ---------------------------------------------
	def verify(self):
		""" End-of-loop verification """
		self.rndcorr = {}
		## resources
		for i in range(self._nSectors):
			es  = self.getComp("EnergySector", i)
			tot = es._r
			for l in range(self._nRegions):
				ea   = self.getComp("EnergyAgent" , l, i)
				tot -= ea._x
			if tot<0: 
				self.rndcorr["v%d"%i] = "dn"
				return False
			if tot>es._rcrit:
				self.rndcorr["v%d"%i] = "up"
				return False
		## all good
		return True

	## verifyStep
	## ---------------------------------------------
	def verifyStep(self):
		""" Verification after every iteration step """
		self.rndcorr = {}
		## capital
		if self._kbar < 0:
			self.rndcorr["c"] = "dn"
			return False
		## consumption
		crit = 0		
		for l in range(self._nRegions):
			r = self.getComp("Region", l)
			crit += r._y
			for i in range(self._nSectors):
				es = self.getComp("EnergySector", i)
				ea = self.getComp("EnergyAgent" , l, i)
				crit += es._ecc*ea._x
		self._cCrit = self._ccrit*crit ## C_t^{crit}
		if self._c < self._cCrit:
			self.rndcorr["c"] = "up"
			return False
		## all good
		return True



