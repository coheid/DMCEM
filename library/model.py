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
		self.tmp       = "%s/tmp/"%self.out
		mkdir(self.out)
		mkdir(self.tmp)
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
		self._cbar     = None ## \bar{C}(t)    == global consumption [cfg if dothrow==False]
		self._cbarexpl = None ## \bar{C}(t)    == global consumption taken from Excel sheet explicitly [cfg]
		self._cbarprev = None ## \bar{C}(t)    == global consumption of previous step
		self._ctry0    = None ##               == mean     for throwing C [cfg]
		self._ctry1    = None ##               == variance for throwing C [cfg]
		self._ccrit    = None ## C^{crit}      == critical value for C [cfg]
		self._deltaa   = None ## \delta_a      == convergence for labor productivity [cfg]
		self._deltan   = None ## \delta_n      == convergence for population size [cfg]
		self._ekk      = None ## k(t)          == sum of k_{i,t}^l over i,l
		self._enn      = None ## n_0^l(t)      == labor cost share of final sector
		self._eqq      = None ## q(t)          == discount factor [cfg]
		self._eqqprev  = None ## q(t)          == discount factor of previous step
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
		setPars(self)

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
		return ## FIXME
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
		mkdir(self.tmp)
		## prepare the giant csv
		csv = {}
		csv[self.name]      = endPars(self, self)
		csv[self.name]["t"] = [self.tstart+x*self.tstep for x in [0]+self.tsteps]
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
		order = list(csv.keys())
		order.sort()
		order = [self.name] + [x for x in order if x!=self.name] ## make sure the model itself is first
		for comp in order:
			vals = csv[comp]
			data = {"%s_%s"%(comp, k): v for k,v in vals.items()}
			df1  = pd.DataFrame.from_dict(data)
			if len(df.index)==0: 
				df = df1.copy()
				continue
			df = df.join(df1)
		df.to_csv("%s/%s.csv"%(self.tmp, name))

	## end
	## ---------------------------------------------
	def end(self, par=None, var="central"):
		""" End sequence, save everything to disk """
		self.dump("final")
		x    = "LF" if self.dolf else "OT"
		name = "%s_sysvar_%s_%s"%(x, par, var) if par and var else "%s_central"%x
		out  = "%s%s/"%(self.out, name)
		mv(self.tmp, out)
		## FIXME: copy logfile too? (but only after execution..)

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
		extractPars(self, self, ["alpha0","beta","cbar","ctry0","ctry1","ccrit","deltaa","deltan","eqq","kbar","nSectors","nRegions","nu0","rho","sigma","vtry1"])
		startPars  (self, self) ## FIXME: not good; should do it with start of others but need nRegions et al BEFORE loading components, so keep it here for now (no syst var though)
		if self.cfg.cache.has("cbarexpl"):
			## this is a hack, but we take cbarexpl from config without making it an observable
			self._cbarexpl = self.cfg.cache.get("cbarexpl").value

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
		self._taubar = 0 if self.dolf else self._taubar ## no tax in laissez-faire
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
		## initial discount parameter and tax
		self._eqq      = 1.0
		self._eqqprev  = 1.0
		self._cbarprev = self._cbar
		self._tau      = 0.0
		#self._tau  = self._taubar*yt ## \tau_t^{eff}, Eqn (28)

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
			if self.dothrow: self.throw() ## throw variables
			## compute pre-iteration step
			self.t = 0
			self.tui.msg("Computing t=%02d"%self.t)
			self.initVals()
			for jt in range(self.nIts):
				self.runStepOne(jt, True)
				self.runStepTwo(jt, True)
			self.store()
			## push to first step (t=1)
			self.runPush(True)
			## main temporal iteration
			v = False
			for self.t in self.tsteps:
				#if any([not c.valid for c in self.objs.values()]): break ## stop if a component failed ## FIXME
				self.tui.msg("Computing t=%02d"%self.t)
				for jt in range(self.nIts):
					self.runStepOne(jt)
					self.runStepTwo(jt)
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

	## runErr
	## ---------------------------------------------
	def runErr(self, idx, isInit=False):
		""" Computes capital and resource prices and pushes Cbar """
		## compute capital price
		errold    = self._err
		r0        = self.getComp("Region", 0)
		self._err = self._alpha0*r0._y/r0._k ## r_t, Eqn (3)
		## new consumption
		if not isInit:
			if self._cbarexpl: self._cbar = self._cbarexpl[self.t] ## take \bar{C}_t values from Excel sheet, explicitly (avoid shooting per iteration)
			else             : self._cbar = math.pow(self._beta*self._err, 1.0/self._sigma)*self._cbarprev ## \bar{C}_t, Eqn (20)
		## new resource prices and consumption
		if not isInit:
			for i in range(self._nSectors):
				es = self.getComp("EnergySector", i)
				es.runEvv()	
		## initial capital price
		if isInit: self._err0 = self._err
		else     : self._eqq  = self._eqqprev*(1.0/self._err)

	## runPush
	## ---------------------------------------------
	def runPush(self, isFirst=False):
		""" Pushes all relevant observables from one timestep to the next """
		## new parameters for labor supply
		for l in range(self._nRegions):
			r = self.getComp("Region", l)
			r.pushEnn()
		## new capital
		yt = 0
		ct = 0
		for l in range(self._nRegions):
			r   = self.getComp("Region", l)
			yt += r._y
			for i in range(self._nSectors):
				es = self.getComp("EnergySector",    i)
				ea = self.getComp("EnergyAgent" , l, i)
				ct += es._ecc*ea._x
		self._kbar = self._y - self._cbar - ct  ## \bar{K}_{t+1}, Eqn (25)
		## push "prev" variables
		self._eqqprev  = self._eqq
		self._cbarprev = self._cbar
		for i in range(self._nSectors):
			es = self.getComp("EnergySector", i)
			es.pushEvv()	
		self.getComp("ClimateModel").pushSs()

	## runStepOne
	## ---------------------------------------------
	def runStepOne(self, idx, isInit=False):
		""" Computations of individual variables within one time step """
		## step Ia: labor allocation
		self._nbar = 0
		for l in range(self._nRegions):
			r           = self.getComp("Region", l)
			r._n        = r._enn*r._eaa  ## N_t^{l,s} = \bar{N}_t^l, page 19
			self._nbar += r._n           ## \bar{N}_t = sum of \bar{N}_t^l over l, page 35 implicit
		self._enn = 1-self._alpha0-self._nu0 ## n_{0,t}^l, page 35
		for l in range(self._nRegions):
			r      = self.getComp("Region", l)
			totenn = self._enn  ## sum of n_{i,t}^l over i, page 35 implicit
			for i in range(self._nSectors):
				es      = self.getComp("EnergySector", i   )
				ea      = self.getComp("EnergyAgent" , l, i)
				ea._enn = (1-es._alpha-es._nu)*self._nu0*ea._eta  ## n_{i,t}^l, page 35
				totenn += ea._enn
			r._n0  = self._enn / totenn * r._n  ## N_{0,t}^l, page 35
			for i in range(self._nSectors):
				ea      = self.getComp("EnergyAgent" , l, i)
				ea._n   = ea._enn / totenn * r._n  ## N_{i,t}^l, page 35
		## step Ib: capital allocation
		totekk = 0  ## sum of k_{i,t}^l over i,l, page 36
		for l in range(self._nRegions):
			r       = self.getComp("Region", l)
			r._ekk0 = self._alpha0*r._y  ## k_{0,t}^l, page 36
			totekk += r._ekk0
			for i in range(self._nSectors):
				es      = self.getComp("EnergySector", i)
				ea      = self.getComp("EnergyAgent" , l, i)
				ea._ekk = es._alpha*self._nu0*ea._eta*r._y  ## k_{i,t}^l, page 36
				totekk += ea._ekk  ## k_t, page 36
		for l in range(self._nRegions):
			r    = self.getComp("Region", l)
			r._k = r._ekk0 / totekk * self._kbar  ## K_{0,t}^l, page 36
			for i in range(self._nSectors):
				ea    = self.getComp("EnergyAgent", l, i)
				ea._k = ea._ekk / totekk * self._kbar  ## K_{i,t}^l, page 36
		## capital and resource prices and consumption
		self.runErr(idx, isInit)
		## regional consumption
		tot = 0
		for l in range(self._nRegions):
			r      = self.getComp("Region", l)
			r._eww = self._enn*r._y/r._n0   ## w_t^l, Eqn (3)
			r._w  += self._eqq*r._eww*r._n  ## W_t^l, page 7
			zx     = 0
			for i in range(self._nSectors):
				es  = self.getComp("EnergySector", i)
				ea  = self.getComp("EnergyAgent" , l, i)
				zx += es._zeta*ea._x
			tf     = self._eqq*self._tau*zx  ## T^l(t) without \theta^l, Eqn (12)
## FIXME: what about \theta^l?
			r._t  += 0 if self.dolf else tf  ## no transfer in laissez-faire
			tot   += self._err0*r._k0 + r._w + r._pi + r._theta*r._t 
		for l in range(self._nRegions):
			r     = self.getComp("Region", l)
			r._mu = (self._err0*r._k0 + r._w + r._pi + r._theta*r._t) / tot  ## \mu^l, Eqn (19)
			r._c  = r._mu*self._cbar  ## C_t^l, Eqn (19)
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
				ea._x    = self._nu0*es._nu*r._y * ea._eta / (es._evv + es._zeta*self._tau*yt)  ## X_{i,t}^l, page 36
				ea._x = round(ea._x, 13)
				ea._z    = es._zeta*ea._x   ## Z_{i,t}^l, Eqn (5)
				self._z += ea._z            ## Z_t

	## runStepTwo
	## ---------------------------------------------
	def runStepTwo(self, idx, isInit=False):
		""" Compute eta and Y """
		## step IIa: climate model
		self.getComp("ClimateModel").run()  ## computing S_t
		## step IIb: energy outputs (total E_t^l determined in IIc below)
		for l in range(self._nRegions):
			for i in range(self._nSectors):
				ea = self.getComp("EnergyAgent", l, i)
				ea.run()  ## computing E_{t,l}^i
		## step IIc: total production (including E_t^l, tax)
		self._y = 0
		yt      = 0
		for l in range(self._nRegions):
			r = self.getComp("Region", l)
			r.run()  ## compute E_t^l and Y_t^l
			self._y += r._y
			yt      += r._gamma*r._y 
			for i in range(self._nSectors):
				es = self.getComp("EnergySector", i)
				ea = self.getComp("EnergyAgent" , l, i)
				ea._eta = es._kappa * (ea._e/r._e)**self._rho  ## \eta_{i,t}^l, Eqn (32)
		self._tau = self._taubar*yt  ## \tau_t^{eff}, Eqn (28)
		## energy prices
		for l in range(self._nRegions):
			r = self.getComp("Region", l)
			for i in range(self._nSectors):
				es      = self.getComp("EnergySector",    i)
				ea      = self.getComp("EnergyAgent" , l, i)
				ea._epp = self._nu0*r._y*(es._kappa/ea._e)*(ea._e/r._e)**self._rho  ## p_{i,t}^l, Eqn (3)

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
		self._cbar = throw(self._c if self._c else self.ctry0, \
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
		if self._cbar < self._cCrit:
			self.rndcorr["c"] = "up"
			return False
		## all good
		return True

