import matplotlib.pyplot as plt
from core.functions import mkdir


## makeCfgKey
## -------------------------------------------------
def makeCfgKey(v, idx1, idx2):
	""" Assembles the key of the variable as in the cfg """
	## note: indices start at 0 in python but at 1 in the cfg!
	if type(idx1) not in [int, type(None)] or type(idx2) not in [int, type(None)]: 
		return None
	if idx1 is not None and idx2 is not None:
		return "%s_%d_%d"%(v, idx1+1, idx2+1)
	if idx1 is not None and idx2 is     None:
		return "%s_%d"   %(v, idx1+1)
	if idx1 is     None and idx2 is not None:
		return "%s_%d"   %(v, idx2+1)
	return v

## endPars
## -------------------------------------------------
def endPars(model, obj):
	""" Save time series for each observable to disk; i.e. both plot
	and return values to csv for each observables """
	mkdir("%s/plots"%model.tmp)
	tsteps = model.tsteps
	csv    = {}
	for k,ts in obj.cache.items():
		res = []
		for t in [0]+tsteps:
			if not t in ts.keys(): res.append(0); continue
			res.append(list(ts.values())[list(ts.keys()).index(t)])
		plt.clf()
		plt.plot(res)
		plt.xlabel("time steps")
		plt.ylabel(k)
		plt.savefig("%s/plots/%s_%s.pdf"%(model.tmp, obj.name, k))
		csv[k] = res
	return csv

## extractPars
## -------------------------------------------------
def extractPars(model, obj, varsToLookFor):
	""" Extracts a given set of parameters from the model's cfg for a given component """
	found = []
	for v in varsToLookFor:
		for key,entry in model.cfg.cache.cache.items():
			k = makeCfgKey(v, obj.idx1, obj.idx2)
			if not k or k!=entry.key: continue
			if type(entry.value)==list and len(entry.value)==3:
				model.setPar(obj.name, entry.value[0], v, "down"   )
				model.setPar(obj.name, entry.value[1], v, "central")
				model.setPar(obj.name, entry.value[2], v, "up"     )
			elif type(entry.value) in [int, float]:
				model.setPar(obj.name, entry.value   , v, "central")
			else:
				continue
			found.append(v)
			break ## go to next v in varsToLookFor
	return found

## setPars
## ---------------------------------------------
def setPars(obj):
	""" Defines the set of observables for this component """
	obj.obs = [x[1:] for x in dir(obj) if x[0:1]=="_" and x[1:2]!="_"]

## startPars
## -------------------------------------------------
def startPars(model, obj, par=None, var="central"):
	""" (Re-)Initialize parameters of this object and reset cache """
	obj.cache = {}
	for p in obj.obs:
		setattr(obj, "_%s"%p, model.getPar(obj.name, p, var if p==par else "central"))
		obj.cache[p] = {} ## reset cache
	## FIXME: review this -- which starting value are we taking? are they OK for new variations??

## storePars
## -------------------------------------------------
def storePars(obj):
	""" Save a evolvable observable in buffer """
	for obs in obj.obs:
		obj.cache[obs][obj.t] = getattr(obj, "_%s"%obs)


## Component
## =================================================
class Component(object):
	""" Parent class for region, sector, model """

	## __init__
	## ---------------------------------------------
	def __init__(self, m):
		self.m     = m     ## pointer to market model
		self.t     = -1    ## time value for sync
		self.l     = None  ## placeholder for region index
		self.i     = None  ## placeholder for energy source index
		self.obs   = []    ## list of observables
		self.cache = {}    ## cache for observable values per timestep
		self.valid = False

	## __str__
	## ---------------------------------------------
	def __str__(self):
		""" Returns a name for this component """
		return self.name()

	## end
	## ---------------------------------------------
	def end(self):
		""" Save time series for each observable to disk; i.e. both plot
		and return values to csv for each observables """
		if not self.m.out: return
		return endPars(self.m, self)

	## idx1
	## ---------------------------------------------
	@property
	def idx1(self):
		""" Returns the first index """
		return self.l if self.l is not None else self.i if self.i is not None else None

	## idx2
	## ---------------------------------------------
	@property
	def idx2(self):
		""" Returns the second index """
		return self.i if self.l is not None and self.i is not None else None

	## init
	## ---------------------------------------------
	def init(self, varsToLookFor):
		""" Extract all variables from cfg """
		self.valid = False
		found      = extractPars(self.m, self, varsToLookFor)
		if len(found)!=len(varsToLookFor):
			n           = makeCfgKey(self.__class__.__name__, self.idx1, self.idx2)
			varsMissing = filter(lambda x: x not in found, varsToLookFor)
			self.m.tui.error("Cannot find all variables for %s in config! %s"%(n, varsMissing))
			return
		self.valid = True

	## name
	## ---------------------------------------------
	@property
	def name(self):
		""" Returns a name for this component """
		return makeCfgKey(self.__class__.__name__, self.idx1, self.idx2)

	## observables
	## ---------------------------------------------
	@property
	def observables(self):
		""" Returns a list of all observables of this component """
		return self.obs

	## start
	## ---------------------------------------------
	def start(self, par=None, var="central"):
		""" (Re-)Initialize parameters of this object and reset cache """
		if not self.valid: return 
		self.t = -1 ## reset time value
		startPars(self.m, self, par, var)

	## store
	## ---------------------------------------------
	def store(self):
		""" Save a evolvable observable in buffer """
		storePars(self)		



