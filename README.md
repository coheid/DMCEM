# DMCEM

## Introduction

`DMCEM` is a software package implementing the simulation of a
Dynamic Multi-Country Equilibrium Model (DMCEM) presented originally
in <a href="https://doi.org/10.1016/j.jet.2018.11.001">https://doi.org/10.1016/j.jet.2018.11.001</a>.
More recent publications involve the April 2020, February 2021, and November 2021 versions of
<a href="http://www.marten-hillebrand.de/research/IMCC_B.pdf">"Who pays the bills?"</a>.


## Program Structure

### Principal Architecture

`DMCEM` is written in python and uses the (private) `WHElight` package.

The program is run through a model class (`Model`, inherits from `Master`) that relies on 
dedicated components, which are implemented as children to their common `Component` class.
A component is an entity that appears standalone within the model and has parameters 
(which have an initial value and evolve in the course of the simulation) as well as a label. 
The components are either the climate model (`ClimateModel`), the different economic 
regions (`Region`), energy sectors (`EnergySector`), or energy agents (`EnergyAgent`). 
An energy agent is the implementation of an energy sector within a given region. 
Thus, components may be equipped with indices corresponding to the indices 
<img src="https://render.githubusercontent.com/render/math?math=i \in \textbf{I}"> and
<img src="https://render.githubusercontent.com/render/math?math=\ell \in \textbf{L}">.
Unlike to the studies citied above, the program does not know dedicated components for the 
final, production, or consumer sectors. Instead, parameters and formulae attributable to these
sectors are absorbed in either the regional components or the model class itself. For example,
the parameter <img src="https://render.githubusercontent.com/render/math?math=\alpha_0">
appears in the computation of the total output, thus, in the 
final sectors of all regions. However, since it is the same parameter for all regions, 
it is made a member of the `Model` class (c.f. next section). In that 
sense, the `Model` class has properties similar to those of the `Component`. 


### Observables

An observable is an economic variable that evolves with time as the simulation progresses, 
e.g. <img src="https://render.githubusercontent.com/render/math?math=Y_t">. An observable 
can belong to a component in case it is logically connected
to it or there are multiple incarnations of the same observable, one for every component
(e.g. one can compute the production output per region, 
<img src="https://render.githubusercontent.com/render/math?math=Y_t^\ell">). 
In case an observable belongs to a component it is indexed. We use the index `i` for the 
energy sectors and `l` for the economic regions. Energy agents (i.e. energy sectors within 
regions) use both indices `l` and `i`. 

Observables are represented as python variables, in fact as members to one of the 
aforementioned components (in case they belong to a component) or the model itself (else). 
The member name of an observable always begins with a single underscore, `_`, e.g. the
observable <img src="https://render.githubusercontent.com/render/math?math=\alpha_0">
is written as `self._alpha0`. Since the Mac filesystem is not case 
sensitive, all python variables including observable members are written in lower-case only. 
However, theory observables come in lower- and upper-cases. To distinguish upper-case observables 
(e.g. <img src="https://render.githubusercontent.com/render/math?math=N_t">)
from lower-case ones
(e.g. <img src="https://render.githubusercontent.com/render/math?math=n_t">), 
we use an `exx` convention:
* upper-case observables (e.g. <img src="https://render.githubusercontent.com/render/math?math=N_t">) 
  are written as is (e.g. `self._n`)
* lower-case observables (e.g. <img src="https://render.githubusercontent.com/render/math?math=n_t">) 
  are written as `exx` (e.g. `self._enn`)
This is motivated by the way of pronouncing the letter n ("enn"); still the convention
is applied to all lower-case obsevables (e.g. `_eqq`, `_eaa`, and so forth).


### Configuration

For each observable, the initial value can be given to the simulation via a configuration (cfg) 
file. The configuration is a text file that obeys the format required by the corresponding 
`WHElight` parser module. In particular, the following rules apply:
* every line specifies an observable (or a non-observable parameter of the program)
* every line is separated into three parts via the delimiter `::`; the first part defines 
  the data type of the parameter or observable value (`bool`, `float`, `int`, `str`, etc.); 
  the second part defines the name of the observable or parameter, while the third part
  defines its values
* the default data type is string (`str`), which can be omitted (then only name and value
  are given, separated by `::`)
* comments can be defined via the hash symbol (`#`), i.e. everything right to the symbol
  will be ignored by the parser
* observable names are given without the preceding underscore, while any indices of the 
  component that the observable belongs to (if any) are given via underscores; note that 
  indexing of the economic observables as well as the naming scheme in the cfg 
  file always start at 1, while indexing within the python program always starts at 0;
  for example:
  - `q_2_1` defines the initial value of the observable 
    <img src="https://render.githubusercontent.com/render/math?math=Q_{1,t}^2">
    which is the sector productivity of the first energy sector `i=1` in region `l=2`;
    if both indices are present, `l` is always considered first and `i` comes second;
  - `alpha_1` defines the initial value of the observable 
    <img src="https://render.githubusercontent.com/render/math?math=\alpha_1"> which is 
    the elasticity of the first energy sector `i=1`
* observable names may contain digits too; in this case one should give the digit _without_
  preceeding underscore; for example: `alpha0` defines the initial value of the observable 
  <img src="https://render.githubusercontent.com/render/math?math=\alpha_0"> which is the 
  elasticity in the final sector (since it is the same for all regions, it is a member of 
  the `Model` class)
* observable values can be passed as scalars or vectors; in the latter case, a list of
  three elements is required corresponding to the central value (second element) and
  the down (first) and up (third) variations; in this way, `DMCEM` will run systematic
  variations to determine the systematic uncertainty of the observable in question [UNTESTED!]


### Methods

The core method is `procedure` in the `Model` class. It first extracts the
obsevables that belong to the model, builds the components, and then extracts their
parameters too. Eventually, the simulation loop is started, consisting of the 
successive calls to the methods `start` (setting initial values of exogenous observables), 
`run` (evolving all observables), and `end` (exporting the values of each observable at
each point in time to a `csv` file). If systematic variations are specified, they are
run therafter once the nominal loop has finished.

Every component has at least the following four methods:
* `init` is responsible for retrieving the observable's initial values from the cfg file
  (it extracts all initial values for all systematic variations)
* `start` picks up the proper initial value for each observable to be used in the 
  iteration therafter; this means selecting among the initial values that `init` has retrieved
* `store` saves the values of all observables after every point in time
* `end` assembles a dictionary of all saved values for all points in time that is used in 
  the export to `csv`


## Installation

`DMCEM` is installed by downloading it through github into a new folder:

```
git clone git@github.com:coheid/DMCEM.git
```

The program requires the (private) `WHElight` core to be installed and configured as well
according to the prescriptions therein. 


## Program Execution

`DMCEM` is executed by passing a configuration file containing all initial values for all
observables:

```
python run.py <path-to-cfg>
```

In addition one can specify an alternative working directory through the option `-d`
(see below).


## Output

In the working directory specified by the user (per default it is the installation 
directory), `DMCEM` produces a sub-folder structure `output/<cfg-name>/` where all
files are dumped. The output produced comprises
* `log.out`: a logfile that allows to a posteriori reconstruct the program verbosity
* a `csv` file for every iteration of the temporal loop; also iterations that were
  not verified as OK are dumped; the output of the good iteration is called `final.csv`
* a series of plots depicting the time series of every observable [UNTESTED!]

