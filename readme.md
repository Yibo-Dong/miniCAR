# MiniCAR

### Organization
```
├── Makefile					// Makefile of this project.
├── bench					// benchmark from HWMCC 15' and 17'
│   └── notsafe					// unsafe / unknown
│       ├── 139442p1.aig
│       ├── ...
│       └── viselevatorp2.aig
├── doc						// documents, use Doxygen to generate. see `make doc` in Makefile
│   └── Doxyfile				// configuration for Doxygen
├── readme.md					// this guide
├── src						// source files
│   ├── checker					// the checkers. They make use of other modules, and do the model checking work. If you want an overview of algorithm, see here.
│   │   ├── bmcChecker.cpp			// a sample BMC checker, which unrolls the transition function `T` and queries the SAT solver. It is not actively maintained.
│   │   ├── bmcChecker.h			 
│   │   ├── carChecker.cpp			// the CAR checker. Since I focus on bug-finding, I use backward-CAR (which is forward search, i.e., use T rather than T^{-1} }).
│   │   └── carChecker.h			// Therefore, something we do not pay attention to: enumerate ~p, use partial states, ... 
│   ├── definition.cpp				// Basic definitions, along with basic utility.
│   ├── definition.h
│   ├── main.cpp				// entrance.
│   ├── solver					// wrappers for SAT solver.
│   │   ├── carsolver.cpp			// the basic wrapper out of SAT solver, which is the basis of others
│   │   ├── carsolver.h
│   │   ├── implysolver.cpp			// a variant to do subsumption test, i.e. UC' \in? set(UC)
│   │   ├── implysolver.h
│   │   ├── invsolver.h				// a variant to check for variants.
│   │   ├── mainsolver.cpp
│   │   ├── mainsolver.h			// the main solver during the search.
│   │   ├── newpartialsolver.cpp		// partial solver, used in forward-CAR to get partial states, i.e., states with partial assignments.
│   │   ├── newpartialsolver.h
│   │   └── startsolver.h			// start solver, used to enumerate states from `\neg P`.
│   ├── thirdParty				// third party sources.
│   │   ├── aiger				// the aiger suit, used to read in AIGER file.
│   │   └── sat					// sat solver that we use.
│   └── utils
│       └── statistics.h			// statistics that we collect.

```
### GET STARTED

#### Build

```shell
# compile it
make -j4 TARGET=<TARGET_NAME>
# or just
make -j4
```

By default, the `TARGET_NAME` is `miniCAR`(as written in `Makefile`).

#### Run

```shell
# make a target directory
mkdir -p test
# run it
# for more options (other than --vb), check `main.cpp` and `definition.h` for more details
miniCAR --vb bench/notsafe/counterp0.aig test
# it shall finish in a few seconds. Afterwards, there will be two files in test:
ls test
# counterp0.log  counterp0.res
# for verify the result, aiger suit (by Armin Biere) is needed.
<pathToAiger>/aigsim -c bench/notsafe/counterp0.aig test/counterp0.res | echo 'passed'
# passed
```

### TIPS

1. The code is messy due to incorporation of **restart** techniques. If you don't care about it, just ignore `ppstoped`,`restart_enabled` and terms like this.
2. Backward-CAR does forward search, i.e., from $I$ towards $\neg P$. We choose to use an approximation of $\neg P$ (represented by a set of UCs). Every time a state reaches this approximation, we do a `finalCheck` to see whether it is in $\neg P$. 