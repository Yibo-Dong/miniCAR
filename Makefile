SOLVER_FLAG = -DMINISAT
SOLVER_SOURCES = src/thirdParty/sat/minisat/core/Solver.cc src/thirdParty/sat/minisat/utils/Options.cc src/thirdParty/sat/minisat/utils/System.cc src/thirdParty/sat/minisat/simp/SimpSolver.cc
SOLVER_OBJECTS = Solver.o Options.o System.o SimpSolver.o
SOLVER_INCLUDE_DIRS = -I./src/thirdParty/sat -I./src/thirdParty/sat/minisat -I./src/thirdParty/sat/minisat/core -I./src/thirdParty/sat/minisat/utils -Isrc/thirdParty/sat/minisat/simp


SEED ?= 0
ifeq ($(SEED),0)
RANDOM_FLAG= 
else
RANDOM_FLAG=-DRANDSEED=$(SEED)
endif

ifeq ($(OPT),NONE)
OPT_LEVEL = 
else
OPT_LEVEL = -Ofast
endif


SOURCES = src/thirdParty/aiger/aiger.c \
			src/checker/carChecker.cpp src/checker/bmcChecker.cpp \
			src/solver/carsolver.cpp src/solver/mainsolver.cpp src/solver/newpartialsolver.cpp src/solver/implysolver.cpp\
			src/definition.cpp src/main.cpp \
			$(SOLVER_SOURCES)
INCLUDE_DIRS = -I./src -I./src/thirdParty/aiger -I./src/checker -I./src/debug -I./src/solver -I./src/utils $(SOLVER_INCLUDE_DIRS)
OBJECTS = carsolver.o implysolver.o newpartialsolver.o mainsolver.o main.o aiger.o $(SOLVER_OBJECTS) carChecker.o bmcChecker.o definition.o

CFLAG = $(INCLUDE_DIRS) -D__STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS -c -g
LFLAG = -g -lz -lpthread 
OPTFLAG = $(OPT_LEVEL) -march=native -frename-registers -funroll-loops -fno-signed-zeros # consider using -fprofile-generate and -fprofile-use

GCC = gcc
GXX = g++

TARGET ?= miniCAR

MODE ?= RELEASE
ifeq ($(MODE),DEBUG)
MODE_MACRO = -DDEBUG
else
MODE_MACRO = -DRELEASE
endif

OPTIONS= ${RANDOM_FLAG} ${SOLVER_FLAG} ${MODE_MACRO}

miniCAR:$(SOURCES)
	$(GCC) ${OPTFLAG} ${OPTIONS} $(CFLAG) $(SOURCES)
	$(GXX) ${OPTFLAG} -g -o ${TARGET} $(OBJECTS) $(LFLAG)

clean: 
	rm -f *.o
	rm -f miniCAR

doc:$(SOURCES) doc/Doxyfile
	cd doc && doxygen Doxyfile

.PHONY: miniCAR clean
