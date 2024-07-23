ifeq ($(SOLVER),MINISAT)
SOLVER_FLAG = -DMINISAT
SOLVER_SOURCES = src/sat/minisat/core/Solver.cc src/sat/minisat/utils/Options.cc src/sat/minisat/utils/System.cc
SOLVER_OBJECTS = Solver.o Options.o System.o 
SOLVER_INCLUDE_DIRS = -I./src/sat -I./src/sat/minisat -I./src/sat/minisat/core -I./src/sat/minisat/utils 
else ifeq ($(SOLVER),GLUCOSE4)
SOLVER_FLAG = -DGLUCOSE
SOLVER_SOURCES = src/sat/glucose-4.2.1/core/Solver.cc src/sat/glucose-4.2.1/utils/Options.cc src/sat/glucose-4.2.1/utils/System.cc src/sat/glucose-4.2.1/core/lcm.cc
SOLVER_OBJECTS = Solver.o Options.o System.o lcm.o
SOLVER_INCLUDE_DIRS = -I./src/sat -I./src/sat/glucose-4.2.1 -I./src/sat/glucose-4.2.1/core -I./src/sat/glucose-4.2.1/utils 
else
SOLVER_FLAG = -DGLUCOSE
SOLVER_SOURCES = src/sat/glucose/core/Solver.cc src/sat/glucose/utils/Options.cc src/sat/glucose/utils/System.cc
SOLVER_OBJECTS = Solver.o Options.o System.o 
SOLVER_INCLUDE_DIRS = -I./src/sat -I./src/sat/glucose -I./src/sat/glucose/core -I./src/sat/glucose/utils
endif

SEED ?= 0
ifeq ($(SEED),0)
RANDOM_FLAG= 
else
RANDOM_FLAG=-DRANDSEED=$(SEED)
endif

SOURCES = src/model/aiger.c \
			src/checker/carChecker.cpp src/checker/bmcChecker.cpp \
			src/solver/carsolver.cpp src/solver/mainsolver.cpp src/solver/newpartialsolver.cpp \
			src/utils/definition.cpp src/solver/implysolver.cpp\
			src/main.cpp \
			$(SOLVER_SOURCES)
INCLUDE_DIRS = -I./src -I./src/model -I./src/checker -I./src/debug -I./src/solver -I./src/utils $(SOLVER_INCLUDE_DIRS)
OBJECTS = carsolver.o implysolver.o newpartialsolver.o mainsolver.o main.o aiger.o $(SOLVER_OBJECTS) carChecker.o bmcChecker.o definition.o

CFLAG = $(INCLUDE_DIRS) -D__STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS -c -g
LFLAG = -g -lz -lpthread 
OPTFLAG = -Ofast -march=native -frename-registers -funroll-loops -fno-signed-zeros # consider using -fprofile-generate and -fprofile-use

GCC = gcc
GXX = g++

TARGET ?= miniCAR

OPTIONS= ${RANDOM_FLAG} ${SOLVER_FLAG}

miniCAR:$(SOURCES)
	$(GCC) ${OPTFLAG} ${OPTIONS} $(CFLAG) $(SOURCES)
	$(GXX) ${OPTFLAG} -g -o ${TARGET} $(OBJECTS) $(LFLAG)

clean: 
	rm *.o
	rm miniCAR

.PHONY: miniCAR clean
