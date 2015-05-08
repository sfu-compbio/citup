
#------------------------------------------------
#  makefile for CITUP
#------------------------------------------------

# change the following two lines according to your system
BOOST_DIRECTORY=/home/andrew/install/include/boost
CPLEX_DIRECTORY=/home/andrew/install/ibm/ILOG/CPLEX_Studio125
CPLEX_BUILD=x86-64_sles10_4.1

CC=g++

all: generatingInput.exe AnalyzeOutputs.exe CplexApplication.exe CitupIter.exe
	mv *.exe ../bin
	rm *.o
	
generatingInput.exe: generatingInput.o generatingFrequencies.o
	$(CC) -ggdb generatingInput.o generatingFrequencies.o -o generatingInput.exe -I$(BOOST_DIRECTORY)
generatingInput.o: generatingInput.cpp generatingFrequencies.h
	$(CC) -ggdb -c generatingInput.cpp -o generatingInput.o -I../External/ -I$(BOOST_DIRECTORY)
generatingFrequencies.o: generatingFrequencies.cpp generatingFrequencies.h
	$(CC) -ggdb -c generatingFrequencies.cpp -o generatingFrequencies.o -I$(BOOST_DIRECTORY)

AnalyzeOutputs.exe: AnalyzeOutputs.cpp
	$(CC) -ggdb AnalyzeOutputs.cpp -o AnalyzeOutputs.exe -O3 -I../External/ -I$(BOOST_DIRECTORY)

CplexApplication.exe: CplexApplication.cpp
	$(CC) -ggdb CplexApplication.cpp -o CplexApplication.exe -O3 -DIL_STD -DILOSTRICTPOD -I../External/ -I$(CPLEX_DIRECTORY)/concert/include -I$(CPLEX_DIRECTORY)/cplex/include -L$(CPLEX_DIRECTORY)/concert/lib/$(CPLEX_BUILD)/static_pic/ -lconcert -L$(CPLEX_DIRECTORY)/cplex/lib/$(CPLEX_BUILD)/static_pic -l ilocplex -l cplex -l m -l pthread

CitupIter.exe: CitupIter.cpp trees.h
	$(CC) -ggdb CitupIter.cpp -o CitupIter.exe -O3 -DIL_STD -DILOSTRICTPOD -I../External/ -I$(BOOST_DIRECTORY) -I$(CPLEX_DIRECTORY)/concert/include -I$(CPLEX_DIRECTORY)/cplex/include -L$(CPLEX_DIRECTORY)/concert/lib/$(CPLEX_BUILD)/static_pic/ -lconcert -L$(CPLEX_DIRECTORY)/cplex/lib/$(CPLEX_BUILD)/static_pic -l ilocplex -l cplex -l m -l pthread

	
