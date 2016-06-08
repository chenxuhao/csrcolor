CXX=g++ 
CXX_FLAGS=-w -O3
EXTRA=-DITERATIONS=10
SRC=main.cc

all:color-omp color-serial

color-omp:
	$(CXX) $(CXX_FLAGS) $(EXTRA) -fopenmp -DENABLE_OPENMP $(SRC) -o $@
	mv $@ ../../bin

color-serial:
	$(CXX) $(CXX_FLAGS) $(EXTRA) $(SRC) -o $@

clean:
	rm color-omp color-serial
