include ../common.mk
all: gc
gc:	graphColoring.cu graphDriver.cpp tree.cpp
	$(NVCC) $(NVFLAGS) graphDriver.cpp graphColoring.cu -o gc
	mv $@ $(BIN)

clean:
	rm gc
