include ../common.mk
INCLUDES += -I../../cub-1.1.1
EXTRA := $(NVFLAGS) $(INCLUDES) -DITERATIONS=10
SRC=main.cu

all: topo_base
topo_base:
	$(NVCC) $(EXTRA) $(SRC) -o topo_base
	mv $@ $(BIN)

clean:
	rm topo_naive topo_ldg
