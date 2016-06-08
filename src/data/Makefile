include ../common.mk
INCLUDES += -I../../cub-1.1.1
B40_DIR=./back40computing-read-only
B40C_INC=-I$(B40_DIR) -I$(B40_DIR)/test
EXTRA := $(NVFLAGS) $(CFLAGS) $(INCLUDES) -DITERATIONS=10 -DBLKSIZE=128
SRC=main.cu

all: data_base data_bitset data_wlc data_ldb data_fusion data_tc data_ldg data_pq data_best

data_base:
	$(NVCC) $(EXTRA) $(SRC) -DVARIANT=0 -o $@
	mv $@ $(BIN)

data_ldg:
	$(NVCC) $(EXTRA) $(SRC) -DVARIANT=1 -o $@
	mv $@ $(BIN)

data_bitset:
	$(NVCC) $(EXTRA) $(SRC) -DVARIANT=2 -o $@
	mv $@ $(BIN)

data_ldb:
	$(NVCC) $(EXTRA) $(B40C_INC) $(SRC) -DVARIANT=6 -o $@
	mv $@ $(BIN)

data_wlc:
	$(NVCC) $(EXTRA) $(SRC) -DVARIANT=5 -o $@
	mv $@ $(BIN)

data_fusion:
	$(NVCC) $(EXTRA) $(SRC) -DVARIANT=4 -o $@
	mv $@ $(BIN)

data_tc:
	$(NVCC) $(EXTRA) $(SRC) -DVARIANT=3 -o $@
	mv $@ $(BIN)

data_pq:
	$(NVCC) $(EXTRA) $(SRC) -DVARIANT=7 -o $@
	mv $@ $(BIN)

data_best:
	$(NVCC) $(EXTRA) $(SRC) -DVARIANT=8 -o $@
	mv $@ $(BIN)

data_comb1:
	$(NVCC) $(EXTRA) $(SRC) -DVARIANT=9 -o $@
	mv $@ $(BIN)

data_comb2:
	$(NVCC) $(EXTRA) $(SRC) -DVARIANT=9 -o $@
	mv $@ $(BIN)

clean:
	rm data_base data_ldg data_bitset data_ldb data_fusion data_tc data_pq

