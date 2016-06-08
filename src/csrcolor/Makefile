include ../common.mk
EXE=csrcolor
all: csrcolor.cu
	$(NVCC) $(NVFLAGS) $(INCLUDES) -lcusparse csrcolor.cu -o $(EXE)
	mv $(EXE) $(BIN)

clean:
	rm csrcolor
