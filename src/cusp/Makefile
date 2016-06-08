all: vc
vc:
	nvcc -w -O3 -I./ -I~/cusplibrary-0.5.1 vertex_coloring.cpp -o vc
	cp $@ $(BIN)
		
clean:
	rm vc
