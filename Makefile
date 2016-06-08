TOPLEVEL := .
APPS := csrcolor serial topo data GM
INPUT_URL := http://
INPUT := csrcolor-inputs.tar.bz2

.PHONY: all clean

all: $(APPS)

$(APPS):
	make -C src/$@

#include common.mk

inputs:
	@echo "Downloading inputs ..."
	@wget $(INPUT_URL) -O $(INPUT)
	@echo "Uncompressing inputs ..."
	@tar xvf $(INPUT)
	@rm $(INPUT)
	@echo "Inputs available at $(TOPLEVEL)/inputs/"

clean:
	for APP in $(APPS); do make -C apps/$$APP clean; done

