#!/bin/bash
CC = ../pochoir
ICC = icpc
OPT_FLAGS = -DCHECK_SHAPE -rdynamic -ldl -O2 -DNDEBUG -Wall -Werror -unroll-aggressive -funroll-loops -xHOST -fno-alias -fno-fnalias -fp-model precise -std=c++0x 
POCHOIR_DEBUG_FLAGS = -O0 -g3 -DDEBUG -debug -Wall -Werror -std=c++0x -I${INTEL_CILK_HEADER} -include cilk_stub.h
ICC_DEBUG_FLAGS = -O0 -g3 -DDEBUG -std=c++0x -I${INTEL_CILK_HEADER} -include cilk_stub.h

# FILTER_OUT(substring, array)
# Filters out all instances of strings containing substring in array
FILTER_OUT = $(foreach v,$(2),$(if $(findstring $(1),$(v)),,$(v)))

# list all tb_*.cpp files and filter out any _pochoir/_kernel_info files
TARGETS = \
	$(call FILTER_OUT,_kernel_info,\
	$(call FILTER_OUT,_pochoir,\
	$(patsubst tb_%.cpp,%,$(wildcard tb_*.cpp)\
)))

FLAG = DEFAULT
ifdef DEBUG_POCHOIR
	TARGETS := $(addsuffix _pochoir, $(TARGETS))
	FLAG = DEBUG_POCHOIR
endif

ifdef DEBUG_ICC
	TARGETS := $(addsuffix _gdb, $(TARGETS))
	FLAG = DEBUG_ICC
endif

all: $(TARGETS)

%: tb_%.cpp
	${CC} -o $@ $(OPT_FLAGS) $<

%_pochoir: tb_%.cpp
	${CC} -o $@ $(POCHOIR_DEBUG_FLAGS) $<

%_gdb: %_pochoir
	${CC} -o $@ $(ICC_DEBUG_FLAGS) tb_$<.cpp

check: all
	python check.py $(FLAG)

clean: 
	rm -f *.o *.i *_pochoir *_gdb *_pochoir.cpp *_kernel_info.cpp *_gen_kernel.cpp *_gen_kernel.so *_color.dat *.out *.csv $(TARGETS)
