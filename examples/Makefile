#!/bin/bash
CC = ./pochoir
ICC = icpc
OPT_FLAGS = -unroll-pointer 
POCHOIR_DEBUG_FLAGS = -debug
POCHOIR_LIB_FILES = "../pochoir.hpp ../pochoir_common.hpp ../pochoir_dloader.hpp ../pochoir_array.hpp ../pochoir_walk.hpp ../pochoir_walk_recursive.hpp ../pochoir_walk_loops.hpp ../pochoir_kernel.hpp ../pochoir_range.hpp"

# FILTER_OUT(substring, array)
# Filters out all instances of strings containing substring in array
FILTER_OUT = $(foreach v,$(2),$(if $(findstring $(1),$(v)),,$(v)))

# list all tb_*.cpp files and filter out any _pochoir/_kernel_info files
TARGETS = \
	$(call FILTER_OUT,_kernel_info,\
	$(call FILTER_OUT,_gen_kernel,\
	$(call FILTER_OUT,_pochoir,\
	$(patsubst tb_%.cpp,%,$(wildcard tb_*.cpp)\
))))

CXXFLAGS = $(OPT_FLAGS)
ifdef DEBUG
	CXXFLAGS = $(POCHOIR_DEBUG_FLAGS)
endif

all: $(TARGETS)

%: tb_%.cpp
	${CC} -o $@ $(CXXFLAGS) $<

check: all
	python check.py

clean: 
	rm -f *.o *.i *_pochoir *_pochoir.cpp *.out *.dat *_kernel_info.cpp *_gen_kernel.cpp *.so *_gdb $(TARGETS)
