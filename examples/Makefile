#!/bin/bash
CC = ./pochoir
ICC = icpc
OPT_FLAGS = -unroll-pointer 
POCHOIR_DEBUG_FLAGS = -debug
POCHOIR_LIB_FILES = "../pochoir.hpp ../pochoir_common.hpp ../pochoir_dloader.hpp ../pochoir_array.hpp ../pochoir_walk.hpp ../pochoir_walk_recursive.hpp ../pochoir_walk_loops.hpp ../pochoir_kernel.hpp ../pochoir_range.hpp"


#${CC} -o fib_gdb ${POCHOIR_DEBUG_FLAGS} fib.cpp
fib : fib.cpp
#   Phase-II compilation
	${CC} -o fib ${OPT_FLAGS} fib.cpp	
#	Phase-I compilation with debugging aid
#	${CC} -o fib_gdb ${POCHOIR_DEBUG_FLAGS} fib.cpp 

overlap_tile_2D : tb_overlap_tile_2D.cpp 
#   Phase-II compilation
#	${CC} -o overlap_tile_2D ${OPT_FLAGS} tb_overlap_tile_2D.cpp
#	Phase-I compilation with debugging aid
	${CC} -o overlap_tile_2D_gdb ${POCHOIR_DEBUG_FLAGS} tb_overlap_tile_2D.cpp 

overlap_tile_3D : tb_overlap_tile_3D.cpp
#   Phase-II compilation
#	${CC} -o overlap_tile_3D ${OPT_FLAGS} tb_overlap_tile_3D.cpp
#	Phase-I compilation with debugging aid
	${CC} -o overlap_tile_3D_gdb ${POCHOIR_DEBUG_FLAGS} tb_overlap_tile_3D.cpp 

clean: 
	rm -f *.o *.i *_pochoir *_pochoir.cpp *.out *.dat *_kernel_info.cpp *_gen_kernel.cpp *.so *_gdb fib overlap_tile_2D overlap_tile_3D
