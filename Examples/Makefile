#!/bin/bash
CC = ./pochoir
ICC = icpc
# it seems that '-ipo' '-O1' optimization flag will cause unexpected segmentation fault
# sometimes.
# OPT_FLAGS = -O3 -DNDEBUG -Wall -Werror -unroll-aggressive -funroll-loops -ipo -xHOST -fno-alias -fno-fnalias -std=c++0x -I${POCHOIR_LIB_PATH}
OPT_FLAGS = -O3 -DNDEBUG -Wall -Werror -unroll-aggressive -funroll-loops -xHOST -fno-alias -fno-fnalias -std=c++0x -I${POCHOIR_LIB_PATH}
# POCHOIR_DEBUG_FLAGS = -O0 -g3 -DDEBUG -debug -Wall -Werror -std=c++0x -I${CILK_HEADER_PATH} -include cilk_stub.h
POCHOIR_DEBUG_FLAGS = -O0 -g3 -debug -Wall -Werror -std=c++0x -I${CILK_HEADER_PATH} -include cilk_stub.h
ICC_DEBUG_FLAGS = -O0 -g3 -DDEBUG -Wall -Werror -std=c++0x -I${CILK_HEADER_PATH} -include cilk_stub.h
# ICC_DEBUG_FLAGS = -O1 -Wall -Werror -std=c++0x 
POCHOIR_LIB_FILES = "../pochoir.hpp ../pochoir_common.hpp ../pochoir_types.hpp ../pochoir_array.hpp ../pochoir_walk.hpp ../pochoir_walk_recursive.hpp ../pochoir_walk_loops.hpp ../pochoir_wrapper.hpp ../pochoir_range.hpp"

# ICC_DEBUG_FLAGS = -O0 -g3 -DDEBUG -Wall -Werror -std=c++0x
example_1D : tb_example_1D.cpp
#   Phase-II compilation
	${CC} -o tb_example_1D_pochoir ${OPT_FLAGS} tb_example_1D.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_example_1D_pochoir ${POCHOIR_DEBUG_FLAGS} tb_example_1D.cpp

heat_1D_NP : tb_heat_1D_NP.cpp
#   Phase-II compilation
	${CC} -o tb_heat_1D_NP_pochoir ${OPT_FLAGS} tb_heat_1D_NP.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_heat_1D_NP_pochoir ${POCHOIR_DEBUG_FLAGS} tb_heat_1D_NP.cpp

heat_NP : tb_heat_2D_NP.cpp
#   Phase-II compilation
	${CC} -o tb_heat_2D_NP_pochoir ${OPT_FLAGS} -split-caching tb_heat_2D_NP.cpp
#	${CC} -o tb_heat_2D_NP_pochoir ${OPT_FLAGS} tb_heat_2D_NP.cpp
#	${CC} -o tb_heat_2D_NP_pochoir ${OPT_FLAGS} -split-c-pointer tb_heat_2D_NP.cpp
#	${ICC} -o heat_2D_NP ${ICC_DEBUG_FLAGS} tb_heat_2D_NP_pochoir.cpp
#	${ICC} -o heat_2D_NP ${OPT_FLAGS}  tb_heat_2D_NP_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_heat_2D_NP_pochoir ${POCHOIR_DEBUG_FLAGS} tb_heat_2D_NP.cpp

heat_NP_zero : tb_heat_2D_NP_zero.cpp
#   Phase-II compilation
	${CC} -o tb_heat_2D_NP_zero_pochoir ${OPT_FLAGS} -split-pointer tb_heat_2D_NP_zero.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_heat_2D_NP_zero_pochoir ${POCHOIR_DEBUG_FLAGS} tb_heat_2D_NP_zero.cpp

heat_P : tb_heat_2D_P.cpp
#   Phase-II compilation
#	${CC} -o heat_2D_P ${OPT_FLAGS} -all-cond-tile-macro tb_heat_2D_P.cpp
	${ICC} -o heat_2D_P_gdb ${ICC_DEBUG_FLAGS} tb_heat_2D_P_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_heat_2D_P_pochoir ${POCHOIR_DEBUG_FLAGS} tb_heat_2D_P.cpp

heat_3D_NP : tb_heat_3D_NP.cpp
#   Phase-II compilation
	${CC} -o tb_heat_3D_NP_pochoir ${OPT_FLAGS} -split-pointer tb_heat_3D_NP.cpp
#	${ICC} -o tb_heat_3D_NP_pochoir ${OPT_FLAGS} tb_heat_3D_NP_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_heat_3D_NP_pochoir ${POCHOIR_DEBUG_FLAGS} tb_heat_3D_NP.cpp

animwave_irregular : tb_animwave_irregular.cpp tb_animwave_irregular_pochoir.cpp 
#   Phase-II compilation
	${CC} -o animwave_irregular ${OPT_FLAGS} tb_animwave_irregular.cpp
#	${CC} -o animwave_irregular ${POCHOIR_DEBUG_FLAGS} -split-macro-shadow tb_animwave_irregular.cpp
#	${ICC} -o animwave_irregular ${OPT_FLAGS} tb_animwave_irregular_pochoir.cpp
#	${ICC} -o animwave_irregular ${ICC_DEBUG_FLAGS} tb_animwave_irregular_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o animwave_irregular ${POCHOIR_DEBUG_FLAGS} tb_animwave_irregular.cpp 

diamond : tb_diamond.cpp # tb_diamond_pochoir.cpp 
#   Phase-II compilation
	${CC} -o diamond ${OPT_FLAGS} tb_diamond.cpp
#	${CC} -o diamond ${POCHOIR_DEBUG_FLAGS} -split-macro-shadow tb_diamond.cpp
#	${ICC} -o diamond ${OPT_FLAGS} tb_diamond_pochoir.cpp
#	${ICC} -o diamond_gdb ${ICC_DEBUG_FLAGS} tb_diamond_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o diamond ${POCHOIR_DEBUG_FLAGS} tb_diamond.cpp 

diamond_1 : tb_diamond_1.cpp # tb_diamond_1_pochoir.cpp 
#   Phase-II compilation
	${CC} -o diamond_1 ${OPT_FLAGS} -unroll-multi-kernel tb_diamond_1.cpp
#	${CC} -o diamond_1 ${OPT_FLAGS} tb_diamond_1.cpp
#	${ICC} -o diamond_1_gdb ${ICC_DEBUG_FLAGS} tb_diamond_1_pochoir.cpp
#	${CC} -o diamond_1 ${POCHOIR_DEBUG_FLAGS} -split-macro-shadow tb_diamond_1.cpp
#	${ICC} -o diamond_1 ${OPT_FLAGS} tb_diamond_1_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o diamond_1 ${POCHOIR_DEBUG_FLAGS} tb_diamond_1.cpp 
#	${ICC} -o diamond_1 ${OPT_FLAGS} tb_diamond_1.cpp 

diamond_2 : tb_diamond_2.cpp # tb_diamond_2_pochoir.cpp 
#   Phase-II compilation
#	${CC} -o diamond_2 ${OPT_FLAGS} tb_diamond_2.cpp
#	${CC} -o diamond_2 ${OPT_FLAGS} tb_diamond_2.cpp
#	${ICC} -o diamond_2_gdb ${ICC_DEBUG_FLAGS} tb_diamond_2_pochoir.cpp
#	${CC} -o diamond_2 ${POCHOIR_DEBUG_FLAGS} -split-macro-shadow tb_diamond_2.cpp
#	${ICC} -o diamond_2 ${OPT_FLAGS} tb_diamond_2_pochoir.cpp
#	Phase-I compilation with debugging aid
	${CC} -o diamond_2 ${POCHOIR_DEBUG_FLAGS} tb_diamond_2.cpp 
#	${ICC} -o diamond_2 ${OPT_FLAGS} tb_diamond_2.cpp 

tile_2D : tb_tile_2D.cpp # tb_tile_2D_pochoir.cpp 
#   Phase-II compilation
#	${CC} -o tile_2D ${OPT_FLAGS} -unroll-t-tile-opt-pointer tb_tile_2D.cpp
#	${ICC} -o tile_2D_gdb ${ICC_DEBUG_FLAGS} tb_tile_2D_pochoir.cpp
#	${CC} -o tile_2D ${POCHOIR_DEBUG_FLAGS} -split-macro-shadow tb_tile_2D.cpp
#	${ICC} -o tile_2D ${OPT_FLAGS} tb_tile_2D_pochoir.cpp
#	Phase-I compilation with debugging aid
	${CC} -o tile_2D ${POCHOIR_DEBUG_FLAGS} tb_tile_2D.cpp 
#	${ICC} -o tile_2D ${OPT_FLAGS} tb_tile_2D.cpp 

tile_3D : tb_tile_3D.cpp # tb_tile_3D_pochoir.cpp 
#   Phase-II compilation
#	${CC} -o tile_3D ${OPT_FLAGS} -unroll-t-tile-pointer tb_tile_3D.cpp
#	${ICC} -o tile_3D_gdb ${ICC_DEBUG_FLAGS} tb_tile_3D_pochoir.cpp
#	${CC} -o tile_3D ${POCHOIR_DEBUG_FLAGS} -split-macro-shadow tb_tile_3D.cpp
#	${ICC} -o tile_3D ${OPT_FLAGS} tb_tile_3D_pochoir.cpp
#	Phase-I compilation with debugging aid
	${CC} -o tile_3D ${POCHOIR_DEBUG_FLAGS} tb_tile_3D.cpp 
#	${ICC} -o tile_3D ${OPT_FLAGS} tb_tile_3D.cpp 

overlap_tile_2D : tb_overlap_tile_2D.cpp # tb_overlap_tile_2D_pochoir.cpp 
#   Phase-II compilation
#	${CC} -o overlap_tile_2D ${OPT_FLAGS} -all-cond-tile-macro-overlap tb_overlap_tile_2D.cpp
#	${ICC} -o overlap_tile_2D_gdb ${ICC_DEBUG_FLAGS} tb_overlap_tile_2D_pochoir.cpp
#	Phase-I compilation with debugging aid
	${CC} -o overlap_tile_2D ${POCHOIR_DEBUG_FLAGS} tb_overlap_tile_2D.cpp 
#	${ICC} -o overlap_tile_2D ${OPT_FLAGS} tb_overlap_tile_2D.cpp 

overlap_tile_3D : tb_overlap_tile_3D.cpp # tb_overlap_tile_3D_pochoir.cpp 
#   Phase-II compilation
#	${CC} -o overlap_tile_3D ${OPT_FLAGS} -unroll-t-tile-pointer-overlap tb_overlap_tile_3D.cpp
#	${ICC} -o overlap_tile_3D_gdb ${ICC_DEBUG_FLAGS} tb_overlap_tile_3D_pochoir.cpp
#	Phase-I compilation with debugging aid
	${CC} -o overlap_tile_3D ${POCHOIR_DEBUG_FLAGS} tb_overlap_tile_3D.cpp 
#	${ICC} -o overlap_tile_3D ${OPT_FLAGS} tb_overlap_tile_3D.cpp 

overlap_stagger_2D : tb_overlap_stagger_2D.cpp # tb_overlap_stagger_2D_pochoir.cpp 
#   Phase-II compilation
	${CC} -o overlap_stagger_2D ${OPT_FLAGS} -unroll-t-tile-pointer-overlap tb_overlap_stagger_2D.cpp
	${ICC} -o overlap_stagger_2D_gdb ${ICC_DEBUG_FLAGS} tb_overlap_stagger_2D_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o overlap_stagger_2D ${POCHOIR_DEBUG_FLAGS} tb_overlap_stagger_2D.cpp 
#	${ICC} -o overlap_stagger_2D ${OPT_FLAGS} tb_overlap_stagger_2D.cpp 

overlap_2D : tb_overlap_2D.cpp # tb_overlap_2D_pochoir.cpp 
#   Phase-II compilation
	${CC} -o overlap_2D ${OPT_FLAGS} -all-cond-tile-pointer-overlap tb_overlap_2D.cpp
	${ICC} -o overlap_2D_gdb ${ICC_DEBUG_FLAGS} tb_overlap_2D_pochoir.cpp
#	${CC} -o overlap_2D ${POCHOIR_DEBUG_FLAGS} -split-macro-shadow tb_overlap_2D.cpp
#	${ICC} -o overlap_2D ${OPT_FLAGS} tb_overlap_2D_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o overlap_2D ${POCHOIR_DEBUG_FLAGS} tb_overlap_2D.cpp 
#	${ICC} -o overlap_2D ${OPT_FLAGS} tb_overlap_2D.cpp 

tile_stagger : tb_tile_stagger.cpp # tb_tile_stagger_pochoir.cpp 
#   Phase-II compilation
	${CC} -o tile_stagger ${OPT_FLAGS} -unroll-t-tile-pointer tb_tile_stagger.cpp
#	${ICC} -o tile_stagger_gdb ${ICC_DEBUG_FLAGS} tb_tile_stagger_pochoir.cpp
#	${CC} -o tile_stagger ${POCHOIR_DEBUG_FLAGS} -split-macro-shadow tb_tile_stagger.cpp
#	${ICC} -o tile_stagger ${OPT_FLAGS} tb_tile_stagger_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tile_stagger ${POCHOIR_DEBUG_FLAGS} tb_tile_stagger.cpp 
#	${ICC} -o tile_stagger ${OPT_FLAGS} tb_tile_stagger.cpp 

tile_simple : tb_tile_simple.cpp # tb_tile_simple_pochoir.cpp 
#   Phase-II compilation
	${CC} -o tile_simple ${OPT_FLAGS} -unroll-t-tile-macro tb_tile_simple.cpp
#	${ICC} -o tile_simple_gdb ${ICC_DEBUG_FLAGS} tb_tile_simple_pochoir.cpp
#	${CC} -o tile_simple ${POCHOIR_DEBUG_FLAGS} -split-macro-shadow tb_tile_simple.cpp
#	${ICC} -o tile_simple ${OPT_FLAGS} tb_tile_simple_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tile_simple ${POCHOIR_DEBUG_FLAGS} tb_tile_simple.cpp 
#	${ICC} -o tile_simple ${OPT_FLAGS} tb_tile_simple.cpp 

meta : tb_meta_meta.cpp
	icpc -o meta tb_meta_meta.cpp -std=c++0x

func : tb_function_object.cpp
	${ICC} -o func tb_function_object.cpp -std=c++0x

animwave : tb_animwave.cpp
#   Phase-II compilation
#	${CC} -o animwave_pochoir ${OPT_FLAGS} tb_animwave.cpp
#	Phase-I compilation with debugging aid
	${CC} -o animwave_pochoir ${POCHOIR_DEBUG_FLAGS} tb_animwave.cpp

apop : apop.cpp
#   Phase-II compilation
	${CC} -o apop_pochoir ${OPT_FLAGS} apop.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o apop_pochoir ${POCHOIR_DEBUG_FLAGS} apop.cpp

life : tb_life.cpp
#   Phase-II compilation
	${CC} -o tb_life_pochoir ${OPT_FLAGS} tb_life.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_life_pochoir ${POCHOIR_DEBUG_FLAGS} tb_life.cpp

psa_struct : tb_psa_struct.cpp
#   Phase-II compilation
#	${CC} -o tb_psa_struct_pochoir ${OPT_FLAGS} tb_psa_struct.cpp
#	Phase-I compilation with debugging aid
	${CC} -o tb_psa_struct_pochoir ${POCHOIR_DEBUG_FLAGS} tb_psa_struct.cpp

rna : tb_rna.cpp
#   Phase-II compilation
	${CC} -o tb_rna_pochoir ${OPT_FLAGS} tb_rna.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_rna_pochoir ${POCHOIR_DEBUG_FLAGS} tb_rna.cpp

lcs : tb_lcs.cpp
#   Phase-II compilation
#	${CC} -o tb_lcs_pochoir ${OPT_FLAGS} tb_lcs.cpp
#	Phase-I compilation with debugging aid
	${CC} -o tb_lcs_pochoir ${POCHOIR_DEBUG_FLAGS} tb_lcs.cpp

berkeley3d7pt : tb_3d7pt.cpp
#   Phase-II compilation
#	${CC} -o tb_3d7pt_pochoir ${OPT_FLAGS} tb_3d7pt.cpp
#	Phase-I compilation with debugging aid
	${CC} -o tb_3d7pt_pochoir ${POCHOIR_DEBUG_FLAGS} tb_3d7pt.cpp

berkeley3d27pt : tb_3d27pt.cpp
#   Phase-II compilation
#	${CC} -o tb_3d27pt_pochoir ${OPT_FLAGS} tb_3d27pt.cpp
#	Phase-I compilation with debugging aid
	${CC} -o tb_3d27pt_pochoir ${POCHOIR_DEBUG_FLAGS} tb_3d27pt.cpp

tb_3dfd : tb_3dfd.cpp
#   Phase-II compilation
	${CC} -o tb_3dfd_pochoir ${OPT_FLAGS} tb_3dfd.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_3dfd_pochoir ${POCHOIR_DEBUG_FLAGS} tb_3dfd.cpp

clean: 
	rm -f *.o *.i *_pochoir *_pochoir.cpp *.out
