#!/bin/bash
CC = ./pochoir
ICC = icpc
OPT_FLAGS = -O3 -DNDEBUG -Wall -Werror -unroll-aggressive -funroll-loops -xHOST -fno-alias -fno-fnalias -fp-model precise -std=c++0x
POCHOIR_DEBUG_FLAGS = -O0 -g3 -DDEBUG -debug -Wall -Werror -std=c++0x
ICC_DEBUG_FLAGS = -O0 -g3 -DDEBUG -Wall -Werror -std=c++0x -I${INTEL_CILK_HEADER} -include cilk_stub.h  
TARGETS = heat_1D_NP heat_2D_NP heat_2D_NP_zero heat_2D_P heat_3D_NP apop life psa_struct rna lcs berkeley3d7pt berkeley3d27pt 3dfd 3dfd_test

ifdef DEBUG
	#	Phase-I compilation with debugging aid
	CFLAGS = ${POCHOIR_DEBUG_FLAGS}
else
	# Phase-II compilation
	CFLAGS = ${OPT_FLAGS}
endif

all: ${TARGETS}

heat_1D_NP : tb_heat_1D_NP.cpp
	${CC} -o tb_heat_1D_NP_pochoir ${CFLAGS} tb_heat_1D_NP.cpp

heat_2D_NP : tb_heat_2D_NP.cpp
	${CC} -o tb_heat_2D_NP_pochoir ${CFLAGS} tb_heat_2D_NP.cpp

heat_2D_NP_zero : tb_heat_2D_NP_zero.cpp
	${CC} -o tb_heat_2D_NP_zero_pochoir ${CFLAGS} tb_heat_2D_NP_zero.cpp

heat_2D_P : tb_heat_2D_P.cpp
	${CC} -o tb_heat_2D_P_pochoir ${CFLAGS} tb_heat_2D_P.cpp
ifdef DEBUG
	${ICC} -o tb_heat_2D_P_pochoir_gdb ${ICC_DEBUG_FLAGS} tb_heat_2D_P_pochoir.cpp
endif

heat_3D_NP : tb_heat_3D_NP.cpp
	${CC} -o tb_heat_3D_NP_pochoir ${CFLAGS} tb_heat_3D_NP.cpp

apop : apop.cpp
	${CC} -o apop_pochoir ${CFLAGS} apop.cpp

life : tb_life.cpp
	${CC} -o tb_life_pochoir ${CFLAGS} tb_life.cpp

psa_struct : tb_psa_struct.cpp
	${CC} -o tb_psa_struct_pochoir ${CFLAGS} tb_psa_struct.cpp

animwave_struct : tb_animwave_struct.cpp
#   Phase-II compilation
	${CC} -o tb_animwave_struct_pochoir -split-c-pointer ${OPT_FLAGS} tb_animwave_struct.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_animwave_struct_pochoir ${POCHOIR_DEBUG_FLAGS} tb_animwave_struct.cpp

parallel_heat_1D_NP : tb_parallel_heat_1D_NP.cpp
#   Phase-II compilation
	${CC} -o tb_parallel_heat_1D_NP_pochoir -split-c-pointer ${OPT_FLAGS} tb_parallel_heat_1D_NP.cpp
	${CC} -o tb_parallel_heat_1D_NP_pochoir_gdb ${POCHOIR_DEBUG_FLAGS} tb_parallel_heat_1D_NP_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_parallel_heat_1D_NP_pochoir ${POCHOIR_DEBUG_FLAGS} tb_parallel_heat_1D_NP.cpp

reaction_diffusion_1D : tb_reaction_diffusion_1D.cpp
#   Phase-II compilation
	${CC} -o tb_reaction_diffusion_1D_pochoir -split-c-pointer ${OPT_FLAGS} tb_reaction_diffusion_1D.cpp
#	${CC} -o tb_reaction_diffusion_1D_pochoir_gdb ${POCHOIR_DEBUG_FLAGS} tb_reaction_diffusion_1D_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_reaction_diffusion_1D_pochoir ${POCHOIR_DEBUG_FLAGS} tb_reaction_diffusion_1D.cpp

rna : tb_rna.cpp
	${CC} -o tb_rna_pochoir ${CFLAGS} tb_rna.cpp

lcs : tb_lcs.cpp
	${CC} -o tb_lcs_pochoir ${CFLAGS} tb_lcs.cpp

berkeley3d7pt : tb_3d7pt.cpp
	${CC} -o tb_berkeley3d7pt_pochoir ${CFLAGS} tb_3d7pt.cpp

berkeley3d27pt : tb_3d27pt.cpp
	${CC} -o tb_berkeley3d27pt_pochoir ${CFLAGS} tb_3d27pt.cpp

3dfd : tb_3dfd.cpp
	${CC} -o tb_3dfd_pochoir ${CFLAGS} tb_3dfd.cpp

3dfd_test : tb_3dfd_test.cpp
	${CC} -o tb_3dfd_test_pochoir ${CFLAGS} tb_3dfd_test.cpp

# Remnants from a time gone by?
cilk_for_1D : tb_cilk_for_1D.cpp
#   Phase-II compilation
	${CC} -o cfor_1D ${OPT_FLAGS} tb_cilk_for_1D.cpp
	${ICC} -o cfor_1D_gdb ${ICC_DEBUG_FLAGS} tb_cilk_for_1D_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_cilk_for_1D_pochoir ${POCHOIR_DEBUG_FLAGS} tb_cilk_for_1D.cpp

cilk_for_2D : tb_cilk_for_2D.cpp
#   Phase-II compilation
	${CC} -o cfor_2D ${OPT_FLAGS} tb_cilk_for_2D.cpp
	${ICC} -o cfor_2D_gdb ${ICC_DEBUG_FLAGS} tb_cilk_for_2D_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_cilk_for_2D_pochoir ${POCHOIR_DEBUG_FLAGS} tb_cilk_for_2D.cpp

check:
	python ./tests.py

clean: 
	rm -f *.o *.i *_pochoir *_pochoir.cpp *.out
