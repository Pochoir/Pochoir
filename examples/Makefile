#!/bin/bash
CC = ../pochoir
ICC = icpc
OPT_FLAGS = -O3 -DNDEBUG -Wall -Werror -funroll-loops -xHOST -fno-alias -fno-fnalias -fp-model precise -std=c++0x
POCHOIR_DEBUG_FLAGS = -O0 -g3 -DDEBUG -debug -Wall -Werror -std=c++0x
#ICC_DEBUG_FLAGS = -O0 -g3 -DDEBUG -Wall -Werror -std=c++0x -I${INTEL_CILK_HEADER} -include cilk_stub.h  
ICC_DEBUG_FLAGS = -O3 -Wall -Werror -std=c++0x -I${INTEL_CILK_HEADER} -include cilk_stub.h  
heat_1D_NP : tb_heat_1D_NP.cpp
#   Phase-II compilation
	${CC} -o heat_1D_NP ${OPT_FLAGS} tb_heat_1D_NP.cpp
#	${ICC} -o heat_1D_NP_gdb ${ICC_DEBUG_FLAGS} tb_heat_1D_NP_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_heat_1D_NP_pochoir ${POCHOIR_DEBUG_FLAGS} tb_heat_1D_NP.cpp

fib : tb_fib.cpp
#   Phase-II compilation
	${CC} -o fib ${OPT_FLAGS} tb_fib.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o fib ${POCHOIR_DEBUG_FLAGS} tb_fib.cpp

heat_NP : tb_heat_2D_NP.cpp
#   Phase-II compilation
	${CC} -o heat_2D_NP ${OPT_FLAGS} tb_heat_2D_NP.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o heat_2D_NP ${POCHOIR_DEBUG_FLAGS} tb_heat_2D_NP.cpp

heat_NP_zero : tb_heat_2D_NP_zero.cpp
#   Phase-II compilation
	${CC} -o heat_2D_NP_zero ${OPT_FLAGS} tb_heat_2D_NP_zero.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o heat_2D_NP_zero ${POCHOIR_DEBUG_FLAGS} tb_heat_2D_NP_zero.cpp

heat_P : tb_heat_2D_P.cpp
#   Phase-II compilation
	${CC} -o heat_2D_P ${OPT_FLAGS} tb_heat_2D_P.cpp
#	${ICC} -o heat_2D_P_gdb ${ICC_DEBUG_FLAGS} tb_heat_2D_P_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o heat_2D_P ${POCHOIR_DEBUG_FLAGS} tb_heat_2D_P.cpp

heat_3D_NP : tb_heat_3D_NP.cpp
#   Phase-II compilation
	${CC} -o heat_3D_NP ${OPT_FLAGS} tb_heat_3D_NP.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o heat_3D_NP ${POCHOIR_DEBUG_FLAGS} tb_heat_3D_NP.cpp

apop : apop.cpp
#   Phase-II compilation
	${CC} -o apop ${OPT_FLAGS} apop.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o apop ${POCHOIR_DEBUG_FLAGS} apop.cpp

life : tb_life.cpp
#   Phase-II compilation
	${CC} -o life ${OPT_FLAGS} tb_life.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o life ${POCHOIR_DEBUG_FLAGS} tb_life.cpp

psa_struct : tb_psa_struct.cpp
#   Phase-II compilation
	${CC} -o psa_struct -split-c-pointer ${OPT_FLAGS} tb_psa_struct.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o psa_struct ${POCHOIR_DEBUG_FLAGS} tb_psa_struct.cpp

parallel_heat_1D_NP : tb_parallel_heat_1D_NP.cpp
#   Phase-II compilation
	${CC} -o parallel_heat_1D_NP -split-c-pointer ${OPT_FLAGS} tb_parallel_heat_1D_NP.cpp
#	${ICC} -o parallel_heat_1D_NP_gdb ${POCHOIR_DEBUG_FLAGS} tb_parallel_heat_1D_NP_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o parallel_heat_1D_NP ${POCHOIR_DEBUG_FLAGS} tb_parallel_heat_1D_NP.cpp

reaction_diffusion_1D : tb_reaction_diffusion_1D.cpp
#   Phase-II compilation
	${CC} -o reaction_diffusion_1D -split-c-pointer ${OPT_FLAGS} tb_reaction_diffusion_1D.cpp
#	${ICC} -o reaction_diffusion_1D ${POCHOIR_DEBUG_FLAGS} tb_reaction_diffusion_1D_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o reaction_diffusion_1D ${POCHOIR_DEBUG_FLAGS} tb_reaction_diffusion_1D.cpp

rna : tb_rna.cpp
#   Phase-II compilation
	${CC} -o rna ${OPT_FLAGS} tb_rna.cpp
	${ICC} -o rna_gdb ${POCHOIR_DEBUG_FLAGS} tb_rna_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o rna ${POCHOIR_DEBUG_FLAGS} tb_rna.cpp

lcs : tb_lcs.cpp
#   Phase-II compilation
	${CC} -o lcs ${OPT_FLAGS} tb_lcs.cpp
	${ICC} -o lcs_gdb ${POCHOIR_DEBUG_FLAGS} tb_lcs_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_lcs_pochoir ${POCHOIR_DEBUG_FLAGS} tb_lcs.cpp

berkeley3d7pt : tb_3d7pt.cpp
#   Phase-II compilation
	${CC} -o 3d7pt ${OPT_FLAGS} tb_3d7pt.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o 3d7pt ${POCHOIR_DEBUG_FLAGS} tb_3d7pt.cpp

berkeley3d27pt : tb_3d27pt.cpp
#   Phase-II compilation
	${CC} -o 3d27pt ${OPT_FLAGS} tb_3d27pt.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o 3d27pt ${POCHOIR_DEBUG_FLAGS} tb_3d27pt.cpp

3dfd : tb_3dfd.cpp
#   Phase-II compilation
	${CC} -o 3dfd ${OPT_FLAGS} tb_3dfd.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_3dfd_pochoir ${POCHOIR_DEBUG_FLAGS} tb_3dfd.cpp

3dfd_test : tb_3dfd_test.cpp
#   Phase-II compilation
	${ICC} -o 3dfd_test_pbl ${OPT_FLAGS} tb_3dfd_test.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o 3dfd_test ${POCHOIR_DEBUG_FLAGS} tb_3dfd_test.cpp


clean: 
	rm -f *.o *.i *_pochoir *_gdb *_pochoir.cpp *.out
