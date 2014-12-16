#!/bin/bash
TARGET = heat_1D_NP heat_2D_NP heat_2D_P 3dfd heat_3D_NP heat_4D_NP life rna psa_struct lcs apop 
CC = ../pochoir
ICC = icpc
#OPT_FLAGS = -O3 -DNDEBUG -Wall -Werror -funroll-loops -xHOST -fno-alias -fno-fnalias -fp-model precise -std=c++0x -lrt -lpfm -L./ -I$(POCHOIR_LIB_PATH)

#options without vectorization
#OPT_FLAGS = -O3 -DNDEBUG -Wall -Werror -funroll-loops -fno-alias -fno-fnalias -fp-model precise -std=c++0x -no-vec -lrt -lpfm -L./ -I$(POCHOIR_LIB_PATH)

#more options for lanka cluster
OPT_FLAGS = -O3 -DNDEBUG -Wall -Werror -funroll-loops -xHOST -fno-alias -fno-fnalias -fp-model precise -std=c++0x -march=corei7-avx -mtune=corei7x-avx -lrt -lpfm -L./ -I$(POCHOIR_LIB_PATH)
#use -vec-report for generating vectorization report

#O0 option for lanka
#OPT_FLAGS = -O0 -DNDEBUG -Wall -Werror -funroll-loops -xHOST -fno-alias -fno-fnalias -fp-model precise -std=c++0x -march=corei7-avx -mtune=corei7x-avx -lrt -lpfm -L./ -I$(POCHOIR_LIB_PATH)

#O0 option for cloud
#OPT_FLAGS = -O0 -DNDEBUG -Wall -Werror -funroll-loops -xHOST -fno-alias -fno-fnalias -fp-model precise -std=c++0x -lrt -lpfm -L./ -I$(POCHOIR_LIB_PATH)
#OPT_FLAGS = -O0 -g -Wall -Werror -std=c++0x -lrt -lpfm -L./ -I$(POCHOIR_LIB_PATH)
POCHOIR_DEBUG_FLAGS = -O0 -g3 -DDEBUG -debug -Wall -Werror -std=c++0x
#ICC_DEBUG_FLAGS = -O0 -g3 -DDEBUG -Wall -Werror -std=c++0x -I${INTEL_CILK_HEADER} -include cilk_stub.h  
ICC_DEBUG_FLAGS = -O3 -Wall -Werror -std=c++0x -I${INTEL_CILK_HEADER} -include cilk_stub.h  
#DEPENDS = ../src/stopwatch.c
DEPENDS = 

# make all targets specified
all : $(TARGET)

heat_1D_NP : tb_heat_1D_NP.cpp
#   Phase-II compilation
	${CC} -o heat_1D_NP ${OPT_FLAGS} tb_heat_1D_NP.cpp ${DEPENDS}
#	${ICC} -o heat_1D_NP_gdb ${ICC_DEBUG_FLAGS} tb_heat_1D_NP_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_heat_1D_NP_pochoir ${POCHOIR_DEBUG_FLAGS} tb_heat_1D_NP.cpp

fib : tb_fib.cpp
#   Phase-II compilation
	${CC} -o fib ${OPT_FLAGS} tb_fib.cpp ${DEPENDS}
#	Phase-I compilation with debugging aid
#	${CC} -o fib ${POCHOIR_DEBUG_FLAGS} tb_fib.cpp

heat_2D_NP : tb_heat_2D_NP.cpp
#   Phase-II compilation
	${CC} -o heat_2D_NP ${OPT_FLAGS} tb_heat_2D_NP.cpp ${DEPENDS}
#	Phase-I compilation with debugging aid
#	${CC} -o heat_2D_NP ${POCHOIR_DEBUG_FLAGS} tb_heat_2D_NP.cpp

heat_NP_zero : tb_heat_2D_NP_zero.cpp
#   Phase-II compilation
	${CC} -o heat_2D_NP_zero ${OPT_FLAGS} tb_heat_2D_NP_zero.cpp ${DEPENDS}
#	Phase-I compilation with debugging aid
#	${CC} -o heat_2D_NP_zero ${POCHOIR_DEBUG_FLAGS} tb_heat_2D_NP_zero.cpp

2d_diffusion_loop : 2d_diffusion_loop.cpp
#   Phase-II compilation
	${CC} -o 2d_diffusion_loop ${OPT_FLAGS} 2d_diffusion_loop.cpp ${DEPENDS}
#	${ICC} -o heat_2D_P_gdb ${ICC_DEBUG_FLAGS} tb_heat_2D_P_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o heat_2D_P ${POCHOIR_DEBUG_FLAGS} tb_heat_2D_P.cpp

2d_diffusion : 2d_diffusion.cpp
#   Phase-II compilation
	${CC} -o 2d_diffusion ${OPT_FLAGS} 2d_diffusion.cpp ${DEPENDS}
#	${ICC} -o heat_2D_P_gdb ${ICC_DEBUG_FLAGS} tb_heat_2D_P_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o heat_2D_P ${POCHOIR_DEBUG_FLAGS} tb_heat_2D_P.cpp

heat_2D_P : tb_heat_2D_P.cpp
#   Phase-II compilation
	${CC} -o heat_2D_P ${OPT_FLAGS} tb_heat_2D_P.cpp
#	${ICC} -o heat_2D_P_gdb ${ICC_DEBUG_FLAGS} tb_heat_2D_P_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o heat_2D_P ${POCHOIR_DEBUG_FLAGS} tb_heat_2D_P.cpp

heat_3D_NP : tb_heat_3D_NP.cpp
#   Phase-II compilation
	${CC} -o heat_3D_NP ${OPT_FLAGS} tb_heat_3D_NP.cpp ${DEPENDS}
#	Phase-I compilation with debugging aid
#	${CC} -o heat_3D_NP ${POCHOIR_DEBUG_FLAGS} tb_heat_3D_NP.cpp

heat_4D_NP : tb_heat_4D_NP.cpp
#   Phase-II compilation
	${CC} -o heat_4D_NP ${OPT_FLAGS} tb_heat_4D_NP.cpp ${DEPENDS}
#	Phase-I compilation with debugging aid
#	${CC} -o heat_4D_NP ${POCHOIR_DEBUG_FLAGS} tb_heat_4D_NP.cpp

apop : apop.cpp
#   Phase-II compilation
	${CC} -o apop ${OPT_FLAGS} apop.cpp ${DEPENDS}
#	Phase-I compilation with debugging aid
#	${CC} -o apop ${POCHOIR_DEBUG_FLAGS} apop.cpp

life : tb_life.cpp
#   Phase-II compilation
	${CC} -o life ${OPT_FLAGS} tb_life.cpp ${DEPENDS}
#	Phase-I compilation with debugging aid
#	${CC} -o life ${POCHOIR_DEBUG_FLAGS} tb_life.cpp

psa_struct : tb_psa_struct.cpp
#   Phase-II compilation
	${CC} -o psa_struct -split-c-pointer ${OPT_FLAGS} tb_psa_struct.cpp ${DEPENDS}
#	Phase-I compilation with debugging aid
#	${CC} -o psa_struct ${POCHOIR_DEBUG_FLAGS} tb_psa_struct.cpp

parallel_heat_1D_NP : tb_parallel_heat_1D_NP.cpp
#   Phase-II compilation
	${CC} -o parallel_heat_1D_NP -split-c-pointer ${OPT_FLAGS} tb_parallel_heat_1D_NP.cpp ${DEPENDS}
#	${ICC} -o parallel_heat_1D_NP_gdb ${POCHOIR_DEBUG_FLAGS} tb_parallel_heat_1D_NP_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o parallel_heat_1D_NP ${POCHOIR_DEBUG_FLAGS} tb_parallel_heat_1D_NP.cpp

reaction_diffusion_1D : tb_reaction_diffusion_1D.cpp
#   Phase-II compilation
	${CC} -o reaction_diffusion_1D -split-c-pointer ${OPT_FLAGS} tb_reaction_diffusion_1D.cpp ${DEPENDS}
#	${ICC} -o reaction_diffusion_1D ${POCHOIR_DEBUG_FLAGS} tb_reaction_diffusion_1D_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o reaction_diffusion_1D ${POCHOIR_DEBUG_FLAGS} tb_reaction_diffusion_1D.cpp

rna : tb_rna.cpp
#   Phase-II compilation
	${CC} -o rna ${OPT_FLAGS} tb_rna.cpp ${DEPENDS}
#	${ICC} -o rna_gdb ${POCHOIR_DEBUG_FLAGS} tb_rna_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o rna ${POCHOIR_DEBUG_FLAGS} tb_rna.cpp

lcs : tb_lcs.cpp
#   Phase-II compilation
	${CC} -o lcs ${OPT_FLAGS} tb_lcs.cpp ${DEPENDS}
#	${ICC} -o lcs_gdb ${POCHOIR_DEBUG_FLAGS} tb_lcs_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_lcs_pochoir ${POCHOIR_DEBUG_FLAGS} tb_lcs.cpp

berkeley3d7pt : tb_3d7pt.cpp
#   Phase-II compilation
	${CC} -o 3d7pt ${OPT_FLAGS} tb_3d7pt.cpp ${DEPENDS}
#	Phase-I compilation with debugging aid
#	${CC} -o 3d7pt ${POCHOIR_DEBUG_FLAGS} tb_3d7pt.cpp

berkeley3d27pt : tb_3d27pt.cpp
#   Phase-II compilation
	${CC} -o 3d27pt ${OPT_FLAGS} tb_3d27pt.cpp ${DEPENDS}
#	Phase-I compilation with debugging aid
#	${CC} -o 3d27pt ${POCHOIR_DEBUG_FLAGS} tb_3d27pt.cpp

3dfd : tb_3dfd.cpp
#   Phase-II compilation
	${CC} -o 3dfd ${OPT_FLAGS} tb_3dfd.cpp ${DEPENDS}
#	Phase-I compilation with debugging aid
#	${CC} -o tb_3dfd_pochoir ${POCHOIR_DEBUG_FLAGS} tb_3dfd.cpp

3dfd_test : tb_3dfd_test.cpp
#   Phase-II compilation
	${ICC} -o 3dfd_test_pbl ${OPT_FLAGS} tb_3dfd_test.cpp ${DEPENDS}
#	Phase-I compilation with debugging aid
#	${CC} -o 3dfd_test ${POCHOIR_DEBUG_FLAGS} tb_3dfd_test.cpp


clean: 
	rm -f $(TARGET) *.o *.i *_pochoir *_gdb *_pochoir.cpp *.out
