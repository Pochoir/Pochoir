#/bin/bash
CC = ./pochoir
ICC = icpc
OPT_FLAGS = -O3 -DNDEBUG -Wall -Werror -unroll-aggressive -funroll-loops -ipo -xHOST -fno-alias -fno-fnalias -fp-model precise -std=c++0x
POCHOIR_DEBUG_FLAGS = -O0 -g3 -DDEBUG -debug -Wall -Werror -std=c++0x -I${INTEL_CILK_HEADER} -include cilk_stub.h
ICC_DEBUG_FLAGS = -O0 -g3 -DDEBUG -Wall -Werror -std=c++0x
heat_1D_NP : tb_heat_1D_NP.cpp
#   Phase-II compilation
	${CC} -o tb_heat_1D_NP_pochoir ${OPT_FLAGS} tb_heat_1D_NP.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_heat_1D_NP_pochoir ${POCHOIR_DEBUG_FLAGS} tb_heat_1D_NP.cpp

heat_NP : tb_heat_2D_NP.cpp
#   Phase-II compilation
	${CC} -o tb_heat_2D_NP_pochoir ${OPT_FLAGS} tb_heat_2D_NP.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_heat_2D_NP_pochoir ${POCHOIR_DEBUG_FLAGS} tb_heat_2D_NP.cpp

heat_NP_zero : tb_heat_2D_NP_zero.cpp
#   Phase-II compilation
#	${CC} -o tb_heat_2D_NP_zero_pochoir ${OPT_FLAGS} tb_heat_2D_NP_zero.cpp
#	Phase-I compilation with debugging aid
	${CC} -o tb_heat_2D_NP_zero_pochoir ${POCHOIR_DEBUG_FLAGS} tb_heat_2D_NP_zero.cpp

heat_P : tb_heat_2D_P.cpp
#   Phase-II compilation
#	${CC} -o tb_heat_2D_P_pochoir ${OPT_FLAGS} tb_heat_2D_P.cpp
#	Phase-I compilation with debugging aid
	${CC} -o tb_heat_2D_P_pochoir ${POCHOIR_DEBUG_FLAGS} tb_heat_2D_P.cpp

heat_3D_NP : tb_heat_3D_NP.cpp
#   Phase-II compilation
	${CC} -o tb_heat_3D_NP_pochoir ${OPT_FLAGS} tb_heat_3D_NP.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o tb_heat_3D_NP_pochoir ${POCHOIR_DEBUG_FLAGS} tb_heat_3D_NP.cpp

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
