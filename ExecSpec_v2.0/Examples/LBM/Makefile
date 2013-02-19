SHELL = /bin/bash

CC = ./pochoir
ICC = icc

CFLAGS = -restrict -unroll-aggressive -funroll-loops -scalar-rep -ipo -xHOST -fno-alias -O3 -no-prec-div -std=c++0x
# CFLAGS = -restrict -no-prec-div 
#debugging
#POFLAGS = -debug
DEBUGFLAGS = -g3 -restrict -O0 -DDEBUG -no-prec-div -std=c++0x -debug
POFLAGS = -O3 -unroll-aggressive -funroll-loops -scalar-rep -ipo -xHOST -fno-alias -restrict -std=c++0x
REPORT_FLAGS = -par-report3 -openmp-report1 -vec-report3
MY_OPTS = -DMYOPT

DEFAULT_OPTIONS = $(CFLAGS) $(REPORT_FLAGS) $(MY_OPTS)
OMP_OPTIONS = -D_OPENMP -openmp
NOVEC_OPTIONS = -no-vec
VEC_OPTIONS = -DVECTORIZE
SOA_OPTIONS = -DSOA

INCLUDE := lbm_1d_array.h config.h lbm.h main.h 
INCLUDE_TANG := lbm_1d_array_tang.h config.h lbm_tang.h main_tang.h 
ORIG_SRC := main.c lbm.c
CILK_CO_SRC := main_co.cpp lbm_co.c
PO_SRC := main_po.cpp lbm_po.c
XO_SRC := main_xo.cpp lbm_xo.c
POCHOIR_SRC := main_tang.cpp lbm_tang.cpp 
POCHOIR_SRC_NAIVE := main_tang_naive.cpp lbm_tang_naive.cpp 

lbm_orig: $(INCLUDE) $(ORIG_SRC)
	$(ICC) $(DEFAULT_OPTIONS) $(NOVEC_OPTIONS) -o BIN/$@ $(ORIG_SRC) 

lbm_orig_soa_novec.exe: $(INCLUDE) $(ORIG_SRC)
	$(ICC) $(DEFAULT_OPTIONS) $(SOA_OPTIONS) $(NOVEC_OPTIONS) -o BIN/$@ $(ORIG_SRC) >& COMPILE_LOG/$@.report

lbm_orig_vec.exe: $(INCLUDE) $(ORIG_SRC)
	$(ICC) $(DEFAULT_OPTIONS) $(VEC_OPTIONS) -o BIN/$@ $(ORIG_SRC) >& COMPILE_LOG/$@.report

lbm_orig_soa_vec.exe: $(INCLUDE) $(ORIG_SRC)
	$(ICC) $(DEFAULT_OPTIONS) $(SOA_OPTIONS) $(VEC_OPTIONS) -o BIN/$@ $(ORIG_SRC) >& COMPILE_LOG/$@.report

lbm_omp_novec.exe: $(INCLUDE) $(ORIG_SRC)
	$(ICC) $(DEFAULT_OPTIONS) $(OMP_OPTIONS) $(NOVEC_OPTIONS) -o BIN/$@ $(ORIG_SRC) >& COMPILE_LOG/$@.report

lbm_omp_soa_novec.exe: $(INCLUDE) $(ORIG_SRC)
	$(ICC) $(DEFAULT_OPTIONS) $(SOA_OPTIONS) $(OMP_OPTIONS) $(NOVEC_OPTIONS) -o BIN/$@ $(ORIG_SRC) >& COMPILE_LOG/$@.report

lbm_omp_vec.exe: $(INCLUDE) $(ORIG_SRC)
	$(ICC) $(DEFAULT_OPTIONS) $(OMP_OPTIONS) $(VEC_OPTIONS) -o BIN/$@ $(ORIG_SRC) >& COMPILE_LOG/$@.report

lbm_omp_soa_vec.exe: $(INCLUDE) $(ORIG_SRC)
	$(ICC) $(DEFAULT_OPTIONS) $(SOA_OPTIONS) $(OMP_OPTIONS) $(VEC_OPTIONS) -o BIN/$@ $(ORIG_SRC) >& COMPILE_LOG/$@.report

lbm_cilk_co_novec.exe: $(INCLUDE) $(CILK_CO_SRC)
#	$(ICC) $(DEBUGFLAGS) $(NOVEC_OPTIONS) -I${CILK_HEADER_PATH} -o BIN/$@ $(CILK_CO_SRC) 
	$(CC) $(DEFAULT_OPTIONS) $(NOVEC_OPTIONS) -I${CILK_HEADER_PATH} -o BIN/$@ $(CILK_CO_SRC) 
#	$(ICC) $(DEFAULT_OPTIONS) $(NOVEC_OPTIONS) -I${CILK_HEADER_PATH} -o BIN/$@ $(CILK_CO_SRC) >& COMPILE_LOG/$@.report

lbm_cilk_co_soa_novec.exe: $(INCLUDE) $(CILK_CO_SRC)
	$(ICC) $(DEFAULT_OPTIONS) $(SOA_OPTIONS) $(NOVEC_OPTIONS) -I${CILK_HEADER_PATH} -o BIN/$@ $(CILK_CO_SRC) 
#	$(ICC) $(DEFAULT_OPTIONS) $(SOA_OPTIONS) $(NOVEC_OPTIONS) -I${CILK_HEADER_PATH} -o BIN/$@ $(CILK_CO_SRC) >& COMPILE_LOG/$@.report

lbm_cilk_co_vec.exe: $(INCLUDE) $(CILK_CO_SRC)
	$(ICC) $(DEFAULT_OPTIONS) $(VEC_OPTIONS) -I${CILK_HEADER_PATH} -o BIN/$@ $(CILK_CO_SRC) 
#	$(ICC) $(DEFAULT_OPTIONS) $(VEC_OPTIONS) -I${CILK_HEADER_PATH} -o BIN/$@ $(CILK_CO_SRC) >& COMPILE_LOG/$@.report

lbm_cilk_co_soa_vec.exe: $(INCLUDE) $(CILK_CO_SRC)
	$(ICC) $(DEFAULT_OPTIONS) $(SOA_OPTIONS) $(VEC_OPTIONS) -o BIN/$@ $(CILK_CO_SRC) >& COMPILE_LOG/$@.report

lbm_tbb_co_novec.exe: $(INCLUDE) $(TBB_CO_SRC)
	$(ICC) $(DEFAULT_OPTIONS) $(TBB_OPTIONS) $(NOVEC_OPTIONS) -o BIN/$@ $(TBB_CO_SRC) >& COMPILE_LOG/$@.report

lbm_tbb_co_soa_novec.exe: $(INCLUDE) $(TBB_CO_SRC)
	$(ICC) $(DEFAULT_OPTIONS) $(SOA_OPTIONS) $(TBB_OPTIONS) $(NOVEC_OPTIONS) -o BIN/$@ $(TBB_CO_SRC) >& COMPILE_LOG/$@.report

lbm_tbb_co_vec.exe: $(INCLUDE) $(TBB_CO_SRC)
	$(ICC) $(DEFAULT_OPTIONS) $(TBB_OPTIONS) $(VEC_OPTIONS) -o BIN/$@ $(TBB_CO_SRC) >& COMPILE_LOG/$@.report

lbm_tbb_co_soa_vec.exe: $(INCLUDE) $(TBB_CO_SRC)
	$(ICC) $(DEFAULT_OPTIONS) $(SOA_OPTIONS) $(TBB_OPTIONS) $(VEC_OPTIONS) -o BIN/$@ $(TBB_CO_SRC) >& COMPILE_LOG/$@.report

lbm_po_novec.exe: $(INCLUDE) $(PO_SRC)
	$(ICC) $(DEFAULT_OPTIONS) $(NOVEC_OPTIONS) -o BIN/$@ $(PO_SRC) >& COMPILE_LOG/$@.report

lbm_xo_novec.exe: $(INCLUDE) $(XO_SRC)
	./pochoir $(POFLAGS) -c lbm_xo.cpp -std=c++0x
	$(ICC) $(DEFAULT_OPTIONS) $(NOVEC_OPTIONS) -DNDEBUG -std=c++0x -Wall -Werror -o BIN/$@ main_xo.cpp lbm_xo_pochoir.cpp >& COMPILE_LOG/$@.report

lbm_tang: $(INCLUDE_TANG) $(POCHOIR_SRC)
	${CC} ${POFLAGS} ${DEFAULT_OPTIONS} ${NOVEC_OPTIONS} -DNDEBUG -Wall -Werror -o BIN/$@ ${POCHOIR_SRC}
#	${CC} ${DEBUGFLAGS} ${NOVEC_OPTIONS} -Wall -Werror -o BIN/$@ ${POCHOIR_SRC}
#	${ICC} ${POFLAGS} ${DEFAULT_OPTIONS} ${NOVEC_OPTIONS} -I${POCHOIR_LIB_PATH} -DNDEBUG -Wall -Werror -o BIN/$@ ${POCHOIR_SRC}
#	${ICC} ${DEBUGFLAGS} ${NOVEC_OPTIONS} -I${POCHOIR_LIB_PATH} -DNDEBUG -Wall -Werror -o BIN/$@ ${POCHOIR_SRC}

lbm_tang_naive: $(INCLUDE) $(POCHOIR_SRC_NAIVE)
	${CC} ${POFLAGS} ${DEFAULT_OPTIONS} ${NOVEC_OPTIONS} -DNDEBUG -Wall -Werror -o BIN/$@ ${POCHOIR_SRC_NAIVE}
clean:
	-rm *.o BIN/* COMPILE_LOG/* OUT/*

## You need to provide at least 3 more targets (xtune_clean, xtune_build, xtune_test) for using XTune
xtune_clean: clean

xtune_build: lbm_cilk_co_vec.exe 

xtune_test: 
	./run_ref.csh lbm_cilk_co_vec.exe 8 $@.out
#	./run_train.csh lbm_cilk_co_vec.exe 8 $@.out
