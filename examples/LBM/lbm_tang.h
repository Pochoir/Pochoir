/* $Id: lbm.h,v 1.1 2004/04/20 14:33:59 pohlt Exp $ */

/*############################################################################*/

#ifndef _LBM_H_
#define _LBM_H_

/*############################################################################*/

#include <pochoir.hpp>
#include "config.h"

/*############################################################################*/

typedef struct
{
    double _C;
    double _N;
    double _S;
    double _E;
    double _W;
    double _T;
    double _B;
    double _NE;
    double _NW;
    double _SE;
    double _SW;
    double _NT;
    double _NB;
    double _ST;
    double _SB;
    double _ET;
    double _EB;
    double _WT;
    double _WB;
    unsigned int _FLAGS;
    // double _FLAGS;
} PoCellEntry;

/*############################################################################*/

#define DFL1 (1.0/ 3.0)
#define DFL2 (1.0/18.0)
#define DFL3 (1.0/36.0)

/*############################################################################*/

typedef enum {C = 0,
              N, S, E, W, T, B,
              NE, NW, SE, SW,
              NT, NB, ST, SB,
              ET, EB, WT, WB,
              FLAGS, N_CELL_ENTRIES} CELL_ENTRIES;
#define N_DISTR_FUNCS FLAGS

typedef enum {OBSTACLE    = 1 << 0,
              ACCEL       = 1 << 1,
              IN_OUT_FLOW = 1 << 2} CELL_FLAGS;

#define MARGIN_Z 2

#include "lbm_1d_array_tang.h"

/*############################################################################*/

void LBM_initializeGrid( Pochoir_Array_3D(PoCellEntry) & pa, int t );
void LBM_initializeSpecialCellsForLDC( Pochoir_Array_3D(PoCellEntry) & pa, int t );
void LBM_initializeSpecialCellsForChannel( Pochoir_Array_3D(PoCellEntry) & pa, int t );
void LBM_loadRandomObstacle( Pochoir_Array_3D(PoCellEntry) & pa, int t );
void LBM_loadObstacleFile( Pochoir_Array_3D(PoCellEntry) & pa, int t, const char* filename );
void LBM_showGridStatistics( Pochoir_Array_3D(PoCellEntry) & pa, int t );
void LBM_handleInOutFlow( Pochoir_Array_3D(PoCellEntry) & pa, int t, int z, int y, int x );
void LBM_performStreamCollide( Pochoir_Array_3D(PoCellEntry) & pa, int t, int z, int y, int x );
void LBM_storeVelocityField( Pochoir_Array_3D(PoCellEntry) & pa, int t,
                             const char* filename, const BOOL binary );
void LBM_compareVelocityField( Pochoir_Array_3D(PoCellEntry) & pa, int t,
                               const char* filename, const BOOL binary );

/*############################################################################*/

#define LBM_performStreamCollide_macro( pa, t, z, y, x ) do { \
	double ux, uy, uz, u2, rho; \
\
/*	if( TEST_FLAG_SWEEP( pa, t-1, z, y, x, OBSTACLE )) { \
		DST_C ( pa, t, z, y, x ) = SRC_C ( pa, t-1, z, y, x ); \
		DST_S ( pa, t, z, y, x ) = SRC_N ( pa, t-1, z, y, x ); \
		DST_N ( pa, t, z, y, x ) = SRC_S ( pa, t-1, z, y, x ); \
		DST_W ( pa, t, z, y, x ) = SRC_E ( pa, t-1, z, y, x ); \
		DST_E ( pa, t, z, y, x ) = SRC_W ( pa, t-1, z, y, x ); \
		DST_B ( pa, t, z, y, x ) = SRC_T ( pa, t-1, z, y, x ); \
		DST_T ( pa, t, z, y, x ) = SRC_B ( pa, t-1, z, y, x ); \
		DST_SW( pa, t, z, y, x ) = SRC_NE( pa, t-1, z, y, x ); \
		DST_SE( pa, t, z, y, x ) = SRC_NW( pa, t-1, z, y, x ); \
		DST_NW( pa, t, z, y, x ) = SRC_SE( pa, t-1, z, y, x ); \
		DST_NE( pa, t, z, y, x ) = SRC_SW( pa, t-1, z, y, x ); \
		DST_SB( pa, t, z, y, x ) = SRC_NT( pa, t-1, z, y, x ); \
		DST_ST( pa, t, z, y, x ) = SRC_NB( pa, t-1, z, y, x ); \
		DST_NB( pa, t, z, y, x ) = SRC_ST( pa, t-1, z, y, x ); \
		DST_NT( pa, t, z, y, x ) = SRC_SB( pa, t-1, z, y, x ); \
		DST_WB( pa, t, z, y, x ) = SRC_ET( pa, t-1, z, y, x ); \
		DST_WT( pa, t, z, y, x ) = SRC_EB( pa, t-1, z, y, x ); \
		DST_EB( pa, t, z, y, x ) = SRC_WT( pa, t-1, z, y, x ); \
		DST_ET( pa, t, z, y, x ) = SRC_WB( pa, t-1, z, y, x ); \
        // return;                                           \
		continue; \
	} \*/ \
\
	rho = + SRC_C ( pa, t-1, z, y, x ) + SRC_N ( pa, t-1, z, y, x ) \
	      + SRC_S ( pa, t-1, z, y, x ) + SRC_E ( pa, t-1, z, y, x ) \
	      + SRC_W ( pa, t-1, z, y, x ) + SRC_T ( pa, t-1, z, y, x ) \
	      + SRC_B ( pa, t-1, z, y, x ) + SRC_NE( pa, t-1, z, y, x ) \
	      + SRC_NW( pa, t-1, z, y, x ) + SRC_SE( pa, t-1, z, y, x ) \
	      + SRC_SW( pa, t-1, z, y, x ) + SRC_NT( pa, t-1, z, y, x ) \
	      + SRC_NB( pa, t-1, z, y, x ) + SRC_ST( pa, t-1, z, y, x ) \
	      + SRC_SB( pa, t-1, z, y, x ) + SRC_ET( pa, t-1, z, y, x ) \
	      + SRC_EB( pa, t-1, z, y, x ) + SRC_WT( pa, t-1, z, y, x ) \
	      + SRC_WB( pa, t-1, z, y, x ); \
\
	ux  = + SRC_E ( pa, t-1, z, y, x ) - SRC_W ( pa, t-1, z, y, x ) \
	      + SRC_NE( pa, t-1, z, y, x ) - SRC_NW( pa, t-1, z, y, x ) \
	      + SRC_SE( pa, t-1, z, y, x ) - SRC_SW( pa, t-1, z, y, x ) \
	      + SRC_ET( pa, t-1, z, y, x ) + SRC_EB( pa, t-1, z, y, x ) \
	      - SRC_WT( pa, t-1, z, y, x ) - SRC_WB( pa, t-1, z, y, x );\
\
	uy  = + SRC_N ( pa, t-1, z, y, x ) - SRC_S ( pa, t-1, z, y, x ) \
	      + SRC_NE( pa, t-1, z, y, x ) + SRC_NW( pa, t-1, z, y, x ) \
	      - SRC_SE( pa, t-1, z, y, x ) - SRC_SW( pa, t-1, z, y, x ) \
	      + SRC_NT( pa, t-1, z, y, x ) + SRC_NB( pa, t-1, z, y, x ) \
	      - SRC_ST( pa, t-1, z, y, x ) - SRC_SB( pa, t-1, z, y, x );\
\
	uz  = + SRC_T ( pa, t-1, z, y, x ) - SRC_B ( pa, t-1, z, y, x ) \
	      + SRC_NT( pa, t-1, z, y, x ) - SRC_NB( pa, t-1, z, y, x ) \
	      + SRC_ST( pa, t-1, z, y, x ) - SRC_SB( pa, t-1, z, y, x ) \
	      + SRC_ET( pa, t-1, z, y, x ) - SRC_EB( pa, t-1, z, y, x ) \
	      + SRC_WT( pa, t-1, z, y, x ) - SRC_WB( pa, t-1, z, y, x );\
\
	ux /= rho;\
	uy /= rho;\
	uz /= rho;\
\
	if( TEST_FLAG_SWEEP( pa, t-1, z, y, x, ACCEL )) {\
		ux = 0.005;\
		uy = 0.002;\
		uz = 0.000;\
	}\
\
	u2 = 1.5 * (ux*ux + uy*uy + uz*uz);\
	DST_C ( pa, t, z, y, x ) = (1.0-OMEGA)*SRC_C ( pa, t-1, z, y, x ) + DFL1*OMEGA*rho*(1.0                                 - u2);\
\
	DST_N ( pa, t, z, y, x ) = (1.0-OMEGA)*SRC_N ( pa, t-1, z, y, x ) + DFL2*OMEGA*rho*(1.0 +       uy*(4.5*uy       + 3.0) - u2);\
	DST_S ( pa, t, z, y, x ) = (1.0-OMEGA)*SRC_S ( pa, t-1, z, y, x ) + DFL2*OMEGA*rho*(1.0 +       uy*(4.5*uy       - 3.0) - u2);\
	DST_E ( pa, t, z, y, x ) = (1.0-OMEGA)*SRC_E ( pa, t-1, z, y, x ) + DFL2*OMEGA*rho*(1.0 +       ux*(4.5*ux       + 3.0) - u2);\
	DST_W ( pa, t, z, y, x ) = (1.0-OMEGA)*SRC_W ( pa, t-1, z, y, x ) + DFL2*OMEGA*rho*(1.0 +       ux*(4.5*ux       - 3.0) - u2);\
	DST_T ( pa, t, z, y, x ) = (1.0-OMEGA)*SRC_T ( pa, t-1, z, y, x ) + DFL2*OMEGA*rho*(1.0 +       uz*(4.5*uz       + 3.0) - u2);\
	DST_B ( pa, t, z, y, x ) = (1.0-OMEGA)*SRC_B ( pa, t-1, z, y, x ) + DFL2*OMEGA*rho*(1.0 +       uz*(4.5*uz       - 3.0) - u2);\
\
	DST_NE( pa, t, z, y, x ) = (1.0-OMEGA)*SRC_NE( pa, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (+ux+uy)*(4.5*(+ux+uy) + 3.0) - u2);\
	DST_NW( pa, t, z, y, x ) = (1.0-OMEGA)*SRC_NW( pa, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (-ux+uy)*(4.5*(-ux+uy) + 3.0) - u2);\
	DST_SE( pa, t, z, y, x ) = (1.0-OMEGA)*SRC_SE( pa, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (+ux-uy)*(4.5*(+ux-uy) + 3.0) - u2);\
	DST_SW( pa, t, z, y, x ) = (1.0-OMEGA)*SRC_SW( pa, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (-ux-uy)*(4.5*(-ux-uy) + 3.0) - u2);\
	DST_NT( pa, t, z, y, x ) = (1.0-OMEGA)*SRC_NT( pa, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (+uy+uz)*(4.5*(+uy+uz) + 3.0) - u2);\
	DST_NB( pa, t, z, y, x ) = (1.0-OMEGA)*SRC_NB( pa, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (+uy-uz)*(4.5*(+uy-uz) + 3.0) - u2);\
	DST_ST( pa, t, z, y, x ) = (1.0-OMEGA)*SRC_ST( pa, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (-uy+uz)*(4.5*(-uy+uz) + 3.0) - u2);\
	DST_SB( pa, t, z, y, x ) = (1.0-OMEGA)*SRC_SB( pa, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (-uy-uz)*(4.5*(-uy-uz) + 3.0) - u2);\
	DST_ET( pa, t, z, y, x ) = (1.0-OMEGA)*SRC_ET( pa, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (+ux+uz)*(4.5*(+ux+uz) + 3.0) - u2);\
	DST_EB( pa, t, z, y, x ) = (1.0-OMEGA)*SRC_EB( pa, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (+ux-uz)*(4.5*(+ux-uz) + 3.0) - u2);\
	DST_WT( pa, t, z, y, x ) = (1.0-OMEGA)*SRC_WT( pa, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (-ux+uz)*(4.5*(-ux+uz) + 3.0) - u2);\
	DST_WB( pa, t, z, y, x ) = (1.0-OMEGA)*SRC_WB( pa, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (-ux-uz)*(4.5*(-ux-uz) + 3.0) - u2);\
} while (0)

/*############################################################################*/

#define LBM_handleInOutFlow_macro( pa, t, z, y, x ) do {\
	double ux , uy , uz , rho ,\
	       ux1, uy1, uz1, rho1,\
	       ux2, uy2, uz2, rho2,\
	       u2, px, py;\
	/* inflow */\
    if (z == 0 + MARGIN_Z) {\
		rho1 = + GRID_ENTRY( pa, t-1, z+1, y, x, C  ) + GRID_ENTRY( pa, t-1, z+1, y, x, N  )\
		       + GRID_ENTRY( pa, t-1, z+1, y, x, S  ) + GRID_ENTRY( pa, t-1, z+1, y, x, E  )\
		       + GRID_ENTRY( pa, t-1, z+1, y, x, W  ) + GRID_ENTRY( pa, t-1, z+1, y, x, T  )\
		       + GRID_ENTRY( pa, t-1, z+1, y, x, B  ) + GRID_ENTRY( pa, t-1, z+1, y, x, NE )\
		       + GRID_ENTRY( pa, t-1, z+1, y, x, NW ) + GRID_ENTRY( pa, t-1, z+1, y, x, SE )\
		       + GRID_ENTRY( pa, t-1, z+1, y, x, SW ) + GRID_ENTRY( pa, t-1, z+1, y, x, NT )\
		       + GRID_ENTRY( pa, t-1, z+1, y, x, NB ) + GRID_ENTRY( pa, t-1, z+1, y, x, ST )\
		       + GRID_ENTRY( pa, t-1, z+1, y, x, SB ) + GRID_ENTRY( pa, t-1, z+1, y, x, ET )\
		       + GRID_ENTRY( pa, t-1, z+1, y, x, EB ) + GRID_ENTRY( pa, t-1, z+1, y, x, WT )\
		       + GRID_ENTRY( pa, t-1, z+1, y, x, WB );                                      \
		rho2 = + GRID_ENTRY( pa, t-1, z+2, y, x, C  ) + GRID_ENTRY( pa, t-1, z+2, y, x, N  )\
		       + GRID_ENTRY( pa, t-1, z+2, y, x, S  ) + GRID_ENTRY( pa, t-1, z+2, y, x, E  )\
		       + GRID_ENTRY( pa, t-1, z+2, y, x, W  ) + GRID_ENTRY( pa, t-1, z+2, y, x, T  )\
		       + GRID_ENTRY( pa, t-1, z+2, y, x, B  ) + GRID_ENTRY( pa, t-1, z+2, y, x, NE )\
		       + GRID_ENTRY( pa, t-1, z+2, y, x, NW ) + GRID_ENTRY( pa, t-1, z+2, y, x, SE )\
		       + GRID_ENTRY( pa, t-1, z+2, y, x, SW ) + GRID_ENTRY( pa, t-1, z+2, y, x, NT )\
		       + GRID_ENTRY( pa, t-1, z+2, y, x, NB ) + GRID_ENTRY( pa, t-1, z+2, y, x, ST )\
		       + GRID_ENTRY( pa, t-1, z+2, y, x, SB ) + GRID_ENTRY( pa, t-1, z+2, y, x, ET )\
		       + GRID_ENTRY( pa, t-1, z+2, y, x, EB ) + GRID_ENTRY( pa, t-1, z+2, y, x, WT )\
		       + GRID_ENTRY( pa, t-1, z+2, y, x, WB );\
\
		rho = 2.0*rho1 - rho2;\
\
		px = (x / (0.5*(SIZE_X-1))) - 1.0;\
		py = (y / (0.5*(SIZE_Y-1))) - 1.0;\
		ux = 0.00;\
		uy = 0.00;\
		uz = 0.01 * (1.0-px*px) * (1.0-py*py);\
\
		u2 = 1.5 * (ux*ux + uy*uy + uz*uz);\
\
		LOCAL( pa, t-1, z, y, x, C ) = DFL1*rho*(1.0                                 - u2);\
\
		LOCAL( pa, t-1, z, y, x, N ) = DFL2*rho*(1.0 +       uy*(4.5*uy       + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, S ) = DFL2*rho*(1.0 +       uy*(4.5*uy       - 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, E ) = DFL2*rho*(1.0 +       ux*(4.5*ux       + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, W ) = DFL2*rho*(1.0 +       ux*(4.5*ux       - 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, T ) = DFL2*rho*(1.0 +       uz*(4.5*uz       + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, B ) = DFL2*rho*(1.0 +       uz*(4.5*uz       - 3.0) - u2);\
\
		LOCAL( pa, t-1, z, y, x, NE) = DFL3*rho*(1.0 + (+ux+uy)*(4.5*(+ux+uy) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, NW) = DFL3*rho*(1.0 + (-ux+uy)*(4.5*(-ux+uy) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, SE) = DFL3*rho*(1.0 + (+ux-uy)*(4.5*(+ux-uy) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, SW) = DFL3*rho*(1.0 + (-ux-uy)*(4.5*(-ux-uy) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, NT) = DFL3*rho*(1.0 + (+uy+uz)*(4.5*(+uy+uz) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, NB) = DFL3*rho*(1.0 + (+uy-uz)*(4.5*(+uy-uz) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, ST) = DFL3*rho*(1.0 + (-uy+uz)*(4.5*(-uy+uz) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, SB) = DFL3*rho*(1.0 + (-uy-uz)*(4.5*(-uy-uz) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, ET) = DFL3*rho*(1.0 + (+ux+uz)*(4.5*(+ux+uz) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, EB) = DFL3*rho*(1.0 + (+ux-uz)*(4.5*(+ux-uz) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, WT) = DFL3*rho*(1.0 + (-ux+uz)*(4.5*(-ux+uz) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, WB) = DFL3*rho*(1.0 + (-ux-uz)*(4.5*(-ux-uz) + 3.0) - u2);\
    } /* end if (z == 0 + MARGIN_Z) */\
\
	/* outflow */\
    if (z == SIZE_Z - 1 + MARGIN_Z) {\
		rho1 = + GRID_ENTRY( pa, t-1, z-1, y, x, C  ) + GRID_ENTRY( pa, t-1, z-1, y, x, N  )\
		       + GRID_ENTRY( pa, t-1, z-1, y, x, S  ) + GRID_ENTRY( pa, t-1, z-1, y, x, E  )\
		       + GRID_ENTRY( pa, t-1, z-1, y, x, W  ) + GRID_ENTRY( pa, t-1, z-1, y, x, T  )\
		       + GRID_ENTRY( pa, t-1, z-1, y, x, B  ) + GRID_ENTRY( pa, t-1, z-1, y, x, NE )\
		       + GRID_ENTRY( pa, t-1, z-1, y, x, NW ) + GRID_ENTRY( pa, t-1, z-1, y, x, SE )\
		       + GRID_ENTRY( pa, t-1, z-1, y, x, SW ) + GRID_ENTRY( pa, t-1, z-1, y, x, NT )\
		       + GRID_ENTRY( pa, t-1, z-1, y, x, NB ) + GRID_ENTRY( pa, t-1, z-1, y, x, ST )\
		       + GRID_ENTRY( pa, t-1, z-1, y, x, SB ) + GRID_ENTRY( pa, t-1, z-1, y, x, ET )\
		       + GRID_ENTRY( pa, t-1, z-1, y, x, EB ) + GRID_ENTRY( pa, t-1, z-1, y, x, WT )\
		       + GRID_ENTRY( pa, t-1, z-1, y, x, WB );                                      \
		ux1  = + GRID_ENTRY( pa, t-1, z-1, y, x, E  ) - GRID_ENTRY( pa, t-1, z-1, y, x, W  )\
		       + GRID_ENTRY( pa, t-1, z-1, y, x, NE ) - GRID_ENTRY( pa, t-1, z-1, y, x, NW )\
		       + GRID_ENTRY( pa, t-1, z-1, y, x, SE ) - GRID_ENTRY( pa, t-1, z-1, y, x, SW )\
		       + GRID_ENTRY( pa, t-1, z-1, y, x, ET ) + GRID_ENTRY( pa, t-1, z-1, y, x, EB )\
		       - GRID_ENTRY( pa, t-1, z-1, y, x, WT ) - GRID_ENTRY( pa, t-1, z-1, y, x, WB );\
		uy1  = + GRID_ENTRY( pa, t-1, z-1, y, x, N  ) - GRID_ENTRY( pa, t-1, z-1, y, x, S  )\
		       + GRID_ENTRY( pa, t-1, z-1, y, x, NE ) + GRID_ENTRY( pa, t-1, z-1, y, x, NW )\
		       - GRID_ENTRY( pa, t-1, z-1, y, x, SE ) - GRID_ENTRY( pa, t-1, z-1, y, x, SW )\
		       + GRID_ENTRY( pa, t-1, z-1, y, x, NT ) + GRID_ENTRY( pa, t-1, z-1, y, x, NB )\
		       - GRID_ENTRY( pa, t-1, z-1, y, x, ST ) - GRID_ENTRY( pa, t-1, z-1, y, x, SB );\
		uz1  = + GRID_ENTRY( pa, t-1, z-1, y, x, T  ) - GRID_ENTRY( pa, t-1, z-1, y, x, B  )\
		       + GRID_ENTRY( pa, t-1, z-1, y, x, NT ) - GRID_ENTRY( pa, t-1, z-1, y, x, NB )\
		       + GRID_ENTRY( pa, t-1, z-1, y, x, ST ) - GRID_ENTRY( pa, t-1, z-1, y, x, SB )\
		       + GRID_ENTRY( pa, t-1, z-1, y, x, ET ) - GRID_ENTRY( pa, t-1, z-1, y, x, EB )\
		       + GRID_ENTRY( pa, t-1, z-1, y, x, WT ) - GRID_ENTRY( pa, t-1, z-1, y, x, WB );\
		ux1 /= rho1;\
		uy1 /= rho1;\
		uz1 /= rho1;\
\
		rho2 = + GRID_ENTRY( pa, t-1, z-2, y, x, C  ) + GRID_ENTRY( pa, t-1, z-2, y, x, N  )\
		       + GRID_ENTRY( pa, t-1, z-2, y, x, S  ) + GRID_ENTRY( pa, t-1, z-2, y, x, E  )\
		       + GRID_ENTRY( pa, t-1, z-2, y, x, W  ) + GRID_ENTRY( pa, t-1, z-2, y, x, T  )\
		       + GRID_ENTRY( pa, t-1, z-2, y, x, B  ) + GRID_ENTRY( pa, t-1, z-2, y, x, NE )\
		       + GRID_ENTRY( pa, t-1, z-2, y, x, NW ) + GRID_ENTRY( pa, t-1, z-2, y, x, SE )\
		       + GRID_ENTRY( pa, t-1, z-2, y, x, SW ) + GRID_ENTRY( pa, t-1, z-2, y, x, NT )\
		       + GRID_ENTRY( pa, t-1, z-2, y, x, NB ) + GRID_ENTRY( pa, t-1, z-2, y, x, ST )\
		       + GRID_ENTRY( pa, t-1, z-2, y, x, SB ) + GRID_ENTRY( pa, t-1, z-2, y, x, ET )\
		       + GRID_ENTRY( pa, t-1, z-2, y, x, EB ) + GRID_ENTRY( pa, t-1, z-2, y, x, WT )\
		       + GRID_ENTRY( pa, t-1, z-2, y, x, WB );                                      \
		ux2  = + GRID_ENTRY( pa, t-1, z-2, y, x, E  ) - GRID_ENTRY( pa, t-1, z-2, y, x, W  )\
		       + GRID_ENTRY( pa, t-1, z-2, y, x, NE ) - GRID_ENTRY( pa, t-1, z-2, y, x, NW )\
		       + GRID_ENTRY( pa, t-1, z-2, y, x, SE ) - GRID_ENTRY( pa, t-1, z-2, y, x, SW )\
		       + GRID_ENTRY( pa, t-1, z-2, y, x, ET ) + GRID_ENTRY( pa, t-1, z-2, y, x, EB )\
		       - GRID_ENTRY( pa, t-1, z-2, y, x, WT ) - GRID_ENTRY( pa, t-1, z-2, y, x, WB );\
		uy2  = + GRID_ENTRY( pa, t-1, z-2, y, x, N  ) - GRID_ENTRY( pa, t-1, z-2, y, x, S  )\
		       + GRID_ENTRY( pa, t-1, z-2, y, x, NE ) + GRID_ENTRY( pa, t-1, z-2, y, x, NW )\
		       - GRID_ENTRY( pa, t-1, z-2, y, x, SE ) - GRID_ENTRY( pa, t-1, z-2, y, x, SW )\
		       + GRID_ENTRY( pa, t-1, z-2, y, x, NT ) + GRID_ENTRY( pa, t-1, z-2, y, x, NB )\
		       - GRID_ENTRY( pa, t-1, z-2, y, x, ST ) - GRID_ENTRY( pa, t-1, z-2, y, x, SB );\
		uz2  = + GRID_ENTRY( pa, t-1, z-2, y, x, T  ) - GRID_ENTRY( pa, t-1, z-2, y, x, B  )\
		       + GRID_ENTRY( pa, t-1, z-2, y, x, NT ) - GRID_ENTRY( pa, t-1, z-2, y, x, NB )\
		       + GRID_ENTRY( pa, t-1, z-2, y, x, ST ) - GRID_ENTRY( pa, t-1, z-2, y, x, SB )\
		       + GRID_ENTRY( pa, t-1, z-2, y, x, ET ) - GRID_ENTRY( pa, t-1, z-2, y, x, EB )\
		       + GRID_ENTRY( pa, t-1, z-2, y, x, WT ) - GRID_ENTRY( pa, t-1, z-2, y, x, WB );\
		ux2 /= rho2;\
		uy2 /= rho2;\
		uz2 /= rho2;\
\
		rho = 1.0;\
\
		ux = 2*ux1 - ux2;\
		uy = 2*uy1 - uy2;\
		uz = 2*uz1 - uz2;\
\
		u2 = 1.5 * (ux*ux + uy*uy + uz*uz);\
\
		LOCAL( pa, t-1, z, y, x, C ) = DFL1*rho*(1.0                                 - u2);\
\
		LOCAL( pa, t-1, z, y, x, N ) = DFL2*rho*(1.0 +       uy*(4.5*uy       + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, S ) = DFL2*rho*(1.0 +       uy*(4.5*uy       - 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, E ) = DFL2*rho*(1.0 +       ux*(4.5*ux       + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, W ) = DFL2*rho*(1.0 +       ux*(4.5*ux       - 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, T ) = DFL2*rho*(1.0 +       uz*(4.5*uz       + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, B ) = DFL2*rho*(1.0 +       uz*(4.5*uz       - 3.0) - u2);\
\
		LOCAL( pa, t-1, z, y, x, NE) = DFL3*rho*(1.0 + (+ux+uy)*(4.5*(+ux+uy) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, NW) = DFL3*rho*(1.0 + (-ux+uy)*(4.5*(-ux+uy) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, SE) = DFL3*rho*(1.0 + (+ux-uy)*(4.5*(+ux-uy) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, SW) = DFL3*rho*(1.0 + (-ux-uy)*(4.5*(-ux-uy) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, NT) = DFL3*rho*(1.0 + (+uy+uz)*(4.5*(+uy+uz) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, NB) = DFL3*rho*(1.0 + (+uy-uz)*(4.5*(+uy-uz) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, ST) = DFL3*rho*(1.0 + (-uy+uz)*(4.5*(-uy+uz) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, SB) = DFL3*rho*(1.0 + (-uy-uz)*(4.5*(-uy-uz) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, ET) = DFL3*rho*(1.0 + (+ux+uz)*(4.5*(+ux+uz) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, EB) = DFL3*rho*(1.0 + (+ux-uz)*(4.5*(+ux-uz) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, WT) = DFL3*rho*(1.0 + (-ux+uz)*(4.5*(-ux+uz) + 3.0) - u2);\
		LOCAL( pa, t-1, z, y, x, WB) = DFL3*rho*(1.0 + (-ux-uz)*(4.5*(-ux-uz) + 3.0) - u2);\
    }\
} while (0)

#endif /* _LBM_H_ */
