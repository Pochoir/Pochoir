/* $Id: lbm.c,v 1.6 2004/05/03 08:23:51 pohlt Exp $ */

/*############################################################################*/

#include <pochoir.hpp>
#include "lbm_tang.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


#if !defined(SPEC_CPU)
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

extern int SIZE_X, SIZE_Y, SIZE_Z;
/*############################################################################*/

void LBM_initializeGrid( Pochoir_Array_3D(PoCellEntry) & pa, int t ) {
	/*voption indep*/
#if !defined(SPEC_CPU)
#ifdef _OPENMP
#pragma omp parallel for
#endif
#endif
    for (int z = 0; z < SIZE_Z + 2 * MARGIN_Z; ++z) {
        for (int y = 0; y < SIZE_Y; ++y) {
    for (int x = 0; x < SIZE_X; ++x) {
		LOCAL( pa, t, z, y, x, C  ) = DFL1;
		LOCAL( pa, t, z, y, x, N  ) = DFL2;
		LOCAL( pa, t, z, y, x, S  ) = DFL2;
		LOCAL( pa, t, z, y, x, E  ) = DFL2;
		LOCAL( pa, t, z, y, x, W  ) = DFL2;
		LOCAL( pa, t, z, y, x, T  ) = DFL2;
		LOCAL( pa, t, z, y, x, B  ) = DFL2;
		LOCAL( pa, t, z, y, x, NE ) = DFL3;
		LOCAL( pa, t, z, y, x, NW ) = DFL3;
		LOCAL( pa, t, z, y, x, SE ) = DFL3;
		LOCAL( pa, t, z, y, x, SW ) = DFL3;
		LOCAL( pa, t, z, y, x, NT ) = DFL3;
		LOCAL( pa, t, z, y, x, NB ) = DFL3;
		LOCAL( pa, t, z, y, x, ST ) = DFL3;
		LOCAL( pa, t, z, y, x, SB ) = DFL3;
		LOCAL( pa, t, z, y, x, ET ) = DFL3;
		LOCAL( pa, t, z, y, x, EB ) = DFL3;
		LOCAL( pa, t, z, y, x, WT ) = DFL3;
		LOCAL( pa, t, z, y, x, WB ) = DFL3;
		CLEAR_ALL_FLAGS_SWEEP( pa, t, z, y, x );
    }
        }
    }
}

/*############################################################################*/

void LBM_loadRandomObstacle( Pochoir_Array_3D(PoCellEntry) & pa, int t ) {
	for( int z = 0 + MARGIN_Z; z < SIZE_Z + MARGIN_Z; ++z ) {
		for( int y = 0; y < SIZE_Y; ++y ) {
	for( int x = 0; x < SIZE_X; ++x ) {
        //if (rand() & 0x1)
        //    SET_FLAG_SWEEP( pa, t, z, y, x, OBSTACLE );
	}
		}
	}
}

/*############################################################################*/


void LBM_loadObstacleFile( Pochoir_Array_3D(PoCellEntry) & pa, int t, const char* filename ) {
	FILE* file = fopen( filename, "rb" );

	for( int z = 0 + MARGIN_Z; z < SIZE_Z + MARGIN_Z; ++z ) {
		for( int y = 0; y < SIZE_Y; ++y ) {
	for( int x = 0; x < SIZE_X; ++x ) {
		if( fgetc( file ) != '.' ) 
            SET_FLAG_SWEEP( pa, t, z, y, x, OBSTACLE );
	}
		fgetc( file );
		}
		fgetc( file );
	}

	fclose( file );
}

/*############################################################################*/

void LBM_initializeSpecialCellsForLDC( Pochoir_Array_3D(PoCellEntry) & pa, int t ) {
	/*voption indep*/
#if !defined(SPEC_CPU)
#ifdef _OPENMP
#pragma omp parallel for private( x, y )
#endif
#endif
	for( int z = 0; z < SIZE_Z+2*MARGIN_Z; ++z ) {
		for( int y = 0; y < SIZE_Y; ++y ) {
	for( int x = 0; x < SIZE_X; ++x ) {
		if( x == 0 || x == SIZE_X-1 ||
		    y == 0 || y == SIZE_Y-1 ||
		    z == 0 + MARGIN_Z || z == SIZE_Z-1+MARGIN_Z ) {
			SET_FLAG_SWEEP( pa, t, z, y, x, OBSTACLE );
		} else {
            if( (z == 1 + MARGIN_Z || z == SIZE_Z-2 + MARGIN_Z) &&
			     x > 1 && x < SIZE_X-2 &&
			     y > 1 && y < SIZE_Y-2 ) {
				SET_FLAG_SWEEP( pa, t, z, y, x, ACCEL );
            }
		}
	}
		}
	}
}

/*############################################################################*/

void LBM_initializeSpecialCellsForChannel( Pochoir_Array_3D(PoCellEntry) & pa, int t ) {
	/*voption indep*/
#if !defined(SPEC_CPU)
#ifdef _OPENMP
#pragma omp parallel for private( x, y )
#endif
#endif
	for( int z = 0; z < SIZE_Z+2*MARGIN_Z; ++z ) {
		for( int y = 0; y < SIZE_Y; ++y ) {
	for( int x = 0; x < SIZE_X; ++x ) {
		if( x == 0 || x == SIZE_X-1 ||
		    y == 0 || y == SIZE_Y-1 ) {
			SET_FLAG_SWEEP( pa, t, z, y, x, OBSTACLE );

            /* is this logic correct? 
             * The following branch probably will never get executed!
             */
			if( (z == 0 + MARGIN_Z || z == SIZE_Z-1+MARGIN_Z) &&
			    ! TEST_FLAG_SWEEP( pa, t, z, y, x, OBSTACLE )) {
                printf("SET_FLAG_SWEEP(pa, %d, %d, %d, %d) = IN_OUT_FLOW\n", t, z, y, x);
				SET_FLAG_SWEEP( pa, t, z, y, x, IN_OUT_FLOW );
            }
		}
	}
		}
	}
}

/*############################################################################*/

void LBM_performStreamCollide( Pochoir_Array_3D(PoCellEntry) & pa, int t, int z, int y, int x ) {
	double ux, uy, uz, u2, rho;

	/*voption indep*/
#if !defined(SPEC_CPU)
#ifdef _OPENMP
#pragma omp parallel for private( ux, uy, uz, u2, rho )
#endif
#endif
	if( TEST_FLAG_SWEEP( pa.interior, t-1, z, y, x, OBSTACLE )) {
		DST_C ( pa.interior, t, z, y, x ) = SRC_C ( pa.interior, t-1, z, y, x );
		DST_S ( pa.interior, t, z, y, x ) = SRC_N ( pa.interior, t-1, z, y, x );
		DST_N ( pa.interior, t, z, y, x ) = SRC_S ( pa.interior, t-1, z, y, x );
		DST_W ( pa.interior, t, z, y, x ) = SRC_E ( pa.interior, t-1, z, y, x );
		DST_E ( pa.interior, t, z, y, x ) = SRC_W ( pa.interior, t-1, z, y, x );
		DST_B ( pa.interior, t, z, y, x ) = SRC_T ( pa.interior, t-1, z, y, x );
		DST_T ( pa.interior, t, z, y, x ) = SRC_B ( pa.interior, t-1, z, y, x );
		DST_SW( pa.interior, t, z, y, x ) = SRC_NE( pa.interior, t-1, z, y, x );
		DST_SE( pa.interior, t, z, y, x ) = SRC_NW( pa.interior, t-1, z, y, x );
		DST_NW( pa.interior, t, z, y, x ) = SRC_SE( pa.interior, t-1, z, y, x );
		DST_NE( pa.interior, t, z, y, x ) = SRC_SW( pa.interior, t-1, z, y, x );
		DST_SB( pa.interior, t, z, y, x ) = SRC_NT( pa.interior, t-1, z, y, x );
		DST_ST( pa.interior, t, z, y, x ) = SRC_NB( pa.interior, t-1, z, y, x );
		DST_NB( pa.interior, t, z, y, x ) = SRC_ST( pa.interior, t-1, z, y, x );
		DST_NT( pa.interior, t, z, y, x ) = SRC_SB( pa.interior, t-1, z, y, x );
		DST_WB( pa.interior, t, z, y, x ) = SRC_ET( pa.interior, t-1, z, y, x );
		DST_WT( pa.interior, t, z, y, x ) = SRC_EB( pa.interior, t-1, z, y, x );
		DST_EB( pa.interior, t, z, y, x ) = SRC_WT( pa.interior, t-1, z, y, x );
		DST_ET( pa.interior, t, z, y, x ) = SRC_WB( pa.interior, t-1, z, y, x );
        return;
//		continue;
	}

	rho = + SRC_C ( pa.interior, t-1, z, y, x ) + SRC_N ( pa.interior, t-1, z, y, x )
	      + SRC_S ( pa.interior, t-1, z, y, x ) + SRC_E ( pa.interior, t-1, z, y, x )
	      + SRC_W ( pa.interior, t-1, z, y, x ) + SRC_T ( pa.interior, t-1, z, y, x )
	      + SRC_B ( pa.interior, t-1, z, y, x ) + SRC_NE( pa.interior, t-1, z, y, x )
	      + SRC_NW( pa.interior, t-1, z, y, x ) + SRC_SE( pa.interior, t-1, z, y, x )
	      + SRC_SW( pa.interior, t-1, z, y, x ) + SRC_NT( pa.interior, t-1, z, y, x )
	      + SRC_NB( pa.interior, t-1, z, y, x ) + SRC_ST( pa.interior, t-1, z, y, x )
	      + SRC_SB( pa.interior, t-1, z, y, x ) + SRC_ET( pa.interior, t-1, z, y, x )
	      + SRC_EB( pa.interior, t-1, z, y, x ) + SRC_WT( pa.interior, t-1, z, y, x )
	      + SRC_WB( pa.interior, t-1, z, y, x );

	ux  = + SRC_E ( pa.interior, t-1, z, y, x ) - SRC_W ( pa.interior, t-1, z, y, x )
	      + SRC_NE( pa.interior, t-1, z, y, x ) - SRC_NW( pa.interior, t-1, z, y, x )
	      + SRC_SE( pa.interior, t-1, z, y, x ) - SRC_SW( pa.interior, t-1, z, y, x )
	      + SRC_ET( pa.interior, t-1, z, y, x ) + SRC_EB( pa.interior, t-1, z, y, x )
	      - SRC_WT( pa.interior, t-1, z, y, x ) - SRC_WB( pa.interior, t-1, z, y, x );

	uy  = + SRC_N ( pa.interior, t-1, z, y, x ) - SRC_S ( pa.interior, t-1, z, y, x )
	      + SRC_NE( pa.interior, t-1, z, y, x ) + SRC_NW( pa.interior, t-1, z, y, x )
	      - SRC_SE( pa.interior, t-1, z, y, x ) - SRC_SW( pa.interior, t-1, z, y, x )
	      + SRC_NT( pa.interior, t-1, z, y, x ) + SRC_NB( pa.interior, t-1, z, y, x )
	      - SRC_ST( pa.interior, t-1, z, y, x ) - SRC_SB( pa.interior, t-1, z, y, x );

	uz  = + SRC_T ( pa.interior, t-1, z, y, x ) - SRC_B ( pa.interior, t-1, z, y, x )
	      + SRC_NT( pa.interior, t-1, z, y, x ) - SRC_NB( pa.interior, t-1, z, y, x )
	      + SRC_ST( pa.interior, t-1, z, y, x ) - SRC_SB( pa.interior, t-1, z, y, x )
	      + SRC_ET( pa.interior, t-1, z, y, x ) - SRC_EB( pa.interior, t-1, z, y, x )
	      + SRC_WT( pa.interior, t-1, z, y, x ) - SRC_WB( pa.interior, t-1, z, y, x );

	ux /= rho;
	uy /= rho;
	uz /= rho;

	if( TEST_FLAG_SWEEP( pa.interior, t-1, z, y, x, ACCEL )) {
		ux = 0.005;
		uy = 0.002;
		uz = 0.000;
	}

	u2 = 1.5 * (ux*ux + uy*uy + uz*uz);
	DST_C ( pa.interior, t, z, y, x ) = (1.0-OMEGA)*SRC_C ( pa.interior, t-1, z, y, x ) + DFL1*OMEGA*rho*(1.0                                 - u2);

	DST_N ( pa.interior, t, z, y, x ) = (1.0-OMEGA)*SRC_N ( pa.interior, t-1, z, y, x ) + DFL2*OMEGA*rho*(1.0 +       uy*(4.5*uy       + 3.0) - u2);
	DST_S ( pa.interior, t, z, y, x ) = (1.0-OMEGA)*SRC_S ( pa.interior, t-1, z, y, x ) + DFL2*OMEGA*rho*(1.0 +       uy*(4.5*uy       - 3.0) - u2);
	DST_E ( pa.interior, t, z, y, x ) = (1.0-OMEGA)*SRC_E ( pa.interior, t-1, z, y, x ) + DFL2*OMEGA*rho*(1.0 +       ux*(4.5*ux       + 3.0) - u2);
	DST_W ( pa.interior, t, z, y, x ) = (1.0-OMEGA)*SRC_W ( pa.interior, t-1, z, y, x ) + DFL2*OMEGA*rho*(1.0 +       ux*(4.5*ux       - 3.0) - u2);
	DST_T ( pa.interior, t, z, y, x ) = (1.0-OMEGA)*SRC_T ( pa.interior, t-1, z, y, x ) + DFL2*OMEGA*rho*(1.0 +       uz*(4.5*uz       + 3.0) - u2);
	DST_B ( pa.interior, t, z, y, x ) = (1.0-OMEGA)*SRC_B ( pa.interior, t-1, z, y, x ) + DFL2*OMEGA*rho*(1.0 +       uz*(4.5*uz       - 3.0) - u2);

	DST_NE( pa.interior, t, z, y, x ) = (1.0-OMEGA)*SRC_NE( pa.interior, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (+ux+uy)*(4.5*(+ux+uy) + 3.0) - u2);
	DST_NW( pa.interior, t, z, y, x ) = (1.0-OMEGA)*SRC_NW( pa.interior, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (-ux+uy)*(4.5*(-ux+uy) + 3.0) - u2);
	DST_SE( pa.interior, t, z, y, x ) = (1.0-OMEGA)*SRC_SE( pa.interior, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (+ux-uy)*(4.5*(+ux-uy) + 3.0) - u2);
	DST_SW( pa.interior, t, z, y, x ) = (1.0-OMEGA)*SRC_SW( pa.interior, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (-ux-uy)*(4.5*(-ux-uy) + 3.0) - u2);
	DST_NT( pa.interior, t, z, y, x ) = (1.0-OMEGA)*SRC_NT( pa.interior, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (+uy+uz)*(4.5*(+uy+uz) + 3.0) - u2);
	DST_NB( pa.interior, t, z, y, x ) = (1.0-OMEGA)*SRC_NB( pa.interior, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (+uy-uz)*(4.5*(+uy-uz) + 3.0) - u2);
	DST_ST( pa.interior, t, z, y, x ) = (1.0-OMEGA)*SRC_ST( pa.interior, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (-uy+uz)*(4.5*(-uy+uz) + 3.0) - u2);
	DST_SB( pa.interior, t, z, y, x ) = (1.0-OMEGA)*SRC_SB( pa.interior, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (-uy-uz)*(4.5*(-uy-uz) + 3.0) - u2);
	DST_ET( pa.interior, t, z, y, x ) = (1.0-OMEGA)*SRC_ET( pa.interior, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (+ux+uz)*(4.5*(+ux+uz) + 3.0) - u2);
	DST_EB( pa.interior, t, z, y, x ) = (1.0-OMEGA)*SRC_EB( pa.interior, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (+ux-uz)*(4.5*(+ux-uz) + 3.0) - u2);
	DST_WT( pa.interior, t, z, y, x ) = (1.0-OMEGA)*SRC_WT( pa.interior, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (-ux+uz)*(4.5*(-ux+uz) + 3.0) - u2);
	DST_WB( pa.interior, t, z, y, x ) = (1.0-OMEGA)*SRC_WB( pa.interior, t-1, z, y, x ) + DFL3*OMEGA*rho*(1.0 + (-ux-uz)*(4.5*(-ux-uz) + 3.0) - u2);
}

/*############################################################################*/

void LBM_handleInOutFlow( Pochoir_Array_3D(PoCellEntry) & pa, int t, int z, int y, int x ) {
	double ux , uy , uz , rho ,
	       ux1, uy1, uz1, rho1,
	       ux2, uy2, uz2, rho2,
	       u2, px, py;
	/* inflow */
	/*voption indep*/
#if !defined(SPEC_CPU)
#ifdef _OPENMP
#pragma omp parallel for private( ux, uy, uz, rho, ux1, uy1, uz1, rho1, \
                                  ux2, uy2, uz2, rho2, u2, px, py )
#endif
#endif
    if (z == 0 + MARGIN_Z) {
		rho1 = + GRID_ENTRY( pa.interior, t-1, z+1, y, x, C  ) + GRID_ENTRY( pa.interior, t-1, z+1, y, x, N  )
		       + GRID_ENTRY( pa.interior, t-1, z+1, y, x, S  ) + GRID_ENTRY( pa.interior, t-1, z+1, y, x, E  )
		       + GRID_ENTRY( pa.interior, t-1, z+1, y, x, W  ) + GRID_ENTRY( pa.interior, t-1, z+1, y, x, T  )
		       + GRID_ENTRY( pa.interior, t-1, z+1, y, x, B  ) + GRID_ENTRY( pa.interior, t-1, z+1, y, x, NE )
		       + GRID_ENTRY( pa.interior, t-1, z+1, y, x, NW ) + GRID_ENTRY( pa.interior, t-1, z+1, y, x, SE )
		       + GRID_ENTRY( pa.interior, t-1, z+1, y, x, SW ) + GRID_ENTRY( pa.interior, t-1, z+1, y, x, NT )
		       + GRID_ENTRY( pa.interior, t-1, z+1, y, x, NB ) + GRID_ENTRY( pa.interior, t-1, z+1, y, x, ST )
		       + GRID_ENTRY( pa.interior, t-1, z+1, y, x, SB ) + GRID_ENTRY( pa.interior, t-1, z+1, y, x, ET )
		       + GRID_ENTRY( pa.interior, t-1, z+1, y, x, EB ) + GRID_ENTRY( pa.interior, t-1, z+1, y, x, WT )
		       + GRID_ENTRY( pa.interior, t-1, z+1, y, x, WB );                            
		rho2 = + GRID_ENTRY( pa.interior, t-1, z+2, y, x, C  ) + GRID_ENTRY( pa.interior, t-1, z+2, y, x, N  )
		       + GRID_ENTRY( pa.interior, t-1, z+2, y, x, S  ) + GRID_ENTRY( pa.interior, t-1, z+2, y, x, E  )
		       + GRID_ENTRY( pa.interior, t-1, z+2, y, x, W  ) + GRID_ENTRY( pa.interior, t-1, z+2, y, x, T  )
		       + GRID_ENTRY( pa.interior, t-1, z+2, y, x, B  ) + GRID_ENTRY( pa.interior, t-1, z+2, y, x, NE )
		       + GRID_ENTRY( pa.interior, t-1, z+2, y, x, NW ) + GRID_ENTRY( pa.interior, t-1, z+2, y, x, SE )
		       + GRID_ENTRY( pa.interior, t-1, z+2, y, x, SW ) + GRID_ENTRY( pa.interior, t-1, z+2, y, x, NT )
		       + GRID_ENTRY( pa.interior, t-1, z+2, y, x, NB ) + GRID_ENTRY( pa.interior, t-1, z+2, y, x, ST )
		       + GRID_ENTRY( pa.interior, t-1, z+2, y, x, SB ) + GRID_ENTRY( pa.interior, t-1, z+2, y, x, ET )
		       + GRID_ENTRY( pa.interior, t-1, z+2, y, x, EB ) + GRID_ENTRY( pa.interior, t-1, z+2, y, x, WT )
		       + GRID_ENTRY( pa.interior, t-1, z+2, y, x, WB );

		rho = 2.0*rho1 - rho2;

		px = (x / (0.5*(SIZE_X-1))) - 1.0;
		py = (y / (0.5*(SIZE_Y-1))) - 1.0;
		ux = 0.00;
		uy = 0.00;
		uz = 0.01 * (1.0-px*px) * (1.0-py*py);

		u2 = 1.5 * (ux*ux + uy*uy + uz*uz);

		LOCAL( pa.interior, t-1, z, y, x, C ) = DFL1*rho*(1.0                                 - u2);

		LOCAL( pa.interior, t-1, z, y, x, N ) = DFL2*rho*(1.0 +       uy*(4.5*uy       + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, S ) = DFL2*rho*(1.0 +       uy*(4.5*uy       - 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, E ) = DFL2*rho*(1.0 +       ux*(4.5*ux       + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, W ) = DFL2*rho*(1.0 +       ux*(4.5*ux       - 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, T ) = DFL2*rho*(1.0 +       uz*(4.5*uz       + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, B ) = DFL2*rho*(1.0 +       uz*(4.5*uz       - 3.0) - u2);

		LOCAL( pa.interior, t-1, z, y, x, NE) = DFL3*rho*(1.0 + (+ux+uy)*(4.5*(+ux+uy) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, NW) = DFL3*rho*(1.0 + (-ux+uy)*(4.5*(-ux+uy) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, SE) = DFL3*rho*(1.0 + (+ux-uy)*(4.5*(+ux-uy) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, SW) = DFL3*rho*(1.0 + (-ux-uy)*(4.5*(-ux-uy) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, NT) = DFL3*rho*(1.0 + (+uy+uz)*(4.5*(+uy+uz) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, NB) = DFL3*rho*(1.0 + (+uy-uz)*(4.5*(+uy-uz) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, ST) = DFL3*rho*(1.0 + (-uy+uz)*(4.5*(-uy+uz) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, SB) = DFL3*rho*(1.0 + (-uy-uz)*(4.5*(-uy-uz) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, ET) = DFL3*rho*(1.0 + (+ux+uz)*(4.5*(+ux+uz) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, EB) = DFL3*rho*(1.0 + (+ux-uz)*(4.5*(+ux-uz) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, WT) = DFL3*rho*(1.0 + (-ux+uz)*(4.5*(-ux+uz) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, WB) = DFL3*rho*(1.0 + (-ux-uz)*(4.5*(-ux-uz) + 3.0) - u2);
    } /* end if (z == 0 + MARGIN_Z) */

	/* outflow */
	/*voption indep*/
#if !defined(SPEC_CPU)
#ifdef _OPENMP
#pragma omp parallel for private( ux, uy, uz, rho, ux1, uy1, uz1, rho1, \
                                  ux2, uy2, uz2, rho2, u2, px, py )
#endif
#endif

    if (z == SIZE_Z - 1 + MARGIN_Z) {
		rho1 = + GRID_ENTRY( pa.interior, t-1, z-1, y, x, C  ) + GRID_ENTRY( pa.interior, t-1, z-1, y, x, N  )
		       + GRID_ENTRY( pa.interior, t-1, z-1, y, x, S  ) + GRID_ENTRY( pa.interior, t-1, z-1, y, x, E  )
		       + GRID_ENTRY( pa.interior, t-1, z-1, y, x, W  ) + GRID_ENTRY( pa.interior, t-1, z-1, y, x, T  )
		       + GRID_ENTRY( pa.interior, t-1, z-1, y, x, B  ) + GRID_ENTRY( pa.interior, t-1, z-1, y, x, NE )
		       + GRID_ENTRY( pa.interior, t-1, z-1, y, x, NW ) + GRID_ENTRY( pa.interior, t-1, z-1, y, x, SE )
		       + GRID_ENTRY( pa.interior, t-1, z-1, y, x, SW ) + GRID_ENTRY( pa.interior, t-1, z-1, y, x, NT )
		       + GRID_ENTRY( pa.interior, t-1, z-1, y, x, NB ) + GRID_ENTRY( pa.interior, t-1, z-1, y, x, ST )
		       + GRID_ENTRY( pa.interior, t-1, z-1, y, x, SB ) + GRID_ENTRY( pa.interior, t-1, z-1, y, x, ET )
		       + GRID_ENTRY( pa.interior, t-1, z-1, y, x, EB ) + GRID_ENTRY( pa.interior, t-1, z-1, y, x, WT )
		       + GRID_ENTRY( pa.interior, t-1, z-1, y, x, WB );                          
		ux1  = + GRID_ENTRY( pa.interior, t-1, z-1, y, x, E  ) - GRID_ENTRY( pa.interior, t-1, z-1, y, x, W  )
		       + GRID_ENTRY( pa.interior, t-1, z-1, y, x, NE ) - GRID_ENTRY( pa.interior, t-1, z-1, y, x, NW )
		       + GRID_ENTRY( pa.interior, t-1, z-1, y, x, SE ) - GRID_ENTRY( pa.interior, t-1, z-1, y, x, SW )
		       + GRID_ENTRY( pa.interior, t-1, z-1, y, x, ET ) + GRID_ENTRY( pa.interior, t-1, z-1, y, x, EB )
		       - GRID_ENTRY( pa.interior, t-1, z-1, y, x, WT ) - GRID_ENTRY( pa.interior, t-1, z-1, y, x, WB );
		uy1  = + GRID_ENTRY( pa.interior, t-1, z-1, y, x, N  ) - GRID_ENTRY( pa.interior, t-1, z-1, y, x, S  )
		       + GRID_ENTRY( pa.interior, t-1, z-1, y, x, NE ) + GRID_ENTRY( pa.interior, t-1, z-1, y, x, NW )
		       - GRID_ENTRY( pa.interior, t-1, z-1, y, x, SE ) - GRID_ENTRY( pa.interior, t-1, z-1, y, x, SW )
		       + GRID_ENTRY( pa.interior, t-1, z-1, y, x, NT ) + GRID_ENTRY( pa.interior, t-1, z-1, y, x, NB )
		       - GRID_ENTRY( pa.interior, t-1, z-1, y, x, ST ) - GRID_ENTRY( pa.interior, t-1, z-1, y, x, SB );
		uz1  = + GRID_ENTRY( pa.interior, t-1, z-1, y, x, T  ) - GRID_ENTRY( pa.interior, t-1, z-1, y, x, B  )
		       + GRID_ENTRY( pa.interior, t-1, z-1, y, x, NT ) - GRID_ENTRY( pa.interior, t-1, z-1, y, x, NB )
		       + GRID_ENTRY( pa.interior, t-1, z-1, y, x, ST ) - GRID_ENTRY( pa.interior, t-1, z-1, y, x, SB )
		       + GRID_ENTRY( pa.interior, t-1, z-1, y, x, ET ) - GRID_ENTRY( pa.interior, t-1, z-1, y, x, EB )
		       + GRID_ENTRY( pa.interior, t-1, z-1, y, x, WT ) - GRID_ENTRY( pa.interior, t-1, z-1, y, x, WB );

		ux1 /= rho1;
		uy1 /= rho1;
		uz1 /= rho1;

		rho2 = + GRID_ENTRY( pa.interior, t-1, z-2, y, x, C  ) + GRID_ENTRY( pa.interior, t-1, z-2, y, x, N  )
		       + GRID_ENTRY( pa.interior, t-1, z-2, y, x, S  ) + GRID_ENTRY( pa.interior, t-1, z-2, y, x, E  )
		       + GRID_ENTRY( pa.interior, t-1, z-2, y, x, W  ) + GRID_ENTRY( pa.interior, t-1, z-2, y, x, T  )
		       + GRID_ENTRY( pa.interior, t-1, z-2, y, x, B  ) + GRID_ENTRY( pa.interior, t-1, z-2, y, x, NE )
		       + GRID_ENTRY( pa.interior, t-1, z-2, y, x, NW ) + GRID_ENTRY( pa.interior, t-1, z-2, y, x, SE )
		       + GRID_ENTRY( pa.interior, t-1, z-2, y, x, SW ) + GRID_ENTRY( pa.interior, t-1, z-2, y, x, NT )
		       + GRID_ENTRY( pa.interior, t-1, z-2, y, x, NB ) + GRID_ENTRY( pa.interior, t-1, z-2, y, x, ST )
		       + GRID_ENTRY( pa.interior, t-1, z-2, y, x, SB ) + GRID_ENTRY( pa.interior, t-1, z-2, y, x, ET )
		       + GRID_ENTRY( pa.interior, t-1, z-2, y, x, EB ) + GRID_ENTRY( pa.interior, t-1, z-2, y, x, WT )
		       + GRID_ENTRY( pa.interior, t-1, z-2, y, x, WB );                           
		ux2  = + GRID_ENTRY( pa.interior, t-1, z-2, y, x, E  ) - GRID_ENTRY( pa.interior, t-1, z-2, y, x, W  )
		       + GRID_ENTRY( pa.interior, t-1, z-2, y, x, NE ) - GRID_ENTRY( pa.interior, t-1, z-2, y, x, NW )
		       + GRID_ENTRY( pa.interior, t-1, z-2, y, x, SE ) - GRID_ENTRY( pa.interior, t-1, z-2, y, x, SW )
		       + GRID_ENTRY( pa.interior, t-1, z-2, y, x, ET ) + GRID_ENTRY( pa.interior, t-1, z-2, y, x, EB )
		       - GRID_ENTRY( pa.interior, t-1, z-2, y, x, WT ) - GRID_ENTRY( pa.interior, t-1, z-2, y, x, WB );
		uy2  = + GRID_ENTRY( pa.interior, t-1, z-2, y, x, N  ) - GRID_ENTRY( pa.interior, t-1, z-2, y, x, S  )
		       + GRID_ENTRY( pa.interior, t-1, z-2, y, x, NE ) + GRID_ENTRY( pa.interior, t-1, z-2, y, x, NW )
		       - GRID_ENTRY( pa.interior, t-1, z-2, y, x, SE ) - GRID_ENTRY( pa.interior, t-1, z-2, y, x, SW )
		       + GRID_ENTRY( pa.interior, t-1, z-2, y, x, NT ) + GRID_ENTRY( pa.interior, t-1, z-2, y, x, NB )
		       - GRID_ENTRY( pa.interior, t-1, z-2, y, x, ST ) - GRID_ENTRY( pa.interior, t-1, z-2, y, x, SB );
		uz2  = + GRID_ENTRY( pa.interior, t-1, z-2, y, x, T  ) - GRID_ENTRY( pa.interior, t-1, z-2, y, x, B  )
		       + GRID_ENTRY( pa.interior, t-1, z-2, y, x, NT ) - GRID_ENTRY( pa.interior, t-1, z-2, y, x, NB )
		       + GRID_ENTRY( pa.interior, t-1, z-2, y, x, ST ) - GRID_ENTRY( pa.interior, t-1, z-2, y, x, SB )
		       + GRID_ENTRY( pa.interior, t-1, z-2, y, x, ET ) - GRID_ENTRY( pa.interior, t-1, z-2, y, x, EB )
		       + GRID_ENTRY( pa.interior, t-1, z-2, y, x, WT ) - GRID_ENTRY( pa.interior, t-1, z-2, y, x, WB );

		ux2 /= rho2;
		uy2 /= rho2;
		uz2 /= rho2;

		rho = 1.0;

		ux = 2*ux1 - ux2;
		uy = 2*uy1 - uy2;
		uz = 2*uz1 - uz2;

		u2 = 1.5 * (ux*ux + uy*uy + uz*uz);

		LOCAL( pa.interior, t-1, z, y, x, C ) = DFL1*rho*(1.0                                 - u2);

		LOCAL( pa.interior, t-1, z, y, x, N ) = DFL2*rho*(1.0 +       uy*(4.5*uy       + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, S ) = DFL2*rho*(1.0 +       uy*(4.5*uy       - 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, E ) = DFL2*rho*(1.0 +       ux*(4.5*ux       + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, W ) = DFL2*rho*(1.0 +       ux*(4.5*ux       - 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, T ) = DFL2*rho*(1.0 +       uz*(4.5*uz       + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, B ) = DFL2*rho*(1.0 +       uz*(4.5*uz       - 3.0) - u2);

		LOCAL( pa.interior, t-1, z, y, x, NE) = DFL3*rho*(1.0 + (+ux+uy)*(4.5*(+ux+uy) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, NW) = DFL3*rho*(1.0 + (-ux+uy)*(4.5*(-ux+uy) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, SE) = DFL3*rho*(1.0 + (+ux-uy)*(4.5*(+ux-uy) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, SW) = DFL3*rho*(1.0 + (-ux-uy)*(4.5*(-ux-uy) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, NT) = DFL3*rho*(1.0 + (+uy+uz)*(4.5*(+uy+uz) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, NB) = DFL3*rho*(1.0 + (+uy-uz)*(4.5*(+uy-uz) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, ST) = DFL3*rho*(1.0 + (-uy+uz)*(4.5*(-uy+uz) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, SB) = DFL3*rho*(1.0 + (-uy-uz)*(4.5*(-uy-uz) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, ET) = DFL3*rho*(1.0 + (+ux+uz)*(4.5*(+ux+uz) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, EB) = DFL3*rho*(1.0 + (+ux-uz)*(4.5*(+ux-uz) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, WT) = DFL3*rho*(1.0 + (-ux+uz)*(4.5*(-ux+uz) + 3.0) - u2);
		LOCAL( pa.interior, t-1, z, y, x, WB) = DFL3*rho*(1.0 + (-ux-uz)*(4.5*(-ux-uz) + 3.0) - u2);
    }
}

/*############################################################################*/

void LBM_showGridStatistics( Pochoir_Array_3D(PoCellEntry) & pa, int t ) {
	int nObstacleCells = 0,
	    nAccelCells    = 0,
	    nFluidCells    = 0;
	double ux, uy, uz;
	double minU2  = 1e+30, maxU2  = -1e+30, u2;
	double minRho = 1e+30, maxRho = -1e+30, rho;
	double mass = 0;

    for (int z = 0 + MARGIN_Z; z < SIZE_Z + MARGIN_Z; ++z) {
        for (int y = 0; y < SIZE_Y; ++y) {
    for (int x = 0; x < SIZE_X; ++x) {
		rho = + LOCAL( pa, t, z, y, x, C  ) + LOCAL( pa, t, z, y, x, N  )
		      + LOCAL( pa, t, z, y, x, S  ) + LOCAL( pa, t, z, y, x, E  )
		      + LOCAL( pa, t, z, y, x, W  ) + LOCAL( pa, t, z, y, x, T  )
		      + LOCAL( pa, t, z, y, x, B  ) + LOCAL( pa, t, z, y, x, NE )
		      + LOCAL( pa, t, z, y, x, NW ) + LOCAL( pa, t, z, y, x, SE )
		      + LOCAL( pa, t, z, y, x, SW ) + LOCAL( pa, t, z, y, x, NT )
		      + LOCAL( pa, t, z, y, x, NB ) + LOCAL( pa, t, z, y, x, ST )
		      + LOCAL( pa, t, z, y, x, SB ) + LOCAL( pa, t, z, y, x, ET )
		      + LOCAL( pa, t, z, y, x, EB ) + LOCAL( pa, t, z, y, x, WT )
		      + LOCAL( pa, t, z, y, x, WB );
		if( rho < minRho ) minRho = rho;
		if( rho > maxRho ) maxRho = rho;
		mass += rho;

		if( TEST_FLAG_SWEEP( pa, t, z, y, x, OBSTACLE )) {
			nObstacleCells++;
		}
		else {
			if( TEST_FLAG_SWEEP( pa, t, z, y, x, ACCEL ))
				nAccelCells++;
			else
				nFluidCells++;

			ux = + LOCAL( pa, t, z, y, x, E  ) - LOCAL( pa, t, z, y, x, W  )
			     + LOCAL( pa, t, z, y, x, NE ) - LOCAL( pa, t, z, y, x, NW )
			     + LOCAL( pa, t, z, y, x, SE ) - LOCAL( pa, t, z, y, x, SW )
			     + LOCAL( pa, t, z, y, x, ET ) + LOCAL( pa, t, z, y, x, EB )
			     - LOCAL( pa, t, z, y, x, WT ) - LOCAL( pa, t, z, y, x, WB );
			uy = + LOCAL( pa, t, z, y, x, N  ) - LOCAL( pa, t, z, y, x, S  )
			     + LOCAL( pa, t, z, y, x, NE ) + LOCAL( pa, t, z, y, x, NW )
			     - LOCAL( pa, t, z, y, x, SE ) - LOCAL( pa, t, z, y, x, SW )
			     + LOCAL( pa, t, z, y, x, NT ) + LOCAL( pa, t, z, y, x, NB )
			     - LOCAL( pa, t, z, y, x, ST ) - LOCAL( pa, t, z, y, x, SB );
			uz = + LOCAL( pa, t, z, y, x, T  ) - LOCAL( pa, t, z, y, x, B  )
			     + LOCAL( pa, t, z, y, x, NT ) - LOCAL( pa, t, z, y, x, NB )
			     + LOCAL( pa, t, z, y, x, ST ) - LOCAL( pa, t, z, y, x, SB )
			     + LOCAL( pa, t, z, y, x, ET ) - LOCAL( pa, t, z, y, x, EB )
			     + LOCAL( pa, t, z, y, x, WT ) - LOCAL( pa, t, z, y, x, WB );
			u2 = (ux*ux + uy*uy + uz*uz) / (rho*rho);
			if( u2 < minU2 ) minU2 = u2;
			if( u2 > maxU2 ) maxU2 = u2;
		}
    }
        }
    }
        printf( "LBM_showGridStatistics:\n"
        "\tnObstacleCells: %7i nAccelCells: %7i nFluidCells: %7i\n"
        "\tminRho: %8.4f maxRho: %8.4f mass: %e\n"
        "\tminU: %e maxU: %e\n\n",
        nObstacleCells, nAccelCells, nFluidCells,
        minRho, maxRho, mass,
        sqrt( minU2 ), sqrt( maxU2 ) );

}

/*############################################################################*/

static void storeValue( FILE* file, OUTPUT_PRECISION* v ) {
	const int litteBigEndianTest = 1;
	if( (*((unsigned char*) &litteBigEndianTest)) == 0 ) {         /* big endian */
		const char* vPtr = (char*) v;
		char buffer[sizeof( OUTPUT_PRECISION )];
		int i;

		for (i = 0; i < sizeof( OUTPUT_PRECISION ); i++)
			buffer[i] = vPtr[sizeof( OUTPUT_PRECISION ) - i - 1];

		fwrite( buffer, sizeof( OUTPUT_PRECISION ), 1, file );
	}
	else {                                                     /* little endian */
		fwrite( v, sizeof( OUTPUT_PRECISION ), 1, file );
	}
}

/*############################################################################*/

static void loadValue( FILE* file, OUTPUT_PRECISION* v ) {
	const int litteBigEndianTest = 1;
	if( (*((unsigned char*) &litteBigEndianTest)) == 0 ) {         /* big endian */
		char* vPtr = (char*) v;
		char buffer[sizeof( OUTPUT_PRECISION )];
		int i;

		fread( buffer, sizeof( OUTPUT_PRECISION ), 1, file );

		for (i = 0; i < sizeof( OUTPUT_PRECISION ); i++)
			vPtr[i] = buffer[sizeof( OUTPUT_PRECISION ) - i - 1];
	}
	else {                                                     /* little endian */
		fread( v, sizeof( OUTPUT_PRECISION ), 1, file );
	}
}

/*############################################################################*/

void LBM_storeVelocityField( Pochoir_Array_3D(PoCellEntry) & pa, int t, 
                             const char* filename, const int binary ) {
	OUTPUT_PRECISION rho, ux, uy, uz;

	FILE* file = fopen( filename, (binary ? "wb" : "w") );

	for( int z = 0 + MARGIN_Z; z < SIZE_Z + MARGIN_Z; ++z ) {
		for( int y = 0; y < SIZE_Y; ++y ) {
	for( int x = 0; x < SIZE_X; ++x ) {
		rho = + GRID_ENTRY( pa, t, z, y, x, C  ) + GRID_ENTRY( pa, t, z, y, x, N  )
		      + GRID_ENTRY( pa, t, z, y, x, S  ) + GRID_ENTRY( pa, t, z, y, x, E  )
		      + GRID_ENTRY( pa, t, z, y, x, W  ) + GRID_ENTRY( pa, t, z, y, x, T  )
		      + GRID_ENTRY( pa, t, z, y, x, B  ) + GRID_ENTRY( pa, t, z, y, x, NE )
		      + GRID_ENTRY( pa, t, z, y, x, NW ) + GRID_ENTRY( pa, t, z, y, x, SE )
		      + GRID_ENTRY( pa, t, z, y, x, SW ) + GRID_ENTRY( pa, t, z, y, x, NT )
		      + GRID_ENTRY( pa, t, z, y, x, NB ) + GRID_ENTRY( pa, t, z, y, x, ST )
		      + GRID_ENTRY( pa, t, z, y, x, SB ) + GRID_ENTRY( pa, t, z, y, x, ET )
		      + GRID_ENTRY( pa, t, z, y, x, EB ) + GRID_ENTRY( pa, t, z, y, x, WT )
		      + GRID_ENTRY( pa, t, z, y, x, WB );                           
		ux  = + GRID_ENTRY( pa, t, z, y, x, E  ) - GRID_ENTRY( pa, t, z, y, x, W  ) 
		      + GRID_ENTRY( pa, t, z, y, x, NE ) - GRID_ENTRY( pa, t, z, y, x, NW ) 
		      + GRID_ENTRY( pa, t, z, y, x, SE ) - GRID_ENTRY( pa, t, z, y, x, SW ) 
		      + GRID_ENTRY( pa, t, z, y, x, ET ) + GRID_ENTRY( pa, t, z, y, x, EB ) 
		      - GRID_ENTRY( pa, t, z, y, x, WT ) - GRID_ENTRY( pa, t, z, y, x, WB );
		uy  = + GRID_ENTRY( pa, t, z, y, x, N  ) - GRID_ENTRY( pa, t, z, y, x, S  ) 
		      + GRID_ENTRY( pa, t, z, y, x, NE ) + GRID_ENTRY( pa, t, z, y, x, NW ) 
		      - GRID_ENTRY( pa, t, z, y, x, SE ) - GRID_ENTRY( pa, t, z, y, x, SW ) 
		      + GRID_ENTRY( pa, t, z, y, x, NT ) + GRID_ENTRY( pa, t, z, y, x, NB ) 
		      - GRID_ENTRY( pa, t, z, y, x, ST ) - GRID_ENTRY( pa, t, z, y, x, SB );
		uz  = + GRID_ENTRY( pa, t, z, y, x, T  ) - GRID_ENTRY( pa, t, z, y, x, B  ) 
		      + GRID_ENTRY( pa, t, z, y, x, NT ) - GRID_ENTRY( pa, t, z, y, x, NB ) 
		      + GRID_ENTRY( pa, t, z, y, x, ST ) - GRID_ENTRY( pa, t, z, y, x, SB ) 
		      + GRID_ENTRY( pa, t, z, y, x, ET ) - GRID_ENTRY( pa, t, z, y, x, EB ) 
		      + GRID_ENTRY( pa, t, z, y, x, WT ) - GRID_ENTRY( pa, t, z, y, x, WB );
		ux /= rho;
		uy /= rho;
		uz /= rho;

		if( binary ) {
			/*
			fwrite( &ux, sizeof( ux ), 1, file );
			fwrite( &uy, sizeof( uy ), 1, file );
			fwrite( &uz, sizeof( uz ), 1, file );
			*/
			storeValue( file, &ux );
			storeValue( file, &uy );
			storeValue( file, &uz );
		} else
			fprintf( file, "%e %e %e\n", ux, uy, uz );
	}
		}
	}

	fclose( file );
}

/*############################################################################*/

void LBM_compareVelocityField( Pochoir_Array_3D(PoCellEntry) & pa, int t, 
                             const char* filename, const int binary ) {
	double rho, ux, uy, uz;
	OUTPUT_PRECISION fileUx, fileUy, fileUz,
	                 dUx, dUy, dUz,
	                 diff2, maxDiff2 = -1e+30;

	FILE* file = fopen( filename, (binary ? "rb" : "r") );

	for( int z = 0+MARGIN_Z; z < SIZE_Z+MARGIN_Z; ++z ) {
		for( int y = 0; y < SIZE_Y; ++y ) {
	for( int x = 0; x < SIZE_X; ++x ) {
		rho = + GRID_ENTRY( pa, t, z, y, x, C  ) + GRID_ENTRY( pa, t, z, y, x, N  )
		      + GRID_ENTRY( pa, t, z, y, x, S  ) + GRID_ENTRY( pa, t, z, y, x, E  )
		      + GRID_ENTRY( pa, t, z, y, x, W  ) + GRID_ENTRY( pa, t, z, y, x, T  )
		      + GRID_ENTRY( pa, t, z, y, x, B  ) + GRID_ENTRY( pa, t, z, y, x, NE )
		      + GRID_ENTRY( pa, t, z, y, x, NW ) + GRID_ENTRY( pa, t, z, y, x, SE )
		      + GRID_ENTRY( pa, t, z, y, x, SW ) + GRID_ENTRY( pa, t, z, y, x, NT )
		      + GRID_ENTRY( pa, t, z, y, x, NB ) + GRID_ENTRY( pa, t, z, y, x, ST )
		      + GRID_ENTRY( pa, t, z, y, x, SB ) + GRID_ENTRY( pa, t, z, y, x, ET )
		      + GRID_ENTRY( pa, t, z, y, x, EB ) + GRID_ENTRY( pa, t, z, y, x, WT )
		      + GRID_ENTRY( pa, t, z, y, x, WB );                           
		ux  = + GRID_ENTRY( pa, t, z, y, x, E  ) - GRID_ENTRY( pa, t, z, y, x, W  ) 
		      + GRID_ENTRY( pa, t, z, y, x, NE ) - GRID_ENTRY( pa, t, z, y, x, NW ) 
		      + GRID_ENTRY( pa, t, z, y, x, SE ) - GRID_ENTRY( pa, t, z, y, x, SW ) 
		      + GRID_ENTRY( pa, t, z, y, x, ET ) + GRID_ENTRY( pa, t, z, y, x, EB ) 
		      - GRID_ENTRY( pa, t, z, y, x, WT ) - GRID_ENTRY( pa, t, z, y, x, WB );
		uy  = + GRID_ENTRY( pa, t, z, y, x, N  ) - GRID_ENTRY( pa, t, z, y, x, S  ) 
		      + GRID_ENTRY( pa, t, z, y, x, NE ) + GRID_ENTRY( pa, t, z, y, x, NW ) 
		      - GRID_ENTRY( pa, t, z, y, x, SE ) - GRID_ENTRY( pa, t, z, y, x, SW ) 
		      + GRID_ENTRY( pa, t, z, y, x, NT ) + GRID_ENTRY( pa, t, z, y, x, NB ) 
		      - GRID_ENTRY( pa, t, z, y, x, ST ) - GRID_ENTRY( pa, t, z, y, x, SB );
		uz  = + GRID_ENTRY( pa, t, z, y, x, T  ) - GRID_ENTRY( pa, t, z, y, x, B  ) 
		      + GRID_ENTRY( pa, t, z, y, x, NT ) - GRID_ENTRY( pa, t, z, y, x, NB ) 
		      + GRID_ENTRY( pa, t, z, y, x, ST ) - GRID_ENTRY( pa, t, z, y, x, SB ) 
		      + GRID_ENTRY( pa, t, z, y, x, ET ) - GRID_ENTRY( pa, t, z, y, x, EB ) 
		      + GRID_ENTRY( pa, t, z, y, x, WT ) - GRID_ENTRY( pa, t, z, y, x, WB );
		ux /= rho;
		uy /= rho;
		uz /= rho;

		if( binary ) {
			loadValue( file, &fileUx );
			loadValue( file, &fileUy );
			loadValue( file, &fileUz );
		}
		else {
			if( sizeof( OUTPUT_PRECISION ) == sizeof( double )) {
				fscanf( file, "%lf %lf %lf\n", &fileUx, &fileUy, &fileUz );
			}
			else {
				fscanf( file, "%f %f %f\n", &fileUx, &fileUy, &fileUz );
			}
		}

		dUx = ux - fileUx;
		dUy = uy - fileUy;
		dUz = uz - fileUz;
		diff2 = dUx*dUx + dUy*dUy + dUz*dUz;
		if( diff2 > maxDiff2 ) maxDiff2 = diff2;
	}
		}
	}

#if defined(SPEC_CPU)
	printf( "LBM_compareVelocityField: maxDiff = %e  \n\n",
	        sqrt( maxDiff2 )  );
#else
	printf( "LBM_compareVelocityField: maxDiff = %e  ==>  %s\n\n",
	        sqrt( maxDiff2 ),
	        sqrt( maxDiff2 ) > 1e-5 ? "##### ERROR #####" : "OK" );
#endif
	fclose( file );
}

