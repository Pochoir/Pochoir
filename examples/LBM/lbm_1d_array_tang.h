/* $Id: lbm_1d_array.h,v 1.1 2004/04/20 14:33:59 pohlt Exp $ */

#ifndef _LBM_MACROS_H_
#define _LBM_MACROS_H_

/*############################################################################*/

#define GRID_ENTRY(pa, t, z, y, x, e)          (pa(t, z, y, x)._##e)
// #define GRID_ENTRY_SWEEP(g,dx,dy,dz,e) ((g)[CALC_INDEX(dx, dy, dz, e)+(i)])

#define LOCAL(pa, t, z, y, x, e)    (GRID_ENTRY( pa, t, z+0, y+0, x+0, e ))
#define NEIGHBOR_C(pa, t, z, y, x)  (GRID_ENTRY( pa, t, z+0, y+0, x+0, C ))
#define NEIGHBOR_N(pa, t, z, y, x)  (GRID_ENTRY( pa, t, z+0, y+1, x+0, N ))
#define NEIGHBOR_S(pa, t, z, y, x)  (GRID_ENTRY( pa, t, z+0, y-1, x+0, S ))
#define NEIGHBOR_E(pa, t, z, y, x)  (GRID_ENTRY( pa, t, z+0, y+0, x+1, E ))
#define NEIGHBOR_W(pa, t, z, y, x)  (GRID_ENTRY( pa, t, z+0, y+0, x-1, W ))
#define NEIGHBOR_T(pa, t, z, y, x)  (GRID_ENTRY( pa, t, z+1, y+0, x+0, T ))
#define NEIGHBOR_B(pa, t, z, y, x)  (GRID_ENTRY( pa, t, z-1, y+0, x+0, B ))
#define NEIGHBOR_NE(pa, t, z, y, x) (GRID_ENTRY( pa, t, z+0, y+1, x+1, NE ))
#define NEIGHBOR_NW(pa, t, z, y, x) (GRID_ENTRY( pa, t, z+0, y+1, x-1, NW ))
#define NEIGHBOR_SE(pa, t, z, y, x) (GRID_ENTRY( pa, t, z+0, y-1, x+1, SE ))
#define NEIGHBOR_SW(pa, t, z, y, x) (GRID_ENTRY( pa, t, z+0, y-1, x-1, SW ))
#define NEIGHBOR_NT(pa, t, z, y, x) (GRID_ENTRY( pa, t, z+1, y+1, x+0, NT ))
#define NEIGHBOR_NB(pa, t, z, y, x) (GRID_ENTRY( pa, t, z-1, y+1, x+0, NB ))
#define NEIGHBOR_ST(pa, t, z, y, x) (GRID_ENTRY( pa, t, z+1, y-1, x+0, ST ))
#define NEIGHBOR_SB(pa, t, z, y, x) (GRID_ENTRY( pa, t, z-1, y-1, x+0, SB ))
#define NEIGHBOR_ET(pa, t, z, y, x) (GRID_ENTRY( pa, t, z+1, y+0, x+1, ET ))
#define NEIGHBOR_EB(pa, t, z, y, x) (GRID_ENTRY( pa, t, z-1, y+0, x+1, EB ))
#define NEIGHBOR_WT(pa, t, z, y, x) (GRID_ENTRY( pa, t, z+1, y+0, x-1, WT ))
#define NEIGHBOR_WB(pa, t, z, y, x) (GRID_ENTRY( pa, t, z-1, y+0, x-1, WB ))

#define COLLIDE_STREAM
// #undef COLLIDE_STREAM
#ifdef COLLIDE_STREAM

#define SRC_C(pa, t, z, y, x)  (LOCAL( pa, t, z, y, x, C  ))
#define SRC_N(pa, t, z, y, x)  (LOCAL( pa, t, z, y, x, N  ))
#define SRC_S(pa, t, z, y, x)  (LOCAL( pa, t, z, y, x, S  ))
#define SRC_E(pa, t, z, y, x)  (LOCAL( pa, t, z, y, x, E  ))
#define SRC_W(pa, t, z, y, x)  (LOCAL( pa, t, z, y, x, W  ))
#define SRC_T(pa, t, z, y, x)  (LOCAL( pa, t, z, y, x, T  ))
#define SRC_B(pa, t, z, y, x)  (LOCAL( pa, t, z, y, x, B  ))
#define SRC_NE(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, NE ))
#define SRC_NW(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, NW ))
#define SRC_SE(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, SE ))
#define SRC_SW(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, SW ))
#define SRC_NT(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, NT ))
#define SRC_NB(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, NB ))
#define SRC_ST(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, ST ))
#define SRC_SB(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, SB ))
#define SRC_ET(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, ET ))
#define SRC_EB(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, EB ))
#define SRC_WT(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, WT ))
#define SRC_WB(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, WB ))

#define DST_C(pa, t, z, y, x)  (NEIGHBOR_C ( pa, t, z, y, x ))
#define DST_N(pa, t, z, y, x)  (NEIGHBOR_N ( pa, t, z, y, x ))
#define DST_S(pa, t, z, y, x)  (NEIGHBOR_S ( pa, t, z, y, x ))
#define DST_E(pa, t, z, y, x)  (NEIGHBOR_E ( pa, t, z, y, x ))
#define DST_W(pa, t, z, y, x)  (NEIGHBOR_W ( pa, t, z, y, x ))
#define DST_T(pa, t, z, y, x)  (NEIGHBOR_T ( pa, t, z, y, x ))
#define DST_B(pa, t, z, y, x)  (NEIGHBOR_B ( pa, t, z, y, x ))
#define DST_NE(pa, t, z, y, x) (NEIGHBOR_NE( pa, t, z, y, x ))
#define DST_NW(pa, t, z, y, x) (NEIGHBOR_NW( pa, t, z, y, x ))
#define DST_SE(pa, t, z, y, x) (NEIGHBOR_SE( pa, t, z, y, x ))
#define DST_SW(pa, t, z, y, x) (NEIGHBOR_SW( pa, t, z, y, x ))
#define DST_NT(pa, t, z, y, x) (NEIGHBOR_NT( pa, t, z, y, x ))
#define DST_NB(pa, t, z, y, x) (NEIGHBOR_NB( pa, t, z, y, x ))
#define DST_ST(pa, t, z, y, x) (NEIGHBOR_ST( pa, t, z, y, x ))
#define DST_SB(pa, t, z, y, x) (NEIGHBOR_SB( pa, t, z, y, x ))
#define DST_ET(pa, t, z, y, x) (NEIGHBOR_ET( pa, t, z, y, x ))
#define DST_EB(pa, t, z, y, x) (NEIGHBOR_EB( pa, t, z, y, x ))
#define DST_WT(pa, t, z, y, x) (NEIGHBOR_WT( pa, t, z, y, x ))
#define DST_WB(pa, t, z, y, x) (NEIGHBOR_WB( pa, t, z, y, x ))

#else /* COLLIDE_STREAM */

/* In following macros, the field in left and the field in right
 * doesn't match up with each other, but it's the case in original
 * code, so we just preserve it!!
 */
#define SRC_C(pa, t, z, y, x)  (NEIGHBOR_C ( pa, t, z, y, x ))
#define SRC_N(pa, t, z, y, x)  (NEIGHBOR_S ( pa, t, z, y, x ))
#define SRC_S(pa, t, z, y, x)  (NEIGHBOR_N ( pa, t, z, y, x ))
#define SRC_E(pa, t, z, y, x)  (NEIGHBOR_W ( pa, t, z, y, x ))
#define SRC_W(pa, t, z, y, x)  (NEIGHBOR_E ( pa, t, z, y, x ))
#define SRC_T(pa, t, z, y, x)  (NEIGHBOR_B ( pa, t, z, y, x ))
#define SRC_B(pa, t, z, y, x)  (NEIGHBOR_T ( pa, t, z, y, x ))
#define SRC_NE(pa, t, z, y, x) (NEIGHBOR_SW( pa, t, z, y, x ))
#define SRC_NW(pa, t, z, y, x) (NEIGHBOR_SE( pa, t, z, y, x ))
#define SRC_SE(pa, t, z, y, x) (NEIGHBOR_NW( pa, t, z, y, x ))
#define SRC_SW(pa, t, z, y, x) (NEIGHBOR_NE( pa, t, z, y, x ))
#define SRC_NT(pa, t, z, y, x) (NEIGHBOR_SB( pa, t, z, y, x ))
#define SRC_NB(pa, t, z, y, x) (NEIGHBOR_ST( pa, t, z, y, x ))
#define SRC_ST(pa, t, z, y, x) (NEIGHBOR_NB( pa, t, z, y, x ))
#define SRC_SB(pa, t, z, y, x) (NEIGHBOR_NT( pa, t, z, y, x ))
#define SRC_ET(pa, t, z, y, x) (NEIGHBOR_WB( pa, t, z, y, x ))
#define SRC_EB(pa, t, z, y, x) (NEIGHBOR_WT( pa, t, z, y, x ))
#define SRC_WT(pa, t, z, y, x) (NEIGHBOR_EB( pa, t, z, y, x ))
#define SRC_WB(pa, t, z, y, x) (NEIGHBOR_ET( pa, t, z, y, x ))

#define DST_C(pa, t, z, y, x)  (LOCAL( pa, t, z, y, x, C  ))
#define DST_N(pa, t, z, y, x)  (LOCAL( pa, t, z, y, x, N  ))
#define DST_S(pa, t, z, y, x)  (LOCAL( pa, t, z, y, x, S  ))
#define DST_E(pa, t, z, y, x)  (LOCAL( pa, t, z, y, x, E  ))
#define DST_W(pa, t, z, y, x)  (LOCAL( pa, t, z, y, x, W  ))
#define DST_T(pa, t, z, y, x)  (LOCAL( pa, t, z, y, x, T  ))
#define DST_B(pa, t, z, y, x)  (LOCAL( pa, t, z, y, x, B  ))
#define DST_NE(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, NE ))
#define DST_NW(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, NW ))
#define DST_SE(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, SE ))
#define DST_SW(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, SW ))
#define DST_NT(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, NT ))
#define DST_NB(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, NB ))
#define DST_ST(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, ST ))
#define DST_SB(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, SB ))
#define DST_ET(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, ET ))
#define DST_EB(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, EB ))
#define DST_WT(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, WT ))
#define DST_WB(pa, t, z, y, x) (LOCAL( pa, t, z, y, x, WB ))

#endif /* COLLIDE_STREAM */

#define MAGIC_CAST(v) ((unsigned int*) ((void*) (&(v))))
#define FLAG_VAR(v) unsigned int* const _aux_ = MAGIC_CAST(v)

#if 1
#define TEST_FLAG(v,f)     ((v) & (f))
#define SET_FLAG(v,f)      {(v) |= (f);}
#define CLEAR_FLAG(v,f)    {(v) &= ~(f);}
#define CLEAR_ALL_FLAGS(v) {(v)  = 0;}

#define TEST_FLAG_SWEEP(pa, t, z, y, x, f)     ((pa(t, z, y, x)._FLAGS) & (f))
#define SET_FLAG_SWEEP(pa, t, z, y, x, f)      {(pa(t, z, y, x)._FLAGS) |= (f);}
#define CLEAR_FLAG_SWEEP(pa, t, z, y, x, f)    {(pa(t, z, y, x)._FLAGS) &= ~(f);}
#define CLEAR_ALL_FLAGS_SWEEP(pa, t, z, y, x) {(pa(t, z, y, x)._FLAGS) = 0;}
#else
#define TEST_FLAG(v,f)     ((*MAGIC_CAST(v)) & (f))
#define SET_FLAG(v,f)      {FLAG_VAR(v); (*_aux_) |=  (f);}
#define CLEAR_FLAG(v,f)    {FLAG_VAR(v); (*_aux_) &= ~(f);}
#define CLEAR_ALL_FLAGS(v) {FLAG_VAR(v); (*_aux_)  =    0;}

#define TEST_FLAG_SWEEP(pa, t, z, y, x, f)     ((*MAGIC_CAST(pa(t, z, y, x)._FLAGS) & (f)))
#define SET_FLAG_SWEEP(pa, t, z, y, x, f)      {FLAG_VAR(pa(t, z, y, x)._FLAGS); (*_aux_) |=  (f);}
#define CLEAR_FLAG_SWEEP(pa, t, z, y, x, f)    {FLAG_VAR(pa(t, z, y, x)._FLAGS); (*_aux_) &= ~(f);}
#define CLEAR_ALL_FLAGS_SWEEP(pa, t, z, y, x) {FLAG_VAR(pa(t, z, y, x)._FLAGS); (*_aux_)  =    0;}
#endif

/*############################################################################*/

#endif /* _LBM_MACROS_H_ */
