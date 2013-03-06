/*
 * ============================================================================
 *
 *       Filename:  PSort.cpp
 *
 *    Description:  sort on pochoir's input
 *
 *        Version:  1.0
 *        Created:  02/06/2012 07:20:12 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Tang (), yuantang@csail.mit.edu
 *        Company:  CSAIL, M.I.T.
 *
 * ============================================================================
 */

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <algorithm>

using namespace std;

struct rec {
    int kernel[2]; /* kernel index */
    int tile[3]; /* tile index */
    long key; /* key for sorting */
    bool operator < (rec const & r) {
        return (key < r.key);
    }
};


int p_cmp (const void * l, const void * r) {
    return (l->key - r->key);
}

#define ARRAY_LENGTH 262710

#define READ_A 0
#define SIZE 55

int main (void) {
    rec p_array[ARRAY_LENGTH];
    int i = 0;
    char fname[10];
#if READ_A
    sprintf(fname, "a_%d.dat", SIZE);
    fprintf(stderr, "READ_A\n");
#else
    sprintf(fname, "b_%d.dat", SIZE);
    fprintf(stderr, "READ_B\n");
#endif
    FILE * fp = fopen(fname, "r");
    if (fp == NULL) {
        fprintf(stderr, "fopen error!\n");
        exit(1);
    }
    while (!feof(fp) && i < ARRAY_LENGTH) {
        rec l_rec;
#if READ_A
        fscanf(fp, "<k_%d_%d> : rec_a(%d, %d, %d)\n", &(l_rec.kernel[0]), &(l_rec.kernel[1]),
                &(l_rec.tile[0]), &(l_rec.tile[1]), &(l_rec.tile[2]));
#else
        fscanf(fp, "<k_%d_%d> : loop_b(%d, %d, %d)\n", &(l_rec.kernel[0]), &(l_rec.kernel[1]),
                &(l_rec.tile[0]), &(l_rec.tile[1]), &(l_rec.tile[2]));
#endif
        // long key = l_rec.tile[0] * 100 * 100 + l_rec.tile[1] * 100 + l_rec.tile[2];
        // p_map.insert(pair<long, rec>(key, l_rec));
        l_rec.key = l_rec.tile[0] * 1000 * 1000 + l_rec.tile[1] * 1000 + l_rec.tile[2];
        p_array[i] = l_rec;
        ++i;
    }
    fclose(fp);

    //showing contents 
    // quick2(p_map, 0, p_map.size()-1);
    // sort(p_map.begin(), p_map.end(), p_cmp);

    qsort(p_array, i, sizeof(rec), p_cmp);
    char l_str[100];
    for (int j = 0; j < i; ++j) {
        sprintf(l_str, "loop_b(%d, %d, %d) : k_%d_%d\n", 
                p_array[j].tile[0], p_array[j].tile[1], p_array[j].tile[2],
                p_array[j].kernel[0], p_array[j].kernel[1]);
        cout << l_str;
    }
    return 0;
}
