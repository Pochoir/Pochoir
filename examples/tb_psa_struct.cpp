/*
 **********************************************************************************
 *  Copyright (C) 2010-2011  Massachusetts Institute of Technology
 *  Copyright (C) 2010-2011  Yuan Tang <yuantang@csail.mit.edu>
 * 		                     Charles E. Leiserson <cel@mit.edu>
 * 	 
 *  Written by: Rezaul Alam Chowdhury <rezaul@mit.edu>
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *   Suggestsions:                  yuantang@csail.mit.edu
 *   Bugs:                          yuantang@csail.mit.edu
 *
 *********************************************************************************
 */

/* Test bench - Gotoh's algorithm for pairwise sequence alignment with affine gap costs */
#include <cstdio>
#include <cstddef>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include <cmath>

#include <pochoir.hpp>

using namespace std;
#define SIMPLE 0
#define N_RANK 1

#ifndef DEFAULT_GAP_OPEN_COST
  #define DEFAULT_GAP_OPEN_COST 2
#endif

#ifndef DEFAULT_GAP_EXTENSION_COST
  #define DEFAULT_GAP_EXTENSION_COST 1
#endif

#ifndef DEFAULT_MISMATCH_COST
  #define DEFAULT_MISMATCH_COST 1
#endif

#define MAX_ID_LENGTH 100

#ifndef MAX_SEQ_LENGTH
  #define MAX_SEQ_LENGTH  100000
#endif


typedef struct
{
   int vD;
   int vI;
   int vG;
} NODE;


enum err_msgs{ SEQUENCE_READ, LENGTH_READ, NO_SEQUENCE, SEQUENCE_TOO_LONG, INVALID_SEQUENCE_LENGTH, FILE_OPEN_ERROR, MEM_ALLOC_ERROR };


int read_next_seq( char *fn, char **sq, int *len )
{
  FILE *fp;

  if ( fn != NULL )
    {
     if ( ( fp = fopen( fn, "rt" ) ) == NULL ) return FILE_OPEN_ERROR;
    }
  else fp = stdin;

  if ( fscanf( fp, "%d", len ) != 1 ) return NO_SEQUENCE;
  
  ( *sq ) = ( char * ) malloc( ( *len + 2 ) * sizeof( char ) );
  
  if ( *sq == NULL ) return MEM_ALLOC_ERROR;

  int i = 1, c;
  
  while ( ( i <= *len ) && ( ( c = fgetc( fp ) ) != NULL ) )
    {
      ( *sq )[ i++ ] = c;
    } 
  
  ( *sq )[ i ] = 0;

  if ( fn != NULL ) fclose( fp );

  return SEQUENCE_READ;
}


int read_sequences( char *fnX, int *nx, char **X, char *fnY, int *ny, char **Y )
{
  int l;

  for ( int i = 0; i < 2; i++ )
    {      
     if ( i == 0 ) l = read_next_seq( fnX, X, nx );
     else l = read_next_seq( fnY, Y, ny );

     if ( l == FILE_OPEN_ERROR )
       {
        printf( "Error: Unable to open input file ( %s )!\n\n", ( i == 0 ) ? fnX : fnY );
        return 0;
       }

     if ( l == MEM_ALLOC_ERROR )
       {
        printf( "Error: Memory allocation error!\n\n" );
        return 0;
       }

     if ( l == NO_SEQUENCE )
       {
        printf( "Error: Failed to read sequence %d!\n\n", i + 1 );
        return 0;
       }
    }

  return 1;
}



int initRandSeq( int nX, char **X, int nY, char **Y )
{
    ( *X ) = ( char * ) malloc( ( nX + 2 ) * sizeof( char ) );
    ( *Y ) = ( char * ) malloc( ( nY + 2 ) * sizeof( char ) );    
  
    if ( ( *X == NULL ) || ( *Y == NULL ) ) 
      {
        printf( "Error: Memory allocation error!\n\n" );      
        return 0;
      }   

    char *SYM = "ACGT";    
    int nS = 4;//sizeof( SYM ) / sizeof( SYM[ 0 ] );

    srand( time( NULL ) );

    for ( int i = 1; i <= nX; ++i )
       ( *X )[ i ] = SYM[ rand( ) % nS ];

    for ( int i = 1; i <= nY; ++i )
       ( *Y )[ i ] = SYM[ rand( ) % nS ];
       
    ( *X )[ nX + 1 ] = ( *Y )[ nY + 1 ] = 0;       
    
    return 1;
}



void print_usage( char *prog )
{
  printf( "Usage: %s [ options ]\n\n", prog );

  printf( "Options:\n" );

  printf( "\t-o value         : gap open cost ( nonnegative integer, default = %d )\n", DEFAULT_GAP_OPEN_COST );
  printf( "\t-e value         : gap extention cost ( nonnegative integer, default = %d )\n", DEFAULT_GAP_EXTENSION_COST );
  printf( "\t-m value         : mismatch cost ( nonnegative integer, default = %d )\n\n", DEFAULT_MISMATCH_COST );

  printf( "\t-s fname1 fname2 : input files containing the two seqeuneces (file format: sequence length in line 1 followed by the sequence in line 2)\n" );
  printf( "\t-r len1 len2     : generate two random sequences of length len1 and len2\n\n" );

  printf( "\t-d               : Run standard DP algorithm\n" );
  printf( "\t-i               : Run iterative stencil\n\n" );  
   
  printf( "\t-h               : print this help screen\n\n" );
}



int read_command_line( int argc, char *argv[ ], int &go_cost, int &ge_cost, int &mm_cost, int &nx, char **fnX, int &ny, char **fnY, 
                       int &Run_standard_dp, int &Run_iter_stencil )
{
  go_cost = DEFAULT_GAP_OPEN_COST;
  ge_cost = DEFAULT_GAP_EXTENSION_COST;
  mm_cost = DEFAULT_MISMATCH_COST;
  nx = ny = 0;
  Run_standard_dp = 0;
  Run_iter_stencil = 0;

  for ( int i = 1; i < argc; )
    {
     int j = i;

     if ( !strcasecmp( argv[ i ], "-d" ) )
       {
        Run_standard_dp = 1;
         
        i++;
        
        if ( i >= argc ) break;
       }
     

     if ( !strcasecmp( argv[ i ], "-i" ) )
       {
        Run_iter_stencil = 1;
         
        i++;
        
        if ( i >= argc ) break;
       }

     
     if ( !strcasecmp( argv[ i ], "-o" ) )
       {
        if ( i + 1 >= argc )
          {
           printf( "Error: Missing gap open cost ( specify -o gap-open-cost )!\n\n" );
           return 0;
          }

        go_cost = atoi( argv[ i + 1 ] );

        if ( go_cost < 0 )
          {
           printf( "Error: Specify non-negative gap introduction cost!\n\n" );
           return 0;
          }

        i += 2;
        
        if ( i >= argc ) break;
       }

       
     if ( !strcasecmp( argv[ i ], "-e" ) )
       {
        if ( i + 1 >= argc )
          {
           printf( "Error: Missing gap extention cost ( specify -e gap-extention-cost )!\n\n" );
           return 0;
          }

        ge_cost = atoi( argv[ i + 1 ] );

        if ( ge_cost < 0 )
          {
           printf( "Error: Specify non-negative gap extention cost!\n\n" );
           return 0;
          }

        i += 2;
        
        if ( i >= argc ) break;
       }
     
     
     if ( !strcasecmp( argv[ i ], "-m" ) )
       {
        if ( i + 1 >= argc )
          {
           printf( "Error: Missing mismatch cost ( specify -m mismatch-cost )!\n\n" );
           return 0;
          }

        mm_cost = atoi( argv[ i + 1 ] );

        if ( mm_cost < 0 )
          {
           printf( "Error: Specify non-negative mismatch cost!\n\n" );
           return 0;
          }

        i += 2;

        if ( i >= argc ) break;
       }

     if ( !strcasecmp( argv[ i ], "-r" ) )
       {
        if ( ( i + 2 >= argc ) || ( argv[ i + 1 ][ 0 ] == '-' ) || ( argv[ i + 2 ][ 0 ] == '-' ) )
          {
           printf( "Error: Missing sequence lengths ( specify -r length-sequence-1 length-sequence-2 )!\n\n" );
           return 0;
          }

        nx = atoi( argv[ i + 1 ] );
        ny = atoi( argv[ i + 2 ] );        

        if ( ( nx < 1 ) || ( ny < 1 ) )
          {
           printf( "Error: Specify positive sequence lengths!\n\n" );
           return 0;
          }
          
        if ( ( nx > MAX_SEQ_LENGTH ) || ( ny > MAX_SEQ_LENGTH ) )
          {
           printf( "Error: Cannot generate sequences longer than %d!\n\n", MAX_SEQ_LENGTH );
           return 0;
          }          

        i += 3;

        if ( i >= argc ) break;
       }
       
     if ( !strcasecmp( argv[ i ], "-s" ) )
       {
        if ( ( i + 2 >= argc ) || ( argv[ i + 1 ][ 0 ] == '-' ) || ( argv[ i + 2 ][ 0 ] == '-' ) )
          {
           printf( "Error: Missing input file names ( specify -s sequence-file-1 sequence-file-2 )!\n\n" );
           return 0;
          }

        ( *fnX ) = strdup( argv[ i + 1 ] );
        ( *fnY ) = strdup( argv[ i + 2 ] );  

        i += 3;
        
        if ( i >= argc ) break;
       }
                     
     if ( !strcasecmp( argv[ i ], "-h" ) || !strcasecmp( argv[ i ], "-help" ) || !strcasecmp( argv[ i ], "--help" ) )
       {
        print_usage( argv[ 0 ] );
        exit( 0 );
        i++;

        if ( i >= argc ) break;
       }


     if ( i == j )
       {
        printf( "Error: Unknown option ( %s )!\n\n", argv[ i ] );
        return 0;
       }
    }

   printf( "Paremeters:\n\n" );
   printf( "\tgap open cost = %d\n", go_cost );
   printf( "\tgap extention cost = %d\n", ge_cost );
   printf( "\tmismatch cost = %d\n\n", mm_cost );
    
   return 1;
}



int stencilPSA4Arrays( int nX, char *X, int nY, char *Y, int goCost, int geCost, int *mmCost )
{
    Pochoir_Shape< N_RANK > pSeq_shape_G[ ] = { { 1, 0 }, { 0, 0 }, { 0, -1 } };    
    Pochoir_Shape< N_RANK > pSeq_shape_G2[ ] = { { 1, 0 }, { 0, 0 } };    
    Pochoir_Shape< N_RANK > pSeq_shape_D[ ] = { { 1, 0 }, { 0, 0 } };    
    Pochoir_Shape< N_RANK > pSeq_shape_I[ ] = { { 1, 0 }, { 0, -1 } };    
    Pochoir< N_RANK > pSeq(pSeq_shape_G);
    Pochoir_Array< int, N_RANK > vG( nY + 1 ), vG2( nY + 1 ), vD( nY + 1 ), vI( nY + 1 );
    Pochoir_Domain J( 0, nY + 1 );
                
    pSeq.Register_Array( vG );
    pSeq.Register_Array( vG2 );
    pSeq.Register_Array( vD );
    pSeq.Register_Array( vI );            
    pSeq.Register_Domain( J );
    
    vG( 0, 0 ) = vG2( 0, 0 ) = 0;
    
    Pochoir_Kernel_1D( pSeq_fn, t, j )

       register int i = t + 1 - j, v;
                  
       if ( ( i >= 0 ) && ( i <= nX ) )
         {             
           if ( ( i > 0 ) && ( j > 0 ) )
             {
                vD( t + 1, j ) = min( vD( t, j ), vG( t, j ) + goCost ) + geCost;
                vI( t + 1, j ) = min( vI( t, j - 1 ), vG( t, j - 1 ) + goCost ) + geCost;                
  
                v = min( vD( t + 1, j ), vI( t + 1, j ) );
  
                vG( t + 1, j ) = min( v, vG2( t, j - 1 ) + mmCost[ X[ i ] == Y[ j ] ] );                                     
             }  
           else
             {
                v = goCost + ( i + j ) * geCost;
                
                vG( t + 1, j ) = v++;                                     
                
                if ( !i ) vD( t + 1, j ) = v;
                if ( !j ) vI( t + 1, j ) = v;                                
             }
             
           vG2( t + 1, j ) = vG( t, j );                    
         }
                      	      
    Pochoir_Kernel_End

    int t = nX + nY;

    pSeq.Run( t, pSeq_fn );
    
    int ug = vG.interior( t, nY ), ud = vD.interior( t, nY ), ui = vI.interior( t, nY );
    int optCost = min( ug, min( ud, ui ) );

    return optCost;    
}


int stencilPSA( int nX, char *X, int nY, char *Y, int goCost, int geCost, int *mmCost )
{
    Pochoir_Shape< N_RANK > pSeq_shape_G[ ] = { { 2, 0 }, { 1, 0 }, { 0, -1 }, { 1, -1 } };
    Pochoir< N_RANK > pSeq(pSeq_shape_G); 
    Pochoir_Array< int, N_RANK > vG( nY + 1 ), vD( nY + 1 ), vI( nY + 1 );
    pSeq.Register_Array( vG );
    pSeq.Register_Array( vD );
    pSeq.Register_Array( vI );            

    vG( 0, 0 ) = vG( 1, 0 ) = 0;
    
    Pochoir_Kernel_1D( pSeq_fn, t, j )
      
       int i = t + 1 - j, v, k;

       if ( ( i >= 0 ) && ( i <= nX ) )
         {                      
           if ( ( i > 0 ) && ( j > 0 ) )
             {
                vD( t + 2, j ) = min( vD( t + 1, j ), vG( t + 1, j ) + goCost ) + geCost;
                vI( t + 2, j ) = min( vI( t + 1, j - 1 ), vG( t + 1, j - 1 ) + goCost ) + geCost;                
                v = min( vD( t + 2, j ), vI( t + 2, j ) );
                vG( t + 2, j ) = min( v, vG( t, j - 1 ) + mmCost[ X[ i ] == Y[ j ] ] );
             }  
           else
             {
                v = goCost + ( i + j ) * geCost;
                
                vG( t + 2, j ) = v++;                                     
                
                if ( !i ) vD( t + 2, j ) = v;                
                if ( !j ) vI( t + 2, j ) = v;                                
             }
         }
                      	      
    Pochoir_Kernel_End

    int t = nX + nY;

    pSeq.Run( t, pSeq_fn );
    
    int optCost = min( vG.interior( t + 1, nY ), min( vD.interior( t + 1, nY ), vI.interior( t + 1, nY ) ) );

    return optCost;    
}

Pochoir_Shape< N_RANK > pSeq_shape[ ] = { { 2, 0 }, { 1, 0 }, { 0, -1 }, { 1, -1 } };    

int stencilPSAStruct( int nX, char *X, int nY, char *Y, int goCost, int geCost, int *mmCost )
{
    Pochoir_Shape< N_RANK > pSeq_shape[ ] = { { 2, 0 }, { 1, 0 }, { 0, -1 }, { 1, -1 } };    
    Pochoir< N_RANK > pSeq(pSeq_shape);    
    Pochoir_Array< NODE, N_RANK > pArray( nY + 1 );    
    pSeq.Register_Array( pArray );
    
    pArray( 0, 0 ).vG = pArray( 1, 0 ).vG = 0;
    
    for ( int i = 0; i < nY; ++i )
      {
//         NODE(pArray( 0, i )).vG = i;
        pArray( 0, i ).vG = i;
//        int v = NODE( pArray( 0, i ) ).vG;
        int v = NODE( pArray( 0, i ) ).vG;
//        printf( "v = %d\n", v );
      }
    
    Pochoir_Kernel_1D( pSeq_fn, t, j )
      
       int i = t + 1 - j, v, k;

       if ( ( i >= 0 ) && ( i <= nX ) )
         {                      
           if ( ( i > 0 ) && ( j > 0 ) )
             {             
                int D = pArray( t + 1, j ).vD, GD = pArray( t + 1, j ).vG;
                int I =  pArray( t + 1, j - 1 ).vI, GI =  pArray( t + 1, j - 1 ).vG;                

                int D1 = min( D, GD + goCost ) + geCost;
                int I1 = min( I, GI + goCost ) + geCost;                

                v = min( D1, I1 );
                int G =  pArray( t, j - 1 ).vG; 
                int G1 = min( v, G + mmCost[ X[ i ] == Y[ j ] ] );

//                printf( "D = %d, GD = %d, I = %d, GI = %d, G = %d, D1 = %d, I1 = %d, G1 = %d\n", D, GD, I, GI, G, D1, I1, G1 );  
                
                 pArray( t + 2, j ).vD = D1;
                 pArray( t + 2, j ).vI = I1;                
                 pArray( t + 2, j ).vG = G1;                

//                 pArray( t + 2, j ) ).D = min(  pArray( t + 1, j ) ).D,  pArray( t + 1, j ) ).G + goCost ) + geCost;
//                 pArray( t + 2, j ) ).I = min(  pArray( t + 1, j - 1 ) ).I,  pArray( t + 1, j - 1 ) ).G + goCost ) + geCost;                
//                v = min(  pArray( t + 2, j ) ).D,  pArray( t + 2, j ) ).I );
//                printf( "t = %d, j = %d, v = %d\n", t, j, v );                
//                 pArray( t + 2, j ) ).G = min( v,  pArray( t, j - 1 ) ).G + mmCost[ X[ i ] == Y[ j ] ] );
             }  
           else
             {
                v = goCost + ( i + j ) * geCost;
  
//                printf( "t = %d, j = %d, v = %d\n", t, j, v );                
                
                 pArray( t + 2, j ).vG = v++;                                     
                
                int w =  pArray( t + 2, j ).vG;
                
//                printf( " pArray( t + 2, j ) ).G = %d\n", w );
                
                if ( !i )  pArray( t + 2, j ).vD = v;                
                if ( !j )  pArray( t + 2, j ).vI = v;                                
             }
         }
                      	      
    Pochoir_Kernel_End

    int t = nX + nY;

    pSeq.Run( t, pSeq_fn );
    
    int optCost = min( NODE( pArray( t + 1, nY ) ).vG, min( NODE( pArray( t + 1, nY ) ).vD, NODE( pArray( t + 1, nY ) ).vI ) );

    return optCost;    
}


int standardDPPSA( int nX, char *X, int nY, char *Y, int goCost, int geCost, int *mmCost )
{
    Pochoir_Array< int, N_RANK > vG( nY + 1 ), vD( nY + 1 ), vI( nY + 1 );
    vG.Register_Shape(pSeq_shape);
    vD.Register_Shape(pSeq_shape);
    vI.Register_Shape(pSeq_shape);

    vG.interior( 0, 0 ) = vD.interior( 0, 0 ) = vI.interior( 0, 0 ) = 0;        
    
    for ( int j = 1; j <= nY; ++j )
      {
        vD.interior( 0, j ) = goCost + ( j + 1 ) * geCost;
        vI.interior( 0, j ) = 0;
        vG.interior( 0, j ) = goCost + j * geCost;
      }
    
    for ( int i = 1; i <= nX; ++i ) 
      {
        vD.interior( i, 0 ) = 0;
        vI.interior( i, 0 ) = goCost + ( i + 1 ) * geCost;
        vG.interior( i, 0 ) = goCost + i * geCost;        

        for ( int j = 1; j <= nY; ++j ) 
          {
            vD.interior( i, j ) = min( vD.interior( i - 1, j ), vG.interior( i - 1, j ) + goCost ) + geCost;
            vI.interior( i, j ) = min( vI.interior( i, j - 1 ), vG.interior( i, j - 1 ) + goCost ) + geCost;

            int mDI = min( vD.interior( i, j ), vI.interior( i, j ) );

            vG.interior( i, j ) = min( mDI, vG.interior( i - 1, j - 1 ) + mmCost[ X[ i ] == Y[ j ] ] );
          }          
      }        

    int optCost = min( vG.interior( nX, nY ), min( vD.interior( nX, nY ), vI.interior( nX, nY ) ) );    
    
    return optCost;    
}



int iterativeStencilPSA( int nX, char *X, int nY, char *Y, int goCost, int geCost, int *mmCost )
{
    Pochoir_Array< int, N_RANK > vG( nY + 1 ), vD( nY + 1 ), vI( nY + 1 );
    vG.Register_Shape(pSeq_shape);
    vD.Register_Shape(pSeq_shape);
    vI.Register_Shape(pSeq_shape);
    
    vG( 0, 0 ) = vG( 1, 0 ) = 0;
    
    for ( int t = 0; t < nX + nY; ++t )
       cilk_for ( int j = 0; j <= nY; ++j )
          {
            int i = t + 1 - j;
            
            if ( ( i >= 0 ) && ( i <= nX ) )
              {            
                if ( ( i > 0 ) && ( j > 0 ) )
                   {
                    vD.interior( t + 2, j ) = min( vD.interior( t + 1, j ), vG.interior( t + 1, j ) + goCost ) + geCost;
                    vI.interior( t + 2, j ) = min( vI.interior( t + 1, j - 1 ), vG.interior( t + 1, j - 1 ) + goCost ) + geCost;                
                    int v = min( vD.interior( t + 2, j ), vI.interior( t + 2, j ) );
                    vG.interior( t + 2, j ) = min( v, vG.interior( t, j - 1 ) + mmCost[ X[ i ] == Y[ j ] ] );
                  }  
                else
                  {
                    int v = goCost + ( i + j ) * geCost;
                    
                    vG.interior( t + 2, j ) = v++;                                     
                    
                    if ( !i ) vD.interior( t + 2, j ) = v;                
                    if ( !j ) vI.interior( t + 2, j ) = v;                                
                  }
             } 
          }
          
    int t = nX + nY + 1;

    int optCost = min( vG.interior( t, nY ), min( vD.interior( t, nY ), vI.interior( t, nY ) ) );

    return optCost;    
}



int main( int argc, char *argv[ ] )
{
    printf( "\nStencil-based DP for optimal pairwise sequence alignment cost with affine gap costs ( Run with option -h for help ).\n\n" );

    int nX, nY;
    char *X = NULL, *Y = NULL;    
    char *fnX = NULL, *fnY = NULL;
    int goCost, geCost, mmCost[ ] = { 0, 0 };
    int RunStandardDP, RunIterativeStencil;

    if ( !read_command_line( argc, argv, goCost, geCost, mmCost[ 0 ], nX, &fnX, nY, &fnY, RunStandardDP, RunIterativeStencil ) )
      {
        print_usage( argv[ 0 ] );
        return 1;
      }

    int gotSeq = 0;

    if ( ( nX < 1 ) || ( nY < 1 ) ) 
      {
       if ( ( fnX == NULL ) || ( fnY == NULL ) )
         {
           printf( "Error: Missing input file names ( specify -s sequence-file-1 sequence-file-2 )!\n\n" );
           print_usage( argv[ 0 ] );
           return 1;         
         }
          
       gotSeq = read_sequences( fnX, &nX, &X, fnY, &nY, &Y );
      }
    else gotSeq = initRandSeq( nX, &X, nY, &Y );

    if ( gotSeq )
      {
        printf( "Sequence lengths = < %d, %d >\n\n", nX, nY );
      
        struct timeval start, end;

        printf( "Running pochoir-based DP ( with struct )..." );
        fflush( stdout );
                   
        gettimeofday( &start, 0 );        
        int optCost = stencilPSAStruct( nX, X, nY, Y, goCost, geCost, mmCost );    
        gettimeofday( &end, 0 );

        double t0 = tdiff( &end, &start );
              
        printf( "\n\nPochoir ( with struct ):\n" );
        printf( "\t alignment cost = %d\n", optCost );    
        //printf( "\t Running time = %.3lf sec\n\n", t0 );    
        printf( "\t Running time = %.3lf ms\n\n", 1.0e3 * t0 );    

        /*printf( "Running pochoir-based DP ( without struct )..." );
        fflush( stdout );
                      
        gettimeofday( &start, 0 );        
        int optCost2 = stencilPSA( nX, X, nY, Y, goCost, geCost, mmCost );    
        gettimeofday( &end, 0 );

        double t1 = tdiff( &end, &start );
      
        printf( "\n\nPochoir ( without struct ):\n" );
        printf( "\t alignment cost = %d\n", optCost2 );    
        if ( t0 > 0 ) printf( "\t Running time = %.3lf sec ( %.3lf x Pochoir-Struct )\n\n", t1, t1 / t0 );    
        else printf( "\t Running time = %.3lf sec\n\n", t1 );   */ 
      
        if ( RunIterativeStencil )
          {
            printf( "Running iterative stencil..." );
            fflush( stdout );
                          
            gettimeofday( &start, 0 );        
            int optCostITST = iterativeStencilPSA( nX, X, nY, Y, goCost, geCost, mmCost );    
            gettimeofday( &end, 0 );

            double t1 = tdiff( &end, &start );
          
            printf( "\n\nIterative Stencil:\n" );
            printf( "\t alignment cost = %d\n", optCostITST );    
            if ( t0 > 0 ) printf( "\t Running time = %.3lf sec ( %.3lf x Pochoir-Struct )\n\n", t1, t1 / t0 );    
            else printf( "\t Running time = %.3lf sec\n\n", t1 );    
          }

        if ( RunStandardDP )
          {
            printf( "Running standard DP..." );
            fflush( stdout );
                          
            gettimeofday( &start, 0 );        
            int optCostSDP = standardDPPSA( nX, X, nY, Y, goCost, geCost, mmCost );    
            gettimeofday( &end, 0 );
            
            double t1 = tdiff( &end, &start );
          
            printf( "\n\nStandard DP:\n" );
            printf( "\t alignment cost = %d\n", optCostSDP );    
            if ( t0 > 0 ) printf( "\t Running time = %.3lf sec ( %.3lf x Pochoir )\n\n", t1, t1 / t0 );    
            else printf( "\t Running time = %.3lf sec\n\n", t1 );    
          }
      }

    if ( X != NULL ) free( X );
    if ( Y != NULL ) free( Y );    
       
    return 0;
}
