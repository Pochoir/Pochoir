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

/* Test bench - length of a longest common subsequence between two given sequences */
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

#ifndef MAX_SEQ_LENGTH
  #define MAX_SEQ_LENGTH  10000000
#endif

enum err_msgs{ SEQUENCE_READ, LENGTH_READ, NO_SEQUENCE, SEQUENCE_TOO_LONG, INVALID_SEQUENCE_LENGTH, FILE_OPEN_ERROR, MEM_ALLOC_ERROR };

Pochoir_Shape< N_RANK > LCS_shape[ ] = { { 2, 0 }, { 1, 0 }, { 0, -1 }, { 1, -1 } };    

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

    char *SYM = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";    
    int nS = 26;//sizeof( SYM ) / sizeof( SYM[ 0 ] );

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

  printf( "\t-s fname1 fname2 : input files containing the two seqeuneces (file format: sequence length in line 1 followed by the sequence in line 2)\n" );
  printf( "\t-r len1 len2     : generate two random sequences of length len1 and len2\n\n" );

  printf( "\t-d               : Run standard DP algorithm\n" );
  printf( "\t-i               : Run iterative stencil\n\n" );  
   
  printf( "\t-h               : print this help screen\n\n" );
}



int read_command_line( int argc, char *argv[ ], int &nx, char **fnX, int &ny, char **fnY, 
                       int &Run_standard_dp, int &Run_iter_stencil )
{
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

   return 1;
}



int stencilLCS( int nX, char *X, int nY, char *Y )
{
    Pochoir< N_RANK > LCS(LCS_shape);    
    Pochoir_Array< int, N_RANK > L( nY + 1 );
    LCS.Register_Array(L);
    Pochoir_Domain J( 0, nY + 1 );
    
    L( 0, 0 ) = L( 1, 0 ) = L( 1, 1 ) = 0;
    
    Pochoir_Kernel_1D( LCS_fn, t, j )
      
       int i = t + 1 - j, v, k;

       if ( ( i >= 0 ) && ( i <= nX ) )
         {                      
           if ( ( i > 0 ) && ( j > 0 ) )
             {
                if ( X[ i ] == Y[ j ] ) L( t + 2, j ) = 1 + L( t, j - 1 );
                else L( t + 2, j ) = max( L( t + 1, j - 1 ), L( t + 1, j ) );
             }  
           else L( t + 2, j ) = 0;
         }
                      	      
    Pochoir_Kernel_End

    LCS.Register_Domain( J );

    int t = nX + nY;

    LCS.Run( t, LCS_fn );
    
    int optLen = L.interior( t + 1, nY );

    return optLen;    
}



int standardDPLCS( int nX, char *X, int nY, char *Y )
{
    Pochoir_Array< int, N_RANK > L( nY + 1 );
    L.Register_Shape(LCS_shape);

    for ( int j = 0; j <= nY; ++j )
        L.interior( 0, j ) = 0;
    
    for ( int i = 1; i <= nX; ++i ) 
      {
        L.interior( i, 0 ) = 0;

        for ( int j = 1; j <= nY; ++j ) 
          {
            if ( X[ i ] == Y[ j ] ) L.interior( i, j ) = 1 + L.interior( i - 1, j - 1 );
            else L.interior( i, j ) = max( L.interior( i, j - 1 ), L.interior( i - 1, j ) );
          }          
      }        

    int optLen = L.interior( nX, nY );    
    
    return optLen;    
}



int iterativeStencilLCS( int nX, char *X, int nY, char *Y )
{
    Pochoir_Array< int, N_RANK > L( nY + 1 );
    L.Register_Shape(LCS_shape);
    
    L.interior( 0, 0 ) = L.interior( 1, 0 ) = L.interior( 1, 1 ) = 0;
    
    for ( int t = 0; t < nX + nY; ++t )
       cilk_for ( int j = max( t + 1 - nX, 0 ); j <= min( t + 1, nY ); ++j )
          {
            int i = t + 1 - j;
            
            if ( ( i > 0 ) && ( j > 0 ) )
               {
                  if ( X[ i ] == Y[ j ] ) L.interior( t + 2, j ) = 1 + L.interior( t, j - 1 );
                  else L.interior( t + 2, j ) = max( L.interior( t + 1, j - 1 ), L.interior( t + 1, j ) );
               }  
            else L.interior( t + 2, j ) = 0;
          }
          
    int t = nX + nY + 1;

    int optLen = L.interior( t, nY );

    return optLen;    
}



int main( int argc, char *argv[ ] )
{
    printf( "\nStencil-based DP for the length of the longest common subsequence ( Run with option -h for help ).\n\n" );

    int nX, nY;
    char *X = NULL, *Y = NULL;    
    char *fnX = NULL, *fnY = NULL;
    int RunStandardDP, RunIterativeStencil;

    if ( !read_command_line( argc, argv, nX, &fnX, nY, &fnY, RunStandardDP, RunIterativeStencil ) )
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

        printf( "Running pochoir-based DP..." );
        fflush( stdout );
               
        gettimeofday( &start, 0 );        
        int optLen = stencilLCS( nX, X, nY, Y );    
        gettimeofday( &end, 0 );

        double t0 = tdiff( &end, &start );
              
        printf( "\n\nPochoir:\n" );
        printf( "\t LCS length = %d\n", optLen );    
        printf( "\t Running time = %.3lf sec\n\n", t0 );    
      
        if ( RunIterativeStencil )
          {
            printf( "Running iterative stencil..." );
            fflush( stdout );
                          
            gettimeofday( &start, 0 );        
            int optLenITST = iterativeStencilLCS( nX, X, nY, Y );    
            gettimeofday( &end, 0 );

            double t1 = tdiff( &end, &start );
          
            printf( "\n\nIterative Stencil:\n" );
            printf( "\t LCS length = %d\n", optLenITST );    
            if ( t0 > 0 ) printf( "\t Running time = %.3lf sec ( %.3lf x Pochoir )\n\n", t1, t1 / t0 );    
            else printf( "\t Running time = %.3lf sec\n\n", t1 );    
          }

        if ( RunStandardDP )
          {
            printf( "Running standard DP..." );
            fflush( stdout );
                          
            gettimeofday( &start, 0 );        
            int optLenSDP = standardDPLCS( nX, X, nY, Y );    
            gettimeofday( &end, 0 );
            
            double t1 = tdiff( &end, &start );
          
            printf( "\n\nStandard DP:\n" );
            printf( "\t LCS length = %d\n", optLenSDP );    
            if ( t0 > 0 ) printf( "\t Running time = %.3lf sec ( %.3lf x Pochoir )\n\n", t1, t1 / t0 );    
            else printf( "\t Running time = %.3lf sec\n\n", t1 );    
          }
      }

    if ( X != NULL ) free( X );
    if ( Y != NULL ) free( Y );    
       
    return 0;
}
