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

/* Test bench - calculate the price of American put option */
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

#define DEFAULT_S 100.00
#define DEFAULT_E 100.00
#define DEFAULT_r   0.10
#define DEFAULT_V   0.25
#define DEFAULT_T   1.00

#define DEFAULT_s 100
#define DEFAULT_t 100

void print_usage( char *prog )
{
  printf( "Usage: %s [ options ]\n\n", prog );

  printf( "Options:\n" );

  printf( "\t-S value : spot price ( default: %0.2lf )\n", DEFAULT_S );
  printf( "\t-E value : exercise price ( default: %0.2lf )\n", DEFAULT_E );
  printf( "\t-r value : interest rate ( default: %0.2lf\% )\n", DEFAULT_r * 100 );    
  printf( "\t-V value : volatility ( default: %0.2lf\% )\n", DEFAULT_V * 100 );  
  printf( "\t-T value : time to mature in years ( default: %0.2lf )\n\n", DEFAULT_T );    

  printf( "\t-s value : steps in space dimension ( default: %d )\n", DEFAULT_s );
  printf( "\t-t value : steps in time dimension ( default: %d )\n\n", DEFAULT_t );

  printf( "\t-i               : Run iterative stencil\n\n" );  
   
  printf( "\t-h               : print this help screen\n\n" );
}

   Pochoir_Shape< N_RANK > APOP_shape[ ] = { { 1, 0 }, { 0, -1 }, { 0, 0 }, { 0, 1 } };    


int read_command_line( int argc, char *argv[ ], 
		       double &S, double &E, double &r, double &V, double &T, 
		       int &ns, int &nt,
                       int &Run_iter_stencil )
{
  S = DEFAULT_S;
  E = DEFAULT_E;
  r = DEFAULT_r;
  V = DEFAULT_V;
  T = DEFAULT_T;
  
  ns = DEFAULT_s;
  nt = DEFAULT_t;
  
  Run_iter_stencil = 0;

  for ( int i = 1; i < argc; )
    {
     int j = i;

     if ( !strcmp( argv[ i ], "-i" ) )
       {
        Run_iter_stencil = 1;
         
        i++;
        
        if ( i >= argc ) break;
       }


     if ( !strcmp( argv[ i ], "-S" ) )
       {
        if ( i + 1 >= argc )
          {
           printf( "Error: Missing spot price ( specify -S spot-price )!\n\n" );
           return 0;
          }

        S = atof( argv[ i + 1 ] );

        if ( S < 0 )
          {
           printf( "Error: Spot price must be non-negative!\n\n" );
           return 0;
          }

        i += 2;

        if ( i >= argc ) break;
       }


     if ( !strcmp( argv[ i ], "-E" ) )
       {
        if ( i + 1 >= argc )
          {
           printf( "Error: Missing exercise price ( specify -E exercise-price )!\n\n" );
           return 0;
          }

        E = atof( argv[ i + 1 ] );

        if ( E < 0 )
          {
           printf( "Error: Exercise price must be non-negative!\n\n" );
           return 0;
          }

        i += 2;

        if ( i >= argc ) break;
       }
       

     if ( !strcmp( argv[ i ], "-r" ) )
       {
        if ( i + 1 >= argc )
          {
           printf( "Error: Missing interest rate ( specify -r interest-rate-in-percentage )!\n\n" );
           return 0;
          }

        r = atof( argv[ i + 1 ] ) / 100.0;

        if ( r < 0 )
          {
           printf( "Error: Interest rate must be non-negative!\n\n" );
           return 0;
          }

        i += 2;

        if ( i >= argc ) break;
       }


     if ( !strcmp( argv[ i ], "-V" ) )
       {
        if ( i + 1 >= argc )
          {
           printf( "Error: Missing volatility ( specify -V volatility-in-percentage )!\n\n" );
           return 0;
          }

        V = atof( argv[ i + 1 ] ) / 100.0;

        if ( V < 0 )
          {
           printf( "Error: Volatility must be non-negative!\n\n" );
           return 0;
          }

        i += 2;

        if ( i >= argc ) break;
       }


     if ( !strcmp( argv[ i ], "-T" ) )
       {
        if ( i + 1 >= argc )
          {
           printf( "Error: Missing maturity time ( specify -T time-to-maturity-in-years )!\n\n" );
           return 0;
          }

        T = atof( argv[ i + 1 ] );

        if ( T <= 0 )
          {
           printf( "Error: Time to maturity must be positive!\n\n" );
           return 0;
          }

        i += 2;

        if ( i >= argc ) break;
       }


     if ( !strcmp( argv[ i ], "-s" ) )
       {
        if ( i + 1 >= argc )
          {
           printf( "Error: Missing steps in space dimension ( specify -s steps-in-space-dimension )!\n\n" );
           return 0;
          }

        ns = atoi( argv[ i + 1 ] );

        if ( ns <= 0 )
          {
           printf( "Error: Number of steps in space dimension must be positive!\n\n" );
           return 0;
          }

        i += 2;

        if ( i >= argc ) break;
       }


     if ( !strcmp( argv[ i ], "-t" ) )
       {
        if ( i + 1 >= argc )
          {
           printf( "Error: Missing steps in time dimension ( specify -t steps-in-time-dimension )!\n\n" );
           return 0;
          }

        nt = atoi( argv[ i + 1 ] );

        if ( nt <= 0 )
          {
           printf( "Error: Number of steps in time dimension must be positive!\n\n" );
           return 0;
          }

        i += 2;

        if ( i >= argc ) break;
       }

                     
     if ( !strcmp( argv[ i ], "-h" ) || !strcmp( argv[ i ], "-help" ) || !strcmp( argv[ i ], "--help" ) )
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
    
   printf( "Parameters:\n\n" );
  
   printf( "\t spot price = %0.2lf\n", S );
   printf( "\t exercise price = %0.2lf\n", E );
   printf( "\t interest rate = %0.2lf\%\n", r * 100 );    
   printf( "\t volatility = %0.2lf\%\n", V * 100 );  
   printf( "\t time to mature ( in years ) = %0.2lf\n\n", T );    
  
   printf( "\t steps in space dimension = %d\n", ns );
   printf( "\t steps in time dimension = %d\n\n", nt );
    
   return 1;
}



void computeCoeffs( double r, double V, double T, int ns, int nt,
		    Pochoir_Array< double, N_RANK > &c )
{
   double V2 = V * V;
   double dt = T / nt;
   double r1 = 1.0 / ( 1.0 + r * dt );
   double r2 = dt / ( 1.0 + r * dt );
   
   cilk_for ( int i = 0; i <= ns; ++i )
     {
       c.interior( 0, i ) = r2 * 0.5 * i * ( - r + V2 * i );       
       c.interior( 1, i ) = r1 * ( 1 - V2 * i * i * dt );
       c.interior( 2, i ) = r2 * 0.5 * i * ( r + V2 * i );       
     }
}


Pochoir_Boundary_1D( apop_bv_1D, arr, t, i )

   if ( ( i < arr.size( 0 ) - 1 ) || ( t == 0 ) ) return arr.get( t, i );
   else return 0;

Pochoir_Boundary_End



double stencilAPOP( double S, double E, double r, double V, double T, 
		 int ns, int nt )
{
   ns = ns + ( ns & 1 );
   double dS = 2.0 * S / ns;
   
   Pochoir< N_RANK > APOP(APOP_shape);
   Pochoir_Array< double, N_RANK > c( ns + 1 );
   Pochoir_Array< double, N_RANK > f( ns + 1 );
   APOP.Register_Array( f );    
   APOP.Register_Array( c );
   
   computeCoeffs( r, V, T, ns, nt, c );   
   
   
   cilk_for ( int i = 0; i <= ns; ++i )
       f.interior( 0, i ) = max( 0.0, E - i * dS );
       
   f.interior( 1, 0 ) = E;    
       
   Pochoir_Domain I( 1, ns );
    
   Pochoir_Kernel_1D( APOP_fn, t, i )
        
       double v = c( 0, i ) * f( t, i - 1 )
                + c( 1, i ) * f( t, i )
       	        + c( 2, i ) * f( t, i + 1 );
        
       f( t + 1, i ) = max( v, E - i * dS );			   
  
   Pochoir_Kernel_End

   APOP.Register_Domain( I );   
   f.Register_Boundary( apop_bv_1D );

   APOP.Run( nt, APOP_fn );
    
   return f.interior( nt, ( ns >> 1 ) );    
}



double iterativeStencilAPOP( double S, double E, double r, double V, double T, 
		             int ns, int nt )
{
   ns = ns + ( ns & 1 );
   double dS = 2.0 * S / ns;
   
   Pochoir_Array< double, N_RANK > c( ns + 1 );
   c.Register_Shape(APOP_shape);
   
   computeCoeffs( r, V, T, ns, nt, c );   
   
   Pochoir_Array< double, N_RANK > f( ns + 1 );
   f.Register_Shape(APOP_shape);
   
   cilk_for ( int i = 0; i <= ns; ++i )
       f.interior( 0, i ) = max( 0.0, E - i * dS );
       
   f.interior( 1, 0 ) = E;    

   for ( int t = 0; t < nt; ++t )
     {        
       f.interior( t + 1, ns ) = 0;
     
       cilk_for ( int i = 1; i < ns; ++i )
         {
           double v = c.interior( 0, i ) * f.interior( t, i - 1 )
                    + c.interior( 1, i ) * f.interior( t, i )
            	    + c.interior( 2, i ) * f.interior( t, i + 1 );
            
           f.interior( t + 1, i ) = max( v, E - i * dS );			   
         }
     }    
      
   return f.interior( nt, ( ns >> 1 ) );  
}



int main( int argc, char *argv[ ] )
{
    printf( "\nStencil-based DP for the price of American put option ( Run with option -h for help ).\n\n" );

    double S, E, r, V, T; 
    int ns, nt;
    
    int RunIterativeStencil;

    if ( !read_command_line( argc, argv, S, E, r, V, T, ns, nt, RunIterativeStencil ) )
      {
        print_usage( argv[ 0 ] );
        return 1;
      }

    struct timeval start, end;

    printf( "Running pochoir-based DP..." );
    fflush( stdout );
           
    gettimeofday( &start, 0 );        
    double price0 = stencilAPOP( S, E, r, V, T, ns, nt );    
    gettimeofday( &end, 0 );

    double t0 = tdiff( &end, &start );
          
    printf( "\n\nPochoir:\n" );
    printf( "\t option price = %.2lf\n", price0 );    
    //printf( "\t Running time = %.3lf sec\n\n", t0 );    
    printf( "\t Running time = %.3lf ms\n\n", 1.0e3 * t0 );    
  
    if ( RunIterativeStencil )
      {
        printf( "Running iterative stencil..." );
        fflush( stdout );
                      
        gettimeofday( &start, 0 );        
        double price1 = iterativeStencilAPOP( S, E, r, V, T, ns, nt );    
        gettimeofday( &end, 0 );

        double t1 = tdiff( &end, &start );
      
        printf( "\n\nIterative Stencil:\n" );
        printf( "\t option price = %.2lf\n", price1 );    
        if ( t0 > 0 ) printf( "\t Running time = %.3lf sec ( %.3lf x Pochoir )\n\n", t1, t1 / t0 );    
        else printf( "\t Running time = %.3lf sec\n\n", t1 );    
      }
       
    return 0;
}
