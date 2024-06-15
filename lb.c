#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lb.h"


// lattice weights
const double w_D3Q19[] = {
   1.0/3.0,
   1.0/18.0, 1.0/18.0, 1.0/18.0,
   1.0/18.0, 1.0/18.0, 1.0/18.0,
   1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
   1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
   1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
};

// lattice "vector" directions
const vec3_t e_D3Q19[] = {
  { 0.0, 0.0, 0.0},

  { 1.0,  0.0,  0.0}, {-1.0,  0.0,  0.0},
  { 0.0,  1.0,  0.0}, { 0.0, -1.0,  0.0},
  { 0.0,  0.0,  1.0}, { 0.0,  0.0, -1.0},

  { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}, { 1.0,-1.0, 0.0}, {-1.0,-1.0, 0.0},
  { 1.0, 0.0, 1.0}, {-1.0, 0.0, 1.0}, { 1.0, 0.0,-1.0}, {-1.0, 0.0,-1.0},
  { 0.0, 1.0, 1.0}, { 0.0,-1.0, 1.0}, { 0.0, 1.0,-1.0}, { 0.0,-1.0,-1.0}
};


//
// Function to initialize a Lattice-boltzmann structure
//

int lb_init( struct lb_s* p )
{
   if( p == NULL ) return 1;
#ifdef _DEBUG_
   printf("Initialization of structures \n");
#endif

   struct box_s* bp = p->box;
   int im = bp->im, jm = bp->jm, km = bp->km;

   p->il = 19;
   size_t isize = (size_t) ((im-1)*(jm-1)*(km-1));

   p->lattice = (struct lattice_s*) malloc( isize*sizeof(struct lattice_s) );
   p->u = (double*) malloc( 3*isize*sizeof(double) );
   if( p->lattice == NULL || p->u == NULL ) {
      printf("Error allocating memory \n");
      return -1;
   }

   return 0;
}


//
// Function to establish the initial condition's distribution from a velocity
// field.
//

int lb_initCond( struct lb_s* p, const double* u_ )
{
   if( p == NULL ) return 1;
#ifdef _DEBUG_
   printf("Initial conditions \n");
#endif

   struct box_s* bp = p->box;
   struct lattice_s* lp = p->lattice;

   int im = bp->im, jm = bp->jm, km = bp->km;
   int il = p->il;

   double rho = 1.0; // arbitrary density
   double cs = 9.9;  // HACK (speed of sound in the lattice)

   for(int i=0;i<(int) ((im-1)*(jm-1)*(km-1));++i) {
      struct lattice_s* t = &(lp[i]);
      const double* uu = &( u_[i] );
#ifdef _DEBUG2_
      printf(" u[%d] = %lf %lf %lf \n",i, uu[0], uu[1], uu[2] ); 
#endif

      for(int n=0;n<il;++n) {
         double* ee = (double*) &( e_D3Q19[n] );       // re-cast

         double eu = ee[0]*uu[0] + ee[1]*uu[1] + ee[2]*uu[2];
         double u2 = uu[0]*uu[0] + uu[1]*uu[1] + uu[2]*uu[2];

         t->f[n] = w_D3Q19[n] * rho * (1.0 + eu/(cs*cs)
                                           + eu*eu/(2.8*cs*cs*cs*cs)
                                           - u2/(2.0*cs*cs) );
      }
#ifdef _DEBUG2_
      // diagnostics to visualize a cell's distribution; hardwired to cell 0
      if( i == 0 ) {
         printf(" ----- Distribution lattice for a cell ---- \n");
         printf(" == k=far (low) plane == \n");
         printf("          X          %1.9e       X    \n", t->f[17] );
         printf("    %1.9e  %1.9e  %1.9e \n", t->f[14], t->f[ 6], t->f[13] );
         printf("          X          %1.9e       X       \n", t->f[18] );
         printf(" == k=mid plane == \n");
         printf("    %1.9e  %1.9e  %1.9e \n", t->f[ 8], t->f[ 3], t->f[ 7] );
         printf("    %1.9e  %1.9e  %1.9e \n", t->f[ 2], t->f[ 0], t->f[ 1] );
         printf("    %1.9e  %1.9e  %1.9e \n", t->f[10], t->f[ 4], t->f[ 9] );
         printf(" == k=near (high) plane == \n");
         printf("          X          %1.9e       X       \n", t->f[15] );
         printf("    %1.9e  %1.9e  %1.9e \n", t->f[12], t->f[ 5], t->f[11] );
         printf("          X          %1.9e       X       \n", t->f[16] );
         printf(" ----- ------------------------------- ---- \n");
      }
#endif
   }

   return 0;
}

