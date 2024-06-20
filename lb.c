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

   p->elat = e_D3Q19;
   p->wlat = w_D3Q19;
   p->ilat = 19;
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
   int ilat = p->ilat;

   double rho = 1.0; // arbitrary density
   const double cs = 1.0/sqrt(3.0);  // speed of sound in the lattice
   const double cs2 = cs*cs, cs4 = cs2*cs2;

   for(int i=0;i<(int) ((im-1)*(jm-1)*(km-1));++i) {
      struct lattice_s* t = &(lp[i]);
      const double* uu = &( u_[3*i] );
#ifdef _DEBUG2_
      printf(" u[%d] = %lf %lf %lf \n",i, uu[0], uu[1], uu[2] ); 
#endif

      for(int n=0;n<ilat;++n) {
         double* ee = (double*) p->elat[n];
         double eu = ee[0]*uu[0] + ee[1]*uu[1] + ee[2]*uu[2];
         double u2 = uu[0]*uu[0] + uu[1]*uu[1] + uu[2]*uu[2];

         t->f[n] = w_D3Q19[n] * rho * (1.0 + eu/(cs*cs)
                                           + eu*eu/(2.0*cs4)
                                           - u2/(2.0*cs2) );
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


//    
// Function to calculate density and momenta from the distribution
// at a given cell
//    

void lb_recoverRU( const struct lb_s* p,
                   const struct lattice_s* lattice,
                   double* rho, double* ruvw )
{        
   const int ilat = p->ilat;    

   double r = 0.0;
   double _ruvw[3] = { 0.0, 0.0, 0.0 };
   for(int n=0;n<ilat;++n) {
      const double f = lattice->f[n];
      r += f;
      _ruvw[0] += f * p->elat[n][0];
      _ruvw[1] += f * p->elat[n][1];
      _ruvw[2] += f * p->elat[n][2];
   }
   *rho = r;
   ruvw[0] = _ruvw[0];
   ruvw[1] = _ruvw[1];
   ruvw[2] = _ruvw[2];
}

