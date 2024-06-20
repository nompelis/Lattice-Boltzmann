#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mesh.h"
#include "lb.h"

int main( int argc, char *argv[])
{
   struct box_s b;
   b.im = 21;
   b.jm = 11;
   b.km = 2;
   b.im = 6;   // HACK
   b.jm = 4;   // HACK
   b.xs[0] = 0.0; b.xs[1] = 0.0; b.xs[2] = 0.0;
   b.xe[0] = 20.0; b.xe[1] = 10.0; b.xe[2] = 1.0/((double) b.km);
   make_box( &b );

   struct lb_s lb1;
   lb1.box = &b;
   lb_init( &lb1 );

   // initial condition (macroscopic)
   for(int i=0;i<(int) ((b.im-1)*(b.jm-1)*(b.km-1));++i) {
      lb1.u[3*i  ] = 1.0;
      lb1.u[3*i+1] = 0.0;
      lb1.u[3*i+2] = 0.0;
   }
   lb_initCond( &lb1, (const double*) lb1.u );

   // test recovery of density and momenta
   for(int i=0;i<(int) ((b.im-1)*(b.jm-1)*(b.km-1));++i) {
      double rho;
      vec3_t ruvw;

      lb_recoverRU( &lb1, &( lb1.lattice[i] ), &rho, ruvw );
      printf(" Rho: %lf  Ruvw: %lf %lf %lf \n", rho, ruvw[0],ruvw[1],ruvw[2] );
   }

   return 0;
}

