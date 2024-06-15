#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lb.h"
#include "mesh.h"


//
// Function to make a notional, generic "box" for any kind of use that needs
// a Cartesian-style discretization defined by "start -> end" points in 3D
// and number of nodes in each direction.
// (This function does not do much, but it may get extendedin the future...)
//

int make_box( struct box_s* p )
{
   if( p == NULL ) return 1;

   p->dx = (p->xe[0] - p->xs[0]) / ((double) (p->im-1));
   p->dy = (p->xe[1] - p->xs[1]) / ((double) (p->jm-1));
   p->dz = (p->xe[2] - p->xs[2]) / ((double) (p->km-1));
   printf("Sizes: %lf,%lf,%lf  %d %d %d \n",
          p->dx, p->dy, p->dz,  p->im, p->jm, p->km );

   return 0;
}

