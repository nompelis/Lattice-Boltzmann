
#ifndef _LB_H_
#define _LB_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mesh.h"


#ifdef __cplusplus
extern "C" {
#endif


typedef double vec3_t[3];

#define MID_CEN     0


struct lattice_s {
   double f[19];
};

struct lb_s {
   struct box_s* box;
   struct lattice_s* lattice;
   int il;
   double* u;

};


// --------- function prototypes --------

int lb_init( struct lb_s* p );
int lb_initCond( struct lb_s* p, const double* u_ );


#ifdef __cplusplus
}
#endif
#endif

