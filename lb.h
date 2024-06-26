
#ifndef _LB_H_
#define _LB_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mesh.h"


#ifdef __cplusplus
extern "C" {
#endif

// a convenient 3-dimensional vector of doubles
typedef double vec3_t[3];

// UNFINISHED definitions to help access the distribution lattice (TODO)
#define MID_CEN     0

// a structure to store the distribution of the lattice (to be extended)
struct lattice_s {
   double f[19];
};

// a structure to store an LB simulation's settings and data
struct lb_s {
   struct box_s* box;
   const vec3_t* elat;
   const double* wlat;
   int ilat;
   struct lattice_s* lattice;
   double* u;

};


// --------- function prototypes --------

// initialization of a Lattice-boltzmann structure
int lb_init( struct lb_s* p );

// a function to store the initial condition distribution from a velocity field
int lb_initCond( struct lb_s* p, const double* u_ );

// get density and momenta at a single cell
void lb_recoverRU( const struct lb_s* p,
                   const struct lattice_s* lattice,
                   double* rho, double* ruvw );


#ifdef __cplusplus
}
#endif
#endif

