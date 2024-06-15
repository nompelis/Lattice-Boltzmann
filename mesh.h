
#ifndef _MESH_H_
#define _MESH_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif


// a structure to store information about a generic Cartesian box
// (has redundancy via the dx,dy,dz for efficiency)
struct box_s {
   int im,jm,km;
   double xs[3],xe[3];
   double dx,dy,dz;
};


// a function to establish a Cartesian box
int make_box( struct box_s* p );


#ifdef __cplusplus
}
#endif
#endif

