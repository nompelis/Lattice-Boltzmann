
#ifndef _MESH_H_
#define _MESH_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif


struct box_s {
   int im,jm,km;
   double xs[3],xe[3];
   double dx,dy,dz;
};



int make_box( struct box_s* p );


#ifdef __cplusplus
}
#endif
#endif

