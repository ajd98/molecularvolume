#ifndef RECURSIVEFILL_H 
#define RECURSIVEFILL_H

#include <math.h>

double dist(double x1, double y1, double z1, double x2, double y2, double z2);

char is_free(double voxx, double voxy, double voxz, double* solute_pos,
             double* solute_rad, int nsolute, double solvent_rad);

void recurse(int ix, int iy, int iz, int nx, int ny, int nz, double voxel_len,
             unsigned char* visited_grid, unsigned char* grid, int* moves,
             double* solute_pos, double* solute_rad, int nsolute, 
             double solvent_rad);

#endif


