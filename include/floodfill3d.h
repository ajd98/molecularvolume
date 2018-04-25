#ifndef FLOODFILL3D_H 
#define FLOODFILL3D_H

double dist(double x1, double y1, double z1, double x2, double y2, double z2);

unsigned char is_free(double voxx, double voxy, double voxz, 
                      const double* restrict solute_pos,
                      const double* restrict solute_rad, int nsolute, 
                      double solvent_rad);

void floodfill(int ix, int iy, int iz, int nx, int ny, int nz, double voxel_len,
               unsigned char* restrict visited_grid, 
               unsigned char* restrict grid, const int* restrict moves,
               const double* restrict solute_pos, 
               const double* restrict solute_rad, int nsolute,
               double solvent_rad);

#endif


