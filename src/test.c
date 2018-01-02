#include "floodfill3d.h"
#include <stdlib.h>
#include <stdio.h>

int main ()
{
    int nsolute = 1;

    int nx = 500;
    int ny = 500;
    int nz = 500;
    unsigned char *grid = (unsigned char *)malloc(nx*ny*nz*sizeof(unsigned char));
    unsigned char *visited_grid = (unsigned char *)malloc(nx*ny*nz*sizeof(unsigned char));
    int ix = 0;
    int iy = 0;
    int iz = 0;
    double voxel_len = 1;
    double* solute_pos = (double*)malloc(nsolute*3*sizeof(double));
    double* solute_rad = (double*)malloc(nsolute*sizeof(double));

    for(int i=0; i<(nx*ny*nz); i++){
        grid[i] = 0;
        visited_grid[i] = 0;
    }

    solute_pos[3*0+0] = 50;
    solute_pos[3*0+1] = 50;
    solute_pos[3*0+2] = 50;

    solute_rad[0] = 10;

    double solvent_rad = 1.4;

    floodfill(ix, iy, iz, nx, ny, nz, voxel_len, visited_grid, grid, 
              solute_pos, solute_rad, nsolute, solvent_rad);

    free(visited_grid);
    free(grid);
    free(solute_pos);
    free(solute_rad);
    return 0;
}
