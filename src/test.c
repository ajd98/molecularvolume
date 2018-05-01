#include "floodfill3d.h"
#include <stdlib.h>
#include <stdio.h>

int 
main (void)
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
    int* moves = (int*)malloc(26*3*sizeof(int));
    double* solute_pos = (double*)malloc(nsolute*3*sizeof(double));
    double* solute_rad = (double*)malloc(nsolute*sizeof(double));

    for(int i=0; i<(nx*ny*nz); i++){
        grid[i] = 0;
        visited_grid[i] = 0;
    }

    moves[3*0+0] = 0;
    moves[3*0+1] = 0;
    moves[3*0+2] = -1;

    moves[3*1+0] = 0;
    moves[3*1+1] = 0;
    moves[3*1+2] = 1;

    moves[3*2+0] = 0;
    moves[3*2+1] = -1;
    moves[3*2+2] = -1;

    moves[3*3+0] = 0;
    moves[3*3+1] = -1;
    moves[3*3+2] = 0;

    moves[3*4+0] = 0;
    moves[3*4+1] = -1;
    moves[3*4+2] = 1;

    moves[3*5+0] = 0;
    moves[3*5+1] = 1;
    moves[3*5+2] = -1;

    moves[3*6+0] = 0;
    moves[3*6+1] = 1;
    moves[3*6+2] = 0;

    moves[3*7+0] = 0;
    moves[3*7+1] = 1;
    moves[3*7+2] = 1;

    moves[3*8+0] = -1;
    moves[3*8+1] = -1;
    moves[3*8+2] = -1;

    moves[3*9+0] = -1;
    moves[3*9+1] = -1;
    moves[3*9+2] = 0;

    moves[3*10+0] = -1;
    moves[3*10+1] = -1;
    moves[3*10+2] = 1;

    moves[3*11+0] = -1;
    moves[3*11+1] = 0;
    moves[3*11+2] = -1;

    moves[3*12+0] = -1;
    moves[3*12+1] = 0;
    moves[3*12+2] = 0;

    moves[3*13+0] = -1;
    moves[3*13+1] = 0;
    moves[3*13+2] = 1;

    moves[3*14+0] = -1;
    moves[3*14+1] = 1;
    moves[3*14+2] = -1;

    moves[3*15+0] = -1;
    moves[3*15+1] = 1;
    moves[3*15+2] = 0;

    moves[3*16+0] = -1;
    moves[3*16+1] = 1;
    moves[3*16+2] = 1;

    moves[3*17+0] = 1;
    moves[3*17+1] = -1;
    moves[3*17+2] = -1;

    moves[3*18+0] = 1;
    moves[3*18+1] = -1;
    moves[3*18+2] = 0;

    moves[3*19+0] = 1;
    moves[3*19+1] = -1;
    moves[3*19+2] = 1;

    moves[3*20+0] = 1;
    moves[3*20+1] = 0;
    moves[3*20+2] = -1;

    moves[3*21+0] = 1;
    moves[3*21+1] = 0;
    moves[3*21+2] = 0;

    moves[3*22+0] = 1;
    moves[3*22+1] = 0;
    moves[3*22+2] = 1;

    moves[3*23+0] = 1;
    moves[3*23+1] = 1;
    moves[3*23+2] = -1;

    moves[3*24+0] = 1;
    moves[3*24+1] = 1;
    moves[3*24+2] = 0;

    moves[3*25+0] = 1;
    moves[3*25+1] = 1;
    moves[3*25+2] = 1;

    solute_pos[3*0+0] = 50;
    solute_pos[3*0+1] = 50;
    solute_pos[3*0+2] = 50;

    solute_rad[0] = 10;

    double solvent_rad = 1.4;

    floodfill(ix, iy, iz, nx, ny, nz, voxel_len, visited_grid, grid, moves, 
              solute_pos, solute_rad, nsolute, solvent_rad);

    free(visited_grid);
    free(grid);
    free(moves);
    free(solute_pos);
    free(solute_rad);
    return 0;
}
