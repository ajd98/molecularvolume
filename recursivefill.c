#include <math.h>
#include "recursivefill.h"

double dist(double x1, double y1, double z1, double x2, double y2, double z2) {
    /*
     * Calculate the distance between (x1, y1, z1) and (x2, y2, z2)
     */ 
    return sqrt(pow(x2-x1, 2) + pow(z2-z1, 2) + pow(z2-z1, 2));
}


unsigned char is_free(double voxx, double voxy, double voxz, double* solute_pos,
                      double* solute_rad, int nsolute, double solvent_rad)
{
    /*
     * ---------------
     * Parameters
     * ---------------
     * voxx: (double) the x-position of the voxel
     * voxy: (double) the y-position of the voxel
     * voxz: (double) the z-position of the voxel
     * solute_pos: shape (n,3) array of positions of solute atoms
     * solute_rad: shape (n,) array of radii of solute atoms
     * nsolute: (int) The value ``n`` in the above two lines, ie, the number of
     *     solute atoms.
     * solvent_rad: (double) the radius of the solvent 

     * --------
     * Returns
     * -------
     * True if a sphere of radius ``solvent_rad`` centered at (voxx, voxy, voxz) 
     * does not intersect any sphere of radius solute_rad[i] centered at 
     * solute_pos[i], for any i=1...n

     * False otherwise
     */
    char check = 1;
    double solx;
    double soly;
    double solz;
    double solrad;

    for (int isolute = 0; isolute < nsolute; isolute++){
        solx = solute_pos[isolute*3];
        soly = solute_pos[isolute*3 + 1];
        solz = solute_pos[isolute*3 + 2];
        solrad = solute_rad[isolute];
        if (dist(voxx, voxy, voxz, solx, soly, solz) < solrad + solvent_rad) {
            check = 0;
            break;
        }
    }
    return check;
}


void recurse(int ix, int iy, int iz, int nx, int ny, int nz, double voxel_len,
             unsigned char* visited_grid, unsigned char* grid, int* moves,
             double* solute_pos, double* solute_rad, int nsolute, 
             double solvent_rad)
{
    /*
     * ----------
     * Parameters
     * ----------
     * ix: x index of current grid position
     * iy: y index of current grid position
     * iz: z index of current grid position
     * nx: length of grid along x-axis, in voxels
     * ny: length of grid along y-axis, in voxels
     * nz: length of grid along z-axis, in voxels
     * voxel_len: edge length of (cubic) voxel
     * visited_grid: shape (nx, ny, nz) array; 1 if visited, 0 otherwise
     * grid: shape (nx, ny, nz) array; 1 if accessible by moves on solvent; 
     *     0 otherwise
     * moves: shape (26,3) array; possible moves
     * solute_pos: shape (nsolute, 3) array; (x,y,z) coordinates of each solute
     *     atom
     * solute_rad: shape (nsolute,) array; radius of each solute atom
     * nsolute: number of solute atoms
     * solvent_rad: radius of solvent (which is approximated as a sphere)
     */
    int i;

    // check if we are in bounds of the grid
    if ((ix < 0) || (iy < 0) || (iz < 0) || (ix >= nx) || (iy >=ny) || (iz >=nz)){
        return;
    } 

    // check if the grid point is visited
    if (visited_grid[ix*ny*nz + iy*nz + iz] == 1) {
        return;
    }

    // now enter the main part of the program
      
    // specify that ix, iy, iz has been visited
    visited_grid[ix*ny*nz + iy*nz + iz] = 1;


    if (is_free(ix*voxel_len, iy*voxel_len, iz*voxel_len, solute_pos, solute_rad, nsolute, solvent_rad)){
        // specify that grid index (ix, iy, iz) is accessible to solvent
        grid[ix*ny*nz + iy*nz + iz] = 1;
        for (i=0;i<26;i++) {
            recurse(ix+moves[3*i], 
                    iy+moves[3*i+1], 
                    iz+moves[3*i+2], 
                    nx, ny, nz, voxel_len, visited_grid,
                    grid, moves, solute_pos, solute_rad, nsolute, solvent_rad);
        }
    }
    return;
}
