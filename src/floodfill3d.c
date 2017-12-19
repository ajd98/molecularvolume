#include "queue.h"
#include <math.h>

double dist(double x1, double y1, double z1, double x2, double y2, double z2) {
    /*
     * Calculate the distance between (x1, y1, z1) and (x2, y2, z2)
     */ 
    double dx = x2-x1;
    double dy = y2-y1;
    double dz = z2-z1;
    return sqrt(dx*dx + dy*dy + dz*dz);
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
    double cut;

    for (int isolute = 0; isolute < nsolute; isolute++){
        solx = solute_pos[isolute*3];
        soly = solute_pos[isolute*3 + 1];
        solz = solute_pos[isolute*3 + 2];
        solrad = solute_rad[isolute];
        cut = solrad + solvent_rad;
        // the distance beetween (solx, soly, solz) and (voxx, voxy, voxz) is
        // bounded below by fabs(solx-voxx) etc.  Check these first to avoid 
        // expensive sqrt calculations.
        if ((fabs(solx - voxx)>cut) || (fabs(soly - voxy) > cut) || (fabs(solz - voxz) > cut)){
            continue;
        }
        else if (dist(voxx, voxy, voxz, solx, soly, solz) < cut) {
            check = 0;
            break;
        }
    }
    return check;
}


void floodfill(int ix, int iy, int iz, int nx, int ny, int nz, double voxel_len,
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

    // Create a *huge* queue
    struct Queue* queue = newQueue(nx*ny*nz); 
    struct Triple coord;
    struct Triple newcoord;
    int i;
    coord.x = ix;
    coord.y = iy;
    coord.z = iz;
    append(coord, queue);
    visited_grid[coord.x*ny*nz + coord.y*nz + coord.z] = 1;

    // The main loop.
    while (!isEmpty(queue)){
        coord = pop(queue);
        
        if (is_free(coord.x*voxel_len,
                   coord.y*voxel_len,
                   coord.z*voxel_len,
                   solute_pos,
                   solute_rad,
                   nsolute,
                   solvent_rad))
        {
            grid[coord.x*ny*nz + coord.y*nz + coord.z] = 1;

            for (i=0;i<26;i++) {
                newcoord.x = coord.x+moves[3*i];
                newcoord.y = coord.y+moves[3*i+1];
                newcoord.z = coord.z+moves[3*i+2];
                if ((visited_grid[newcoord.x*ny*nz + newcoord.y*nz + newcoord.z]==0) 
                     && (newcoord.x >= 0) && (newcoord.x < nx) 
                     && (newcoord.y >= 0) && (newcoord.y < ny) 
                     && (newcoord.z >= 0) && (newcoord.z < nz))
                {
                    visited_grid[newcoord.x*ny*nz + newcoord.y*nz + newcoord.z] = 1;
                    append(newcoord, queue);
                } 
            }
        }  
    }
}
