#!python
#cython: boundscheck=False, wraparound=False, profile=True, linetrace=True
#distutils: define_macros=CYTHON_TRACE_NOGIL=1
cimport cython
cimport numpy
import numpy
from libc.math cimport floor, ceil, sqrt
from libc.stdlib cimport malloc, free
'''
Make sure that the solute is whole (rather than split over a periodic boundary)
before running this.
'''

cdef extern from "floodfill3d.h":
    double dist(double x1, double y1, double z1, double x2, double y2, double z2) nogil
    
    unsigned char is_free(double voxx, double voxy, double voxz, double* solute_pos,
                          double* solute_rad, int nsolute, double solvent_rad) nogil
    
    void floodfill(int ix, int iy, int iz, int nx, int ny, int nz, double voxel_len,
                   unsigned char* visited_grid, unsigned char* grid, int* moves,
                   double* solute_pos, double* solute_rad, int nsolute, 
                   double solvent_rad) nogil
    

#cpdef double volume(numpy.ndarray[numpy.float64_t, ndim=2] _solute_pos, 
def volume(numpy.ndarray[numpy.float64_t, ndim=2] _solute_pos, 
                    numpy.ndarray[numpy.float64_t, ndim=1] _solute_rad, 
                    numpy.ndarray[numpy.float64_t, ndim=2] _solvent_pos, 
                    double _solvent_rad,
                    double _voxel_len):
    '''
    -----------
    Parameters
    -----------
    _solute_pos: shape (n,3) array of (x,y,z) coordinates of solute atoms
    _solute_rad: shape (n,) array of solute atom radii
    _solvent_pos: shape (m,3) array of (x,y,z) coordinates of solvent atoms
    _solvent_rad: solvent molecule radius; solvent is approximated as sphere
    _voxel_len: edge length of voxels

    -------------
    Returns
    -------------
    The volume that is accessible to the centers of the water molecules
    '''
    cdef:
        # Have Cython convert the arrays with a leading underscore from numpy
        # arrays to typed memoryviews
        double *solute_pos = <double *>malloc(_solute_pos.shape[0]*3*sizeof(double))
        double *solute_rad = <double *>malloc(_solute_rad.shape[0]*sizeof(double))
        double *solvent_pos = <double *>malloc(_solvent_pos.shape[0]*3*sizeof(double))
        double *solvent_floor = <double *>malloc(_solute_pos.shape[0]*3*sizeof(double))
        double *solvent_ceil = <double *>malloc(_solute_pos.shape[0]*3*sizeof(double))

        double box_buffer 
        double x_min, y_min, z_min
        double x_max, y_max, z_max
        int nx, ny, nz
        double x, y, z, X, Y, Z
        int ix, iy, iz, iX, iY, iZ
        double solx, soly, solz
        int nsolvent = _solute_pos.shape[0]
        int nsolute = _solvent_pos.shape[0]
        int i, j, k
        double voxel_len = _voxel_len
        double solvent_rad = _solvent_rad

    # assign solute_pos, solute_rad, and solvent_pos
    with nogil:
        for i in range(nsolute):
            solute_pos[3*i]   = _solute_pos[i,0]
            solute_pos[3*i+1] = _solute_pos[i,1]
            solute_pos[3*i+2] = _solute_pos[i,2]
            solute_rad[i] = _solute_rad[i]
        for i in range(nsolvent):
            solvent_pos[3*i]   = _solvent_pos[i,0]
            solvent_pos[3*i+1] = _solvent_pos[i,1]
            solvent_pos[3*i+2] = _solvent_pos[i,2]
            solvent_floor[3*i]   = 0
            solvent_floor[3*i+1] = 0
            solvent_floor[3*i+2] = 0
            solvent_ceil[3*i]   = 0
            solvent_ceil[3*i+1] = 0
            solvent_ceil[3*i+2] = 0


    # We want the grids to be large enough that solvent can completely surround
    # the solute; calculate a buffer size to do this.
    box_buffer = 2*(_solute_rad.max() + solvent_rad)
    x_min = min(_solute_pos[:,0].min(), _solvent_pos[:,0].min())
    y_min = min(_solute_pos[:,1].min(), _solvent_pos[:,1].min())
    z_min = min(_solute_pos[:,2].min(), _solvent_pos[:,2].min())

    x_max = max(_solute_pos[:,0].max(), _solvent_pos[:,0].max())
    y_max = max(_solute_pos[:,1].max(), _solvent_pos[:,1].max())
    z_max = max(_solute_pos[:,2].max(), _solvent_pos[:,2].max())

    # shift the solute positions by the appropriate amount
    with nogil:
        for i in range(nsolute):
            solute_pos[3*i]   += (-1*x_min + box_buffer)
            solute_pos[3*i+1] += (-1*y_min + box_buffer)
            solute_pos[3*i+2] += (-1*z_min + box_buffer)
        for i in range(nsolvent):
            solvent_pos[3*i]   += (-1*x_min + box_buffer)
            solvent_pos[3*i+1] += (-1*y_min + box_buffer)
            solvent_pos[3*i+2] += (-1*z_min + box_buffer)
    

    # the extreme values for the (rectangular) grid
    x_min -= box_buffer
    y_min -= box_buffer
    z_min -= box_buffer
    x_max += box_buffer
    y_max += box_buffer
    z_max += box_buffer

    # Find the number of voxels to use; avoid fencepost error
    #
    #        * -- * -- * -- * -- * -- * -- *
    #        1    2    3    4    5    6    7 
    # 
    nx = numpy.ceil((x_max - x_min)/voxel_len) + 1
    ny = numpy.ceil((y_max - y_min)/voxel_len) + 1
    nz = numpy.ceil((z_max - z_min)/voxel_len) + 1

    # unsigned char is 8-bit; I use it as a bool 
    cdef unsigned char *grid = <unsigned char *>malloc(nx*ny*nz*sizeof(unsigned char))
    cdef unsigned char *visited_grid = <unsigned char *>malloc(nx*ny*nz*sizeof(unsigned char))
    with nogil:
        for i in range(nx*ny*nz):
            grid[i] = 0
            visited_grid[i] = 0



    # Find the positions of the voxels surround each solvent 
    with nogil:
        for i in range(nsolute):
            solx = solvent_pos[3*i+0]
            soly = solvent_pos[3*i+1]
            solz = solvent_pos[3*i+2]

            solvent_floor[3*i+0] = floor(solx/voxel_len)*voxel_len
            solvent_ceil[3*i+0] = ceil(solx/voxel_len)*voxel_len

            solvent_floor[3*i+1] = floor(soly/voxel_len)*voxel_len
            solvent_ceil[3*i+1] = ceil(soly/voxel_len)*voxel_len

            solvent_floor[3*i+2] = floor(solz/voxel_len)*voxel_len
            solvent_ceil[3*i+2] = ceil(solz/voxel_len)*voxel_len


    # possible directions to move; first we make _moves, which is a memoryview
    # then we convert the memoryview to a true C-style array
    cdef int[:,:] _moves = numpy.array([[0,0,-1],
                                        [0,0,1],
                                                   
                                        [0,-1,-1],
                                        [0,-1,0],
                                        [0,-1,1],
                                                   
                                        [0,1,-1],
                                        [0,1,0],
                                        [0,1,1],
                                                   
                                        [-1,-1,-1],
                                        [-1,-1,0],
                                        [-1,-1,1],
                                                   
                                        [-1,0,-1],
                                        [-1,0,0],
                                        [-1,0,1],
                                                   
                                        [-1,1,-1],
                                        [-1,1,0],
                                        [-1,1,1],
                                                   
                                        [1,-1,-1],
                                        [1,-1,0],
                                        [1,-1,1],
                                                   
                                        [1,0,-1],
                                        [1,0,0],
                                        [1,0,1],
                                                   
                                        [1,1,-1],
                                        [1,1,0],
                                        [1,1,1]], dtype=numpy.int32)
    # Start making the true c-style array
    cdef int *moves = <int *>malloc(3*26*sizeof(int))
    with nogil:
        for i in range(26):
            moves[3*i+0] = _moves[i,0]
            moves[3*i+1] = _moves[i,1]
            moves[3*i+2] = _moves[i,2]

    cdef int ptcnt = 0
    cdef double vol

    with nogil:
        # 
        # The algorithm:
        #   For each water molecule, consider which grid points are possible by
        #   moving that molecule around. The molecule is approximated by a
        #   sphere of size ``solvent_rad``, and it cannot move to voxel centers
        #   where the sphere of radius ``solvent_rad`` centered at that point
        #   would overlap with any of the spheres centered at a coordinate in
        #   ``solute_pos``, with the corresponding radius from ``solute_rad``.
        #   This process is repeated for every water molecule.
        # 
        for i in range(nsolvent):
            # The water falls between 8 voxel coordinates
            #
            #
            #          * -------- *
            #          |\         |\
            #          | \     |  | \
            #          |  \    |  |  \ 
            #          |   * -------- *
            #          * --|---w--*   |
            #           \  |   |   \  |
            #            \ |   |    \ |
            #             \|         \|
            #              * -------- *
            #
            x = solvent_floor[3*i+0]
            X = solvent_ceil[3*i+0]
            y = solvent_floor[3*i+1]
            Y = solvent_ceil[3*i+1]
            z = solvent_floor[3*i+2]
            Z = solvent_ceil[3*i+2]

            # the index of ``grid`` corresponding to the given voxel
            ix = int(x/voxel_len)
            iX = int(X/voxel_len)
            iy = int(y/voxel_len)
            iY = int(Y/voxel_len)
            iz = int(z/voxel_len)
            iZ = int(Z/voxel_len)
            # check if any solute molecule would overlap a solvent placed at
            # any of the following:
            #  (x,y,z)
            #  (x,y,Z)
            #  (x,Y,z)
            #  (x,Y,Z)
            #  (X,y,z)
            #  (X,y,Z)
            #  (X,Y,z)
            #  (X,Y,Z)
            # 

            if is_free(x, y, z, solute_pos, solute_rad, nsolute, solvent_rad):
                floodfill(ix, iy, iz, nx, ny, nz, voxel_len, visited_grid, grid, 
                        moves, solute_pos, solute_rad, nsolute, solvent_rad)

            if is_free(x, y, Z, solute_pos, solute_rad, nsolute, solvent_rad):
                floodfill(ix, iy, iZ, nx, ny, nz, voxel_len, visited_grid, grid, 
                        moves, solute_pos, solute_rad, nsolute, solvent_rad)

            if is_free(x, Y, z, solute_pos, solute_rad, nsolute, solvent_rad):
                floodfill(ix, iY, iz, nx, ny, nz, voxel_len, visited_grid, grid, 
                        moves, solute_pos, solute_rad, nsolute, solvent_rad)

            if is_free(x, Y, Z, solute_pos, solute_rad, nsolute, solvent_rad):
                floodfill(ix, iY, iZ, nx, ny, nz, voxel_len, visited_grid, grid, 
                        moves, solute_pos, solute_rad, nsolute, solvent_rad)

            if is_free(X, y, z, solute_pos, solute_rad, nsolute, solvent_rad):
                floodfill(iX, iy, iz, nx, ny, nz, voxel_len, visited_grid, grid, 
                        moves, solute_pos, solute_rad, nsolute, solvent_rad)

            if is_free(X, y, Z, solute_pos, solute_rad, nsolute, solvent_rad):
                floodfill(iX, iy, iZ, nx, ny, nz, voxel_len, visited_grid, grid, 
                        moves, solute_pos, solute_rad, nsolute, solvent_rad)

            if is_free(X, Y, z, solute_pos, solute_rad, nsolute, solvent_rad):
                floodfill(iX, iY, iz, nx, ny, nz, voxel_len, visited_grid, grid, 
                        moves, solute_pos, solute_rad, nsolute, solvent_rad)

            if is_free(X, Y, Z, solute_pos, solute_rad, nsolute, solvent_rad):
                floodfill(iX, iY, iZ, nx, ny, nz, voxel_len, visited_grid, grid, 
                        moves, solute_pos, solute_rad, nsolute, solvent_rad)

        #find the volume
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    ptcnt += grid[i*ny*nz+j*nz+k]

        vol = (1-float(ptcnt)/float(nx*ny*nz))*(x_max-x_min)*(y_max-y_min)*(z_max-z_min)
        free(grid)
        free(visited_grid)
        free(solute_pos)
        free(solute_rad)
        free(solvent_pos)
        free(solvent_floor)
        free(solvent_ceil)

    #return vol
    ### for debugging ###
    _grid = numpy.zeros((int(nx),int(ny),int(nz)))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                _grid[i,j,k] = grid[i*ny*nz+j*nz+k]
    return vol, _grid
