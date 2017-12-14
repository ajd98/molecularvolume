cimport cython
cimport numpy
import numpy
from libc.math cimport floor, ceil, sqrt
'''
Make sure that the solute is whole (rather than split over a periodic boundary)
before running this.
'''

cdef bool is_free(double voxx, double voxy, double voxz,
        double[:,:] solute_pos, double[:] solute_rad,
        int nsolute, double solvent_rad):
    '''
    ---------------
    Parameters
    ---------------
    voxx: (double) the x-position of the voxel
    voxy: (double) the y-position of the voxel
    voxz: (double) the z-position of the voxel
    solute_pos: shape (n,3) array of positions of solute atoms
    solute_rad: shape (n,) array of radii of solute atoms
    nsolute: (int) The value ``n`` in the above two lines, ie, the number of
        solute atoms.
    solvent_rad: (double) the radius of the solvent 

    --------
    Returns
    -------
    True if a sphere of radius ``solvent_rad`` centered at (voxx, voxy, voxz) 
    does not intersect any sphere of radius solute_rad[i] centered at 
    solute_pos[i], for any i=1...n

    False otherwise
    '''

    cdef bool check = True

    for isolute in range(nsolute):
        solx = solute_pos[isolute,0]
        soly = solute_pos[isolute,1]
        solz = solute_pos[isolute,2]
        solrad = solute_rad[isolute]
        if dist(x,y,z,solx,soly,solz) < solrad + solvent_rad:
            check = False
            break
    return check


cdef double dist(double x1, double y1, double z1, 
                 double x2, double y2, double z2):
    '''
    Calculate the distance between points (x1,y1,z1) and (x2,y2,z2)
    '''
    return sqrt((x2-x1)**2, (y2-y1)**2, (z2-z1)**2)


cpdef double volume(numpy.ndarray[numpy.float64_t, ndim=2] _solute_pos, 
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
        # arrays to C-style arrays
        double[:,:] solute_pos = _solute_pos
        double[:] solute_rad = _solute_rad
        double[:,:] solvent_pos = _solvent_pos
        double box_buffer 
        double x_min, y_min, z_min
        double x_max, y_max, z_max
        int nx, ny, nz
        double x, y, z, X, Y, Z
        double solx, soly, solz
        double solvent_floor[:,:]
        double solvent_ceil[:,:]
        bool check


    # We want the grids to be large enough that solvent can completely surround
    # the solute; calculate a buffer size to do this.
    box_buffer = 2*(_solute_rad.max() + solvent_rad)
    x_min = min(_solute_pos[:,0].min(), _solvent_pos[:,0].min())
    y_min = min(_solute_pos[:,1].min(), _solvent_pos[:,1].min())
    z_min = min(_solute_pos[:,2].min(), _solvent_pos[:,2].min())

    x_max = max(_solute_pos[:,0].max(), _solvent_pos[:,0].max())
    y_max = max(_solute_pos[:,1].max(), _solvent_pos[:,1].max())
    z_max = max(_solute_pos[:,2].max(), _solvent_pos[:,2].max())
    

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
    ny = numpy.ceil((x_max - x_min)/voxel_len) + 1
    nz = numpy.ceil((x_max - x_min)/voxel_len) + 1

    cdef bool[nx,ny,nz] grid = numpy.zeros((nx,ny,nz), dtype=bool)
    cdef bool[nx,ny,nz] visited_grid = numpy.zeros((nx,ny,nz), dtype=bool)

    # shift the solute positions by the appropriate amount
    # (will this work with C-style arrays?)
    solute_pos[:,0] -= xmin
    solute_pos[:,1] -= ymin
    solute_pos[:,2] -= zmin

    solvent_pos[:,0] -= xmin
    solvent_pos[:,1] -= ymin
    solvent_pos[:,2] -= zmin

    solvent_floor = numpy.floor(solvent_pos)
    solvent_ceil = numpy.floor(solvent_pos)
    nsolvent = solvent_pos.shape[0]
    nsolute = solute_pos.shape[0]

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
        for iwater in range(nsolvent):
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
            x = solvent_floor[iwater, 0]
            X = solvent_ciel[iwater, 0]
            y = solvent_floor[iwater, 1]
            Y = solvent_ciel[iwater, 1]
            z = solvent_floor[iwater, 2]
            Z = solvent_ciel[iwater, 2]

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
            '''
            open voxels surrounding water
            assign ``visited_grid`` to one for each of the open voxels

            while there are open voxels:
                for each open voxel:
                    if voxel is free:
                        assign corresponding value in ``grid`` to one
                        for each adjacent voxel:
                            if adjacent voxel is not visited:
                                open voxel
                                assign corresponding value in ``visited_grid`` to one 
                    close voxel
                    
            '''
            #                         
            #                 
            #                 
            #             
            check = True
            if is_free(x, y, z, solute_pos, solute_rad, nsolute):
                open_voxel(
                visited_grid[ix,iy,iz] = 1

            if is_free(x, y, Z, solute_pos, solute_rad, nsolute):
                visited_grid[ix,iy,iZ] = 1

            if is_free(x, Y, z, solute_pos, solute_rad, nsolute):
                visited_grid[ix,iY,iz] = 1

            if is_free(x, Y, Z, solute_pos, solute_rad, nsolute):
                visited_grid[ix,iY,iZ] = 1

            if is_free(X, y, z, solute_pos, solute_rad, nsolute):
                visited_grid[iX,iy,iz] = 1

            if is_free(X, y, Z, solute_pos, solute_rad, nsolute):
                visited_grid[iX,iy,iZ] = 1

            if is_free(X, Y, z, solute_pos, solute_rad, nsolute):
                visited_grid[iX,iY,iz] = 1

            if is_free(X, Y, Z, solute_pos, solute_rad, nsolute):
                visited_grid[iX,iY,iZ] = 1



























cdef bool open_arr_is_empty(int[:,:] *open_arr, ngridpts):
    bool empty = True
    for i in range(ngridpts):
        if *open_arr[i,0] != -1 or *open_arr[i,1] != -1 or *open_arr[i,2] != -1:
            empty = False
            break
    return empty

cdef void flood(bool[:,:,:] *grid, ix, iy, iz, nx, ny, nz,
                double[:,:] solute_pos, double[:] solute_rad,
                int nsolute, double voxel_len, int ngridpts):
    '''
    flood the grid from the voxel indexed by (ix, iy, iz)
    

    '''
    cdef double x 
    cdef double y 
    cdef double z 
    cdef int[ngridpts,3] open_arr

    # initialize open_arr to -1
    for i in range(ngridpts):
        open_arr[i,0] = -1
        open_arr[i,1] = -1
        open_arr[i,2] = -1

    open_arr[0,0] = ix 
    open_arr[0,0] = iy 
    open_arr[0,0] = iz 

    # Find open returns true if it finds any voxel open
    # and returns False otherwise.  It also changes the 
    # values of ix, iy, and iz to the indices of the first
    # open voxel it finds.
    # Finally, find_open removes (ix, iy, iz) from the open
    # array
    # Before closing an open voxel, grid is assigned a value of True/False
    while find_open(&open_arr, &ix, &iy, &iz, ngridpts):
        x = ix*voxel_len
        y = iy*voxel_len
        z = iz*voxel_len
        if is_free(x, y, z, solute_pos, solute_rad, nsolute):
            add_to_open(

        



cdef bool find_open(int[:,:] *open_arr, int *ix, int *iy, int *iz, int ngridpts):
    cdef int i
    for i in range(ngridpts):
        if *open_arr[i,0] != -1:
            *ix = open_arr[i,0]
            *iy = open_arr[i,1]
            *iz = open_arr[i,2]
            open_arr[i,0] = -1
            open_arr[i,1] = -1
            open_arr[i,2] = -1
            return True
    return False

