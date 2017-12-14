cimport cython
cimport numpy
import numpy
from libc.math cimport floor, ceil, sqrt
'''
Make sure that the solute is whole (rather than split over a periodic boundary)
before running this.
'''

cdef bint is_free(double voxx, double voxy, double voxz,
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

    cdef bint check = True

    for isolute in range(nsolute):
        solx = solute_pos[isolute,0]
        soly = solute_pos[isolute,1]
        solz = solute_pos[isolute,2]
        solrad = solute_rad[isolute]
        if dist(voxx,voxy,voxz,solx,soly,solz) < solrad + solvent_rad:
            check = False
            break
    return check


cdef double dist(double x1, double y1, double z1, 
                 double x2, double y2, double z2):
    '''
    Calculate the distance between points (x1,y1,z1) and (x2,y2,z2)
    '''
    return sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)


cdef bint find_open(int[:,:]* open_arr, int* ix, int* iy, int* iz, int ngridpts):
    cdef int i
    for i in range(ngridpts):
        if open_arr[0][i,0] != -1:
            ix[0] = open_arr[i,0]
            iy[0] = open_arr[i,1]
            iz[0] = open_arr[i,2]
            open_arr[i,0] = -1
            open_arr[i,1] = -1
            open_arr[i,2] = -1
            return True
    return False


cdef void open_voxel(int[:,:]* open_arr, int arraylen, int ix, int iy, int iz):
    '''
    Open the voxel (ix, iy, iz)
    '''
    for i in range(arraylen):
        if open_arr[0][i,0] == -1:
            open_arr[0][i,0] = ix
            open_arr[0][i,1] = iy
            open_arr[0][i,2] = iz
            break
    return 

cdef void close_voxel(int[:,:]* open_arr, int arraylen, int ix, int iy, int iz):
    '''
    close the voxel (ix, iy, iz)
    '''
    for i in range(arraylen):
        if open_arr[0][i,0] == ix and open_arr[0][i,1] == iy and open_arr[0][i,2] == iz:
            open_arr[0][i,0] = -1
            open_arr[0][i,1] = -1
            open_arr[0][i,2] = -1
            break
    return 


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
        # arrays to typed memoryvieews
        double[:,:] solute_pos = _solute_pos
        double[:] solute_rad = _solute_rad
        double[:,:] solvent_pos = _solvent_pos
        double box_buffer 
        double x_min, y_min, z_min
        double x_max, y_max, z_max
        int nx, ny, nz
        double x, y, z, X, Y, Z
        int ix, iy, iz, iX, iY, iZ
        double solx, soly, solz
        double[:,:] solvent_floor = numpy.zeros(_solute_pos.shape)
        double[:,:] solvent_ceil = numpy.zeros(_solute_pos.shape)
        int nsolvent
        int solute
        bint check
        int i, j, k
        double voxel_len = _voxel_len


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

    cdef bint[:,:,:] grid = numpy.zeros((nx, ny, nz), dtype=numpy.bool_)
    cdef bint[:,:,:] visited_grid = numpy.zeros((nx, ny, nz), dtype=numpy.bool_)
    # initialize grid to 0

    cdef int ngridpts = nx*ny*nz
    cdef int[:,:,:] open_arr = numpy.ones((ngridpts,3), dtype=numpy.int32)

    with nogil:
        for i in range(ngridpts):
            for j in range(3):
                open_arr[i,j] = -1 

    # shift the solute positions by the appropriate amount
    # (will this work with C-style arrays?)
    solute_pos[:,0] -= xmin
    solute_pos[:,1] -= ymin
    solute_pos[:,2] -= zmin

    solvent_pos[:,0] -= xmin
    solvent_pos[:,1] -= ymin
    solvent_pos[:,2] -= zmin

    # These will be used to find the voxels surrounding each solvent
    nsolvent = _solvent_pos.shape[0]
    nsolute = _solute_pos.shape[0]

    # Find the positions of the voxels surround each solvent 
    with nogil:
        solvent_floor = numpy.zeros(_solvent_pos)
        solvent_ceil = numpy.zeros(_solvent_pos)
        for isolute in range(nsolute):
            solx = solvent_pos[isolute,0]
            soly = solvent_pos[isolute,1]
            solz = solvent_pos[isolute,2]

            solvent_floor[isolute,0] = int(solx/voxel_len)*voxel_len
            solvent_ceil[isolute,0] = (int(solx/voxel_len)+1)*voxel_len

            solvent_floor[isolute,1] = int(soly/voxel_len)*voxel_len
            solvent_ceil[isolute,1] = (int(soly/voxel_len)+1)*voxel_len

            solvent_floor[isolute,2] = int(solz/voxel_len)*voxel_len
            solvent_ceil[isolute,2] = (int(solz/voxel_len)+1)*voxel_len


    cdef int[:,:] movelist = numpy.zeros((26,3), dtype=numpy.int32)
    movelist[0,:]  = [0,0,-1]
    movelist[1,:]  = [0,0,1]

    movelist[2,:]  = [0,-1,-1]
    movelist[3,:]  = [0,-1,0]
    movelist[4,:]  = [0,-1,1]

    movelist[5,:]  = [0,1,-1]
    movelist[6,:]  = [0,1,0]
    movelist[7,:]  = [0,1,1]

    movelist[8,:]  = [-1,-1,-1]
    movelist[9,:]  = [-1,-1,0]
    movelist[10,:] = [-1,-1,1]

    movelist[11,:] = [-1,0,-1]
    movelist[12,:] = [-1,0,0]
    movelist[13,:] = [-1,0,1]

    movelist[14,:] = [-1,1,-1]
    movelist[15,:] = [-1,1,0]
    movelist[16,:] = [-1,1,1]

    movelist[17,:] = [1,-1,-1]
    movelist[18,:] = [1,-1,0]
    movelist[19,:] = [1,-1,1]

    movelist[20,:] = [1,0,-1]
    movelist[21,:] = [1,0,0]
    movelist[22,:] = [1,0,1]

    movelist[23,:] = [1,1,-1]
    movelist[24,:] = [1,1,0]
    movelist[25,:] = [1,1,1]

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
            if is_free(x, y, z, solute_pos, solute_rad, nsolute, solvent_rad):
                open_voxel(&open_arr, ngridpts, ix, iy, iz)
                visited_grid[ix,iy,iz] = 1

            if is_free(x, y, Z, solute_pos, solute_rad, nsolute, solvent_rad):
                open_voxel(&open_arr, ngridpts, ix, iy, iZ)
                visited_grid[ix,iy,iZ] = 1

            if is_free(x, Y, z, solute_pos, solute_rad, nsolute, solvent_rad):
                open_voxel(&open_arr, ngridpts, ix, iY, iz)
                visited_grid[ix,iY,iz] = 1

            if is_free(x, Y, Z, solute_pos, solute_rad, nsolute, solvent_rad):
                open_voxel(&open_arr, ngridpts, ix, iY, iZ)
                visited_grid[ix,iY,iZ] = 1

            if is_free(X, y, z, solute_pos, solute_rad, nsolute, solvent_rad):
                open_voxel(&open_arr, ngridpts, iX, iy, iz)
                visited_grid[iX,iy,iz] = 1

            if is_free(X, y, Z, solute_pos, solute_rad, nsolute, solvent_rad):
                open_voxel(&open_arr, ngridpts, iX, iy, iZ)
                visited_grid[iX,iy,iZ] = 1

            if is_free(X, Y, z, solute_pos, solute_rad, nsolute, solvent_rad):
                open_voxel(&open_arr, ngridpts, iX, iY, iz)
                visited_grid[iX,iY,iz] = 1

            if is_free(X, Y, Z, solute_pos, solute_rad, nsolute, solvent_rad):
                open_voxel(&open_arr, ngridpts, iX, iY, iZ)
                visited_grid[iX,iY,iZ] = 1

            # now we use ix, iy, and iz for general indexing of the grids
            # find_open returns false when nothing is open
            # it also sticks in new values for ix, iy, and iz, so that they 
            # correspond to open voxels
            while find_open(&open_arr, &ix, &iy, &iz, ngridpts):

                x = voxel_len*ix
                y = voxel_len*iy
                z = voxel_len*iz

                if is_free(x, y, z, solute_pos, solute_rad, nsolute, solvent_rad):
                    grid[ix, iy, iz] = True
                    # iterate over adjacent voxels
                    for imove in range(26):
                        iX = movelist[imove][0] + ix
                        iY = movelist[imove][1] + iy
                        iZ = movelist[imove][2] + iz
                        # Check if the move is outside the grid
                        if iX >= 0 and iY >= 0 and iZ >= 0 and iX < nx and iY < nY and iZ < nZ:
                            if not visited_grid[iX, iY, iZ]:
                                open_voxel(&open_arr, ngridpts, iX, iY, iZ)
                                visted_grid[iX, iY, iZ] = True
                close_voxel(&open_arr, ngridpts, ix, iy, iz)
        #find the volume
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    ptcnt += grid[i,j,k]
        vol = float(ptcnt)/float(nx*ny*nz)*(x_max-x_min)*(y_max-y_min)*(z_max-z_min)
        return vol

