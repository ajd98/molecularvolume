#!/usr/bin/env python
import numpy
import matplotlib.pyplot as pyplot
from mpl_toolkits.mplot3d import Axes3D
import sys
sys.path.append('../')
import volume
import pdb2volume

def test_simple():
    solute_pos = numpy.array(((0,0,0),), dtype=numpy.float64)
    solute_rad = numpy.array((3,), dtype=numpy.float64)
    solvent_pos = numpy.array(((8,0,0),), dtype=numpy.float64)
    solvent_rad = 1.4
    voxel_len = 0.05

    print("Test: single atom of radius 3")
    vol = volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, 
                              voxel_len)
    
    analytical_vol = numpy.pi*4./3*(solute_rad[0]+solvent_rad)**3

    print("  calculated volume: {:f}".format(vol))
    print("  analytical volume: {:f}".format(analytical_vol))
    perr = 100*abs(vol-analytical_vol)/analytical_vol
    print("  percent error:     {:f}%".format(perr))
    if perr < 1:
        print("  TEST PASSED")
        return 0
    else:
        print("  TEST FAILED")
        return 1

def test_simple_overlap():
    solute_pos = numpy.array(((0,0,0),
                              (0,0,0)), dtype=numpy.float64)
    solute_rad = numpy.array((3,3), dtype=numpy.float64)
    solvent_pos = numpy.array(((8,0,0),), dtype=numpy.float64)
    solvent_rad = 1.4
    voxel_len = 0.05

    print("Test: two atoms of radius 3, completely overlapping")
    vol = volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, 
                              voxel_len)
    
    analytical_vol = numpy.pi*4./3*(solute_rad[0]+solvent_rad)**3

    print("  calculated volume: {:f}".format(vol))
    print("  analytical volume: {:f}".format(analytical_vol))
    perr = 100*abs(vol-analytical_vol)/analytical_vol
    print("  percent error:     {:f}%".format(perr))
    if perr < 1:
        print("  TEST PASSED")
        return 0
    else:
        print("  TEST FAILED")
        return 1

def test_2sphere():
    solute_pos = numpy.array(((0,0,0),
                              (8,0,0)), dtype=numpy.float64)
    sr = 3
    solute_rad = numpy.array(((sr),
                              (sr)), dtype=numpy.float64)

    solvent_pos = numpy.array(((16,0,0),), dtype=numpy.float64)
    solvent_rad = 1.4
    voxel_len = 0.02

    print("Test: two non-overlapping atoms of radius 3")
    vol = volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, 
                              voxel_len)
    
    r = solute_rad[0]+solvent_rad
    analytical_vol = numpy.pi*4./3*r**3 \
                   + numpy.pi*4./3*r**3 

    print("  calculated volume: {:f}".format(vol))
    print("  analytical volume: {:f}".format(analytical_vol))
    perr = 100*abs(vol-analytical_vol)/analytical_vol
    print("  percent error:     {:f}%".format(perr))
    if perr < 1:
        print("  TEST PASSED")
        return 0
    else:
        print("  TEST FAILED")
        return 1

def test_2sphere_overlapping():
    solute_pos = numpy.array(((0,0,0),
                              (3,0,0)), dtype=numpy.float64)
    sr = 3
    solute_rad = numpy.array(((sr),
                              (sr)), dtype=numpy.float64)

    solvent_pos = numpy.array(((11,0,0),), dtype=numpy.float64)
    solvent_rad = 1.4
    voxel_len = 0.05

    print("Test: two partially overlapping atoms of radius 3")
    vol = volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, 
                              voxel_len)
    
    # volume of sphere 0 + volume of sphere 1 - volume of intersection
    h = (sr+solvent_rad) - (solute_pos[1,0] - solute_pos[0,0])/2
    r = sr+solvent_rad
    analytical_vol = numpy.pi*4./3*r**3 \
                   + numpy.pi*4./3*r**3 \
                   - 2*1./3*numpy.pi*h**2*(3*r-h)

    print("  calculated volume: {:f}".format(vol))
    print("  analytical volume: {:f}".format(analytical_vol))
    perr = 100*abs(vol-analytical_vol)/analytical_vol
    print("  percent error:     {:f}%".format(perr))
    if perr < 1:
        print("  TEST PASSED")
        return 0
    else:
        print("  TEST FAILED")
        return 1
    

def test_protein():
    print("Test: full protein in explicit waters")

    try:
        vol = pdb2volume.PDBVolume('villin.pdb', 
                                   'radii.lib', voxel_len=0.5).run()
        print("  calculated volume: {:f}".format(vol))
        print("  analytical volume unknown.")
        print("  TEST PASSED")
        return 0
    except:
        print("  TEST FAILED")
        return 1

def test_void():
    print("Test: capsid structure with void")
    s = 1.4


    # Make a cube of spheres.  the edge length of the cube is 4*s, where s is
    # the radius of cube.
    #
    #
    #        *      *      *
    #    
    #
    #        *      *      *
    #             
    #          
    #        *      *      *
    # The space between surfaces of spheres on the diagonal is 2*(sqrt(2) - 1)*s
    # Since 2*(sqrt(2)-1) < 1, a water of radius s will not fit through the spaces.
    #

    solute_pos = numpy.array(((0  ,0  ,0  ),
                              (2*s,0  ,0  ),
                              (4*s,0  ,0  ),
                              (6*s,0  ,0  ),
                              (0  ,2*s,0  ),
                              (2*s,2*s,0  ),
                              (4*s,2*s,0  ),
                              (6*s,2*s,0  ),
                              (0  ,4*s,0  ),
                              (2*s,4*s,0  ),
                              (4*s,4*s,0  ),
                              (6*s,4*s,0  ),
                              (0  ,6*s,0  ),
                              (2*s,6*s,0  ),
                              (4*s,6*s,0  ),
                              (6*s,6*s,0  ),

                              (0  ,0  ,2*s),
                              (2*s,0  ,2*s),
                              (4*s,0  ,2*s),
                              (6*s,0  ,2*s),
                              (0  ,2*s,2*s),
                             #(2*s,2*s,2*s), # skip this one, since it's in the middle of the cube
                             #(4*s,2*s,2*s),
                              (6*s,2*s,2*s),
                              (0  ,4*s,2*s),
                             #(2*s,4*s,2*s),
                             #(4*s,4*s,2*s),
                              (6*s,4*s,2*s),
                              (0  ,6*s,2*s),
                              (2*s,6*s,2*s),
                              (4*s,6*s,2*s),
                              (6*s,6*s,2*s),

                              (0  ,0  ,4*s),
                              (2*s,0  ,4*s),
                              (4*s,0  ,4*s),
                              (6*s,0  ,4*s),
                              (0  ,2*s,4*s),
                             #(2*s,2*s,4*s),
                             #(4*s,2*s,4*s),
                              (6*s,2*s,4*s),
                              (0  ,4*s,4*s),
                             #(2*s,4*s,4*s),
                             #(4*s,4*s,4*s),
                              (6*s,4*s,4*s),
                              (0  ,6*s,4*s),
                              (2*s,6*s,4*s),
                              (4*s,6*s,4*s),
                              (6*s,6*s,4*s),

                              (0  ,0  ,6*s),
                              (2*s,0  ,6*s),
                              (4*s,0  ,6*s),
                              (6*s,0  ,6*s),
                              (0  ,2*s,6*s),
                              (2*s,2*s,6*s),
                              (4*s,2*s,6*s),
                              (6*s,2*s,6*s),
                              (0  ,4*s,6*s),
                              (2*s,4*s,6*s),
                              (4*s,4*s,6*s),
                              (6*s,4*s,6*s),
                              (0  ,6*s,6*s),
                              (2*s,6*s,6*s),
                              (4*s,6*s,6*s),
                              (6*s,6*s,6*s)), dtype=numpy.float64)


    solute_rad = numpy.ones(solute_pos.shape[0])*s

    # one solvent outside cube
    solvent_pos = numpy.array(((10*s,0,0),), dtype=numpy.float64)
    solvent_rad = s

    voxel_len = 0.1

    vol_empty = volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, 
                              voxel_len)
    print("  Calculated volume with solvent only outside cage: {:f}".format(vol_empty))

    # Now include a solvent inside the cube.
    solvent_pos = numpy.array(((10*s,0,0),
                               (3*s,3*s,3*s)), dtype=numpy.float64)


    vol_filled = volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, 
                               voxel_len)
    print("  Calculated volume with solvent also inside cage:  {:f}".format(vol_filled))
    
    dv_actual = vol_empty - vol_filled
    # the difference should be between:
    dv_lb = (2*s)**3
    # based on close packing of spheres in fcc latice; the pitch of tetrahedral 
    # arangement of sphere centers is sqrt(6)*d/3, where d=2*s
    dv_ub = (2*s+2*(2*s-6.0**(0.5)*2*s/3))**3
    print("  Difference in volumes should be between {:f} and {:f}".format(dv_lb, dv_ub))
    print("  Calculated difference in volumes: {:f}".format(vol_empty-vol_filled))
    if dv_lb < dv_actual and dv_actual < dv_ub:
        print("  TEST PASSED")
    else:
        print("  TEST FAILED")


def show_protein_surface():
    vol, grid = pdb2volume.PDBVolume('villin.pdb', 
                                     'radii.lib', voxel_len=0.5).run()
    print(vol)
    fig = pyplot.figure()
    ax = fig.add_subplot(111, projection='3d')
    w = numpy.where(grid==0)
    ax.scatter(w[0], w[1], w[2], marker='o')
    ax.set_xlim(0,grid.shape[0])
    ax.set_ylim(0,grid.shape[1])
    ax.set_zlim(0,grid.shape[2])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    pyplot.show()

if __name__ == "__main__":
    test_simple()
    test_simple_overlap()
    test_2sphere()
    test_2sphere_overlapping()
    test_protein()
    test_void()
