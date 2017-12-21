#!/usr/bin/env python
import pdb2volume
import numpy
import matplotlib.pyplot as pyplot
from mpl_toolkits.mplot3d import Axes3D
import sys
sys.path.append('../')
import volume

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
