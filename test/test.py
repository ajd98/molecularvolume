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
    solute_rad = numpy.array(((3),), dtype=numpy.float64)
    solvent_pos = numpy.array(((8,0,0),), dtype=numpy.float64)
    solvent_rad = 1.4
    voxel_len = 0.05

    vol = volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, 
                              voxel_len)
    
    analytical_vol = numpy.pi*4./3*(solute_rad[0]+solvent_rad)**3

    print("Test: single atom of radius 3")
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
                              (3,0,0)), dtype=numpy.float64)
    sr = 3
    solute_rad = numpy.array(((sr),
                              (sr)), dtype=numpy.float64)

    solvent_pos = numpy.array(((11,0,0),), dtype=numpy.float64)
    solvent_rad = 1.4
    voxel_len = 0.5

    vol = volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, 
                              voxel_len)
    
    # volume of sphere 0 + volume of sphere 1 - volume of intersection
    h = (solute_rad[0]+solvent_rad)/2
    r = solute_rad[0]+solvent_rad
    analytical_vol = numpy.pi*4./3*r**3 \
                   + numpy.pi*4./3*r**3 \
                   - 2*1./3*numpy.pi*h**2*(3*r-h)

    print("Test: two overlapping atoms of radius 3")
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
    vol = pdb2volume.PDBVolume('villin.pdb', 
                                     'cavity.lib.autogen', voxel_len=0.5).run()
    print(vol)

def show_protein_surface():
    vol, grid = pdb2volume.PDBVolume('villin.pdb', 
                                     'cavity.lib.autogen', voxel_len=0.5).run()
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
    #test_simple()
    test_2sphere()
