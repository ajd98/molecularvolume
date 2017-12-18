#!/usr/bin/env python
import pdb2volume
import numpy
import matplotlib.pyplot as pyplot
from mpl_toolkits.mplot3d import Axes3D

if __name__ == "__main__":
    vol, grid = pdb2volume.PDBVolume('villin.pdb', 
                                     'cavity.lib.autogen', voxel_len=0.2).run()
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
