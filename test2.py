import volume
import numpy
import matplotlib.pyplot as pyplot
from mpl_toolkits.mplot3d import Axes3D

def main():
    #print("testing...")
    #solute_pos = numpy.array(((0,0,0),), dtype=numpy.float64)
    #solute_rad = numpy.array((3,), dtype=numpy.float64)
    #solvent_pos = numpy.array(((0,10,0),), dtype=numpy.float64)

    #solvent_rad = 1.4
    #voxel_len = 0.5

    #vol, grid = volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, voxel_len)

    #fig = pyplot.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #w = numpy.where(grid==0)
    #ax.scatter(w[0], w[1], w[2], marker='o')
    #ax.set_xlim(0,40)
    #ax.set_ylim(0,40)
    #ax.set_zlim(0,40)
    #ax.set_xlabel('x')
    #ax.set_ylabel('y')
    #ax.set_zlabel('z')
    #print("volume: ", vol)
    #pyplot.show()

    #pyplot.clf()

    solute_pos = numpy.array(((0,0,0),
                              (12,0,0),
                              (0,-12,0)), dtype=numpy.float64)

    solute_rad = numpy.array((3,
                              3,
                              2), dtype=numpy.float64)

    solvent_pos = numpy.array(((0,10,0),), dtype=numpy.float64)
    solvent_rad = 1.4
    voxel_len = 0.5

    vol, grid = volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, voxel_len)
    print("volume: ", vol)
    print("actual: ", 2*4./3*numpy.pi*(4.4**3)+4./3*numpy.pi*(3.4**3))

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
    main()
