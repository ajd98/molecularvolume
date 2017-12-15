import volume
import numpy
import matplotlib.pyplot as pyplot
from mpl_toolkits.mplot3d import Axes3D

def main():
    print("testing...")
    solute_pos = numpy.array(((0,0,0),), dtype=numpy.float64)
    solute_rad = numpy.array((3,), dtype=numpy.float64)
    solvent_pos = numpy.array(((0,10,0),), dtype=numpy.float64)

    solvent_rad = 1.4
    voxel_len = 0.5

    vol, grid = volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, voxel_len)

    fig = pyplot.figure()
    ax = fig.add_subplot(111, projection='3d')
    w = numpy.where(grid==0)
    ax.scatter(w[0], w[1], w[2], marker='o')
    ax.set_xlim(0,40)
    ax.set_ylim(0,40)
    ax.set_zlim(0,40)
    pyplot.show()

    #pyplot.savefig('test.pdf')

if __name__ == "__main__":
    main()
