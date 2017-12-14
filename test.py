import volume
import numpy

def main():
    print("testing...")
    solute_pos = numpy.array(((0,0,0),), dtype=numpy.float64)
    solute_rad = numpy.array((1,), dtype=numpy.float64)
    solvent_pos = numpy.array(((5,0,0),), dtype=numpy.float64)
    solvent_rad = 1.4
    voxel_len = 0.1
    vol = volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, voxel_len)
    print(vol)

if __name__ == "__main__":
    main()
