import volume
import numpy
import cProfile
import pstats

def main():
    print("testing...")
    solute_pos = numpy.array(((0,0,0),), dtype=numpy.float64)
    solute_rad = numpy.array((1,), dtype=numpy.float64)
    solvent_pos = numpy.array(((5,0,0),), dtype=numpy.float64)

    solvent_rad = 1.4
    voxel_len = 1.0
    vol = volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, voxel_len)
    print(vol)
    
    voxel_len = 0.9
    vol = volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, voxel_len)
    print(vol)

    voxel_len = 0.8
    vol = volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, voxel_len)
    print(vol)

    voxel_len = 0.5
    vol = volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, voxel_len)
    print(vol)

    cProfile.runctx("volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, voxel_len)",
                    globals(),
                    locals(),
                    "Profile.prof")
    s = pstats.Stats("Profile.prof")
    s.strip_dirs().sort_stats("time").print_stats()

if __name__ == "__main__":
    main()
