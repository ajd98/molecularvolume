import volume
import numpy
import cProfile
import pstats
import matplotlib.pyplot as pyplot

def main():
    print("testing...")
    solute_pos = numpy.array(((0,0,0),), dtype=numpy.float64)
    solute_rad = numpy.array((3,), dtype=numpy.float64)
    solvent_pos = numpy.array(((0,10,0),), dtype=numpy.float64)

    solvent_rad = 1.4
    voxel_len = 0.5
    voxel_sizes = [2,1.9,1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1,1.0,0.9,0.8,0.7,0.6,0.5]
    vols = []
    for voxel_len in voxel_sizes:
        vols.append(volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, voxel_len))

    pyplot.plot(voxel_sizes, vols, color='black')
    actual = sum([numpy.pi*4./3*(solrad+solvent_rad)**3 for solrad in solute_rad])
    pyplot.plot((min(voxel_sizes), max(voxel_sizes)), (actual, actual), color='gray')
    pyplot.savefig('test.pdf')


    #solute_pos = numpy.array(((0,0,0),
    #                          (1,0,0),
    #                          (2,0,0),
    #                          (3,0,0),
    #                          (4,0,0),
    #                          (5,0,0)), dtype=numpy.float64)
    #solute_rad = numpy.array((1,
    #                          1,
    #                          1,
    #                          1,
    #                          1,
    #                          1), dtype=numpy.float64)
    #solvent_pos = numpy.array(((0,5,0),), dtype=numpy.float64)

    #solvent_rad = 1.4
    #voxel_len = 1.0
    #vol = volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, voxel_len)
    #print(vol)
    
    #voxel_len = 0.9
    #vol = volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, voxel_len)
    #print(vol)

    #voxel_len = 0.8
    #vol = volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, voxel_len)
    #print(vol)

    #voxel_len = 0.5
    #vol = volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, voxel_len)
    #print(vol)

    #cProfile.runctx("volume.volume(solute_pos, solute_rad, solvent_pos, solvent_rad, voxel_len)",
    #                globals(),
    #                locals(),
    #                "Profile.prof")
    #s = pstats.Stats("Profile.prof")
    #s.strip_dirs().sort_stats("time").print_stats()

if __name__ == "__main__":
    main()
