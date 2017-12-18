#!/usr/bin/env python
import numpy
import sys
sys.path.append('../')
import volume

class PDBStructure(object):
    def __init__(self):
        pass

    def from_file(self, pdbpath):
        self.atoms = []
        with open(pdbpath, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    self.atoms.append({'atomname': line[12:16],
                                       'resname': line[17:20],
                                       'position': numpy.fromstring(line[32:54], 
                                                           dtype=numpy.float64, 
                                                           sep=' ')
                                       }
                                      )

class PDBVolume(object):
    def __init__(self, pdbpath, radiipath, solventname='WAT', solventrad=1.4, voxel_len=0.5):
        self.solventname = solventname
        self.solventrad = 1.4
        self.voxel_len = 0.5
        self.radii = {}
        self.loadradii(radiipath)
        self.pdb = PDBStructure()
        self.pdb.from_file(pdbpath)
        self.run()

    def loadradii(self, radiipath):
        with open(radiipath, 'r') as f:
            for line in f:
                if line.startswith("SPAT"):
                    # 6:10 (if wildcard is included)
                    resfield = line[6:9]

                    # 10:14
                    atomfield = line[10:14]

                    # 16:21
                    radius = float(line[16:21])
                    try:
                        self.radii[resfield][atomfield] = radius
                    except KeyError:
                        self.radii[resfield] = {}
                        self.radii[resfield][atomfield] = radius

    def run(self):
        # make numpy arrays for the solvent and solute
        solute = [] 
        solvent = []
        solute_rad = []
        for atom in self.pdb.atoms:
            if atom['resname'] == self.solventname:
                solvent.append(atom['position'])
            else:
                try:
                    solute_rad.append(self.radii[atom['resname']][atom['atomname']])
                    solute.append(atom['position'])
                except KeyError:
                    print("no radius specified for residue, atom name ({:s}, {:s}). "\
                          "Excluding from volume calculation".format(atom['resname'],
                                                                     atom['atomname']))

        solute = numpy.array(solute, dtype=numpy.float64)
        solvent = numpy.array(solvent, dtype=numpy.float64)
        solute_rad = numpy.array(solute_rad, dtype=numpy.float64)
        

        vol, grid = volume.volume(solute, solute_rad, solvent, self.solventrad, self.voxel_len)
        print(vol)

if __name__ == "__main__":
    PDBVolume('villin.pdb', 
              'cavity.lib.autogen')
