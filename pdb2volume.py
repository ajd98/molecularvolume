#!/usr/bin/env python
import argparse
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
    def __init__(self, pdbpath, radiipath, explicitsolvent=False, solventname='WAT', solventrad=1.4, voxel_len=0.5):
        self.use_explicit_solvent = explicitsolvent
        self.solventname = solventname
        self.solventrad = 1.4
        self.voxel_len = voxel_len
        self.radii = {}
        self.loadradii(radiipath)
        self.pdb = PDBStructure()
        self.pdb.from_file(pdbpath)

    def loadradii(self, radiipath):
        with open(radiipath, 'r') as f:
            for line in f:
                # 0:3 
                resfield = line[0:3]

                # 4:7
                atomfield = line[4:8]

                # 9:15
                radius = float(line[8:14])
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
        

        if self.use_explicit_solvent:
            return volume.volume_explicit_sol(solute, solute_rad, solvent, 
                                              self.solventrad, self.voxel_len)
        else:
            return volume.volume(solute, solute_rad, self.solventrad, 
                                 self.voxel_len)

class PDB2VolumeTool(PDBVolume):
    def __init__(self):
        self._parse_args()
        super(PDB2VolumeTool, self).__init__(self.args.pdbpath, 
                                             self.args.radiipath,
                                             explicitsolvent=self.args.explicitsolvent,
                                             solventname=self.args.solventname,
                                             solventrad=self.args.solventrad,
                                             voxel_len=self.args.voxel_len)
        v = self.run()
        print(u"Volume: {:.01f} \u00c5^3")


    def _parse_args(self):
        parser = argparse.ArgumentParser()

        parser.add_argument("-p", "--pdbpath", dest='pdbpath', 
                            type=str,
                            help="Path to the input PDB file")

        parser.add_argument("-r", "--radiipath", dest="radiipath",
                            type=str,
                            default="radii.lib",
                            help="Path to the radii library file")

        parser.add_argument("-e", "--explicit-solvent", dest="explicitsolvent", 
                            action="store_true",
                            help="Use explicit water positions to account for "
                            "void spaces filled with water. Void spaces with "
                            "internal waters will not count toward the "
                            "molecular volume")

        parser.add_argument("-sr", "--solvent-radius", dest="solventrad", default=1.4,
                            type=float,
                            help="The radius of the solvent molecules (which "
                            "are assummed to be spherical)")

        parser.add_argument("-sn", "--solvent-name", dest="solventname", default="WAT",
                            type=str,
                            help="The residue name of the solvent (used only if "
                            "-e/--explicit-solvent is specified)")

        parser.add_argument("-v", "--voxel-len", dest="voxel_len", default=0.1, 
                            type=float,
                            help="The edge length of the (cubic) voxels.")
        self.args = parser.parse_args()

if __name__ == "__main__":
    PDB2VolumeTool()
