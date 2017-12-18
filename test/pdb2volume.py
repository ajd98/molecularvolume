#!/usr/bin/env python
import numpy
import volume

class PDBVolume(object):
    def __init__(self, pdbpath, radiipath):
        self.loadradii(radiipath)

    def loadradii(radiipath):
        with open(radiipath, 'r') as f:
            for line in f:
                if line.startswith("SPAT"):
                    # 6:10
                    resfield = line[6:10]

                    # 10:14
                    atomfield = line[10:14]

                    # 16:21
                    radius = float(line[16:21])


