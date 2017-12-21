# molecularvolume

## Purpose
This software implements voxel-based volume calculations for molecular systems,
via a Python interface.

The volume is defined as the region of space that is not accessible to solvent 
molecules; more exactly, the volume is calculated as the sum of the volumes of
those voxels such that a (spherical) solvent molecule of radius `r` positioned
at the center of the voxel does not intersect any solute atom `i`, which is
represented as a sphere of radius `s<sub>i</sub>` and position `p<sub>i</sub>`. 
The solvent radius `r`, and the radius `s<sub>i</sub>` and position 
`p<sub>i</sub>` of each solute atom is specified by the user. A flood-fill 
algorithm (think of the paint bucket in MS Paint, except in three dimensions) 
is used to determine which voxels are accessible to solvent, such that void 
space in the interior of the solute is included as part of the volume.


## Requirements: 
Software requirements:
* gcc (including libm)
* Python (version 2.7 or 3.6)

Required Python packages:
* Numpy 
* Cython
* argparse
* distutils

Hardware requirements:
* 12+ GB ram

## Installation
After cloning this git repository, navigate to the `./src` directory via a
command line, and run `make install`.

Following installation, navigate to the `./test` directory, make sure that
`python` is in your `$PATH`, and run `python test.py`.  All of the tests should
pass; note that a large amount of memory (8+ GB) is needed for some of the tests.

## Use

Multiple python functions are available:

`volume_explicit_sol`: Positions of solvent molecules must be explicitely 
provided, and the flood-fill algorithm extends from the position of every
solvent molecule.  Voids in the solute with clathrate solvent are thus not
included in the volume, while voids not occupied by solvent are included in the
volume.

`volume`: Assumes all provided atomic coordinates correspond to solute; the 
flood-fill algorithm extends from the region exterior to the provided atoms.
Since the mesh of voxels need not be large enough to include solvent molecules,
this function will be faster than `volume_explicit_sol`.
