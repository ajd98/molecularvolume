# load_radii.tcl
#
# Written 18.05.05 by Alex DeGrave
#
# This tcl script for VMD updates atomic radii to be consistent with a radii.lib
# file. 
#
# To use, open the tcl console in VMD and type "source load_radii.tcl". After
# loading your molecule in VMD, type the following in the tcl console:
#
#     updateRadii "<path to radii.lib>" <molecule id>
#
# where <path to radii.lib> is the file path to the radii.lib file, and 
# <molecule id> is the molid of the molecule for which the van der Waals
# radii should be updated. For example, you might type:
#
#     updateRadii "radii.lib" 0
#

proc updateRadii {radiusfile molid} {
    set f [open $radiusfile r]
    set contents [read $f]
    close $f
    set lines [split $contents "\n"]
    foreach line $lines {

        if {[regexp {[a-zA-Z]} $line]} {
            set fields [regexp -all -inline {\S+} $line]
            set resname [lindex $fields 0]
            set atomname [lindex $fields 1]
            set radius [lindex $fields 2]

            # Update resname and atomname to to play nicely with + and - signs
            set resname [regsub {([a-zA-Z0-9]+)([+-])} $resname {"\1\\\2"}]
            set atomname [regsub {([a-zA-Z0-9]+)([+-])} $atomname {"\1\\\2"}]

            puts "Working on residue $resname"
            puts "  Setting selection 'resname $resname and name $atomname' to radius $radius"
            set sel [atomselect $molid "resname $resname and name $atomname"]
            $sel set radius $radius
            }
        }
    }
}
