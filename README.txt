===dedalusBatch.job===
A batch script used on Bridges. Calls SSX_model_A.py to simulate spheromak evolution, then calls merge1.py and merge2.py to merge data together into an h5 file  named merged.h5.

===mergeBatch.job===
A batch script used on Bridges. Calls merge1.py and merge2.py to merge data into a merged file. Use in case dedalusBatch.job fails to merge (out of time). --cleanup is true by default.

===SSX_model_A.py===
The Dedalus script used to set up and simulate spheromak evolution. Requires spheromak.py to be in same directory.
Parameters--
Eta: Magnetic diffusivity
Mu: Viscosity*rho0
Kappa: Temperature diffusivity*rho0
Gamma: Ideal gas adiabatic index

===Spheromak.py===
Dependency of SSX_model_A.py. Sets spheromak initial current density. On top of the regular spheromak functions there's a tanh mask, because the functions are periodic, and we only want one spheromak.

===merge1.py===
Merges data from seperate processessors into coherent, time ordered file.
Usage:
    merge1.py <base_path> [--cleanup]
Options:
    --cleanup   Delete distributed files after merging

base_path is a directory containing the files to be merged

===merge2.py===
Merges data from seperate time-ordered h5 files (output of merge1.py), and creates one large, merged file.
Usage:
    merge2.py <joint_path> <set_path> [--cleanup]

Options:
    --cleanup   Delete distributed files after merging

joint_path is file name of finished, merged file
set_path is a directory containing files to be merged

===toVapor.py===
Script to create a VAPOR-readable .vdf file from an h5 file produced by Dedalus. Uses netCDF as an intermediate format.
Pulls all variables output by Dedalus and captures time evolution if present.
Usage:
    toVapor.py <fileIn> <fileOutName> [<dimensionRatio>]

fileIn is the relevant h5 file
fileOutName is the name (no extension) of the resultant file (eg. "ssx")
dimensionRatio is an optional argument specifying the desired visualized dimensions. Syntax: "minZ:minY:minX:maxZ:maxY:maxX"

FOR HELP regarding getting data into VAPOR more generally:
see Eric Hester's comment here https://groups.google.com/forum/#!topic/dedalus-users/2tS6PS-zKLM

===particlePuncture.py===
Makes a puncture plot of a particle trajectory given information about where to take the slice, and with what precision.

===particleStaticFieldSpheromak.py===
Produces netCDF file with data about particle orbit and fields. Use ncdfvdfcreate and ncdf2vdf to make VAPOR-readable vdf files.
For more info, see the above link.
