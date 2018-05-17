Building PULMsim:
=================


	========================


	A simple makefile is provided that can be edited to compile on any linux/mac machine with a C++ compiler.  First edit the makefile and change the inputs:
	====================================

	Before compilation, ensure that the Eigen include directory $EIGENPATH is added to the include path. In VS2013 this can be found from the toolbar in  �Project -> [project name] Properties -> Configuration Properties -> VC++ Directories� and enter the path for the �Eigen� directory into �include Directories�. This may need to be repeated for Debug and Release versions of the code.
	==================
===============


Outputs
=======

The simulation has four outputs. First, the "flux_*.csv" file(s) contain information for each tree in the simulation. Global properties (total gas volume, concentration at mouth etc.) are outputted into its respective flux file at each time step. Perturbed trees have the suffix "*ltree_[perturbation type]_[lobe region]_[generation of perturbation].csv". These contain just the linearised difference of each quantity compared to the baseline simulation (contained in the "stree" file). 

Second, the "conc*" outputs can either be in .csv or .vtk format depending on the input argument "OUTPUT". These print out a snapshot of the properties on the whole tree. The frequency of these outputs is set by the input "printerval", which is the number of breaths (counting in and out separately) between each output. These files are suffixed with an index corresponding to their chronological order.

Third, the "lung_unit_map*" files related the lung units in the tree flux files to the lung units in the tree flux files.

Finally, the "summary*" file contains summary information about the simulation.