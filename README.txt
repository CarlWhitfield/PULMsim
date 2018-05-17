Building PULMsim:
=================
Executable files are available in the directory Exec for both linux/mac and windows systems, so if you are only interested in running the code, you can skip straight to Running PULMsim. Otherwise, PULMsim is built like any other C++ program, as follows.
	Install Eigen libraries:
	========================
	The code uses Eigen libraries for linear algebra (v3.3.4). These can be downloaded here. To avoid lookup errors, extract into a directory called Eigen, so that the include (.h) files are contained in $EIGENPATH/Eigen/Eigen/.	Building in linux/mac command line:	===================================

	A simple makefile is provided that can be edited to compile on any linux/mac machine with a C++ compiler.  First edit the makefile and change the inputs:	- CC: The c++ compiler command (e.g. g++).	- CPATH: the path to your C++ compiler (if this is already in $PATH, can be left blank).	- EIGENPATH: The path to eigen libraries	- COPTIONS: other compiler options (e.g. optimisation, warning level etc.)	Compile from the command-line by moving to make directory and entering make.	Building in Microsoft Visual Studio:
	====================================
	The code has also been compiled in Microsoft VS 2013. First, start a new project selecting the type C++ -> Win32 Console Application. Remove any default .hpp or .cpp files created from the solution explorer. Right-click on Header files from solution explorer and add an Existing Item. Select all the files contained in the PULMsim include directory.  Next, right-click on Source files from solution explorer and add an Existing Item. Select all the files contained in the PULMsim src directory.
	Before compilation, ensure that the Eigen include directory $EIGENPATH is added to the include path. In VS2013 this can be found from the toolbar in  Project -> [project name] Properties -> Configuration Properties -> VC++ Directories and enter the path for the Eigen directory into include directories. This may need to be repeated for Debug and Release versions of the code.	Build and run commands can now be operated through the GUI.	Building in Xcode:
	==================	This code has also been compiled in Xcode 9.2. First, start a new project selecting the type Console Application and set the language to C++. Remove any default .hpp or .cpp files created from the finder panel on the left. Right-click on the project folder in the finder panel and click Add files to [project name]. Select and all the files contained in the PULMsim include and src directories. Before compilation, ensure that the Eigen include directory $EIGENPATH is added to the include path. In Xcode this can be found from the finder panel by clicking on project (called [project name], should be the top item in the panel with a little Xcode symbol beside it) and then clicking on Build Settings in the main panel. Scroll down to Search Paths and add the Eigen include path as a new entry in Header Search Paths.	Build and run commands can now be operated through the GUI.Running PULMsim
===============
Executables (PULMsim_linux for linux or mac, and PULMsim.exe for windows) are run from the command line. The command takes a single argument which is the path to an input file. Several examples of input files, details of the options available are given in example files in the Input directory. The executable will run directly in the directory it is called from and will output files there too.
 	Running multiple simulations on linux/mac
	=========================================	
	Several simulations can be run on unix command lines using the bash script repeatrun.sh provided. This script searches for any files in the current directory with the suffix .lung and assumes these are input files. It will then sequentially run a simulation for each input file in a new sub-directory of the current directory. If the executable filename changes from PULMsim_linux then this needs to be updated in repeatrun.sh by changing the argument exec_name.
	Running multiple simulations on Windows
	=======================================
	A windows batch script run_all.bat is also provided which performs the same actions as repeatrun.sh. If the executable filename changes from PULMsim.exe then this needs to be updated in run_all.bat by changing the argument exec_name.Warning: do not use . In input file names as this may cause problems with overwriting data when running multiple jobs.

Outputs
=======

The simulation has four outputs. First, the flux_*.csv file(s) contain information for each tree in the simulation. Global properties (total gas volume, concentration at mouth etc.) are outputted into its respective flux file at each time step. Perturbed trees have the suffix *ltree_[perturbation type]_[lobe region]_[generation of perturbation].csv. These contain just the linearised difference of each quantity compared to the baseline simulation (contained in the stree file). 

Second, the conc* outputs can either be in .csv or .vtk format depending on the input argument OUTPUT. These print out a snapshot of the properties on the whole tree. The frequency of these outputs is set by the input printerval, which is the number of breaths (counting in and out separately) between each output. These files are suffixed with an index corresponding to their chronological order.

Third, the lung_unit_map* files related the lung units in the tree flux files to the lung units in the tree flux files.

Finally, the summary* file contains summary information about the simulation.