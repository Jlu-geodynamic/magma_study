# Files used in the study: 
Lithospheric Mantle Heterogeneity Drives Delayed Magmatism and Wide Continent-Ocean Transitions in Rifted Margins.

Yuan Wang, Zhonglan Liu


Model runs in this study used dealii 9.2.0.

The software version used for this study is found at: 

ASPECT 2.3

        https://github.com/Djneu/aspect/tree/fault_analysis
	
FastScape commit 18f2588

(elimated)

        https://github.com/fastscape-lem/fastscapelib-fortran 
	
(backup)

        https://github.com/Jlu-geodynamic/fastscape_backup_commit18f2588

# Additional ASPECT plugins statement:

1."comp" and "rift" were modified from the plugins used in the study: 


Evolution of rift systems and their fault networks in response to surface processes

Derek Neuharth, Sascha Brune, Thilo Wrona, Anne Glerum, Jean Braun, Xiaoping Yuan 


We change some initial setting to fit our study.
Note: The initial topography function of this part is under development, it is not used in this study. 


2."grain" were modified from the sorce code. In this study, we mainly added mantle melting function(including hydrous and unhydrous melting), parameterized melt extraction and crustal accretion.

TODO:Notes in the plugin will be updated in the future.

# To install ASPECT with FastScape

1.Get ASPECT and FastScape from the sites above.

        git clone https://github.com/Djneu/aspect/tree/fault_analysis
        git clone https://github.com/Jlu-geodynamic/fastscape_backup_commit18f2588

2.Create a build directory for fastscape and compile it with an added flag for creating a shared library.

        cmake -DBUILD_FASTSCAPELIB_SHARED=ON /path/to/fastscapemake
        make

3.Compile ASPECT with a flag FASTSCAPE_DIR pointing to the fastscape build folder with the shared library

        cmake -DFASTSCAPE_DIR=/path/to/fastscape/build -DDEAL_II_DIR=/path/to/dealii -DCMAKE_BUILD_TYPE=Release /path/to/aspect
        make

