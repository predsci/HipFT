<img width=500 src="doc/hipft_logo.png" alt="HipFT" />  
  
# HipFT: High-performance Flux Transport
    
[Predictive Science Inc.](https://www.predsci.com)  
 
--------------------------------  
  
## OVERVIEW ##
<img align="right" src="doc/hipft_example.gif" alt="HipFT Example">

HipFT is a flux transport model written in modern Fortran that is used as the computational core of the [Open-source Flux Transport (OFT)](https://github.com/predsci/oft) software suite.  OFT is a complete system for generating full-Sun magnetograms through acquiring & processing observational data, generating realistic convective flows, and running the flux transport model.  
  
HipFT implements advection, diffusion, and data assimilation on the solar surface on a logically rectangular nonuniform spherical grid.  It is parallelized for use with multi-core CPUs and GPUs using a combination of Fortran's standard parallel `do concurrent` (DC), OpenMP Target data directives, and MPI.  It uses high-order numerical methods such as SSPRK(4,3), Strang splitting, WENO3-CS(h), and the super time-stepping scheme RKG2.  The code is designed to be modular, incorporating various differential rotation, meridional flow, super granular convective flow, and data assimilation models.  It can also compute multiple realizations in a single run spanning multiple choices of parameters.  
  
HipFT can be used with MacOS, Linux, and Windows (through WSL) on CPUs and GPUs (NVIDIA or Intel).  
  
--------------------------------  
   
## HOW TO BUILD HIPFT ##
  
HipFT has been tested to work using GCC's `gfortran` (>8), Intel's `ifx` (>21), or NVIDIA's `nvfortran` (>21.5) compilers.  Note that it is NOT compatible with the older Intel `ifort` compiler.  
It is recommended to use the latest compiler version available.  

HipFT requires the [HDF5](https://www.hdfgroup.org/solutions/hdf5) library.  
When building for GPUs, the library must be compiled by the same compiler HipFT is using (e.g. nvfortran).  
Often, this requires building the library from source. 
  
The included `build.sh` script will take a configuration file, generate a Makefile, and build the code.  
The folder `conf` contains example configuration files for various compilers and systems.  
We recommend copying the configuration file closest to your setup and then modifying it to conform to your compiler and system (such as `HDF5` library paths/flags, compiler flags, etc.).  
  
Given a configure script `conf/my_custom_build.conf`, the build script is invoked as:  
```
> ./build.sh ./conf/my_custom_build.conf
```  

GCC:  If you are using the GCC compiler (`gfortran`), you must set the  environment variable `OMP_NUM_THREADS` to the derised number of CPU threads to use per MPI process before building the code.
  
### RUN THE HIPFT TESTSUITE ###
  
To test if the installation is working, we recommend running the testsuite after installation.  
To do this, enter the `testsuite/` directory and run:  
  
`./run_test_suite.sh`  
  
This will run the tests using 1 MPI rank.  
  
IMPORTANT:  If you are building/running HipFT on a multi-core CPU, you will most likely need to  
use the `-mpicall` option to the `run_test_suite.sh` script to set the proper thread affinity.  
For example:  For OpenMPI, one would likely want to use `-mpicall="mpirun --bind-to socket -np"`.
  
--------------------------------  
  
## HOW TO RUN HIPFT ##
  
### Setting Input Options  
  
HipFT uses a namelist in an input text file.  The default name for the input is `hipft.in`  
  
A full working input file with all the default parameter options is provided in the file:  
  
`doc/hipft.in.documentation` 
   
A detailed description of each parameter is also given in that file, and (in addition to this readme) is the current main documentation of the code.  
  
We have also provided example input files for use cases in the `examples/` folder as well as in the testsuite.  
  
### Launching the Code ###
    
To run HipFT, set the desired run parameters in a file (e.g. `hipft.in`), then copy or link the `hipft` executable into the same directory as the input file and run the command:  
  
`<MPI_LAUNCHER> <MPI_OPTIONS> -np <N> ./hipft <input_file>`  
  
where `<N>` is the total number of MPI ranks to use, `<MPI_LAUNCHER>` is your MPI run command (e.g. `mpiexec`,`mpirun`, `ibrun`, `srun`, etc), and `<MPI_OPTIONS>` are additional MPI options that may be needed (such as `--bind-to socket` or `--bind-to numa` for CPUs running with OpenMPI).  

For example:  `mpirun -np 1 ./hipft hipft.in`  
  
The MPI ranks split up the number of realizations to work on.  
Therefore, if you are only running 1 realization, you must use only 1 MPI rank.  
The number of ranks cannot be larger than the number of realizations.  
  
The code is parallelized with Fortran `do concurrent` and OpenMP Target data directives within each MPI rank.  
  
### Running HipFT on CPUs ###
  
On CPUs, the code is multi-threaded for each MPI rank.  This can require proper setting of the `OMP_NUM_THREADS` and `ACC_NUM_CORES` environment variables (and for GCC, setting them before compilation).  
It also requires properly setting the thread affinity in the launch of MPI as shown above.  
For example, running HipFT compiled for GCC and OpenMPI on 4 compute nodes with dual-socket 64-core EPYC CPUs (setup in the BIOS as 4 NUMA domains per socket and no hyperthreading) with more than 16 realizations could be compiled with `OMP_NUM_THREADS=16` and launched with:  
  
`mpirun --bind-to numa --ntasks-per-node 8 ./hipft hipft.in 1>hipft.log 2<hipft.err`  
  
A simpler example of running on a single desktop (1 socket), with `OMP_NUM_THREADS` having been set to the total number of threads (before compilation for GCC, at run time with NV and IFX):
  
`mpirun --bind-to socket -np 1 ./hipft hipft.in 1>hipft.log 2<hipft.err` 
  
Depending on the system setup, it may be difficult to actualize the full possibly performance on CPU nodes.  We therefore highly recommend running HipFT on GPUs.  
  
### Running HipFT on GPUs ###
  
For standard cases, the code should be launched with the number of MPI ranks per node being equal to the number of GPUs per node (assuming you are running at least that many realizations). 
e.g.  
`mpiexec --ntasks-per-node 4 ./hipft 1>hipft.log 2<hipft.err`  
or  
`mpiexec --npersocket 2 ./hipft 1>hipft.log 2<hipft.err`  
  
Having more data on a GPU typically makes it more efficient.  Therefore for a given number of realizations, it is recommended to try different combinations of numbers of GPUs (e.g. 4 realizations on 2 GPUs (2 per GPU) versus 4 realizations on 1 GPU).  
  
### Solution Output ###
  
The output of HipFT are HDF5 map files in phi-theta coordinates.  
When running with multiple realizations, the output becomes a 3D file with one dimension along realization number.
In the `bin/` folder, we have provided python scripts for reading in the data and plotting it.
  
### Processing Scripts ###
  
The `/bin` folder contains several python and bash scripts that can be used to post process and plot results of a HipFT run.  Full documentation on their use is pending, however each script has some level of documentation within it.  Check back here for a list of common commands to process the runs.

--------------------------------

## Sample Data for Convective Flows and Data Assimilation
  
To use HipFT with data assimilation and convective flows uses data generated from the [MagMAP](https://github.com/predsci/magmap) and [ConFlow](https://github.com/predsci/conflow) code packages respectively.  
A sample data set from each code is provided in the zenodo data set:  
  
[HipFT Sample Input Dataset for Convective Flows and Data Assimilation](https://zenodo.org/doi/10.5281/zenodo.10271120)
  
This data package contains a Carrington rotation of convective flows at 15 minute cadence as well as a year of HMI-derived data assimilation maps for 2022.  
We have also provided an example HipFT input file for running a full year simulation using this data in this repo's `examples/flux_transport_1yr_flowCAa_diff175_data_assim/` folder.  
Note that the flow files are auto-repeated in HipFT when run for longer than a Carrington rotation, so they can be used for arbitrary length runs of HipFT.  

--------------------------------
  
