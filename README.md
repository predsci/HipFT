![HipFT](hipft_logo.png)

# HipFT: High-performance Flux Transport 
Predictive Science Inc.  
www.predsci.com

## OVERVIEW ##
  
HipFT is the computational core of the Open-source Flux Transport (OFT) software suite, which is a data-assimilation flux transport model used to generate an ensemble of synchronic radial magnetic field maps for use as boundary conditions for the coronal field models.   
  
HipFT implements advection, diffusion, and data assimilation for the solar surface on a logically rectangular non-uniform spherical grid.  It is written in Fortran and parallelized for use with multi-core CPUs and GPUs using a combination of OpenACC/MP directives and Fortran's standard parallel `do concurrent`.  To alleviate the strict time-step stability criteria for the diffusion equation, we use an extended stability Runge-Kutta super time-stepping algorithm.  The code is designed to be modular, incorporating various differential rotation, meridianal flow, super granular convective flow, and data assimilation models.
  
--------------------------------  
   
## HOW TO BUILD HIPFT ##  

Copy a build script from the `build_examples` folder that is closest to your setup to the base directory.  
Modify the script to set the `HDF5` library paths/flags and compiler flags compatible with your system environment.   
Then, run the script to build `HIPFT` (for example, `./my_build.sh`).  
  
See the multiple build example scripts in the `build_examples` folder for more details.  
  
To use the python post processing scripts, add the `HIPFT` `bin` folder to your `PATH`.  
  
--------------------------------  
  
## HOW TO USE HIPFT ##  
  
### Setting Input Options  
  
`HIPFT` uses a namelist in an input text file called `hipft.dat`.  
  
### Launching the Code ###  
    
To run `HIPFT`, set the desired run parameters in a file (e.g. `hipft.in`), then copy or link the `hipft` executable into the same directory as the input file and run the command:  
  `<MPI_LAUNCHER> -np <N> ./hipft <input_file>`  
where `<N>` is the total number of MPI ranks to use and `<MPI_LAUNCHER>` is your MPI run command (e.g. `mpiexec`,`mpirun`, `ibrun`, `srun`, etc).  
For example:  `mpiexec -np 4 ./hipft hipft.in`  
  
The MPI ranks split up the number of realizations to work on.  
Therefore, if you are only running 1 realization, you should use 1 MPI rank.  
The number of ranks cannot be larger than the number of realizations.  
  
The code is parallelized with OpenACC/OpenMP/StdPar across each MPI rank.  
  
### Running HIPFT on GPUs ###
  
For standard cases, one should launch the code such that the number of MPI ranks per node is equal to the number of GPUs per node (assuming you are running at least that many reaslizations). 
e.g.  
`mpiexec -np <N> --ntasks-per-node 4 ./hipft`  
or  
`mpiexec -np <N> --npersocket 2 ./hipft`  
  
### Solution Output ###  
  
[This section to be filled in later]  
  
### Helpful Scripts ###  
  
Some useful python scripts for reading, extracting, and plotting the HIPFT input/output data can be found in the  `bin` folder.  


