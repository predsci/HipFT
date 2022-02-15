![HipFT](hipft_logo.png)

# HipFT: High-performance Flux Transport 
Predictive Science Inc.  
www.predsci.com

## OVERVIEW ##
  
HipFT is the computational core of the Open-source Flux Transport (OFT) software suite, which is a data-assimilation flux transport model used to generate an ensemble of synchronic radial magnetic field maps for use as boundary conditions for the coronal field models.   
  
HipFT implements advection, diffusion, and data assimilation for the solar surface on a logically rectangular non-uniform spherical grid.  It is written in Fortran and parallelized for use with multi-core CPUs and GPUs using a combination of OpenACC/MP directives and Fortran's standard parallel `do concurrent`.  To alleviate the strict time-step stability criteria for the diffusion equation, we use an extended stability Runge-Kutta super time-stepping algorithm.  The code is designed to be modular, incorporating various differential rotation, meridianal flow, super granular convective flow, and data assimilation models.  Multiple realizations of the evolving flux will be computed in parallel using MPI in order to produce an ensemble of model outputs for uncertainty quantification.  
