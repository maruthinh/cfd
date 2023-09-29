**Introduction**

Solves 1D heat diffusion equation with and without source term (by default)
1. Finite Volume Method (FVM) 
2. PETSc for solving a linear system
3. Can run in parallel
4. Analytical solution is also calculated and printed to the Result.csv file

 
**How to compile and run the code** 

 1. Need Fortran compiler and PETSc library 
 2. to run the executable without source term (heat flux)
    
    `mpiexec.hydra -n 2 ./heat_diff -pc_type gamg`
 4. To run with the source term (heat flux) use the following command line options. By default, it runs without heat flux 

    `mpiexec.hydra -n 2 ./heat_diff -Nx 100 -xl 0.0 -xr 0.02 -q 1000000.0 -area 1.0 -k 0.5 -Tl 100.0 -Tr 200.0`

 
**Results**

Results can be plotted using the jupyter-notebook present in this folder

