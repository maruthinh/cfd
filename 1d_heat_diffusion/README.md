Solves 1D heat diffusion equation 
1. Finite Volume Method (FVM) 
2. PETSc for solving linear system
3. Can run in parallel
 
 **How to compile and run the code** 
 1. Need fortran compiler and PETSc library 
 2. to run the executable 
    a. mpiexec.hydra -n 2 ./heat_diff -pc_type gamg
