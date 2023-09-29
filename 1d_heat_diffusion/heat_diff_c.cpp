/*To solve heat diffusion equation using FVM. AS described in Versteeg H and 
!Malalasekra 

!To set constant temperature boundary condition on both ends*/
#include <petsc.h>
#include <vector> 
#include <iostream> 

//void subroutine bc_const_temp_both_sides(A, b, temp_l, temp_r, k, area, dx, Nx, q){

//}


PetscErrorCode bc_const_temp_both_sides(Mat A, Vec b, PetscReal temp_l, PetscReal temp_r,
    PetscReal k, PetscReal area, PetscReal dx, PetscInt Nx, PetscReal q);

static char help[] = "\
                      Solves the 1D heat diffusion equation.\n \
                      Usage: \n";


int main(int argc, char **argv){

  PetscInt Nx, Istart, Iend, col_idx[3], rank, num_procs;
  PetscReal dx, xl, xr, start_time, anal_term1, anal_term2, start_tot_time; 
  std::vector<PetscReal> grid;
  //PetscScalar, pointer :: soln(:);
  PetscScalar k, ap, ae, aw, su, sp, area, ka_dx, Tlb, Trb, val[3], zero, q;
  Vec x, b, u, vout;
  Mat A;
  KSP ksp; 
  PC pc ;
  PetscViewer viewer ;
  VecScatter ctx;
  PetscBool flg ;
  const PetscScalar *soln;
  
  PetscCall(PetscInitialize(&argc, &argv, (char *)0, help));
  //PetscCall(InitializeOptions(&user));

  //start_tot_time = MPI_Wtime();
  //Default simulation parameters (source term zero)
  xl=0.0;
  xr=0.5;
  Nx=100; //total number of grid points 
  k=1000.0; //thermal conductivity  (W/(m.k))
  area=0.01; //(m^2)
  Tlb = 100.0; //temperature on the left boundary Tlb.....Trb (k)
  Trb = 500.0; //temperature on the right boundary (k)
  q=0.0; //heat flux (W/m^3) 
  
 // Default simulation parameters (with source term)
// xl=0.0
// xr=0.02
// Nx=5 !total number of grid points 
// k=0.5 !thermal conductivity  (W/(m.k))
// area=1.0 !(m^2)
// Tlb = 100.0 !temperature on the left boundary Tlb.....Trb (k)
// Trb = 200.0 !temperature on the right boundary (k)
// q=1.0e6 ! heat flux (W/m^3) 

  //get command line inputs
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-Nx", &Nx, &flg));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-xl", &xl, &flg));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-xr", &xr, &flg));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-k", &k, &flg));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-area", &area, &flg));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-Tl", &Tlb, &flg));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-Tr", &Trb, &flg));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-q", &q, &flg));

  MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  /*std::cout << "Size and ranks are: " << num_procs << "\t" << rank 
    << std::endl;*/

  dx = (xr-xl)/(Nx-1);

  for(int i=0; i<Nx; i++){
    grid.push_back(xl+i*dx);
  }

  /*for(auto elem : grid) {
    std::cout << elem << "\n";
  }*/

  //create solution vector, rhs and vector to store analytical soln
  PetscCall(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, Nx, &x));
  PetscCall(VecSetFromOptions(x));
  PetscCall(VecDuplicate(x,&b));
  PetscCall(VecDuplicate(x,&u));
  
  //create a Matrix to store coefficients
  PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
  PetscCall(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, Nx, Nx));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatSetUp(A));

  //Get a global indices of rows of matrix for each processor
  PetscCall(MatGetOwnershipRange(A, &Istart, &Iend));

  std::cout << "The start and end indices of mat for each rank: " << Istart 
    << "\t" << Iend << "\n";
  
  //coefficients for the interior of the domain
  ka_dx = k*area/dx;
  sp = 0.0;
  aw = ka_dx;
  ae = ka_dx; 
  ap = aw+ae-sp;

  //Insert the coefficients in to the matrix for the interior of the domain
  val[0]=-aw; val[1]=ap; val[2]=-ae;
  for(PetscInt i=Istart+1; i<Iend+1; i++){
    col_idx[0]=i-1; col_idx[1]=i; col_idx[2]=i+1;
    if(i>0 && i<Nx-1) {
      PetscCall(MatSetValues(A, 1, &i, 3, col_idx, val, INSERT_VALUES));
    }
  }

  //initialize rhs with source term 
  su = q*area*dx;        
  PetscCall(VecSet(b, su));
  
  //Apply boundary conditions
  PetscCall(bc_const_temp_both_sides(A, b, Tlb, Trb, k, area, dx, Nx, q));

  //assemble matrix and rhs 
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscCall(VecAssemblyBegin(b));
  PetscCall(VecAssemblyEnd(b));
  
  //to view (print to terminal) assembled matrix 
  //PetscCall(MatView(A, PETSC_VIEWER_STDOUT_WORLD));
  //to view assembled vector 
  //PetscCall(VecView(b, PETSC_VIEWER_STDOUT_WORLD));
  

  //exact solution 
  //T=c1*x+c2, c1=(Trb-Tlb)/(xr-xl), c2=Tlb
  //call VecGetOwnershipRange(u, Istart, Iend, ierr)
  for(PetscInt i=Istart; i<Iend; i++){
    anal_term1 = (Trb-Tlb)/(xr-xl);
    anal_term2 = (q/(2.0*k))*((xr-xl)-grid.at(i));
    PetscCall(VecSetValue(u, i, grid.at(i)*(anal_term1+anal_term2)+Tlb, 
          INSERT_VALUES));
  }
  //create linear solver context
  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp)); 
  PetscCall(KSPSetOperators(ksp, A, A));
  PetscCall(KSPSetFromOptions(ksp));
  
  start_time = MPI_Wtime();

  PetscCall(KSPSolve(ksp, b, x));

  std::cout << "Total time taken for KSP solve: rank: " << rank << "\t" 
    << MPI_Wtime()-start_time << "\n";

  //assemble solution matrix 
  PetscCall(VecAssemblyBegin(x));
  PetscCall(VecAssemblyEnd(x));
 

  //Scatter the solution from all the procs to root rank Vec (vout) 
  PetscCall(VecScatterCreateToZero(x, &ctx, &vout));
  PetscCall(VecScatterBegin(ctx, x, vout, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecScatterEnd(ctx, x, vout, INSERT_VALUES, SCATTER_FORWARD));

  //Print solution 
  /*if(rank==0){
    PetscCall(VecGetArrayRead(vout, &soln));
    for(int i=0; i<Nx; i++){
      std::cout << soln[i] << "\n";
    }
  }*/

  PetscCall(PetscFinalize());
  return 0;
}

PetscErrorCode bc_const_temp_both_sides(Mat A, Vec b, PetscReal temp_l, 
    PetscReal temp_r, PetscReal k, PetscReal area, PetscReal dx, PetscInt Nx, 
    PetscReal q){
  
  PetscReal val[2], ka_dx, sp, su, ae, aw, ap;
  PetscMPIInt rank, num_procs;
  PetscInt col_idx[2];
  
  MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  //BC on the left side of the domain
  ka_dx=k*area/dx;
  sp=-2.0*ka_dx;
  if(rank==0){  
    su=q*area*dx+2.0*ka_dx*temp_l; 
    aw=0.0; ae=ka_dx; ap=ae+aw-sp; 
    val[0]=ap; val[1]=-ae;
    col_idx[0]=0; col_idx[1]=1;
    PetscInt i=0;
    PetscCall(MatSetValues(A, 1, &i, 2, col_idx, val, INSERT_VALUES));
    PetscCall(VecSetValue(b, 0, su, INSERT_VALUES));
  }

  if(rank==num_procs-1){ //right boundary 
    su=q*area*dx+2.0*ka_dx*temp_r; 
    aw=ka_dx; ae=0.0; ap=ae+aw-sp; 
    val[0]=-aw; val[1]=ap;
    col_idx[0]=Nx-2; col_idx[1]=Nx-1;
    auto temp = Nx-1;
    PetscCall(MatSetValues(A, 1, &temp, 2, col_idx, val, INSERT_VALUES));
    PetscCall(VecSetValue(b, Nx-1, su, INSERT_VALUES));
  }
}


