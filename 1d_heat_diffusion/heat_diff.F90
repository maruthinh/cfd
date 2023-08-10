!To solve heat diffusion equation using FVM. AS described in Versteeg H and 
!Malalasekra 

!To set constant temperature boundary condition on both ends 
subroutine bc_const_temp_both_sides(A, b, temp_l, temp_r, k, area, dx, Nx, q)
#include <petsc/finclude/petsc.h>
  use petsc 
  implicit none
 
  Mat, intent(inout) :: A
  Vec, intent(inout) :: b
  PetscInt, intent(in) :: Nx 
  PetscReal, intent(in) :: temp_l, temp_r, k, area, dx, q
  PetscReal val(2), ka_dx, sp, su, ae, aw, ap
  PetscMPIInt rank, num_procs
  PetscInt col_idx(2)
  PetscErrorCode ierr 

  call MPI_Comm_size(PETSC_COMM_WORLD, num_procs, ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
  

  !Bc on the left side of the domain 
  ka_dx=k*area/dx
  sp=-2.0*ka_dx
  if(rank==0) then 
    su=q*area*dx+2.0*ka_dx*temp_l 
    aw=0.0; ae=ka_dx; ap=ae+aw-sp 
    val(1)=ap; val(2)=-ae
    col_idx(1)=0; col_idx(2)=1
    call MatSetValues(A, 1, 0, 2, col_idx, val, INSERT_VALUES, ierr)
    call VecSetValue(b, 0, su, INSERT_VALUES, ierr)
  endif
   
  !BC on the right side of the domain  
  if(rank==num_procs-1) then !right boundary 
    su=q*area*dx+2.0*ka_dx*temp_r 
    aw=ka_dx; ae=0.0; ap=ae+aw-sp 
    val(1)=-aw; val(2)=ap
    col_idx(1)=Nx-2; col_idx(2)=Nx-1
    call MatSetValues(A, 1, Nx-1, 2, col_idx, val, INSERT_VALUES, ierr)
    call VecSetValue(b, Nx-1, su, INSERT_VALUES, ierr)
  endif
end subroutine bc_const_temp_both_sides 

program main

#include <petsc/finclude/petsc.h>

  use petsc 
  implicit none

  PetscInt Nx, i, Istart, Iend, col_idx(3), rank, num_procs
  PetscErrorCode ierr 
  PetscReal dx, xl, xr, start_time, anal_term1, anal_term2 
  PetscReal, allocatable, dimension(:) :: grid
  PetscScalar, pointer :: soln(:)
  PetscScalar k, ap, ae, aw, su, sp, area, ka_dx, Tlb, Trb, val(3), zero, q
  Vec x, b, u, vout
  Mat A
  KSP ksp 
  PC pc 
  PetscViewer viewer 
  VecScatter ctx
  PetscBool flg 

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  
  !Default simulation parameters (source term zero)
  xl=0.0
  xr=0.5
  Nx=100 !total number of grid points 
  k=1000.0 !thermal conductivity  (W/(m.k))
  area=0.01 !(m^2)
  Tlb = 100.0 !temperature on the left boundary Tlb.....Trb (k)
  Trb = 500.0 !temperature on the right boundary (k)
  q=0.0 ! heat flux (W/m^3) 
  
  !Default simulation parameters (with source term)
 ! xl=0.0
 ! xr=0.02
 ! Nx=5 !total number of grid points 
 ! k=0.5 !thermal conductivity  (W/(m.k))
 ! area=1.0 !(m^2)
 ! Tlb = 100.0 !temperature on the left boundary Tlb.....Trb (k)
 ! Trb = 200.0 !temperature on the right boundary (k)
 ! q=1.0e6 ! heat flux (W/m^3) 
 
  !get command line inputs
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-Nx', Nx, flg, ierr)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-xl', xl, flg, ierr)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-xr', xr, flg, ierr)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-k', k, flg, ierr)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-area', area, flg, ierr)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-Tl', Tlb, flg, ierr)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-Tr', Trb, flg, ierr)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-q', q, flg, ierr)
  
  !initialize petsc 
  call MPI_Comm_size(PETSC_COMM_WORLD, num_procs, ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)

  allocate(grid(Nx))
  dx=(xr-xl)/(Nx-1)

  !generate grid
  do i=1,Nx
   grid(i) = xl+(i-1)*dx
  enddo
  
  !if(rank==0) then 
  !  print*, grid
  !endif 
  
  start_time = MPI_WTime()
  !create solution vector, rhs and vector to store analytical soln
  call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, NX, x, ierr)
  call VecSetFromOptions(x, ierr)
  call VecDuplicate(x,b,ierr)
  call VecDuplicate(x,u,ierr)
 
  !create a Matrix to store coefficients
  call MatCreate(PETSC_COMM_WORLD, A, ierr)
  call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, Nx, Nx, ierr)
  call MatSetFromOptions(A, ierr)
  call MatSetUp(A, ierr)

  !Get a global indices of rows of matrix for each processor
  call MatGetOwnershipRange(A, Istart, Iend, ierr)

  !coefficients for the interior of the domain
  ka_dx = k*area/dx
  sp = 0.0
  aw = ka_dx 
  ae = ka_dx 
  ap = aw+ae-sp
  
  !Insert the coefficients in to the matrix for the interior of the domain
  val(1)=-aw; val(2)=ap; val(3)=-ae
  do i=Istart+1, Iend
    col_idx(1)=i-1; col_idx(2)=i; col_idx(3)=i+1
    if(i>0 .and. i<Nx-1) then
      call MatSetValues(A, 1, i, 3, col_idx, val, INSERT_VALUES, ierr)
    endif
  enddo

  !initialize rhs with source term 
  su = q*area*dx        
  call VecSet(b, su, ierr)
  
  !Boundary conditions. Set matrix coefficients and set boundary values for the 
  !rhs vector
  call bc_const_temp_both_sides(A, b, Tlb, Trb, k, area, dx, Nx, q)
 
  !assemble matrix and rhs 
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
  call VecAssemblyBegin(b, ierr)
  call VecAssemblyEnd(b, ierr)

  !to view (print to terminal) assembled matrix 
  !call MatView(A, PETSC_ViEWER_STDOUT_WORLD, ierr)
  !to view assembled vector 
  !call VecView(b, PETSC_ViEWER_STDOUT_WORLD, ierr)

  !exact solution 
  !T=c1*x+c2, c1=(Trb-Tlb)/(xr-xl), c2=Tlb
  call VecGetOwnershipRange(u, Istart, Iend, ierr)
  do i=Istart, Iend-1
    call VecSetValue(u, i, (grid(i+1)*((Trb-Tlb)/(xr-xl)))+Tlb, ierr)
  enddo

  !create linear solver context
  call KSPCreate(PETSC_COMM_WORLD, ksp, ierr) 
  call KSPSetOperators(ksp, A, A, ierr)
  call KSPSetFromOptions(ksp, ierr)
  call KSPSolve(ksp, b, x, ierr)

  !assemble solution matrix 
  call VecAssemblyBegin(x, ierr)
  call VecAssemblyEnd(x, ierr)

  !Scatter the solution from all the procs to root rank Vec (vout) 
  call VecScatterCreateToZero(x, ctx, vout, ierr)
  call VecScatterBegin(ctx, x, vout, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEND(ctx, x, vout, INSERT_VALUES, SCATTER_FORWARD, ierr)

  !get the solution from Petsc vector to fortran pointer 
  if(rank==0) then 
    allocate(soln(0:Nx-1))
    call VecGetArrayReadF90(vout, soln, ierr)
  endif
  
  !To view solution, analytical soln and soln (on root rank) vectors 
  !call VecView(x, PETSC_ViEWER_STDOUT_WORLD, ierr)
  ! call VecView(u, PETSC_ViEWER_STDOUT_WORLD, ierr)
  !if(rank==0) then
    ! call VecView(vout, PETSC_ViEWER_STDOUT_WORLD, ierr)
  !  print*, soln
  !endif

  !writing to file 
  if(rank==0) then
    open(unit=10, file='Results.csv', status='unknown', iostat=ierr)
    write(10, '(A1, A1, A6, A1, A10)') 'x', ',', 'result', ',', 'analytical'

    do i=1,Nx 
      anal_term1 = (Trb-Tlb)/(xr-xl)
      anal_term2 = (q/(2.0*k))*((xr-xl)-grid(i))
      write(10, '(ES12.4, A1, ES12.4, A1, ES12.4)') grid(i), ',', soln(i), ',', &
        grid(i)*(anal_term1+anal_term2)+Tlb
    enddo
    close(unit=10)
  endif

  !destroy petsc vars 
  call KSPDestroy(ksp, ierr)
  call MatDestroy(A, ierr)
  call VecDestroy(b, ierr)
  call VecDestroy(u, ierr)
  call VecDestroy(x, ierr)
  call VecDestroy(vout, ierr)
  call VecScatterDestroy(ctx, ierr)

  print*, "Total time taken for matrix and vector creation as well as KSP &
    solve: ", MPI_WTime()-start_time

  call PetscFinalize(ierr)
end program main 
