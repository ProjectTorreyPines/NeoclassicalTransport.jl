!-----------------------------------------------------------------
! neo_serial.f90
!
! Serial (no-MPI) variant of gacode's neo/src/neo.f90 main program,
! linked out-of-tree against the unmodified neo_lib.a. It never calls
! MPI_INIT, so it runs anywhere a plain process can run: on login
! nodes (where the `neo -e` wrapper may dispatch via srun/mpirun and
! fail outside an allocation) and as a child of a Slurm step task
! (where the stock MPI binary joins the parent step's PMI world and
! blocks in MPI_INIT). Same role as the "serial-modified NEO
! executable" used on NERSC.
!
! Build: make       (in this directory, with the GACODE env sourced)
! Use  : ENV["NEO_EXECUTABLE"] = ".../neo_serial" makes run_neo
!        invoke it directly instead of the `neo -e` wrapper.
!-----------------------------------------------------------------
program neo

  use mpi
  use neo_globals

  implicit none

  integer, external :: omp_get_max_threads

  !----------------------------------------------------------------
  ! Query OpenMP for threads
  !
  n_omp = omp_get_max_threads()

  ! Strictly serial: rank 0 of a world of 1; NEO_COMM_WORLD is set to
  ! the MPI_COMM_WORLD constant for completeness but is never used
  ! because no MPI routine is ever called (n_radial=1 local runs).
  i_proc = 0
  n_proc = 1
  NEO_COMM_WORLD = MPI_COMM_WORLD

  ! Path is cwd:
  path= './'

  call neo_read_input
  call neo_do

end program neo
