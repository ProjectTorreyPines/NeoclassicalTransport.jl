!-----------------------------------------------------------------
! neo_dump.f90
!
! Reference-data driver for the native Julia NEO port
! (src/neo/ in NeoclassicalTransport.jl).
!
! Runs a full serial NEO simulation (exactly like neo_serial, so all the
! usual out.neo.* files are produced; out.neo.f is written at full
! precision because neo_do_dump.f90 -- the stock neo_do.f90 with the
! g(:) output format widened -- is linked in front of neo_lib.a), then
! rebuilds the per-radius setup (profiles, energy basis, equilibrium,
! rotation, collision matrices) and writes every intermediate array the
! Julia port reproduces to out.neo.dump at full double precision.
!
! Format of out.neo.dump: blocks of
!   # <name> <n1> [<n2> ...]
!   one value per line, Fortran (column-major) element order
!
! Build: make neo_dump   (in this directory, GACODE env sourced)
! Use  : run in a directory holding input.neo.gen (from neo_parse.py)
!-----------------------------------------------------------------

program neo_dump

  use mpi
  use neo_globals
  use neo_energy_grid
  use neo_equilibrium
  use neo_rotation

  implicit none

  integer, external :: omp_get_max_threads
  integer, parameter :: io = 77
  integer, parameter :: m0_fcoll = 16
  integer, parameter :: n_lam = 7
  character(len=*), parameter :: fmt = '(es24.16e3)'
  real, dimension(n_lam) :: lam_set
  real, dimension(:,:), allocatable :: f, fb
  integer :: it, is, i

  n_omp = omp_get_max_threads()
  i_proc = 0
  n_proc = 1
  NEO_COMM_WORLD = MPI_COMM_WORLD
  path = './'

  call neo_read_input

  ! 1) full run: all standard outputs, g(:) -> out.neo.f (full precision)
  call neo_do
  if (error_status > 0) then
     print '(a)', trim(error_message)
     stop 1
  endif

  ! 2) neo_do deallocated its state on exit; rebuild the setup pieces
  call neo_make_profiles
  if (error_status > 0) stop 2
  call neo_check
  if (error_status > 0) stop 3

  allocate(thcyc(1-n_theta:2*n_theta))
  do it=1,n_theta
     thcyc(it-n_theta) = it
     thcyc(it) = it
     thcyc(it+n_theta) = it
  enddo
  cderiv(-2) =  1
  cderiv(-1) = -8
  cderiv(0)  =  0
  cderiv(1)  =  8
  cderiv(2)  = -1

  call ENERGY_basis_ints_alloc(1)
  call ENERGY_basis_ints
  call ENERGY_coll_ints_alloc(1)
  call EQUIL_alloc(1)
  call ROT_alloc(1)
  call EQUIL_do(1)
  call ROT_solve_phi(1)
  if (error_status > 0) stop 4
  call ENERGY_coll_ints(1)

  open(unit=io, file=trim(path)//'out.neo.dump', status='replace')

  ! ---- switches and sizes
  call wint('n_species', n_species)
  call wint('n_energy', n_energy)
  call wint('n_xi', n_xi)
  call wint('n_theta', n_theta)
  call wint('collision_model', collision_model)
  call wint('equilibrium_model', equilibrium_model)
  call wint('rotation_model', rotation_model)
  call wint('ae_flag', ae_flag)
  call wint('laguerre_method', laguerre_method)
  call wint('is_ele', is_ele)
  call wint('n_row', n_species*(n_energy+1)*(n_xi+1)*n_theta)

  ! ---- profiles (as used by the kinetic solve, after make_profiles/rotation)
  call w1('z', z, n_species)
  call w1('mass', mass, n_species)
  call w1('dens', dens(:,1), n_species)
  call w1('temp', temp(:,1), n_species)
  call w1('dlnndr', dlnndr(:,1), n_species)
  call w1('dlntdr', dlntdr(:,1), n_species)
  call w1('nu', nu(:,1), n_species)
  call w1('vth', vth(:,1), n_species)
  call w0('r', r(1))
  call w0('rmaj', rmaj(1))
  call w0('q', q(1))
  call w0('rho', rho(1))
  call w0('sign_q', sign_q)
  call w0('sign_bunit', sign_bunit)
  call w0('dphi0dr', dphi0dr(1))
  call w0('epar0', epar0(1))
  call w0('omega_rot', omega_rot(1))
  call w0('omega_rot_deriv', omega_rot_deriv(1))
  call w0('dens_ae', dens_ae(1))
  call w0('temp_ae', temp_ae(1))
  call w0('dlnndr_ae', dlnndr_ae(1))
  call w0('dlntdr_ae', dlntdr_ae(1))

  ! ---- energy/xi basis
  call wi1('e_lag', e_lag, n_xi+1)
  call wi1('xi_beta_l', xi_beta_l, n_xi+1)
  call w1('mygamma2', mygamma2, size(mygamma2))
  call w2('evec_e0', evec_e0, n_energy+1, n_xi+1)
  call w2('evec_e1', evec_e1, n_energy+1, n_xi+1)
  call w2('evec_e2', evec_e2, n_energy+1, n_xi+1)
  call w2('evec_e05', evec_e05, n_energy+1, n_xi+1)
  call w2('evec_e105', evec_e105, n_energy+1, n_xi+1)
  call w4('emat_e05', emat_e05, n_energy+1, n_energy+1, n_xi+1, 2)
  call w4('emat_en05', emat_en05, n_energy+1, n_energy+1, n_xi+1, 2)
  call w4('emat_e05de', emat_e05de, n_energy+1, n_energy+1, n_xi+1, 2)
  call w4('emat_e0', emat_e0, n_energy+1, n_energy+1, n_xi+1, 1)
  call w4('emat_e1', emat_e1, n_energy+1, n_energy+1, n_xi+1, 3)

  ! ---- collision matrices (species,species,energy,energy,xi)
  call w5('emat_coll_test', emat_coll_test, n_species, n_species, &
       n_energy+1, n_energy+1, n_xi+1)
  call w5('emat_coll_field', emat_coll_field, n_species, n_species, &
       n_energy+1, n_energy+1, n_xi+1)

  ! ---- equilibrium
  call w1('theta', theta, n_theta)
  call w0('d_theta', d_theta)
  call w1('k_par', k_par, n_theta)
  call w1('v_drift_x', v_drift_x, n_theta)
  call w1('v_drift_th', v_drift_th, n_theta)
  call w1('gradr', gradr, n_theta)
  call w1('gradpar_gradr', gradpar_gradr, n_theta)
  call w1('w_theta', w_theta, n_theta)
  call w1('Btor', Btor, n_theta)
  call w1('Bpol', Bpol, n_theta)
  call w1('Bmag', Bmag, n_theta)
  call w1('Bmag_rderiv', Bmag_rderiv, n_theta)
  call w1('gradpar_Bmag', gradpar_Bmag, n_theta)
  call w1('bigR', bigR, n_theta)
  call w1('bigR_rderiv', bigR_rderiv, n_theta)
  call w1('gradpar_bigR', gradpar_bigR, n_theta)
  call w0('bigR_th0', bigR_th0)
  call w0('bigR_th0_rderiv', bigR_th0_rderiv)
  call w0('gradr_th0', gradr_th0)
  call w0('Btor_th0', Btor_th0)
  call w0('Bpol_th0', Bpol_th0)
  call w0('Bmag_th0', Bmag_th0)
  call w0('Bmag_th0_rderiv', Bmag_th0_rderiv)
  call w0('ftrap', ftrap)
  call w0('I_div_psip', I_div_psip)
  call w0('Bmag2_avg', Bmag2_avg)
  call w0('Bmag2inv_avg', Bmag2inv_avg)
  call w0('Btor2_avg', Btor2_avg)
  call w0('bigRinv_avg', bigRinv_avg)
  call w0('gradpar_Bmag2_avg', gradpar_Bmag2_avg)

  ! ---- rotation
  call w1('phi_rot', phi_rot, n_theta)
  call w1('phi_rot_deriv', phi_rot_deriv, n_theta)
  call w1('phi_rot_rderiv', phi_rot_rderiv, n_theta)
  call w0('phi_rot_avg', phi_rot_avg)
  call w2('dens_fac', dens_fac, n_species, n_theta)

  ! ---- transport results of the full run (kept in neo_globals)
  ! neo_dke_out(is,:) = pflux, eflux, mflux, eflux-omega*mflux, vpol_th0,
  !                     vtor_th0 + vtor_0order_th0
  ! neo_gv_out(is,:)  = pflux_gv, eflux_gv, mflux_gv, eflux_gv-omega*mflux_gv
  ! neo_dke_1d_out    = jpar, jtor
  call w2('neo_dke_out', neo_dke_out(1:n_species,:), n_species, 6)
  call w2('neo_gv_out', neo_gv_out(1:n_species,:), n_species, 4)
  call w1('neo_dke_1d_out', neo_dke_1d_out, 2)

  close(io)

  ! ---- fcoll tables on a fixed lambda set (separate file; large)
  lam_set(1) = 1.0
  lam_set(2) = 0.05
  lam_set(3) = 0.1
  lam_set(4) = 6.0
  lam_set(5) = 10.0
  lam_set(6) = 3670.0
  lam_set(7) = 1.0/3670.0
  allocate(f(-m0_fcoll:m0_fcoll,-m0_fcoll:m0_fcoll))
  allocate(fb(-m0_fcoll:m0_fcoll,-m0_fcoll:m0_fcoll))
  open(unit=io, file=trim(path)//'out.neo.dump_fcoll', status='replace')
  call wint('m0', m0_fcoll)
  call wint('n_lambda', n_lam)
  call w1('lambda', lam_set, n_lam)
  do i=1,n_lam
     call neo_compute_fcoll(m0_fcoll, lam_set(i), f, fb)
     call w2('fcoll', f, 2*m0_fcoll+1, 2*m0_fcoll+1)
     call w2('fcoll_bar', fb, 2*m0_fcoll+1, 2*m0_fcoll+1)
  enddo
  close(io)
  deallocate(f, fb)

  call ENERGY_basis_ints_alloc(0)
  call ENERGY_coll_ints_alloc(0)
  call EQUIL_alloc(0)
  call ROT_alloc(0)
  deallocate(thcyc)

contains

  subroutine wint(name, v)
    character(len=*), intent(in) :: name
    integer, intent(in) :: v
    write(io,'(a,a,1x,i0)') '# ', name, 1
    write(io,'(i0)') v
  end subroutine wint

  subroutine w0(name, v)
    character(len=*), intent(in) :: name
    real, intent(in) :: v
    write(io,'(a,a,1x,i0)') '# ', name, 1
    write(io,fmt) v
  end subroutine w0

  subroutine w1(name, v, n1)
    character(len=*), intent(in) :: name
    integer, intent(in) :: n1
    real, dimension(n1), intent(in) :: v
    write(io,'(a,a,1x,i0)') '# ', name, n1
    write(io,fmt) v
  end subroutine w1

  subroutine wi1(name, v, n1)
    character(len=*), intent(in) :: name
    integer, intent(in) :: n1
    integer, dimension(n1), intent(in) :: v
    write(io,'(a,a,1x,i0)') '# ', name, n1
    write(io,'(i0)') v
  end subroutine wi1

  subroutine w2(name, v, n1, n2)
    character(len=*), intent(in) :: name
    integer, intent(in) :: n1, n2
    real, dimension(n1,n2), intent(in) :: v
    write(io,'(a,a,2(1x,i0))') '# ', name, n1, n2
    write(io,fmt) v
  end subroutine w2

  subroutine w4(name, v, n1, n2, n3, n4)
    character(len=*), intent(in) :: name
    integer, intent(in) :: n1, n2, n3, n4
    real, dimension(n1,n2,n3,n4), intent(in) :: v
    write(io,'(a,a,4(1x,i0))') '# ', name, n1, n2, n3, n4
    write(io,fmt) v
  end subroutine w4

  subroutine w5(name, v, n1, n2, n3, n4, n5)
    character(len=*), intent(in) :: name
    integer, intent(in) :: n1, n2, n3, n4, n5
    real, dimension(n1,n2,n3,n4,n5), intent(in) :: v
    write(io,'(a,a,5(1x,i0))') '# ', name, n1, n2, n3, n4, n5
    write(io,fmt) v
  end subroutine w5

end program neo_dump
