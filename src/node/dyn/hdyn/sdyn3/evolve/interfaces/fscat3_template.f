
      program fscat_template

c     Template F77 program using scatter3.

c     The "states" defining initial, intermediate, and final
c     configurations of the system are defined in f_scatter3.h.

c     The "f_" and "c_" routines are FORTRAN-callable links to
c     Starlab C++ functions.

      include 'f_scatter3.h'
      include 'f_interface.h'

      iseed  = 0	                ! Seed for random number generator
					! (0 ==> system chooses number)
      n_rand = 0        	        ! Number of times to invoke generator
                                	! before starting "for real" 

      cpu_time_check = 3600     	! Interval for checking CPU time (sec) 
      dt_out = VERY_LARGE_NUMBER 	! Output time interval 
      dt_snap = VERY_LARGE_NUMBER 	! Snapshot time interval 
      snap_cube_size = 10       	! Size of cube for snapshot output 

c     Note: all operations modify the /init/, /inter/, and /final/
c     common arrays (which are kept conveniently hidden from view).

c     Initialization:

      call f_srandinter(iseed, n_rand)
      call f_make_standard_init()
      call f_initialize_angles(0, 0, 0.0d0)  ! Some FORTRANs insist on the "d0"

c     e.g:

      mean_anomaly = f_randinter(0.0d0, 2*PI)
      v_inf = 0.1

c     Perform the scattering calculation. 

      call f_cpu_init()
      call f_scatter3(cpu_time_check, dt_out, dt_snap, snap_cube_size)

c     Print all available information on the initial state. 

      write(6,'(''Initial state:'')')
      write(6,'(''    m1         = '',f10.6)')1-m2 ! convention: m1 + m2 = 1 
      write(6,'(''    m2         = '',f10.6)')m2
      write(6,'(''    m3         = '',f10.6)')m3
      write(6,'(''    r1         = '',f10.6)')r1
      write(6,'(''    r2         = '',f10.6)')r2
      write(6,'(''    r3         = '',f10.6)')r3
      write(6,'(''    sma        = '',f10.6)')init_sma
      write(6,'(''    ecc        = '',f10.6)')init_ecc
      write(6,'(''    v_inf      = '',f10.6)')v_inf
      write(6,'(''    rho        = '',f10.6)')rho
      write(6,'(''    r_init_min = '',f10.6)')r_init_min
      write(6,'(''    r_init_max = '',1p,e14.6)')r_init_max
      write(6,'(''    r_init     = '',f10.6)')r_init
      write(6,'(''    r_stop     = '',1p,e14.6)')r_stop
      write(6,'(''    tidal_tol  = '',1p,e14.6)')tidal_tol
      write(6,'(''    eta        = '',f10.6)')eta
      write(6,'(''    phase:'')')
      write(6,'(''        cos_theta    = '',f10.6)')phase(1)
      write(6,'(''        phi          = '',f10.6)')phase(2)
      write(6,'(''        psi          = '',f10.6)')phase(3)
      write(6,'(''        mean_anomaly = '',f10.6)')phase(4)
      call print_system(
     $        init_index1, init_mass1, init_pos1, init_vel1,
     $        init_index2, init_mass2, init_pos2, init_vel2,
     $        init_index3, init_mass3, init_pos3, init_vel3)

c     Print all available information on the intermediate state. 

      write(6,'(/''Intermediate state:'')')
      write(6,'(''    n_osc      = '',i6)')n_osc
      write(6,'(''    n_kepler   = '',i6)')n_kepler
      write(6,'(''    n_stars    = '',i6)')n_stars
      write(6,'(''    index      = '',3i3)')(index(k),k=1,3)
      write(6,'(''    r_min      = '',1p,3e14.6)')(r_min(k),k=1,3)
      write(6,'(''    r_min_min  = '',1p,e14.6)')r_min_min
      call f_print_intermediate_descriptor(inter_descr)
      call print_system(
     $        inter_index1, inter_mass1, inter_pos1, inter_vel1,
     $        inter_index2, inter_mass2, inter_pos2, inter_vel2,
     $        inter_index3, inter_mass3, inter_pos3, inter_vel3)

c     Print all available information on the final state. 

      write(6,'(/''Final state:'')')
      write(6,'(''    sma              = '',f10.6)')final_sma
      write(6,'(''    ecc              = '',f10.6)')final_ecc
      write(6,'(''    outer_separation = '',1p,e14.6)')outer_sep
      write(6,'(''    escaper          = '',i10)')escaper
      write(6,'(''    error            = '',f10.6)')error
      write(6,'(''    time             = '',1p,e14.6)')time
      write(6,'(''    n_steps          = '',i10)')n_steps
      write(6,'(''    virial_ratio     = '',f10.6)')virial_ratio
      call f_print_final_descriptor(final_descr)
      call print_system(
     $        final_index1, final_mass1, final_pos1, final_vel1,
     $        final_index2, final_mass2, final_pos2, final_vel2,
     $        final_index3, final_mass3, final_pos3, final_vel3)

      write(6,'(/''CPU time = '',f10.6)')f_cpu_time()

      end


      subroutine print_system(index1, mass1, pos1, vel1,
     $        		      index2, mass2, pos2, vel2,
     $        		      index3, mass3, pos3, vel3)

      integer index1, index2, index3
      real*8  mass1, mass2, mass3
      real*8  pos1(3),pos2(3),pos3(3)
      real*8  vel1(3),vel2(3),vel3(3)

      write(6,'(''    system:'')')

      call print_particle(1, index1, mass1, pos1, vel1)
      call print_particle(2, index2, mass2, pos2, vel2)
      call print_particle(3, index3, mass3, pos3, vel3)

      end


      subroutine print_particle(i, index, mass, pos, vel)
      integer i, index
      real*8  mass, pos(3), vel(3)

      write(6,'(''      '',i3,'' ('',i3,'')    mass:'',f14.6)')
     $        i, index, mass
      write(6,'(''                   pos: '',3f14.6)')(pos(k),k=1,3)
      write(6,'(''                   vel: '',3f14.6)')(vel(k),k=1,3)

      end
