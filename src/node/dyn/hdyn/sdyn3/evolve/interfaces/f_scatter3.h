
c     Definitions to allow FORTRAN programs to use Starlab scatter3.

      implicit real*8 (a-h, o-z)

      parameter (VERY_LARGE_NUMBER = 1.e30,
     $           PI = 3.14159265358979323846)

c     Structures defining the initial, intermediate, and final states of a
c     single scattering:

c     Initial state structure:

      real*8  m2		! mass of secondary in target binary (m1+m2=1) 
      real*8  m3		! mass of incoming projectile 
      real*8  r1		! radius of primary 
      real*8  r2		! radius of secondary 
      real*8  r3		! radius of third star 
      real*8  init_sma          ! inner binary semi-major axis 
      real*8  init_ecc          ! inner binary eccentricity 
      real*8  v_inf		! projectile velocity at infinity  (v_crit) 
      real*8  rho		! projectile impact parameter 
      real*8  r_init_min	! min initial distance projectile to target 
      real*8  r_init_max	! max initial distance projectile to target 
      real*8  r_init		! actual initial separation 
      real*8  r_stop		! max permitted distance escaper to binary 
      real*8  tidal_tol         ! tidal perturbation at start/stop 
      real*8  phase(4)          ! phase angles specifying the initial state:
				!	cos_theta, phi, psi, mean_anomaly
      real*8  eta		! accuracy parameter 

c     This is really inelegant, but it is most convenient for the
c     interface routines to assume that the FORTRAN common blocks are
c     precisely equivalent to the corresponding C/C++ structures.
c     In particular, while it would probably (certainly) be better
c     for the "system" parameters below to be expressed as arrays of
c     indices, mass, etc., the Starlab software treats them as arrays
c     of structures.  We could try to remap between the FORTRAN and
c     and C(++) descriptions using "memcpy," but the alignment of data
c     in memory is important.  Since there is no way to know in advance
c     the conventions used on any given machine, we take the simplest
c     way out and mimic the C structure here.  Placing the indices at
c     the end of the structure here (and in scatter3.h) should make
c     things work right on machines that like to align real*8 and
c     larger structures on 8-byte boundaries.

c     This is what you get for insisting on FORTRAN in the 1990s!
c     (Of course, we could use FORTRAN-90...)

      real*8  init_mass1, init_mass2, init_mass3
      real*8  init_pos1(3), init_pos2(3), init_pos3(3)
      real*8  init_vel1(3), init_vel2(3), init_vel3(3)
      integer init_index1, init_index2, init_index3

      common  /init/ m2, m3, r1, r2, r3, init_sma, init_ecc,
     $               v_inf, rho, r_init_min, r_init_max, r_init,
     $               r_stop, tidal_tol, phase, eta,
     $               init_mass1, init_pos1, init_vel1, init_index1,
     $               init_mass2, init_pos2, init_vel2, init_index2,
     $               init_mass3, init_pos3, init_vel3, init_index3

c     Intermediate state structure:

      integer n_osc		! number of "oscillations" in min(r) 
      integer n_kepler          ! number of analytic continuations 
      integer n_stars		! final number of stars 
      integer index(3)          ! final labels 
      real*8  r_min(3)          ! minimum interparticle separations 
      real*8  r_min_min         ! absolute minimum separation 
      integer inter_descr	! intermediate descriptor

c     (Same comments as above.)

      real*8  inter_mass1, inter_mass2, inter_mass3
      real*8  inter_pos1(3), inter_pos2(3), inter_pos3(3)
      real*8  inter_vel1(3), inter_vel2(3), inter_vel3(3)
      integer inter_index1, inter_index2, inter_index3

      common /inter/ n_osc, n_kepler, n_stars, index,
     $               r_min, r_min_min, inter_descr,
     $               inter_mass1, inter_pos1,
     $        			inter_vel1, inter_index1,
     $               inter_mass2, inter_pos2,
     $        			inter_vel2, inter_index2,
     $               inter_mass3, inter_pos3,
     $        			inter_vel3, inter_index3

c     Final state structure:

      real*8  final_sma         ! final binary semi-major axis 
      real*8  final_ecc         ! final binary eccentricity 
      real*8  outer_sep         ! final separation between third star  
				! and binary (if any)   
      integer final_descr	! final descriptor
      integer escaper		! ID of escaping star (0 if none exists) 
      real*8  error		! rel. energy error (unit=binary energy) 
      real*8  time		! termination time 
      integer n_steps		! number of integration steps 
      real*8  virial_ratio	! final ratio K.E./|P.E.| of outer orbit 

c     (Same comments as above.)

      real*8  final_mass1, final_mass2, final_mass3
      real*8  final_pos1(3), final_pos2(3), final_pos3(3)
      real*8  final_vel1(3), final_vel2(3), final_vel3(3)
      integer final_index1, final_index2, final_index3

      common /final/ final_sma, final_ecc, outer_sep, final_descr,
     $               escaper, error, time, n_steps, virial_ratio,
     $               final_mass1, final_pos1, 
     $               		final_vel1, final_index1,
     $               final_mass2, final_pos2,
     $               		final_vel2, final_index2,
     $               final_mass3, final_pos3,
     $               		final_vel3, final_index3
