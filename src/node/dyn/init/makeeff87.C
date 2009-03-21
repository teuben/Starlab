
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Construct an unscaled EFF87 model, with a desired cut-off
//// radius. The new model system is written to standard output.  The
//// model system is shifted to place its center of mass at rest at
//// the origin of coordinates.  Unscaled systems will be in
//// approximate virial equilibrium, based on the continuum limit.
////
//// Usage:   makeeff87 [OPTIONS]
////
//// Options:
////          -c    add a comment to the output snapshot [false]
////          -C    output data in 'col' format [no]
////          -i    number the particles sequentially [don't number]
////          -n    specify number of particles [no default]
////          -o    echo value of random seed [don't echo]
////          -r    specify radius cutoff in Rc [default 33.76]
////          -s    specify random seed [random from system clock]
////          -g    set gamma3D value [no default]
////          -T    dump to stdout: f(E), 
////                                a(E)   = C1 \int_0^E f(e)     de 
////                                b(E)   = C1 \int_0^E f(e)  e  de 
////                                G(E,r) = C2 \int_0^E f(e) v^2 de = a(E)*phi(r) - b(E) 
////
//// Written by Evghenii Gaburov.
////
//// Report makeeff87 bugs to egaburov@strw.leidenuniv.nl
//// Report  Starlab  bugs to starlab@sns.ias.edu.

//............................................................................
//   version 1:    March 2009  Evghenii Gaburov               
//   version 1.1:  March 2009  Evghenii Gaburov               
//   version 1.2:  March 2009  Evghenii Gaburov               
//                           email: egaburov@strw.leideuniv.nl
// ---
//  Reference: Binney & Tremaine, p. 236
//             Elson, Freeman & Fall '86
//............................................................................

#include "dyn.h"
#include <cassert>

#ifdef TOOLBOX

#define  SEED_STRING_LENGTH  60


struct spline {
  vector<double> x, y, y2;
  void init(vector<double> &xvec, 
	    vector<double> &yvec, 
	    bool flag = false, 
	    double yp0 = 0, 
	    double ypN = 0) {
    vector<double> u;
    
    x = xvec;
    y = yvec;
    int n = x.size();
   
    y2.resize(n);
    u.resize(n);
    if (!flag) {
      y2[0] = u[0] = 0.0;
    } else {
      y2[0] = -0.5;
      u [0] = 3.0/(x[1]-x[0])*((y[1]-y[0])/(x[1]-x[0]) - yp0);
    }
    for (int i = 1; i < n-1; i++) {
      assert(x[i  ] - x[i-1] > 0);
      assert(x[i+1] - x[i-1] > 0);
      double sig = (x[i] - x[i-1])/(x[i+1] - x[i-1]);
      double p   = sig*y2[i-1] + 2.0;
      y2[i]      = (sig - 1.0)/p;
      u [i]      = (y[i+1] - y[i])/(x[i+1] - x[i]) - (y[i] - y[i-1])/(x[i] - x[i-1]);
      u [i]      = (6.0*u[i]/(x[i+1] - x[i-1]) - sig*u[i-1])/p;
    }
    double qn, un;
    if (!flag) {
      qn = un = 0.0;
    } else {
      qn = 0.5;
      un = (3.0/(x[n-1] - x[n-2]))*(ypN - (y[n-1] - y[n-2])/(x[n-1] - x[n-2]));
    }
    y2[n-1]     = (un - qn*u[n-2])/(qn*y2[n-2] + 1.0);
    for (int k = n-2; k >= 0; k--)
      y2[k] = y2[k]*y2[k+1] + u[k];
  };

  double eval(double xval) {
    int n = x.size();
    int klo = 0;
    int khi = n-1;
    while (khi - klo > 1) {
      int k = (khi + klo) >> 1;
      if (x[k] > xval) khi = k;
      else             klo = k;
    }
    double h = x[khi] - x[klo];
    assert(h > 0.0);
    double a = (x[khi] - xval)/h;
    double b = (xval - x[klo])/h;
    double v = a*y[klo] + b*y[khi] + (a*(a*a - 1.0)*y2[klo] + b*(b*b - 1.0)*y2[khi])*(h*h)/6.0;
    return v;
  }

  double integrate(double a, double b) {
    double c  = (a + b)/2.0;
    double h3 = (b - a)/6.0;
    return h3*(eval(a) + eval(b) + 4.0*eval(c));
  }
  
  spline integrate() {
    int n = x.size();
    vector<double> yint(n);
    yint[0] = 0.0;
    for (int i = 1; i < n; i++)
      yint[i] = yint[i-1] + integrate(x[i-1], x[i]);
    spline spl;
    spl.init(x, yint);
    return spl;
  }

};


struct eff87 {
  double gamma, rc, rt, rcut, Rhm;
  int nbins;
  double phi_infty;
  spline spl_tmp, spl_mass, spl_phi;
  spline spl_r_phi, spl_r_mass;
  spline spl_fE, spl_aE, spl_bE;
  

  vector<double> r, tmp;
  eff87(double g, double Rcut = 100, int n = 1000) {init(g, Rcut, n); }

  double rho(double x)    { return pow(1+x*x, -0.5*gamma) - pow(1+rcut*rcut, -0.5*gamma); }
  double drho(double x)   { return  -x*gamma*pow(1+x*x,-0.5*gamma - 1); }
  double d2rho(double x)  { return -gamma*pow(1+x*x,-0.5*gamma - 2)*(1-x*x*(gamma+1)); }
  double mass(double x)   { return spl_mass.eval(x); }
  double radius(double m) { return spl_r_mass.eval(m); }
  double phi(double x)    { return spl_phi.eval(x); }
  double dphi(double x)   { return spl_mass.eval(x)/(x*x);}
  double f(double E)      { return exp(spl_fE.eval(E)); }
  double a(double E)      { return spl_aE.eval(E); }
  double b(double E)      { return spl_bE.eval(E); }

  double G(double E, double r) { return phi(r)*a(E) - b(E); }

  void compute_fgE(double E,
		   double &f,
		   double &a, double &b) {
    int n_int = 100;
    double e0 = 0;
    double e1 = E;

    double de = (e1 - e0)/n_int;
    e0 += 0.5*de;
    
    vector<double> elist(n_int), fElist(n_int), aElist(n_int), bElist(n_int);
    double e = e0;
    for (int i = 0; i < n_int; i++) {
      elist [i] = e;
      double x = spl_r_phi.eval(e);
      double g = 2.0 - 4*M_PI*x*x*x*rho(x)/mass(x);
      g = (d2rho(x) + g/x*drho(x))/(dphi(x)*dphi(x));
      fElist[i] =   g/sqrt(E - phi(x));
      aElist[i] = 3*g*sqrt(E - phi(x));
      bElist[i] =   g*sqrt(E - phi(x))*(E + 2*e);
      e += de;
    }
    spline tmp;
    tmp.init(elist, fElist); 
    tmp = tmp.integrate();
    f = tmp.eval(E);

    tmp.init(elist, aElist); 
    tmp = tmp.integrate();
    a = tmp.eval(E);

    tmp.init(elist, bElist);
    tmp = tmp.integrate();
    b = tmp.eval(E);
  }

  void compute_fg() {
    int n = 1000;

    double e0 = 0;
    double e1 = phi_infty;
    
    e0 = log(1.0e-3);
    e1 = log(phi_infty);
    
    double de = (e1 - e0)/n;
    e0 += de*0.5;

    vector<double> elist(n), flist(n), alist(n), blist(n);
    double e = e0; 
    for (int i = 0; i < n; i++) {
      elist[i] = exp(e);
      compute_fgE(exp(e), flist[i], alist[i], blist[i]);
      e += de;
      flist[i] = log(flist[i]);
    }
    spl_fE.init(elist, flist);
    spl_aE.init(elist, alist);
    spl_bE.init(elist, blist);
  }

  double compute_enclosed_mass(double Rcut) {
    
    int n = nbins/2;
    r.resize(nbins);
    tmp.resize(nbins);
    double dr = rt/n;
    for (int i = 0; i < n; i++)
      r[i] = 0.5*rc/n + i*dr;
    
    double r1 = log(rt);
    double r2 = log(Rcut);
    dr = (r2 - r1)/n;
    r1 += 0.5*dr;
    for (int i = 0; i < n; i++)
      r[i + n] = exp(r1 + i*dr);

    
    // compute enclosed mass, m(r)
    for (int i = 0; i < nbins; i++)
      tmp[i] = 4*M_PI*r[i]*r[i]*rho(r[i]);
    spl_tmp.init(r, tmp);
    spl_mass = spl_tmp.integrate();
    spl_r_mass.init(spl_mass.y, spl_mass.x);

    return spl_mass.y[spl_mass.y.size()-1];
  }

  void init(double g, double Rcut, int n) {
    gamma = g;
    nbins = 2*n;

    fprintf(stderr, "Generating EFF87 profile:\n");
    fprintf(stderr, "  gamma3D= %lg\n", gamma);
    
    // build r-array

    rc     = sqrt(pow(2.0, 2.0/gamma) - 1);
    rt     = 2*rc;
    Rcut  *= rc;

    rcut  = Rcut;

    fprintf(stderr, " ** Computing enclosed mass, m(r) & r(m)\n");
    double Mtot = compute_enclosed_mass(Rcut);

    // compute potential, phi(r)
    
    fprintf(stderr, " ** Computing potential, phi(r) & r(phi)\n");
    for (int i = 0; i < nbins; i++)
      tmp[i] = dphi(r[i]);
    spl_tmp.init(r, tmp);
    spl_phi   = spl_tmp.integrate();
    phi_infty = spl_phi.y[nbins-1]; //eval(rcut);
    phi_infty = spl_phi.eval(rcut);
    for (int i = 0; i < nbins; i++) {
      tmp[i] = phi_infty - spl_phi.y[i];
    }
    spl_phi.init(r, tmp);
    // compute r(phi)

    vector<double> rinv(nbins);
    for (int i = 0; i < nbins; i++) {
      int k = nbins - i - 1;
      rinv[i] = r[k];
      tmp [i] = spl_phi.y[k];
    }
    spl_r_phi.init(tmp, rinv);


    fprintf(stderr, " ** Computing distribution functions, f(E), a(E) & b(E)\n");
    compute_fg();
    Rhm = radius(Mtot*0.5);
    fprintf(stderr, "     Rc/Rh  = %lg \n",  rc/Rhm);
    fprintf(stderr, "     Rt/Rc  = %lg \n",  Rcut/rc);
    fprintf(stderr, "     Rt/Rh  = %lg \n",  Rcut/Rhm);
    fprintf(stderr, "     Mc/M    = %lg\n",  mass(rc)/Mtot);
    fprintf(stderr, "     M(Rt)  = %lg\n",  Mtot);
    fprintf(stderr, "     phi(0) = %lg\n", phi_infty);
  };
 
};


local void makeeff87(dyn *b, int n, double Rcut, real gamma3D, bool dump_fg) {
  eff87 m(gamma3D, Rcut, 1000);
  Rcut = m.rcut;
  
  if (dump_fg) {
    fprintf(stderr, " ** Dumping f(E) & g(E) into stdout \n");
    int n = m.spl_fE.x.size();
    for (int i = 0; i < n; i++)
      fprintf(stdout, "i= %d:  E= %lg  f(E)= %lg a(E)= %lg  b(E)= %lg\n",
	      i, m.spl_fE.x[i], m.f(m.spl_fE.x[i]),
	      m.spl_aE.y[i], m.spl_bE.y[i]);
    return;
  }
  fprintf(stderr, " ** Generating positions & velocities \n");
  int  i;
  
  real  partmass;		   // equal mass of each particle	     
  real  radius;		   // absolute value of position vector      
  real  velocity;		   // absolute value of velocity vector      
  real  theta, phi;		   // direction angles of above vectors      
  real  x, y;		           // for use in rejection technique         
  real  scalefactor;             // for converting between different units 
  real  inv_scalefactor;         // inverse scale factor                   
  real  sqrt_scalefactor;        // sqare root of scale factor             
  real  mrfrac;                  // m( rfrac )                             
  real  m_min, m_max;            // mass shell limits for quiet start      
  dyn * bi;
  
  b->set_mass(1.0);
  partmass = 1.0 / ((real) n);

  double Mmax = m.mass(Rcut);
  
  for (i = 0, bi = b->get_oldest_daughter(); i < n;
       i++, bi = bi->get_younger_sister()) {
    
    bi->set_mass(partmass);
    
    // The position coordinates are determined by inverting the cumulative
    // mass-radius relation, with the cumulative mass drawn randomly from
    // [0, mfrac]; cf. Aarseth et al. (1974), eq. (A2).
    
    real rrrr = randinter(0, Mmax);
    radius = m.radius(rrrr);
    
    theta = acos(randinter(-1.0, 1.0));
    phi = randinter(0.0, TWO_PI);
    
    bi->set_pos(vec(radius * sin( theta ) * cos( phi ),
		    radius * sin( theta ) * sin( phi ),
		    radius * cos( theta )));

    // The velocity coordinates are determined using von Neumann's
    // rejection technique, cf. Aarseth et al. (1974), eq. (A4,5).
    // First we take initial values for x, the ratio of velocity and
    // escape velocity (q in Aarseth et al.), and y, as a trick to
    // enter the body of the while loop.
    
#if 0
    // wrong one
    double phi  = m.phi(0);

    double Gmax = phi*m.a(phi) - m.b(phi);
    double G    = randinter(0.0, Gmax);
    double Elo = 0;
    double Ehi = phi;
    while (Ehi - Elo > 0.0001*Ehi) {
      double Etest = 0.5*(Ehi + Elo);
      double Gtest = phi*m.a(Etest) - m.b(Etest);
      if (Gtest > G) Ehi = Etest;
      else           Elo = Etest;
    }
    double E = 0.5*(Elo + Ehi);
    velocity = sqrt(2.0*(phi - E));
#else
    double phi  = m.phi(radius);
    double fmax = m.f(phi);
    x = 0.0;
    y = 1.0;
    while (y*fmax > x*x*m.f(phi*(1-x*x))) {
      x = randinter(0.0, 1.0);
      y = randinter(0.0, 1.0);
    }
    velocity = x*sqrt(2.0*phi);
#endif

      
    theta = acos(randinter(-1.0, 1.0));
    phi = randinter(0.0,TWO_PI);
    
    bi->set_vel(vec(velocity * sin( theta ) * cos( phi ),
		    velocity * sin( theta ) * sin( phi ),
 		    velocity * cos( theta )));
  }

  real xfac = 0.8/m.Rhm;
  real vfac = 1.0/sqrt(xfac);
  
  for_all_daughters(dyn, b, bi) {
    bi->set_pos(xfac*bi->get_pos());
    bi->set_vel(vfac*bi->get_vel());
  }


  b->to_com();
  
  putrq(b->get_log_story(), "initial_mass", 1.0);
    
}

local void swap(dyn* bi, dyn* bj)
{
    if (bi != bj) {

	real mass  = bi->get_mass();
	vec pos = bi->get_pos();
	vec vel = bi->get_vel();

	bi->set_mass(bj->get_mass());
	bi->set_pos(bj->get_pos());
	bi->set_vel(bj->get_vel());

	bj->set_mass(mass);
	bj->set_pos(pos);
	bj->set_vel(vel);
    }
}

// Original reshuffling code was O(N^2)!  Fixed by Steve, July 2000.

local void reshuffle_all(dyn* b, int n)
{
    // Reorder the tree by swapping mass, pos, and vel of each node
    // with a randomly chosen node.

    // First make a list of nodes.

    int i = 0;
    dynptr *list = new dynptr[n];
    for_all_daughters(dyn, b, bi)
	list[i++] = bi;

    // Then go down the list and randomize.

    for (i = 0; i < n; i++)
	swap(list[i], list[(int)randinter(0, n)]);

    delete [] list;
}

//-----------------------------------------------------------------------------
//  main  --
//      usage:
//                makeeff87 -n# -g# [options]  ,
// 
//            where -n# is the number of bodies in the Nbody system.
//                  -g# is the value of 3D gamma of EFF87
//    options:
//            The following options are allowed:
//       seed:
//            -s #   where # stands for an integer either in the range
//                   [1,2147483647], or # = 0 which causes a seed to be
//                   chosen as the value of the UNIX clock in seconds; this
//                   guarantees that no two calls will give the same value for
//                   the seed if they are more than 2.0 seconds apart.
//     cutoff:
//            -r #   radius fraction of an (infinitely extended) EFF87 model
//-----------------------------------------------------------------------------

main(int argc, char ** argv) {
  int  i;
  int  n;
  int  input_seed, actual_seed;
  real  rfrac;
  real  gamma3D = -1;
  
  bool  g_flag = false;
  bool  c_flag = false;
  bool  C_flag = false;
  bool  i_flag = false;
  bool  m_flag = false;
  bool  n_flag = false;
  bool  o_flag = false;
  bool  r_flag = false;
  bool  R_flag = true;
  bool  s_flag = false;
  bool  u_flag = false;
  bool  T_flag = false;
  
  char  *comment;
  char  seedlog[SEED_STRING_LENGTH];
  
  check_help();
  
  extern char *poptarg;
  int c;
  const char *param_string = "c:Cin:os:r:Rug:T";
  
  rfrac = 33.76;
  
  while ((c = pgetopt(argc, argv, param_string,
		      "$Revision$", _SRC_)) != -1)
    switch(c) {
      
    case 'c': c_flag = true;
      comment = poptarg;
      break;
    case 'C': C_flag = true;
      break;
    case 'i': i_flag = true;
      break;
    case 'n': n_flag = true;
      n = atoi(poptarg);
      break;
    case 'o': o_flag = true;
      break;
    case 'r': r_flag = true;
      rfrac = atof(poptarg);
      break;
    case 'R': R_flag = !R_flag;
      break;
    case 's': s_flag = true;
      input_seed = atoi(poptarg);
      break;
    case 'u': u_flag = true;
      break;
    case 'T': T_flag = true;
      break;
    case 'g': gamma3D = atof(poptarg);
      g_flag = true;
      break;
    case '?': params_to_usage(cerr, argv[0], param_string);
      get_help();
      exit(1);
    }            
  
  if (!n_flag) {
    cerr << "makeeff87: specify the number # of";
    cerr << " particles with -n#\n";
    exit(1);
  }
  if (!g_flag) {
    cerr << "makeeff87: specify gamma3D with -g#\n";
    exit(1);
  }
    
  dyn *b, *by, *bo;
  b = new dyn();
  b->set_root();
  
  if (C_flag) b->set_col_output(true);

  if (n > 0) {
    bo = new dyn();
    if (i_flag)
      bo->set_label(1);
    b->set_oldest_daughter(bo);
    bo->set_parent(b);
  }
  
  for (i = 1; i < n; i++) {
    by = new dyn();
    if (i_flag)
      by->set_label(i+1);
    by->set_parent(b);
    bo->set_younger_sister(by);
    by->set_elder_sister(bo);
    bo = by;
  }
  
  if (c_flag)
    b->log_comment(comment);
  b->log_history(argc, argv);
  
  if (!s_flag)
    input_seed = 0;                         	// default
  actual_seed = srandinter(input_seed);
  
  if (o_flag) cerr << "makeeff87: random seed = " << actual_seed << endl;
  
  sprintf(seedlog, "       random number generator seed = %d",actual_seed);
  b->log_comment(seedlog);
  
  if (n != 0) {
    makeeff87(b, n,  rfrac,  gamma3D, T_flag);
  }
  
  if (!T_flag)  put_dyn(b);

  rmtree(b);
}

#endif

