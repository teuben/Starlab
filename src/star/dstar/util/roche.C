
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// roche:  Takes binary parameters (semi-major axis, eccentricity, primary
////         and secondary stellar parameters) and constructs a 
////         ps-file of a Roche-lobe with the stars inside.
////        
////         The Roche-lobe is drawn in normalized units (a=1)
////         but scaling is provided in Rsun.
////         The primary (most massive star) is *always* printed on the left
////         When the primary is Roche-lobe filling, its roche-lobe is
////         hashed, otherwise the star is represented with a hashed circle
////         The size of the hashed circle equals the starllar radius, 
////         with respect to the orbital separation.
////         The stellar core is drawn as a solid circle, the size of which 
////         equals the ratio of the stellar mass to the core mass:
////         So the size of the represents its mass, not its actual radius.
////         The color coding is as follows:
////         I tried to make hotter stars bluer and colder star red, but at the
////         moment that is paramterized rather simplistic.
////         The center of mass of the binary system is idicated with
////         a vertical dotted line.
////
////        Note: that when no command line options are provided input
////              via pipe is expected.
////              The format in which this input is given is controlled via
////              the -I option (see below).
////
//// Options:     -t    time in Myr [0]
////              -a    semi-major axis in Rsun [none]
////              -e    eccentricity [0]
////              -f    output option [1]
////                   0: series of ps files with name SeBa_????.ps
////                   1: series of gif files with name SeBa_????.gif
////                   2: single ps file
////                   3: single gif file
////                   4: xwindows screen
////              -M    primary mass in Msun [none]
////              -m    seconday mass in Msun [none]
////              -R    primary radius in Rsun [none]
////              -r    secondary radius in Rsun [none]
////              -C    primary core mass in Msun [none]
////              -c    secondary core mass in Msun [none]
////              -P/p  primary stellar type (see SeBa) [main-sequence]
////              -S/s  secondary stellar type (see SeBa) [main-sequence]
////              -I    designation for the stellar coding used
////                    0 [default]: SeBa (Portegies Zwart \& Verbunt 1997)
////                    1          : BSE  (Hurley et al 2000)
////              -v    verbose output [true]
////
//// Example:     % roche -t 100 -a 100 -r 0.1 -P 7 -M 10 -R 30 -C 3 \   ~
////                      -S 3 -m 4 -r 3 -c 0.4 
////       
//// Todo   :     add command line possibilties for stellar temperature
////              and stellar types.
//++
//++ Some run-time parameters:
//++
//++ The initial virial radius is read from the input snapshot.  If no
//++ initial virial radius is found there, the value specified by the
//++ "-V" command-line option is used.  If no "-V' option is set, a
//++ default value of 1 is assumed.
//++
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//      26 February 2003
//      Version 1.0       S.F. Portegies Zwart <spz@science.uva.nl>
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <math.h>
#include "stdinc.h"       
#include "cpgplot.h"
#include "roche.h"
#include "double_support.h"

#define MIN_RSTAR 2.e-3
#define MAXIT 100

int get_pgcolor_index(stellar_type_summary stp) {
    int color = 1;
    switch(stp) {
    case ZAMS: color = 3;
      break;
    case Early_Giant: color = 13;
      break;
    case Late_Giant:  color = 2;
      break;
    case  Helium_Remnant:  color = 11;
      break;
    case White_Dwarf:  color = 1;
      break;
    case  Neutron_Remnant:  color = 4;
      break;
    case Inert_Remnant:  color = 14;
      break;
    case Unspecified: color = 1;
    };

    return color;
}

int get_pgcolor_index(stellar_type stype) {
    stellar_type_summary stp = summarize_stellar_type((stellar_type)stype);
    return get_pgcolor_index(stp);
}


// calculates outer limit of roche-lobe
local void rlimit(real q, real L, real x, real *f,real *df, real dummy=0) {

  real qi = 1./q;
  real q11=1./(1.+qi);
  real cnst =qi/L+1./(1.-L)+0.5*(1.+qi)*pow(L-q11, 2);

  real r1=abs(x);
  real r2=abs(1-x);
  real r3=abs(x-q11);

  *f=qi/r1+1./r2+0.5*(1.+qi)*pow(r3, 2)-cnst;
  *df=-qi*x/pow(r1, 3)+(1-x)/pow(r2, 3)+(1.+qi)*(x-q11);

  //  PRC(r1);PRC(r2);PRC(r3);PRC(*f);PRL(*df);

  return;
}

local real rtsafe(void (*funcd)(real, real, real, real *, real *, real), 
	     real q, real L, real x1, real x2, real xll, real xacc)
{
 int j;
 real df,dx,dxold,f,fh,fl;
 real temp,xh,xl,rts;

 (*funcd)(q, L, x1,&fl,&df, xll);
 // PRC(q);PRC(x1);PRC(fl);PRL(df);
 (*funcd)(q, L, x2,&fh,&df, xll);
 // PRC(q);PRC(x2);PRC(fh);PRL(df);
 if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
   cerr << "Error occured in rtsafe, exit" << endl;
   exit(-1);
 }
 if (fl == 0.0) return x1;
 if (fh == 0.0) return x2;
 if (fl < 0.0) {
  xl=x1;
  xh=x2;
 } else {
  xh=x1;
  xl=x2;
 }
 rts=0.5*(x1+x2);
 dxold=fabs(x2-x1);
 dx=dxold;
 (*funcd)(q, L, rts,&f,&df, xll);
 for (j=1;j<=MAXIT;j++) {
  if ((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)
   || (fabs(2.0*f) > fabs(dxold*df))) {
   dxold=dx;
   dx=0.5*(xh-xl);
   rts=xl+dx;
   if (xl == rts) return rts;
  } else {
   dxold=dx;
   dx=f/df;
   temp=rts;
   rts -= dx;
   if (temp == rts) return rts;
  }
  if (fabs(dx) < xacc) return rts;
  (*funcd)(q, L, rts,&f,&df, xll);
  if (f < 0.0)
   xl=rts;
  else
   xh=rts;
 }
 exit(-1);
 return 0.0;
}
#undef MAXIT


local real left_limit(real q, real L) {

  real dummy = 0;
  void (*funcd)(real, real, real, real *, real *, real); 
  funcd = rlimit;

  real x1=-0.5*L;
  real x2=-L;
  real rlef = rtsafe(rlimit,q,L,x1,x2,dummy, ROCHE_ACCURACY);

  return rlef;
}

local real right_limit(real q, real L) {

  real dummy = 0;
  void (*funcd)(real, real, real, real *, real *, real); 
  funcd = rlimit;

  //  const=q/x+1./(1.-x)+0.5*(1.+q)*(x-q11)**2

  //  real qi = 1./q;
  real x1=1.5-0.5*L;
  real x2=2.0-L;
  real rrig = rtsafe(rlimit,q,L,x1,x2,dummy, ROCHE_ACCURACY);

  return rrig;
}

// return roche radius as fraction of the semi-major axis.
// So to obtain the true roche_radius call:
// real Rl = semi_major_axis * roche_radius(m1, m2);
// Eggleton PP., ApJ, 1983, 268, 368.
local real roche_radius(real mthis, real mother) {

  real mr = mthis/mother;
  real q1_3 = pow(mr, ONE_THIRD);
  real q2_3 = pow(q1_3, 2);                  //pow(mr, TWO_THIRD);
  
  return 0.49*q2_3/(0.6*q2_3 + log(1 + q1_3));

}

local real first_Largrangian_point(real qinv) {

  real fL, dfL, dL, L, q11;
  q11=1./(1.+qinv);
  L=0.5+0.2222222*log10(qinv);
  do {
    fL = qinv/pow(L, 2)- 1/pow(1-L, 2) - (1.+qinv)*L+1;
    dfL=-2*qinv/pow(L, 3) - 2./pow(1-L, 3) - (1.+qinv);
    dL = -fL/(dfL*L);
    L *= 1.+dL;
  }
  while(abs(dL)>1.e-6);

  return L;
}

local void rline(real q, real L, real y, real *f,real *df, real xl) {
  //real rline(y,f) {
  //c calculates value of y^2 for x^2 value
  //       common/roche/ q,q11,const,const2,xsq,onexsq

  real xsq=pow(xl,2);
  real onexsq=pow(1.-xl, 2);

  real qi=q;
  real q11=1./(1.+qi);
  real cnst =qi/L+1./(1.-L)+0.5*(1.+qi)*pow(L-q11, 2);
  real cnst2=0.5*(1.+qi)*pow(xl-q11,2)-cnst;

  real r1=sqrt(xsq+y);
  real r2=sqrt(onexsq+y);
  *f=qi/r1+1./r2+cnst2;
*df =-0.5*qi/pow(r1, 3)-0.5/pow(r2, 3);

  return;
}

local void compute_lobes(real q, real L, int npl, real xpl[], real ypl[]) {

  real qi=1/q;
  real q11=1./(1.+qi);
  real cnst=qi/L+1./(1.-L)+0.5*(1.+qi)*pow(L-q11,2);

  real lrl = left_limit(q, L);
  xpl[0]= lrl;
  ypl[0]= 0;

  xpl[npl] = L;
  ypl[npl] = 0;

  real rrl = right_limit(q, L);
  xpl[2*npl-1] = rrl;
  ypl[2*npl-1] = 0;

  real y1 = 0;
  real y2 = pow(L, 2);
  real ysq;
  // left lobe
  real dxl = (xpl[npl]-xpl[0])/npl;
  for(int i=1; i<npl; i++) {
    xpl[i] = xpl[0] + i*dxl;
    ysq = rtsafe(rline,qi,L,y1,y2,xpl[i],ROCHE_ACCURACY);
    ypl[i]=sqrt(ysq);
  }

  //right lobe
  real dxr = (xpl[2*npl-1]-xpl[npl])/npl;
  for(int i=1; i<npl-1; i++) {
    xpl[npl+i] = xpl[npl] + i*dxr;
    ysq = rtsafe(rline,qi,L,y1,y2,xpl[npl+i],ROCHE_ACCURACY);
    ypl[npl+i]=sqrt(ysq);
  }
}

void draw_disc(real xc, real yc, real dx, real dy) {

  float x[4],y[4];

  x[0] = xc - dx;
  x[1] = xc + dx;
  x[2] = xc + dx;
  x[3] = xc - dx;

  y[0] = yc - dy;
  y[1] = yc + dy;
  y[2] = yc - dy;
  y[3] = yc + dy;

  cpgpoly(4, x, y);
}

void draw_jet(real xc, real yc, real l, real theta, real dt) {
  
  float x[4],y[4];

  x[0] = xc +  l*cos(theta+dt);
  x[2] = xc -  l*cos(theta-dt);
  x[1] = xc +  l*cos(theta-dt);
  x[3] = xc -  l*cos(theta+dt);

  y[0] = yc +  l*sin(theta+dt);
  y[2] = yc -  l*sin(theta-dt);
  y[1] = yc +  l*sin(theta-dt);
  y[3] = yc -  l*sin(theta+dt);

  cpgpoly(4, x, y);
}

local bool cpgpulsar(real x, real y, stellar_type type) {

  float dx= 0.02;
  float xpos = x - dx+0.005;
  float ypos = y - dx;

  char pulsar[255];
  bool plt_star = true;

  int ocolor;
  cpgqci(&ocolor);
  cpgsci(1);

  // compact objects are black
  if(type==Xray_Pulsar) {
    plt_star = false;
    draw_jet(x, y, 0.15, 60.*PI/180, 2.*PI/180.);
  }
  else  if(type==Radio_Pulsar) {
    plt_star = false;
    sprintf(pulsar, "\\(2270)"); 
    cpgsch(4.);
    cpgptxt(xpos-0.01, ypos-0.02, -30, 0.5, pulsar);
    cpgsch(2.);
    cpgptxt(xpos, ypos, -30., 0.5, pulsar);
    cpgsch(1.);
  }
  else if(type == Neutron_Star) {
    plt_star = false;
    sprintf(pulsar, "\\(2270)"); 
    cpgsch(2.);
    cpgptxt(xpos, ypos, -30., 0.5, pulsar);
    cpgsch(1.);
  }
  else if(type == Black_Hole) {
    plt_star = false;
    cpgsch(2.);
    cpgcirc(xpos+dx, ypos+dx, 4*MIN_RSTAR);
    cpgsci(0);
    cpgcirc(xpos+dx, ypos+dx, MIN_RSTAR);
  }
  else {
   plt_star = true;
  }

  cpgsch(1.);
  cpgsci(ocolor);

  return plt_star;
}


void print_time(real time) {

  char T[] = "T";
  if(time>0) {
    char tm[255];
    sprintf(tm, "T = %7.1lf Myr", time); 
    cpgmtxt(T, 5.0, 0.0, 0.0, tm);
  }
  else {
    cpgmtxt(T, 5.0, 0.0, 0.0, "Zero age");
  }
}

void print_mass(real mass, float pgy, bool primary=true) {

  if(mass>0) {
    char m[255];
    char T[] = "T";
    if(primary)
      sprintf(m, "M = %7.2lf M\\d\\(2281)\\u", mass); 
    else
      sprintf(m, "m = %7.2lf M\\d\\(2281)\\u", mass); 
    cpgmtxt(T, pgy, 0.85, 0.0, m);
  }
}

void initialize_pgplot(roche r, int pgformat,
		       real xmin, real xmax, real ymin, real ymax) {

  char filename[255];
  char ftype[255];
  if(pgformat==1)
    sprintf(ftype, "gif/gif");
  else if(pgformat==0)
    sprintf(ftype, "ps/cps");
    
  int fid = r.get_step();
  if(fid<10)
    sprintf(filename, "SeBa_000%1d.%s", fid, ftype);
  else if(fid<100)
    sprintf(filename, "SeBa_00%2d.%s", fid, ftype);
  else if(fid<1000)
    sprintf(filename, "SeBa_0%3d.%s", fid, ftype);
  else
    sprintf(filename, "SeBa_%4d.%s", fid, ftype);

  if(pgformat<=1)
    if(cpgbeg(0, filename, 1, 1) != 1)
      exit(-1);

  char *T = "T";

  // initialize temperature color scheme
  float red_hls = 120; 
  float h, dh = (360-red_hls)/100.;
  float l = 0.5;
  float s = 1;
  for(int i=0; i<100; i++) {
    h = red_hls + i*dh;
    cpgshls(17+i, h, l, s);
  }
  cpgscir(17, 117);
  float di = 2./100.;
  for(int i=0; i<100; i++) {
    cpgsci(17+i);
  }
  float dx = 0.1*(xmax-xmin)/100.;
  float dy = ymax/100;

  // startup pgplot
  cpgask(0);
  cpgslw(2);
  cpgscf(2);
  cpgsch(1.);
  cpgsci(1);
  cpgenv(xmin, xmax, -ymax, ymax, 1, -2);
  //  cpgenv(xmin, xmax, -ymax, ymax, 1, 1);
  //  cpgvstd();
  // This works identical to cpgenv
  //cpgsvp(0., 1., 0., 1.);
  //cpgwnad(xmin, xmax, -ymax, ymax);
  //cpgswin(xmin, xmax, -0.4*ymax, ymax);
  //cpgbox("BNTS",1., 4,"BNTSV",1.,4);

  // Allow plotting outside viewport
  cpgsclp(0);

  // plot temperature color scheme
  char Te[255];
  sprintf(Te, "log\\d10\\u T\\deff\\u"); 
  cpgptxt(xmax+0.10, -0.5*ymax, -90., 0.5, Te);
  char T3[255];
  sprintf(T3, "10\\u3\\dK"); 
  char T5[255];
  sprintf(T5, "10\\u5\\dK"); 
  cpgtext(xmax+0.07, 0., T5);
  cpgtext(xmax+0.07, -ymax, T3);

  cpgsfs(1);
  for(int i=0; i<100; i++) {
    cpgsci(17+i);
    cpgrect(xmax + 0.02, xmax+0.07, -ymax+i*dy,  -ymax+(i+1)*dy);
  }
  cpgsci(1);

}

void draw_scale_bar(real size, real xmax, real ymax) {

    float a = size;
    bool inRsun = true;
    if(size<1) {
      inRsun = false;
      a *= 6.96e+5; //km
    }
    int nsub;
    real scale = (real)cpgrnd(a, &nsub);
    if(scale/a>1) {
      scale /= 10;
    }
    char rsun[255];
    if(inRsun) {
      sprintf(rsun, "%d R\\d\\(2281)\\u", (int)scale); 
    }
    else {
      sprintf(rsun, "%d km", (int)scale); 
    }

    float xleng[1] = {scale/a};
    float max[1] = {xmax-xleng[0]};
    float min[1] = {xmax};
    float yval[1] = {ymax+0.02};

    cpgerrx(1, min, max, yval, 2.);
    cpgtext(xmax-0.15, ymax-0.03, rsun);
}

void draw_disrupted_binary(roche r, int pgformat) {
  cerr << "Binary Disrupted" << endl;

  char *T = "T";

  real time   = r.get_time();
  real semi   = r.get_semi_major_axis();
  real ecc    = r.get_eccentricity();
  real mprim  = r.get_primary().mass;
  real msec   = r.get_secondary().mass;   
  real rprim  = r.get_primary().radius;
  real rsec   = r.get_secondary().radius;
  real mpcore = r.get_primary().mcore;
  real mscore = r.get_secondary().mcore;

  real dr = 2*(rprim+rsec);
  float xmin = -0.5;
  float xmax =  1.5;
  float ymin = -0.5;
  float ymax = 0.5;

  initialize_pgplot(r, pgformat, xmin, xmax, ymin, ymax);

  float xp = 0;
  float xs = 1;
  if(r.get_primary().mass<r.get_secondary().mass) {
    float xp = 1;
    float xs = 0;
  }
  real rp = rprim/dr;
  real rs = rsec/dr;

  int pcolor = (int)(Starlab::max(1., 
		     Starlab::min(100., 
				  100*(r.get_primary().temperature-3.4)/2.)));
  cpgsci(1);
  cpgsfs(1);
  float dy = ymax/100;
  cpgrect(xmax + 0.02, xmax+0.07, -ymax+pcolor*dy,  -ymax+(pcolor+1)*dy);
  cpgsci(16+pcolor);  // add 16 for skipping predefined colors in pgplot

  if(cpgpulsar(xp, 0., r.get_primary().type)) {
    cpgcirc(xp, 0., rp);
  }
  //  cpgcirc(xp, 0., 0.2);
  cpgsfs(2);

  // set secondary color
  int scolor = (int)(Starlab::max(1., 
                     Starlab::min(100., 
				  100*(r.get_secondary().temperature-3.4)/2.)));
  cpgsci(1);
  cpgsfs(1);
  cpgrect(xmax + 0.02, xmax+0.07, -ymax+scolor*dy,  -ymax+(scolor+1)*dy);
  cpgsci(16+scolor);

  cpgsfs(1);
  if (cpgpulsar(xs, 0., r.get_secondary().type)) {
    cpgcirc(xs, 0., rs);
  }
  cpgsfs(2);

  cpgsfs(4);

  // back to black
  cpgsci(1);

  //  draw core as sphere with radius mc/m and color for stellar core
  cpgsci(get_pgcolor_index(Neutron_Remnant));

  cpgsfs(1);
  if(mpcore>0)
    //    cpgcirc(0., 0., rp * mpcore/mprim);
        cpgcirc(xp, 0., rp * Starlab::min(0.95, mpcore/mprim));
  if(mscore>0)
    //    cpgcirc(1., 0., rs * Starlab::min(0.5, mscore/msec));
    cpgcirc(xs, 0., rs * Starlab::min(0.95, mscore/msec));
  cpgsfs(1);

  // back to black
  cpgsci(1);
  draw_scale_bar(Starlab::max(rprim, rsec), xmax, ymax);

  char state[255];
  sprintf(state, "%s, %s", type_string(r.get_primary().type), 
	                   type_string(r.get_secondary().type));
  cpgmtxt(T, 3.5, 0.0, 0.0, state);

  char bt[255];
  sprintf(bt, "%s", type_string(r.get_binary_type())); 
  cpgmtxt(T, 2.0, 0.0, 0.0, bt);

  print_time(time);

  // print stellar masses
  print_mass(mprim, 3.75);
  print_mass(msec, 2.5, false);

 if(pgformat<=1) {
   cpgend();
 }
 else {
   cpgpage();
 }

}

void draw_single_merger(roche r, int pgformat) {
  cerr << "Binary is merged single star" << endl;

  char *T = "T";

  real time   = r.get_time();
  real mprim  = r.get_primary().mass;
  real rprim  = r.get_primary().radius;
  real mpcore = r.get_primary().mcore;

  real dr = 2*rprim;
  float xmin = -1.0;
  float xmax =  1.0;
  float ymin = -0.5;
  float ymax =  0.5;

  initialize_pgplot(r, pgformat, xmin, xmax, ymin, ymax);

  float xp = 0;
  real rp = rprim/dr;

  int pcolor = (int)(Starlab::max(1., 
		     Starlab::min(100., 
				  100*(r.get_primary().temperature-3.4)/2.)));
  cpgsci(1);
  cpgsfs(1);
  float dy = ymax/100;
  cpgrect(xmax + 0.02, xmax+0.07, -ymax+pcolor*dy,  -ymax+(pcolor+1)*dy);
  cpgsci(16+pcolor);  // add 16 for skipping predefined colors in pgplot

  if(cpgpulsar(xp, 0., r.get_primary().type)) {
    cpgcirc(xp, 0., rp);
  }
  //  cpgcirc(xp, 0., 0.2);
  cpgsfs(2);

  // back to black
  cpgsci(1);

  //  draw core as sphere with radius mc/m and color for stellar core
  cpgsci(get_pgcolor_index(Neutron_Remnant));

  cpgsfs(1);
  if(mpcore>0)
    //    cpgcirc(0., 0., rp * mpcore/mprim);
    cpgcirc(xp, 0., rp * Starlab::min(0.95, mpcore/mprim));
  cpgsci(1);

  // back to black
  draw_scale_bar(rprim, xmax, ymax);

  char state[255];
  sprintf(state, "%s", type_string(r.get_primary().type)); 
  cpgmtxt(T, 3.5, 0.0, 0.0, state);

  char bt[255];
  sprintf(bt, "%s", type_string(r.get_binary_type())); 
  cpgmtxt(T, 2.0, 0.0, 0.0, bt);

  print_time(time);

  // print stellar masses
  print_mass(mprim, 3.75);

 if(pgformat<=1) {
   cpgend();
 }
 else {
   cpgpage();
 }

}

void draw_strong_encounter(roche r, int pgformat) {
  cerr << "Binary exeriences strong encounter" << endl;

}

//void draw_roche(real semi, real mprim, real msec,  
//		real rprim=0, real rsec=0,
//		real mpcore = 0, real mscore = 0) {
void draw_roche(roche r,
		int pgformat) {

  cerr << "Binary, draw Roche lobes" << endl; 

  real time   = r.get_time();
  real semi   = r.get_semi_major_axis();
  real ecc    = r.get_eccentricity();
  real mprim  = r.get_primary().mass;
  real msec   = r.get_secondary().mass;   
  real rprim  = Starlab::max(MIN_RSTAR, r.get_primary().radius);
  real rsec   = Starlab::max(MIN_RSTAR, r.get_secondary().radius);
  real mpcore = r.get_primary().mcore;
  real mscore = r.get_secondary().mcore;

  stellar_type ptype   = (stellar_type)r.get_primary().type;
  stellar_type stype   = (stellar_type)r.get_secondary().type;

  real Rlp = roche_radius(mprim, msec);
  real Rls = roche_radius(msec, mprim);

  real q = msec/mprim;
  real L1 = first_Largrangian_point(1./q);

  // Assure that we print the primary onthe left.
  if(r.get_primary().mass<r.get_secondary().mass) {
    q = 1/q;
    L1 = first_Largrangian_point(1./q);
    }

  PRC(L1);PRC(q);PRC(Rlp);PRL(Rls);

  real rlim = right_limit(q, L1);
  real llim = left_limit(q, L1);

  int npl = 250;
  real xpl[2*npl];
  real ypl[2*npl];

  compute_lobes(q, L1, npl, xpl, ypl);

  //find extrema
  real xmin = xpl[0];
  real xmax = xpl[2*npl-1];

  if(r.get_primary().mass<r.get_secondary().mass) {
    real xm = xmin;
    xmin = 1-xmax;
    xmax = 1-xm;
  }

  real ymax = -1000;
  for(int i=0; i<2*npl; i++) {
    ymax = Starlab::max(ymax, ypl[i]);
  }

  initialize_pgplot(r, pgformat, xmin, xmax, -ymax, ymax);

  char *T = "T";

  //Draw Roche lobes and stars
  float rxpl[2*npl];
  float ryplt[2*npl];
  float ryplb[2*npl];

  // some tideous fiddeling for q>1 case
  if(r.get_primary().mass<r.get_secondary().mass) {
    for(int i=0; i<2*npl; i++) {
      rxpl[i] = 1-xpl[i];
    }
  }
  else {
    for(int i=0; i<2*npl; i++) {
      rxpl[i] = xpl[i];
    }
  }
  for(int i=0; i<2*npl; i++) {
    ryplt[i] = ypl[i];
    ryplb[i] = -ypl[i];
  }
  cpgline(npl+1, &rxpl[0], &ryplt[0]);
  cpgline(npl+1, &rxpl[0], &ryplb[0]);

  real rp = Rlp/rprim;
  real rs = Rls/rsec;

  bool prim_RLOF = false;
  bool sec_RLOF = false;
  if(rprim>=Rlp-ROCHE_ACCURACY) prim_RLOF = true;
  if(rsec>=Rls-ROCHE_ACCURACY)  sec_RLOF = true;

  // set primary color
  int pcolor = (int)(Starlab::max(1., 
		     Starlab::min(100., 
				  100*(r.get_primary().temperature-3.4)/2.)));
  cpgsci(1);
  cpgsfs(1);
  float dy = ymax/100;
  cpgrect(xmax + 0.02, xmax+0.07, -ymax+pcolor*dy,  -ymax+(pcolor+1)*dy);
  cpgsci(16+pcolor);  // add 16 for skipping predefined colors in pgplot

  if(!prim_RLOF) {
    cpgsfs(1);

    // underfilling Roche lobe
    if(cpgpulsar(0., 0., r.get_primary().type)) {
	cpgcirc(0., 0., rprim);
    }
    cpgsfs(2);
  }
  else {
    // primary fills its Roche-lobe see if secondary has disc
    real rcs = (1+q)*pow(0.5 - 0.227 * log10(q), 4);
    if(rsec<rcs) {                 // draw disc around secondary
      cpgsfs(1);
      draw_disc(1., 0., 0.7*Rls, 0.1*Rls);
    }

    // fill Roche lobe
    cpgsfs(1);
    if(r.get_primary().mass>r.get_secondary().mass) {
      cpgpoly(npl+1, &rxpl[0], &ryplt[0]);
      cpgpoly(npl+1, &rxpl[0], &ryplb[0]);
    }
    else {
      cpgpoly(npl, &rxpl[npl], &ryplt[npl]);
      cpgpoly(npl, &rxpl[npl], &ryplb[npl]);
    }
  }

  // back to black
  cpgsci(1);

  cpgline(npl, &rxpl[npl], &ryplt[npl]);
  cpgline(npl, &rxpl[npl], &ryplb[npl]);

  // set secondary color
  int scolor = (int)(Starlab::max(1., 
                     Starlab::min(100., 
				  100*(r.get_secondary().temperature-3.4)/2.)));
  cpgsci(1);
  cpgsfs(1);
  cpgrect(xmax + 0.02, xmax+0.07, -ymax+scolor*dy,  -ymax+(scolor+1)*dy);
  cpgsci(16+scolor);

  if(!sec_RLOF) {
    cpgsfs(1);
    if (cpgpulsar(1., 0., r.get_secondary().type)) {
      cpgcirc(1., 0., rsec);
    }
    cpgsfs(2);
  }
  else {

    // secondary fills roche-lobe, see if primary has a disc
    real rcp = (1+q)*pow(0.5 - 0.227 * log10(q), 4);
    if(rprim<rcp) {  // draw disc around secondary
      cpgsfs(1);
      draw_disc(0., 0., 0.7*Rlp, 0.1*Rlp);
    }

    cpgsfs(1);
    // Fill secondary Roche lobe
    if(r.get_primary().mass>r.get_secondary().mass) {
      cpgpoly(npl, &rxpl[npl], &ryplt[npl]);
      cpgpoly(npl, &rxpl[npl], &ryplb[npl]);
    }
    else {
      cpgpoly(npl+1, &rxpl[0], &ryplt[0]);
      cpgpoly(npl+1, &rxpl[0], &ryplb[0]);
    }
  }

  cpgsfs(4);

  // back to black
  cpgsci(1);

  //  draw core as sphere with radius mc/m and color for stellar core
  cpgsci(get_pgcolor_index(Neutron_Remnant));

  cpgsfs(1);
  if(mpcore>0)
    cpgcirc(0., 0., Starlab::min(Rlp,rprim) * mpcore/mprim);
  if(mscore>0)
    cpgcirc(1., 0., Starlab::min(Rls, rsec) * mscore/msec);
  cpgsfs(1);

  // back to black
  cpgsci(1);

  // add scale bar
  if(semi>0) {
    draw_scale_bar(semi, xmax, ymax);
    // add orbital period
    real P = semi_to_period(semi, mprim, msec);
    char Pd[255];
    sprintf(Pd, "P = %7.2lf days", P); 
    cpgmtxt(T, 5.0, 0.85, 0.0, Pd);

  } 
  else {
    float min[1] = {L1};
    float yval[1] = {-0.7*ymax};
    char *rsun = "No size scaling";
    cpgtext(xmax-0.15, ymax, rsun);
  }

  // add eccentricity information, if e>0
  char e[255];
  if(ecc>0) 
    sprintf(e, "e = %7.2lf", ecc); 
  else
    sprintf(e, "Circular"); 
  cpgmtxt(T, 1.25, 0.85, 0.0, e);

  // print stellar masses
  print_mass(mprim, 3.75);
  print_mass(msec, 2.5, false);

  // print time if t>0
  print_time(time);

  // initialize binary type and mass transfer type, but only print the 
  // latter if there is RLOF.
  char btt[255];
  char bt[255];
  sprintf(bt, "%s", type_string(r.get_binary_type())); 
  sprintf(btt, "(%s)", type_string(r.get_mass_transfer_type()));
  cpgmtxt(T, 2.0, 0.0, 0.0, bt);

  // illustrate binary configuration in SPZ notation
  char state[255];
  if(prim_RLOF) {
    if(sec_RLOF) {
      sprintf(state, "[%s, %s]", type_string(r.get_primary().type), 
	      type_string(r.get_secondary().type));
      cpgmtxt(T, 2.0, 0.2, 0.0, btt);
    }
    else {
      sprintf(state, "[%s, %s)", type_string(r.get_primary().type), 
	      type_string(r.get_secondary().type));
      cpgmtxt(T, 2.0, 0.2, 0.0, btt);
    }
  }
  else {
    if(sec_RLOF) {
      sprintf(state, "(%s, %s]", type_string(r.get_primary().type), 
	      type_string(r.get_secondary().type));
      cpgmtxt(T, 2.0, 0.2, 0.0, btt);
    }
    else {
      sprintf(state, "(%s, %s)", type_string(r.get_primary().type), 
	      type_string(r.get_secondary().type));
    }
  }
  cpgmtxt(T, 3.5, 0.0, 0.0, state);

  // plot center-of-mass dotted line
  real com = msec/(mprim + msec);
  //  if(q>1) com = q-1;
  float xj[2] = {com, com};
  float yj[2] = {-0.75*ymax, 0.75*ymax};

  cpgsls(4);
  cpgline(2, xj, yj);
  //  cpgptxt(com, -ymax, 0., 0.5, "com");
  cpgsls(1);

 if(pgformat<=1) {
   cpgend();
 }
 else {
   cpgpage();
 }

}


//------------------------------------------------------------------
//  main  -- roche
//           draws roche lobe for an evolving binary
//      usage:
//              roche
//
//      options: many
//
//----------------------------------------------------------------
//
int main(int argc, char ** argv)   {

  bool verbose = false;
  int c;

  int format = 1; // series of gif files
  
  input_type inpttype = SeBa;
  roche *r = new roche();

  char *comment;
  extern char *poptarg;
  int  pgetopt(int, char **, char *);

 while ((c = pgetopt(argc, argv, "a:C:c:e:f:I:M:m:P:p:R:r:S:s:t:v")) != -1)

    switch(c)
      {
      case 'a': r->set_semi_major_axis(atof(poptarg));
	        break;
      case 'e': r->set_eccentricity(atof(poptarg));
	        break;
      case 'f': format = atoi(poptarg);
	        break;
      case 'C': r->set_primary_mcore(atof(poptarg));
	        break;
      case 'c': r->set_secondary_mcore(atof(poptarg));
	        break;
      case 'I': inpttype = (input_type)atoi(poptarg);
	        break;
      case 'M': r->set_primary_mass(atof(poptarg));
	        break;
      case 'm': r->set_secondary_mass(atof(poptarg));
	        break;
      case 'P': 
      case 'p': r->set_primary_type(atoi(poptarg));
	        break;
      case 'R': r->set_primary_radius(atof(poptarg));
	        break;
      case 'r': r->set_secondary_radius(atof(poptarg));
	        break;
      case 'S': 
      case 's': r->set_secondary_type(atoi(poptarg));
	        break;
      case 't': r->set_time(atof(poptarg));
	        break;
      case 'v': verbose = true;
	        break;
      case '?': cerr << "\nusage: roche << infile";
	exit(1);
      }

// format ==0 or 1 are handeled in draw_roche
 if(format==2) {
   cpgbeg(0, "roche.ps/cps", 1, 1);
 }
 else  if(format==3) {
   cpgbeg(0, "roche.gif/gif", 1, 1);
 }
 else if(format==4) {
   cpgbeg(0, "/xw", 1, 1);
 }
 
 int n=0;
 r->set_input_type(inpttype);
 if(r->check_roche_input()) {
   r->set_step(n++);
   draw_roche(*r, format);
 }
 else {
   delete r;
   while((r = read_roche(cin, inpttype)) != NULL) {

     r->set_input_type(inpttype);
     if(r->check_roche_input()) {
       r->set_step(n++);
       cerr << "Draw roche lobes or otherwise" << endl;

       if(r->get_binary_type()==Disrupted) {
	 draw_disrupted_binary(*r, format);
       } 
       //       else if(r->get_semi_major_axis()<=0) {
       else if(r->get_binary_type()==Merged) {
	 draw_single_merger(*r, format);
       }
       else if(r->get_binary_type() == Strong_Encounter) {
	 draw_strong_encounter(*r, format);
       }
       else {
	 draw_roche(*r, format);
       }
       cerr << "Summary: " << endl;
       cerr << *r; 
     } else {
       cerr << "Did not pass if statement" << endl;
       return 1;
     }
     delete r;
   };
	 //   while(!cin.eof());   // eof() does not work propery with cin??
 }

 if(format>1) {
   cpgend();
 }
}




