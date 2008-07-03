//  imf_evolve.C: SPZ:May 1998
//                Semi-analytic evolution of the mass function 
//                of a collisiosn dominated system.
//                The assumptions are rather symplistic:

#include "dyn.h"
#include "star/stdfunc.h"

#define Rsun_pc 2.255e-8                 // R_sun/parsec = 6.960e+10/3.086e+18

#ifdef TOOLBOX

local real encounter_rate(int im, real ni, real mi, real ri, real vi, 
			  real rci,
			  int jm, real nj, real mj, real rj, real vj, 
			  real rcj) 
{
	  
    real to_volume = 4*PI/3.;

    //PRC(im);PRC(ni);PRC(mi);PRC(ri);PRC(vi);PRL(rci);
    //PRC(jm);PRC(nj);PRC(mj);PRC(rj);PRC(vj);PRL(rcj);
  
    real n_rhoi = 1./(to_volume * pow(rci, 3));
    real n_rhoj = 1./(to_volume * pow(rcj, 3));
    real r_tot  = ri + rj;
    real m_tot  = mi + mj;
    real v_rel  = sqrt(pow(vi, 2) + pow(vj, 2));
    real rcore = Starlab::min(rci, rcj);

    real geometric = 6.31e-15 * v_rel * pow(r_tot, 2);
    real grav_foc = 3.61e-9 * m_tot * r_tot / v_rel;

    real rate = 0;
	
    if (im==jm) {

      if (ni>1)
	rate = 0.5 * ((ni-1) * n_rhoi)
                   *  (nj * n_rhoj)
                   * pow(rcore, 3)
                   * (grav_foc + geometric);
    }
    else
      rate = (ni * n_rhoi)
           * (nj * n_rhoj)
           * pow(rcore, 3)
	   * (grav_foc + geometric);
    
#if 0
    if(rate>0)
      cerr << jm <<" "<< im <<" "
	   << mj <<" "<< mi <<" "
	   << nj <<" "<< ni <<" "
	   << rate << endl;
#endif

    return rate;
}



//#else

main(int argc, char **argv)
{
    int nzones = 10000;
    int nstep  = 1;
    real mmin = 0.1;
    real mmax = 100;
    real x_imf= 2.35;
    real rcore= 0.1; //pc
    real vsun = 10;  //km/s
    real dt = 1;     //Myr
    int n = 4;

    extern char *poptarg;
    int c;
    const char *param_string = "N:n:M:m:T:x:r:v:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c) {

	    case 'n': nzones = atoi(poptarg);
		      break;
	    case 'N': nstep = atoi(poptarg);
		      break;
	    case 'M': mmax = atof(poptarg);
		      break;
	    case 'm': mmin = atof(poptarg);
		      break;
	    case 'T': dt = atof(poptarg);
		      break;
	    case 'x': x_imf = atof(poptarg);
		      break;
	    case 'r': rcore = atof(poptarg);
		      break;
	    case 'v': vsun = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    // Loop over input until no more data remain.

    real *imf  = new real[n*nzones+1];
    real *cmf  = new real[n*nzones+1];
    real *fmf  = new real[n*nzones+1];
    real *mass = new real[n*nzones+1];
    real *radi = new real[n*nzones+1];
    real *vdsp = new real[n*nzones+1];
    real *loss = new real[n*nzones+1];
    real *gain = new real[n*nzones+1];

    real dmass  = (mmax-mmin)/(1.*nzones);

    for(int im=1; im<n*nzones; im++) {
      mass[im] = mmin + dmass * im;
      //      radi[im] = 1./pow(mass[im], 2.7);
      radi[im] = 1./pow(mass[im], 1.7);
      vdsp[im] = vsun/sqrt(mass[im]);
      imf[im]  = 0;
      cmf[im]  = 0;
      fmf[im]  = 0;
      loss[im] = 0;
      gain[im] = 0;
    }

    for(int im=1; im<nzones; im++)
      imf[im] = 1./pow(mass[im], x_imf);
    for(int im=1; im<nzones; im++)
      imf[im] = imf[im]/imf[nzones-1];

    //real imtot=0;
    for(int im=1; im<nzones; im++) {
      //imtot += imf[im] * mass[im];
      cmf[im]  = imf[im];
    }

    real rate, total_rate = 0;
    real ncoll, nloss=0;
    int km;
    real total_loss, total_gain, nstar_loss;
    total_loss=total_gain=nstar_loss=0;
    
    for(int is=0; is<nstep; is++) {
      for (int im=1; im<n*nzones; im++) {
	for (int jm=im; jm<n*nzones; jm++)
	  if(cmf[im]>0 && cmf[jm]>0) {

	  rate = encounter_rate(im,  cmf[im], mass[im], radi[im], 
				vdsp[im], rcore, 
				jm, cmf[jm], mass[jm], radi[jm],
				vdsp[jm], rcore);
	  total_rate += rate;

	  ncoll = Starlab::max(0., Starlab::min(rate * dt,
			  Starlab::min(cmf[im]-loss[im],
				  cmf[jm]-loss[jm])));

	  km = Starlab::min(n*nzones, im + jm); // Conserve mass
	  loss[im] += ncoll;
	  loss[jm] += ncoll;
	  gain[km] += ncoll;

	  //PRC(im);PRC(jm);PRL(km);
	  //PRC(loss[im]);PRC(loss[jm]);PRL(gain[km]);
	}
      }

      cerr << "     Total encounter rate = "
	   << total_rate << " [per Myr]" << endl;

      for (int im=1; im<n*nzones; im++) {
	nloss += gain[im] - loss[im];
	total_loss += loss[im] * mass[im];
	total_gain += gain[im] * mass[im];
	fmf[im] = cmf[im] - loss[im] + gain[im];
      }
      for (int im=1; im<=n*nzones; im++) {
	loss[im]=gain[im]=0;
	cmf[im] = fmf[im];
      }
    }

    // Normalizing to imf N=1;
    real imf_tot = 0;
    for (int im=1; im<n*nzones; im++) 
      imf_tot += imf[im];
    real fmtot=0;
    real imtot=0;
    for (int im=1; im<n*nzones; im++) {
      imf[im] = imf[im]/imf_tot;
      fmf[im] = fmf[im]/imf_tot;
      imtot += imf[im]*mass[im];
      fmtot += fmf[im]*mass[im];
    }

    cerr << "Old and new mass function: "<<endl;
    for (km=1; km<n*nzones; km++) {
      if(imf[km]>0 || fmf[km]>0)
	cerr << km <<" "<< mass[km] <<" "
	       << imf[km] << " " << fmf[km] << endl;
    }

    cerr << " N star lost = " << nloss<< endl;
    cerr << "Total mass goes from: "<< imtot << " to " << fmtot <<endl;

}

#endif
