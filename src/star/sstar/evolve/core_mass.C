#include "node.h"
#include "single_star.h"
#include "main_sequence.h"


void dump_core_mass(single_star* str, char * filename) {

  ofstream s(filename, ios::app|ios::out);
  if (!s) cerr << "error: couldn't create file "<<filename<<endl;
  
  s << str->get_current_time() 
    <<" "<< str->get_core_mass()
    <<" "<< str->temperature()
    <<" "<< str->get_luminosity()
    <<" "<< str->get_radius()
    <<" "<< str->get_element_type()
    <<endl;      
  
  s.close();
}


#ifdef TOOLBOX

void main(int argc, char ** argv) {
    int id=1,n=4;
    int  c;
    real mf=1,rf=1, tf=1;
    real m=1.;
    stellar_type type = Main_Sequence;
    real t_start = 0;
    real t_end   = 10000;
    real time = 0;

    extern char *poptarg;

    while ((c = pgetopt(argc, argv, "m:n:",
		    "$Revision$", _SRC_)) != -1)
	switch(c)
	    {
            case 'm': m = atof(poptarg);
                      break;
            case 'n': n = atoi(poptarg);
                      break;
	    }

    mf = m;
    node *b  = mknode(1);
    b->get_starbase()->set_stellar_evolution_scaling(mf, rf, tf);
    addstar(b, t_start, type);

    single_star* str = (single_star*)b->get_oldest_daughter()->get_starbase();
    t_end = str->nucleair_evolution_time();
    cerr<<"t_end"<<t_end<<endl;
    do {
      int i=0;
      str = (single_star*)b->get_oldest_daughter()->get_starbase();
      real t_next= str->get_next_update_age();
//      cerr << "t_next,time"<<t_next<<" "<<time<<endl;
//      real dt=(t_next-time)/n;
      real dt= 0.1 * str->get_evolve_timestep();
      cerr << "dt = " << dt << endl;
      if (dt<=0.) cerr<<"dt <0 !!!!!!!!!!!!!!!!!!!!!!!"<<endl;
//      cerr<<"dt"<<dt<<endl;
//      for(i=0;i<n;i++) {
	time+=dt;
//	cerr<<"time"<<time<<endl;
	b->get_oldest_daughter()->get_starbase()->evolve_element(time);
	str = (single_star*)b->get_oldest_daughter()->get_starbase();
	dump_core_mass(str, "m_core.dat");
//      }
//      dt=.001;
//      time+=dt;
//      b->get_oldest_daughter()->get_starbase()->evolve_element(time);
//      str = (single_star*)b->get_oldest_daughter()->get_starbase();
//      dump_core_mass(str, "m_core.dat");
    } while (time<=t_end);
}

#endif
