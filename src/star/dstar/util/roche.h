
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/*
 *  roche.h  --  basic structures for drawing Roche-lobes
 *           
 *.....................................................................
 *    version 1:  Feb 2003   Simon F. Portegies Zwart
 *    version 2:
 *...................................................................
 *     This file includes:
 *  1) definition of structure loca_star
 *  2) definition of class roche
 *
 *....................................................................
 */
#ifndef  _ROCHE_
#  define  _ROCHE_

#include "double_support.h"
#include "main_sequence.h"
#include "single_star.h"

#define ROCHE_ACCURACY 1.e-4

enum input_type {SeBa = 0, BSE=1};

int get_pgcolor_index(stellar_type_summary stp);
int get_pgcolor_index(stellar_type stype);

struct local_star {
  stellar_type type;
  int id;
  real temperature;
  real mass;
  real radius;
  real mcore;

  local_star() {
    type = NAS;
    temperature = 0;
    id = 0;
    mass = radius = mcore = 0;
  }
};

/*-----------------------------------------------------------------------------
 *  roche  --  basic structures for drawing Roche-lobes
 *-----------------------------------------------------------------------------
 */
class roche
    {
    protected:

      int step;
      binary_type btype;
      mass_transfer_type mttype;
      input_type input_code;
      real time;
      real semi;
      real ecc;
      local_star primary;
      local_star secondary;
      
    private:
    public:

      roche() {
	input_code = SeBa;
	btype = Detached;
	mttype = Unknown;
	step = 1;
	semi = -1;
	ecc = 0;
	primary.mass = -1;
	primary.radius = 0;
	primary.mcore = 0;
	secondary.mass = -1;
	secondary.radius = 0;
	secondary.mcore = 0;
      }
      roche(real mp, real ms) {
	btype = Detached;
	mttype = Unknown;
	step = 1;
	input_code = SeBa;
	time = 0;
	roche(-1.0, 0.0, mp, ms);
	primary.radius = 0;
	primary.mcore = 0;
	secondary.radius = 0;
	secondary.mcore = 0;
      }
      roche(real a, real e, real mp, real ms) {

	roche(mp, ms);
	semi = a;
	ecc = e;
      }
      roche(real a, real e, real mp, real ms, real rp, real rs, 
	    real mcp, real mcs) {
	roche(a, e, mp, ms);
	primary.radius = rp;
	primary.mcore = mcp;
	secondary.radius = rp;
	secondary.mcore = mcp;
      }

      ~roche() {}

      int  get_step() {return step;}
      binary_type  get_binary_type() {return btype;}
      mass_transfer_type get_mass_transfer_type() {return mttype;}
      real get_time() {return time;}
      real get_semi_major_axis() {return semi;}
      real get_eccentricity() {return ecc;}

      local_star get_secondary(){return  secondary;}
      local_star get_primary(){return  primary;}

      void set_step(int s) {step = s;}
      void set_binary_type(binary_type b) {btype = b;}
      void set_mass_transfer_type(mass_transfer_type m) {mttype = m;}
      void set_time(real x) {time = x;}
      void set_semi_major_axis(real x) {semi = x;}
      void set_eccentricity(real x) {ecc = x;}
      void set_primary_id(int i) {primary.id = i;}
      void set_primary_type(int i) {primary.type = (stellar_type)i;}
      void set_primary_mass(real x) {primary.mass = x;}
      void set_primary_radius(real x){primary.radius = x;}
      void set_primary_mcore(real x){primary.mcore = x;}
      void set_secondary_id(int i) {secondary.id = i;}
      void set_secondary_type(int i) {secondary.type = (stellar_type)i;}
      void set_secondary_mass(real x) {secondary.mass = x;}
      void set_secondary_radius(real x) {secondary.radius = x;}
      void set_secondary_mcore(real x){secondary.mcore = x;}

      void set_primary_temperature(real t) { primary.temperature = t;}
      void set_secondary_temperature(real t) {secondary.temperature = t;}

      void set_input_type(input_type i) {input_code = i;}

      bool check_roche_input() {

	if(primary.mcore>primary.mass) primary.mcore = primary.mass;
	if(secondary.mcore>secondary.mass) secondary.mcore = secondary.mass;

	if(semi>0) {
	  primary.radius /= semi;
	  secondary.radius /= semi;
	}

	//	return true;

	//#if 0
	if(time<0 || time>1000000) {
	  return false;
	}

	PRL(ecc);
	if(ecc>1) {
	  return false;
	}

	if(primary.mass<=0 && secondary.mass<=0) {
	  return false;
	}

	if(semi<0) 
	  return false;
	  
	return true;
	//#endif
      }

      friend ostream& operator<<(ostream& s, roche& r) {
	s << r.get_time() <<" "<< r.get_semi_major_axis() 
	  <<" "<< r.get_eccentricity()
	  <<"\t "<< r.get_primary().type 
	  <<" "<<  r.get_primary().mass 
	  <<" "<< r.get_primary().radius <<" "<< r.get_primary().mcore 
	  <<"\t "<< r.get_secondary().type 
	  <<" "<< r.get_secondary().mass 
	  <<" "<< r.get_secondary().radius <<" "<< r.get_secondary().mcore 
	  << endl;

	return s;
      }

    };

roche *read_roche(istream& s, input_type input_code = SeBa) {

  roche *r = new roche();
  if(s.eof() || r == NULL) 
    return NULL;

  int id, pid, sid, ptype, stype;
  int btpe, mttpe;
  real time, semi, ecc, pmass, pradius, pmcore, smass, sradius, smcore;
  real ptemp, stemp;

  if(input_code == SeBa) {
    cerr << "\nInput from SeBa" << endl;
    // Typical SeBa input
    //0 0 138 0.4    0 3 13 4.36074 0.01    1 3 9 3.53452 0.01    
    //0 14.3679 140.447 0.4    0 5 12.6182 12.183 3.04833    1 3 8.99844 4.76709 0.01    

    s >> id >> btpe >> mttpe >> time >> semi >> ecc 
      >> pid >> ptype >> pmass >> pradius >> ptemp >> pmcore 
      >> sid >> stype >> smass >> sradius >> stemp >> smcore; 
    //    if(btpe<Strong_Encounter) {
    //    }

  }
  else if(input_code == BSE) {
    cerr << "Input from BSE" << endl;
    // Typical BSE output
    //     TIME    M1     M2  KW1 KW2    SEP    ECC  R1/ROL1 R2/ROL2  TYPE
    //    0.0000  1.936  0.361  1  0   130.022  0.00  0.023   0.011  INITIAL
    // 1279.5332  1.936  0.361  2  0   130.022  0.00  0.051   0.011  KW CHNGE
    // 1290.1372  1.936  0.361  3  0   130.030  0.00  0.080   0.011  KW CHNGE
    // 1321.0033  1.931  0.361  4  0    93.276  0.00  0.204   0.015  KW CHNGE
    // 1515.2379  1.910  0.361  5  0    94.103  0.00  0.455   0.015  KW CHNGE
    // 1517.2539  1.908  0.361  5  0    70.398  0.00  1.001   0.020  BEG RCHE
    // 1517.2539  0.535  0.361  8  0     4.951  0.00  1.001   0.020  COMENV
    // 1517.2539  0.535  0.361  8  0     4.951  0.00  0.080   0.202  END RCHE
    // 1520.2010  0.532  0.361 11  0     4.969  0.00  0.007   0.201  KW CHNGE
    //10181.4150  0.532  0.361 11  0     1.002  0.00  0.033   1.001  BEG RCHE
    //15000.0000  0.532  0.056 11  0     0.705  0.00  0.034   1.098  MAX TIME

    char state[8];

    id = pid = sid = -1;
    pmcore = smcore = 0;

    cerr << "Reading input" << endl;
    s >> time >> pmass >> smass >> ptype >> stype >> semi >> ecc 
      >> pradius >> sradius;
    //      >> pradius >> sradius >> state;

    cerr << "\nNew binary entry." << endl;
    
  PRC(time); PRC(pmass);  PRC(smass);  PRC(ptype);  PRC(stype);  PRC(semi);  PRC(ecc); 
      PRC(pradius);  PRC(sradius); 

    ptype = convert_BSE_to_SeBa_stellar_type(ptype);
    stype = convert_BSE_to_SeBa_stellar_type(stype);
    pradius *= roche_radius(semi, pmass, smass);
    sradius *= roche_radius(semi, smass, pmass);
  }

  PRC(time);PRC(semi);PRC(ecc); 
  PRC(btpe);PRL(mttpe);
  PRI(4);PRC(pid);PRC(ptype);PRC(pmass);PRC(pradius);PRC(ptemp);PRL(pmcore);
  PRI(4);PRC(sid);PRC(stype);PRC(smass);PRC(sradius);PRC(stemp);PRL(smcore);

  r->set_time(time);
  r->set_binary_type((binary_type)btpe);
  r->set_mass_transfer_type((mass_transfer_type)mttpe);
  r->set_semi_major_axis(semi);
  r->set_eccentricity(ecc);
  r->set_primary_id(pid);
  r->set_primary_type(ptype);
  r->set_primary_mass(pmass);
  r->set_primary_radius(pradius);
  r->set_primary_temperature(ptemp);
  r->set_primary_mcore(pmcore);
  r->set_secondary_id(sid);
  r->set_secondary_type(stype);
  r->set_secondary_mass(smass);
  r->set_secondary_radius(sradius);
  r->set_secondary_temperature(stemp);
  r->set_secondary_mcore(smcore);
    
  return r;
}
      
	
#endif // _ROCHE_
