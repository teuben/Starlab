
//// scatter_hist: keeps track of scatter experiments
////          
//// Options:     none
//-----------------------------------------------------------------------------
//   version 1:  Aug 2000   Simon Portegies Zwart   spz@mit.edu
//                                                  mit
//.............................................................................
//   non-local functions: 
//-----------------------------------------------------------------------------

#include "stdinc.h"
#include "sdyn.h"
#include "scatter_exp.h"
//
// scatter_hist.h
//

#ifndef    _SCATTER_HIST
#  define  _SCATTER_HIST

/*-----------------------------------------------------------------------------
 * scatter_hist  --  a linked list of scatter histories.
 *-----------------------------------------------------------------------------
 */
class scatter_hist : public scatter_exp {
    protected:

  //    int *stellar_identities;      //  {1,2,3,4,5}
  //    int n_relations;              //  4
  //    int *stellar_relations;       //  {1,2, 1,3, 2,3 4,5}

    scatter_hist * past;
    scatter_hist * future;

    public:
       scatter_hist(istream &s) : scatter_exp(s) {

	 past=NULL;
	 future=NULL;
       }

       scatter_hist(scatter_exp &hi) : scatter_exp(hi) {

	 past=NULL;
	 future=NULL;
       }

       scatter_hist() : scatter_exp() {

	 //	 for(int i=0; i<N_RHO_ZONE_MAX; i++) 
	 //	   n_hits[i] = 0;

	 //	   if (s) {
	 //  	      past=s;
	 //	      past->future=this;
	 //	      future = NULL;
	 //	   }
	 //	   else {
	      past=NULL;
	      future=NULL;
	      //	  }
       }

       ~scatter_hist(){
	 
         if (future!=NULL) {
	     scatter_hist *tmp = future;
	     future = NULL;
	     delete tmp;
	 }

	 if (past)
	    past->future = NULL;
       }
  
       scatter_hist* get_past() {return past;}
       void set_past(scatter_hist *sb) {past = sb;}
       scatter_hist* get_future() {return future;}
       scatter_hist* get_first() {
           if (past!=NULL)
              return past->get_first();
           else
              return this;
       }
       scatter_hist* get_last() {
           if (future!=NULL)
              return future->get_last();
           else 
              return this;
       }

       void set_future(scatter_hist* f) {future = f;}
       void set_last(scatter_hist* f) {
            if(future!=NULL) 
              future->set_last(f);
            else
              future=f;
       }

       scatter_hist* get_scatter_hist(sdyn* b=NULL);
       void add_scatter_hist(istream &s);
       void add_scatter_hist(scatter_exp he, int zone);
       scatter_hist* scatter_hist::contains(scatter_hist*);

  scatter_hist* get_identical_scatter_hist();
  void inc_id_scenario();

  void put_scatter_hist(ostream&, bool verbose = false);

  //  bool read_scatter_hist(istream& s);

  void put_state(ostream&);

};

#define for_all_scatter_hist(scatter_hist, base, scatter_hist_next)   \
        for (scatter_hist* scatter_hist_next = base;                  \
	     scatter_hist_next != NULL;                               \
	     scatter_hist_next = scatter_hist_next->get_future())

scatter_hist* get_history(scatter_hist *hi, istream& is);

//bool scenarios_identical(scatter_hist* hi, scatter_hist* ha);
void put_state(scatter_hist * hi, ostream & s);
scatter_hist* initialize_scatter_hist(sdyn* b);

#endif
