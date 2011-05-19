
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/// @file grape6.h  C/C++ interface definitions for using GRAPE-6.
//
//  version 1:  July 2000	Steve McMillan, Jun Makino

#ifndef  STARLAB_GRAPE6_H
#  define  STARLAB_GRAPE6_H

extern "C" void g6_set_tunit_(int *new_tunit);
extern "C" void g6_set_xunit_(int *new_xunit);

extern "C" int  g6_open_(int *cluster_id);
extern "C" int  g6_close_(int *cluster_id);


// **************** For use with sapporo16_steve ***************

extern "C" int g6_open_special(int ngpu, int *list);
extern "C" int g6_open_special2(int ngpu, int *list, int nmax);

// *************************************************************


extern "C" int  g6_reset_(int *cluster_id);
extern "C" int  g6_reset_fofpga_(int *cluster_id);

extern "C" int  g6_npipes_();

extern "C" void g6_set_ti_(int *cluster_id, double *ti);
extern "C" int  g6_set_j_particle_(int *cluster_id,
				   int *address,
				   int *index,
				   double *tj,
				   double *dtj,
				   double *mass,
				   double k18[3],	// k/18
				   double j6[3],	// j/6
				   double a2[3],	// a/2
				   double v[3],
				   double x[3]);

// *** Compiler doesn't like passing vec* as real** in the following
// *** functions.  The GRAPE-4 fix was to declare the arguments as vec*.
// *** Not tested...

extern "C" void g6calc_firsthalf_(int *cluster_id,
				  int *nj,
				  int *ni,
				  int index[],
				  vec xi[],		// double xi[][3],
				  vec vi[],		// double vi[][3],
				  vec aold[],	// double acc[][3],
				  vec j6old[],	// double jerk[][3],
				  double phiold[],
				  double *eps2,
				  double h2[]);

extern "C" int g6calc_lasthalf_(int *cluster_id,
				int *nj,
				int *ni,
				int index[],
				vec xi[],		// double xi[][3],
				vec vi[],		// double vi[][3],
				double *eps2,
				double h2[],
				vec acc[],		// double acc[][3],
				vec jerk[],		// double jerk[][3],
				double pot[]);

extern "C" int g6calc_lasthalf2_(int *cluster_id,
				 int *nj,
				 int *ni,
				 int index[],
				 vec xi[],		// double xi[][3],
				 vec vi[],		// double vi[][3],
				 double *eps2,
				 double h2[],
				 vec acc[],		// double acc[][3],
				 vec jerk[],		// double jerk[][3],
				 double pot[],
				 int nnbindex[]);

extern "C" int g6_read_neighbour_list_(int *cluster_id);
extern "C" int g6_read_neighbour_list_old_(int *cluster_id);

extern "C" int g6_get_neighbour_list_(int *cluster_id,
				      int *pipe,
				      int *max_length,
				      int *nblen,
				      int nbl[]);

extern "C" int g6_initialize_jp_buffer_(int *clusterid, int *size);
extern "C" int g6_flush_jp_buffer_(int *clusterid);

#endif

//=======================================================================//
//  +---------------+        _\|/_        +------------------------------\\ ~
//  |  the end of:  |         /|\         |  inc/grape6.h
//  +---------------+                     +------------------------------//
//========================= STARLAB =====================================\\ ~

