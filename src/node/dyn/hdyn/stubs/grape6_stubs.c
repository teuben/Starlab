
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~


// grape6_stubs.C: GRAPE-6 stubs, based on the real GRAPE-6 library.
//		   Strategy: Void functions have null bodies;
//		   int functions return 0.
//
// Steve McMillan, June 2004


void g6_set_tunit_(int *new_tunit) {};
void g6_set_xunit_(int *new_xunit) {};

int  g6_open_(int *cluster_id) {return 0;};
int  g6_close_(int *cluster_id) {return 0;};

int  g6_reset_(int *cluster_id) {return 0;};
int  g6_reset_fofpga_(int *cluster_id) {return 0;};

int  g6_npipes_() {return 0;};

void g6_set_ti_(int *cluster_id, double *ti) {};
int  g6_set_j_particle_(int *cluster_id,
			int *address,
			int *index,
			double *tj,
			double *dtj,
			double *mass,
			double k18[3],	// k/18
			double j6[3],	// j/6
			double a2[3],	// a/2
			double v[3],
			double x[3]) {return 0;};

void g6calc_firsthalf_(int *cluster_id,
		       int *nj,
		       int *ni,
		       int index[],
		       double xi[][3],
		       double vi[][3],
		       double acc[][3],
		       double jerk[][3],
		       double phiold[],
		       double *eps2,
		       double h2[]) {};

int g6calc_lasthalf_(int *cluster_id,
		     int *nj,
		     int *ni,
		     int index[],
		     double xi[][3],
		     double vi[][3],
		     double *eps2,
		     double h2[],
		     double acc[][3],
		     double jerk[][3],
		     double pot[]) {return 0;};

int g6calc_lasthalf2_(int *cluster_id,
		      int *nj,
		      int *ni,
		      int index[],
		      double xi[][3],
		      double vi[][3],
		      double *eps2,
		      double h2[],
		      double acc[][3],
		      double jerk[][3],
		      double pot[],
		      int nnbindex[]) {return 0;};

int g6_read_neighbour_list_(int *cluster_id) {return 0;};

int g6_get_neighbour_list_(int *cluster_id,
			   int *pipe,
			   int *max_length,
			   int *nblen,
			   int nbl[]) {return 0;};
