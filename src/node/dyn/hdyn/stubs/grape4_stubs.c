
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~


// grape4_stubs.C: GRAPE-4 stubs, based on the real GRAPE-4 library.
//		   Strategy: Void functions have null bodies;
//		   int functions return 0.
//
// Steve McMillan, June 2004


void h3open_() {};

void h3close_() {};

int h3npipe_() {return 0;};

void h3setnboards_(int *i) {};

int h3getnboards_() {return 0;};

unsigned int h3wait_() {return 0;};

void h3setti_(double *ti) {};

void h3setmode_(int *mode, int *iboard) {};

void h3setled_(int *mode, int *iboard) {};

int h3jpmax_() {return 0;};

void h3jpdma_indirect_(int *nj, int *hostindex, double **xj, 
		       double **vj, double **aj, double **jj, double *mj,
		       double *tj, int *mode) {};

void h3mjpdma_indirect_(int *nj, int *hostindex, double **xj, 
			double **vj, double **aj, double **jj, 
			double *mj, double *tj, int *mode,
			int *buff_id) {};

void h3mjpdma_start_(int *buff_id) {};

void h3mjpdma_flush_() {};

void h3calc_(int *nj, int *ni, double **xi, double **vi,
	     double *eps2, double* h2, 
	     double **acc, double **jerk, double *pot) {};

void h3calc_firsthalf_(int *nj, int *ni, double **xi, 
		       double **vi, double *eps2, double *h2) {};

void h3calc_lasthalf_(int *ni, double **acc, double **jerk,
		      double *pot) {};

void h3nbread_(int *nboards) {};

int h3nblist_(int *board, int *chip, int *nblist) {return 0;};

int count_nblist_low(int board, int chip) {return 0;};

void h3setdebuglevel_(int *debug_level) {};

void h3checkdmabuffers_() {};

void set_time_check_mode(int mode) {};
