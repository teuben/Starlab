
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//  grape4.h: C interface definition for GRAPE-4 (HARP-3)
//.............................................................................
//    version 1:  Feb 1995   Piet Hut, Steve McMillan, Jun Makino
//    version 2:  August 1996 Jun Makino
//.............................................................................
//     This file includes:
//  1) definition of C interface to GRAPE-4
//.............................................................................

#ifndef  STARLAB_HARP3_H
#  define  STARLAB_HARP3_H

extern "C" void h3open_();
extern "C" void h3close_();
extern "C" int h3npipe_();
extern "C" void h3setnboards_(int * );
extern "C" int h3getnboards_();
extern "C" unsigned int h3wait_();
extern "C" void h3setti_(double * ti);
extern "C" void h3setmode_(int * mode, int * iboard);
extern "C" void h3setled_(int * mode, int * iboard);
extern "C" int h3jpmax_();
extern "C" void h3jpdma_indirect_(int * nj,int * hostindex, vec * xj, 			     vec * vj,vec * aj,vec * jj,real * mj,real * tj,int *mode);
extern "C" void h3mjpdma_indirect_(int * nj,int * hostindex, vec * xj, 			      vec * vj,vec * aj,vec * jj,              			      real * mj,real * tj,int *mode, int * buff_id);
extern "C" void h3mjpdma_start_(int * buff_id);
extern "C" void h3mjpdma_flush_();
extern "C" void h3calc_(int * nj,int * ni,vec * xi,vec * vi,real * eps2,real* h2,
		   vec * acc,vec * jerk,real * pot);
extern "C" void h3calc_firsthalf_(int * nj,int *ni,vec * xi,vec * vi,real *eps2,real *h2);
extern "C" void h3calc_lasthalf_(int *ni,vec * acc,vec * jerk,real *pot);
extern "C" void h3nbread_(int *nboards);
extern "C" int h3nblist_(int *board, int *chip, int *nblist);
extern "C" int count_nblist_low(int board, int chip);
extern "C" void h3setdebuglevel_(int * debug_level);
extern "C" void h3checkdmabuffers_();
extern "C" void set_time_check_mode(int mode);

#endif

//=======================================================================//
//  +---------------+        _\|/_        +------------------------------\\ ~
//  |  the end of:  |         /|\         |  inc/grape4.h
//  +---------------+                     +------------------------------//
//========================= STARLAB =====================================\\ ~

