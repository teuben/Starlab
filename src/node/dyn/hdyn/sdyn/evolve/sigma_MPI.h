
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/*
 *  sigma_MPI.h: definitions for MPI implementation
 *.............................................................................
 *    version 1:  November 2000   Simon Portegies Zwart
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  1) definition of state structures
 *       ....
 *.............................................................................
 */

#ifndef  STARLAB_SIGMA_MPI_H
#  define  STARLAB_SIGMA_MPI_H

#include "sigma.h"
#ifdef USE_MPI
#include "mpi++.h"
#else
#include "localmpi++.h"
#endif

void slave_process(int master, 
		   sigma_input &input, MPI_Datatype inputtype,
		   scatter_exp &experiment, MPI_Datatype scatter_exp_type);
void execute_sigma_experiment(sigma_input &input);
void terminate_all_processors();
int master_process(sigma_out & out,
		   sigma_input &input, MPI_Datatype inputtype,
		   scatter_exp &experiment, MPI_Datatype scatter_exp_type);

#endif



