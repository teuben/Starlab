
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// worldlines: returns n worldlines from worldbundles 
////
////    options:
////          -B starting  worldbundle [0]
////          -b ending worldbundle [all]
////          -d output timestep [1/64]
////          -t output worldline at time t [not specified]
////
////     input: worldbundle output file from kira
////
////    output: worldlines of selected worldbundles 
////
////    Problem: The output worldline can presently not be used
///              as an input for kira...
////
////    SPZ at MIT in Dec 2000
////
//
//----------------------------------------------------------------------

#include "worldline.h"

#ifdef TOOLBOX

local void create_and_print_interpolated_tree(ostream & s, 
					      worldbundle *wb,
					      real t_world, 
					      bool vel = true,
					      bool brief = false) {

  dyn* b = create_interpolated_tree(wb, t_world, vel);
  if (b) {
    put_dyn(s, *b, brief);
    rmtree(b);
  }
}

local void convert_relative_to_absolute(dyn* b)
{
    if (b->get_parent()) b->inc_pos(b->get_parent()->get_pos());
    for_all_daughters(dyn, b, bb) convert_relative_to_absolute(bb);
}

main(int argc, char *argv[]) {

  bool c_flag = false;

  real dt = 0.015625;	// powers of 2 are preferred, but not essential  
  int b_start = 0;      // starting bundle
  int b_end   = 32767;  // all files in the worldbundle
  bool t_flag = false;
  real t_world;         // output worldline at time = t


  extern char *poptarg;
  char *comment;
  int c;
  char* param_string = "B:b:d:f:t:c:";
  
  check_help();
  
  while ((c = pgetopt(argc, argv, param_string)) != -1)
    switch(c) {
      
	case 'B': b_start = atoi(poptarg);
		break;
	case 'b': b_end = atoi(poptarg);
		break;
	case 'd': dt = atof(poptarg);
		break;
	case 't': t_flag = true;
		  t_world = atof(poptarg);
		break;
	case 'c': c_flag = TRUE;
		comment = poptarg;
		break;
	case '?': params_to_usage(cerr, argv[0], param_string);
		get_help();
		exit(1);
	      }

  if(t_flag) {
    b_start = 0;     
    b_end   = 32767; 
  }
  else if(b_end<=b_start)
    err_exit("begin worldbunlde <= end worldbundle");

  
  real time, t_min, t_max;
  dyn *b;
  worldbundle *wb;
  int itot=0, iw=0;
  if(b_start>0) {
    for(int i=0; i<b_start; i++)  {
      
      cerr << "Skip worldbundle " << i << endl;
      if(!(wb = read_bundle(cin)))
	err_exit("Not that many worldbundles");
    }
  }

  int ib = b_start;
  while (ib<b_end  && (wb = read_bundle(cin))) {

    time = wb->get_t_min();
    t_max = wb->get_t_max();

    ib++;
    cerr << "Time= " << time << " N_bundle= " << ib << endl; 

    if(t_flag) {

      if(t_world>=time && t_world<=t_max) {
	cerr << "Time= " << t_world
	  << " N (b, w, tot)= " << ib << " " << iw <<" " << itot << endl; 

	create_and_print_interpolated_tree(cout, wb, t_world);

	delete wb;
	exit(1);
      }
    }
    else {

      iw = 0;
      do {
	
	itot++;
	iw++;
	create_and_print_interpolated_tree(cout, wb, time);
	cerr << "Time= " << time
	  << " N (b, w, tot)= " << ib << " " << iw <<" " << itot << endl; 
#if 0
	create_and_print_interpolated_tree(cout, wb, t_world);
	b = create_interpolated_tree(wb, time, true);
	//    b = create_interpolated_tree(wb, time);
	if (b) {
	  //      convert_relative_to_absolute(b);
	  cerr << "Time= " << time 
	    << " N (b, w, tot)= " << ib << " " << iw <<" " << itot << endl; 
	  put_dyn(cout, *b, false);
	  rmtree(b);
	}
#endif	
	time += dt;
      }
      while(time<=t_max);
    }

    delete wb;
  };
}

#endif
