
/* Like c_interface.h, but everything is explicitly extern "C". */

extern "C" void c_cpu_init();
extern "C" real c_cpu_time();

extern "C" int  c_srandinter(int, int);
extern "C" real c_randinter(real, real);
extern "C" void c_print_initial_random_parameters();

extern "C" int  c_pgetopt(int, char**, char*);

extern "C" void c_make_standard_init(initial_state3*);
extern "C" void c_initialize_angles(initial_state3*, int, bool, real);

extern "C" void c_print_scatter3_info(initial_state3*,
				      intermediate_state3*,
				      final_state3*,
				      bool, bool, bool,
				      real);
extern "C" void c_print_intermediate_descriptor(enum intermediate_descriptor3);
extern "C" void c_print_final_descriptor(enum final_descriptor3);

extern "C" void  c_scatter3(initial_state3*,
			    intermediate_state3*,
			    final_state3*,
			    real, real, real, real);
