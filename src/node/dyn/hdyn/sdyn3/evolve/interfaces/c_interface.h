
void c_cpu_init();
real c_cpu_time();

int  c_srandinter(int, int);
real c_randinter(real, real);
void c_print_initial_random_parameters();

int  c_pgetopt(int, char**, char*);

void c_make_standard_init(initial_state3*);
void c_initialize_angles(initial_state3*, int, bool, real);

void c_print_scatter3_info(initial_state3*,
			   intermediate_state3*,
			   final_state3*,
			   bool, bool, bool,
			   real);
void c_print_intermediate_descriptor(enum intermediate_descriptor3);
void c_print_final_descriptor(enum final_descriptor3);

void  c_scatter3(initial_state3*,
		 intermediate_state3*,
		 final_state3*,
		 real, real, real, real);
