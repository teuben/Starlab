
c       Declarations of interface routines (no trailing underscores!)

	external f_cpu_init		!$pragma C (f_cpu_init)
	external f_cpu_time		!$pragma C (f_cpu_time)

	external f_srandinter		!$pragma C (f_srandinter)
	external f_randinter		!$pragma C (f_randinter)

	external f_print_initial_random_parameters
				!$ pragma C (f_print_initial_random_parameters)

	external f_make_standard_init	!$pragma C (f_make_standard_init)
	external f_initialize_angles	!$pragma C (f_initialize_angles)

	external f_print_scatter3_info	!$pragma C (f_print_scatter3_info)
	external f_print_intermediate_descriptor
	1			!$pragma C (f_print_intermediate_descriptor)
	external f_print_final_descriptor
	1			!$pragma C (f_print_final_descriptor)

c       Note that some brain-dead FORTRANs demand that the $pragma
c       descriptor be on the SAME line as the external declaration.

	external f_scatter3		!$pragma C (f_scatter3)
	
