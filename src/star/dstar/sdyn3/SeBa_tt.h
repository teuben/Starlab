#include "stdinc.h"
#include "scatter3.h"
#include "dstar_to_dyn.h"

void merge_binary_components(dyn* bi);
void ionize_binary(dyn* target);
bool repair_tree(dyn * root);
void update_binary_parameters(dyn *bi, 
			      final_state3 &final3);
void preserve_initial_state(dyn* bi, dyn *si,
			    final_state3& final3);
void exchange_secondary(dyn *bi,
			dyn *si,
			final_state3& final3);
void exchange_primary(dyn *bi,
		      dyn *si,
                      final_state3& final3);
void merge_two_stars(dyn *pi, dyn *si);
void merger_binary_primary(dyn* bi,
			   dyn *si,
			   final_state3& final3);
void merger_binary_secondary(dyn* bi,
			     dyn *si,
			     final_state3& final3);
void merger_binary_bullet(dyn* bi,
			  dyn* si,
                          final_state3& final3);
void merger_escape_primary(dyn *bi, dyn *si, 
			   final_state3 &final3);
void merger_escape_secondary(dyn *bi, dyn *si, 
			     final_state3 &final3);
void merger_escape_bullet(dyn *bi, dyn *si, 
			  final_state3 &final3);
void merger_triple(dyn *bi, dyn *si, 
		   final_state3 &final3);
void stopped_scatter3();
void unknown_final_scatter3();
void error_scatter3();
final_descriptor3 form_scatter_product(dyn * bi, 
				       dyn * si,
				       initial_state3& init3,
				       intermediate_state3& inter3,
				       final_state3 & final3); 

