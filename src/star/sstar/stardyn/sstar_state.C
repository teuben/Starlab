#include "sstar_to_dyn.h"
#include "star_state.h"
#include "single_star.h"
#include "dyn.h"

#ifndef TOOLBOX

//Declaration can be cound in: sstar_to_dyn.h
star_state make_star_state(dyn* bi) {

  star_state str;
  if (bi->get_oldest_daughter() != NULL)
    return str;

  starbase *sb = bi->get_starbase();
  if (bi->get_use_sstar() && has_sstar(bi)) 
    str.make_star_state((star*) sb);
  else if (sb->get_star_story()!=NULL) {
    
    real t_curr, t_rel, m_rel, m_env, m_core, co_core, T_eff, L_eff, R_eff;
    real p_psr, b_fld;
    stellar_type type=NAS;
    story *st = sb->get_star_story();
    extract_story_chapter(type, t_curr, t_rel, m_rel, m_env, m_core, co_core,
			  T_eff, L_eff, p_psr, b_fld, *st);
    R_eff = sqrt(max(0.0, (1130.*L_eff)/pow(T_eff, 4)));
    str.lclass    = get_luminosity_class(T_eff, L_eff);
    str.class_tpe = get_spectral_class(T_eff);
    str.type      = type;
    str.mass      = m_env + m_core;
    str.radius    = R_eff;
    //str.magnetic_field = b_fld;
    //str.rotation_period = p_psr;
  }
  else {
    bi->pretty_print_node(cerr);
    return str;
  }
  return str;
}

#else


main() {
  cerr <<"Hallo"<<endl;
}


#endif
