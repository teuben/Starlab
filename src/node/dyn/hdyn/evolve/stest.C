
#include "hdyn.h"

main()
{
    hdyn *b = get_hdyn(cin);
    if (b == NULL) err_exit("Can't read input snapshot");

    while(1) {
	int i = 0;

	for_all_daughters(hdyn, b, bi)
	    putiq(bi->get_dyn_story(), "nb_check_counter", i++);

	for_all_daughters(hdyn, b, bi) {
	    if (find_qmatch(bi->get_dyn_story(), "nb_check_counter")) {
		i = getiq(bi->get_dyn_story(), "nb_check_counter");
		rmq(bi->get_dyn_story(), "nb_check_counter");
	    }
	}
    }
}
