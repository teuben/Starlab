#include "node.h"

void main(int argc, char ** argv)
{
    node* b;
    b = get_node(cin);

    if (find_qmatch(b->get_dyn_story(), "test"))
	cerr << "hello!\n";

    putrq(b->get_dyn_story(), "test", 42.0);

    put_node(cout, *b);
}
