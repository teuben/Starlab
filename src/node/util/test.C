#include "node.h"

void main(int argc, char ** argv)
{
    node *b;
    b = get_node();

    putrq(b->get_dyn_story(), "test", 42.0);

    put_node(b);
}
