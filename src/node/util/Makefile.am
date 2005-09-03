# 
#        #=======================================================#     _\|/_
#       #  __  _____           ___                    ___       #       /|\
#      #  /      |      ^     |   \  |         ^     |   \     #           _\|/_
#     #   \__    |     / \    |___/  |        / \    |___/    #             /|\
#    #       \   |    /___\   |  \   |       /___\   |   \   #  _\|/_
#   #     ___/   |   /     \  |   \  |____  /     \  |___/  #    /|\
#  #                                                       #             _\|/_
# #=======================================================#               /|\
#
##
##	src/node/util/Makefile.am
##..............................................................................
##	Dec 1994		Piet Hut & Steve McMillan
##	Jun 2004		Ernest Mamikonyan
##..............................................................................
##

noinst_LTLIBRARIES = libnode-1.la
bin_PROGRAMS = node_io display_subtrees rmq snapprune make_single_node \
    add_daughter_node pretty_print_tree print_normal makenode makemass \
    makesecondary mass_dist m_bar renumber setmass std_rename_nodes add_index

libnode_1_la_SOURCES = node_io.C hydrobase_io.C starbase_io.C tree_io.C \
    node_tt.C rmtree.C print_normal.C starbase_init.C makenode.C makemass.C \
    renumber.C

node_io_CPPFLAGS = -DTOOLBOX
display_subtrees_CPPFLAGS = -DTOOLBOX
rmq_CPPFLAGS = -DTOOLBOX
snapprune_CPPFLAGS = -DTOOLBOX
make_single_node_CPPFLAGS = -DTOOLBOX
add_daughter_node_CPPFLAGS = -DTOOLBOX
pretty_print_tree_CPPFLAGS = -DTOOLBOX
print_normal_CPPFLAGS = -DTOOLBOX
makenode_CPPFLAGS = -DTOOLBOX
makemass_CPPFLAGS = -DTOOLBOX
makesecondary_CPPFLAGS = -DTOOLBOX
mass_dist_CPPFLAGS = -DTOOLBOX
m_bar_CPPFLAGS = -DTOOLBOX
renumber_CPPFLAGS = -DTOOLBOX
setmass_CPPFLAGS = -DTOOLBOX
std_rename_nodes_CPPFLAGS = -DTOOLBOX

node_io_LDADD = ../libnode.la
display_subtrees_LDADD = ../libnode.la
rmq_LDADD = ../libnode.la
snapprune_LDADD = ../libnode.la
make_single_node_LDADD = ../libnode.la
add_daughter_node_LDADD = ../libnode.la
pretty_print_tree_LDADD = ../libnode.la
print_normal_LDADD = ../libnode.la
makenode_LDADD = ../libnode.la
makemass_LDADD = ../libnode.la
makesecondary_LDADD = ../libnode.la
mass_dist_LDADD = ../libnode.la
m_bar_LDADD = ../libnode.la
renumber_LDADD = ../libnode.la
setmass_LDADD = ../libnode.la
std_rename_nodes_LDADD = ../libnode.la

node_io_SOURCES = node_io.C
display_subtrees_SOURCES = display_subtrees.C
rmq_SOURCES = rmq.C
snapprune_SOURCES = snapprune.C
make_single_node_SOURCES = make_single_node.C
add_daughter_node_SOURCES = add_daughter_node.C
pretty_print_tree_SOURCES = pretty_print_tree.C
print_normal_SOURCES = print_normal.C
makenode_SOURCES = makenode.C
makemass_SOURCES = makemass.C
makesecondary_SOURCES = makesecondary.C
mass_dist_SOURCES = mass_dist.C
m_bar_SOURCES = m_bar.C
renumber_SOURCES = renumber.C
setmass_SOURCES = setmass.C
std_rename_nodes_SOURCES = std_rename_nodes.C
