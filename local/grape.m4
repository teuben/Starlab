dnl This M4 file defines the GRAPE software interface.
dnl Please don't change the names of the variables, just
dnl replace the second part of the call with your values.

dnl ***** FOR EXAMPLES, SEE grape.m4.examples IN THIS DIRECTORY *****

dnl Where to find GRAPE libraries:

m4_define(GRAPE_LDFLAGS_, [-L/usr2/makino/disk3src/harplibs/x86_64])

dnl GRAPE libraries that must be linked in:

m4_define(GRAPE_LIBS_,    [-lg6lx -lg6lxsim2])
