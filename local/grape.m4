dnl this M4 file defines the GRAPE software interface
dnl please, don't change the names of the variables
dnl just replace the second part of the call with your values

dnl where to find GRAPE libraries
dnl m4_define(GRAPE_LDFLAGS_, [-L/usr/local/grape/lib])
m4_define(GRAPE_LDFLAGS_, [-L/usr2/makino/disk3src/harplibs/x86_64])

dnl list a GRAPE libraries that must be linked in
dnl m4_define(GRAPE_LIBS_, [-lgl6 -lgrape -lblah])
m4_define(GRAPE_LIBS_, [-lg6lx -lg6lxsim2])
