dnl this M4 file defines the GRAPE software interface
dnl please, don't change the names of the variables
dnl just replace the second part of the call with your values

dnl where to find GRAPE include files
m4_define(GRAPE_CPPFLAGS_, [-I/usr/local/grape/include])

dnl where to find GRAPE libraries
m4_define(GRAPE_LDFLAGS_, [-L/usr/local/grape/lib])

dnl list a GRAPE libraries that must be linked in
m4_define(GRAPE_LIBS_, [-lgl6 -lgrape -lblah])
