
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Qt_pp3:  Qt-enabled snapshot browsing tool.  Standalone tool applies
////          pp3 to specific nodes, under user control.  The compiled
////          version operates from within kira, activated by a "QT" file.
////
//// Options:    -o    specify test option [0]
////                       0: normal use
////                       1: test QApplication
////                       2: NULL QApplication
////                       3: test DISPLAY
////                       4: NULL DISPLAY
////
//// Steve McMillan, 6/04

#include "hdyn.h"

#ifdef HAVE_LIBQT
#  include "pp3_widget.h"
#  include <qapplication.h>
#endif

#ifndef TOOLBOX

#ifndef HAVE_LIBQT

int Qt_pp3(hdyn *b,
	   const char *disp)
{
    cerr << "Qt_pp3: Qt support unavailable." << endl;
}

#else

// Note that QApplication will fail and TERMINATE the calling program
// if the application window cannot be opened.  We don't want such a
// failure to terminate a kira run!  Thus, we start cautiously by
// determining the X display to use, and cautiously attepmting to open
// it.  If we fail, we exit the function.  If we succeed, we pass the
// Display pointer to QApplication.  Note the overloaded definition.

int Qt_pp3(hdyn *b,
	   QApplication *app)
{
    // Return status:
    //
    //		0	"Next" pressed or Qt window closed.
    //		1	"Quit" pressed.
    //
    // It is up to the calling function to decide what to do with
    // this information.

    bool del = false;
    pp3_widget *Qtpp3;

    if (app)

	Qtpp3 = (pp3_widget*)app->mainWidget();

    else {

	// Shouldn't happen... Set up a fake argv.

	int argc = 1;

	char argv1[argc][1024];
	strcpy(argv1[0], "Qt_pp3");

	char *argv[argc];
	argv[0] = argv1[0];

	app = new QApplication(argc, argv);
	Qtpp3 = new pp3_widget(1);
	app->setMainWidget(Qtpp3);
	del = true;
    }

    Qtpp3->set_node(b);
    Qtpp3->show();
    int ret = app->exec();
    if (ret == 2) cerr << "(timed out)" << endl;

    if (del) {
	delete Qtpp3;
	delete app;
    }
    return ret;
}

#include <X11/Xlib.h>

int Qt_pp3(hdyn *b,
	   const char *disp)		// default = NULL ==> use $DISPLAY
{
    // Choose an X display name.

    char display[1024];
    if (disp)
	strncpy(display, disp, 1023);
    else
	strncpy(display, getenv("DISPLAY"), 1023);
    display[1023] = '\0';

    // QApplication will terminate the program if we can't open the
    // X display.  Avoid this by attempting to open the display before
    // launching the pp3_widget,

    Display *Xdisp = XOpenDisplay(display);
    if (!Xdisp) {
	cerr << "Can't open X display " << display << endl;
	return 1;
    } else
	cerr << "Using X display " << display << endl;

    // However, using this X display directly doesn't seem to work.
    // Remove it and make up a fake command-line string to let Qt do it.

    XCloseDisplay(Xdisp);
    int argc = 3;

    char argv1[argc][1024];
    strcpy(argv1[0], "Qt_pp3");
    strcpy(argv1[1], "-display");
    strcpy(argv1[2], display);

    char *argv[argc];
    for (int k = 0; k < argc; k++)
	argv[k] = argv1[k];

    QApplication *app = new QApplication(argc, argv);
    pp3_widget *Qtpp3 = new pp3_widget(1);
    app->setMainWidget(Qtpp3);

    Qtpp3->set_node(b);
    Qtpp3->show();
    int ret = app->exec();
    if (ret == 2) cerr << "(timed out)" << endl;

    delete Qtpp3;			// necessary if timers are running?
    delete app;

    return ret;
}

#endif

#else

int main(int argc, char **argv)
{
    check_help();

#ifdef HAVE_LIBQT

    extern char *poptarg;
    int c;
    char* param_string = "o:";

    int option = 0;

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'o': option = atoi(poptarg);
		      break;

            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	}

    hdyn *b = NULL;
    QApplication *app = NULL;

    if (option == 0) {

	// Normal use for the toolbox version is to create the application
	// here, rather than creating multiple applications if Qt_pp3() is
	// called multiple times.  To see the alternative, comment out 
	// the next three lines.

	app = new QApplication(argc, argv);
	pp3_widget *Qtpp3 = new pp3_widget(0);
	app->setMainWidget(Qtpp3);

	while ((b = get_hdyn())) {
	    if (Qt_pp3(b, app) == 1) exit(0);
	    delete b;
	}

    } else {

	// For testing the various argument options.

	if (b = get_hdyn()) {
	    if (option == 1) {
		app = new QApplication(argc, argv);
		pp3_widget *Qtpp3 = new pp3_widget(1);
		app->setMainWidget(Qtpp3);
		cerr << option << ": sending existing QApplication" << endl;
		Qt_pp3(b, app);
	    } else if (option == 2) {
		cerr << option << ": sending NULL QApplication" << endl;
		Qt_pp3(b, app);
	    } else if (option == 3) {
		cerr << option << ": sending X display name" << endl;
		Qt_pp3(b, ":0.0");
	    } else if (option == 4) {
		cerr << option << ": sending NULL X display" << endl;
		Qt_pp3(b, (char*)NULL);
	    }
	}
    }

#else

    cerr << "Qt support unavailable." << endl;

#endif
}

#endif
