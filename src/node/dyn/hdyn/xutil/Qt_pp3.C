
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

	// Shouldn't happen...
	// (We need a better way to set up a fake argv...)

	char argv1[2][1024];
	strcpy(argv1[0], "Qt_pp3");
	strcpy(argv1[1], "");

	char *argv[2];
	argv[0] = argv1[0];
	argv[1] = argv1[1];

	int argc = 1;

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

    // Attempt to open the display.

    Display *Xdisp = XOpenDisplay(display);
    if (!Xdisp) {
	cerr << "Can't open X display " << display << endl;
	return 1;
    } else
	cerr << "Using X display " << display << endl;

    // Display is apparently usable.

    QApplication *app = new QApplication(Xdisp);
    pp3_widget *Qtpp3 = new pp3_widget(1);
    app->setMainWidget(Qtpp3);

    Qtpp3->set_node(b);
    Qtpp3->show();
    int ret = app->exec();
    if (ret == 2) cerr << "(timed out)" << endl;

    delete Qtpp3;			// necessary if timers are running?
    delete app;
    XCloseDisplay(Xdisp);		// Qt won't do this...
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
    char* param_string = "c:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	}

    QApplication *app = NULL;

    // Normal use for the toolbox version is to create the application
    // here, rather than creating multiple applications if Qt_pp3() is
    // called multiple times.  To see the alternative, comment out 
    // the next three lines.

    app = new QApplication(argc, argv);
    pp3_widget *Qtpp3 = new pp3_widget(0);
    app->setMainWidget(Qtpp3);

    hdyn *b;

    while ((b = get_hdyn())) {
	if (Qt_pp3(b, app) == 1) exit(0);
	delete b;
    }

#else

    cerr << "Qt support unavailable." << endl;

#endif
}

#endif
