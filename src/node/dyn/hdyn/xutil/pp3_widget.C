
#include "pp3_widget.h"
#include "hdyn.h"

#include <qapplication.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qslider.h>
#include <qpushbutton.h>
#include <qtextedit.h>
#include <qlayout.h>
#include <qfont.h>
#include <qtimer.h>

#include <sstream>

#define FONT_SIZE	18			// unit = pixels
#define FONT2_SIZE	15

#define INIT_TIMEOUT	10
#define STD_TIMEOUT	120

pp3_widget::pp3_widget(int config,		// default = 0
		       QWidget *parent,		// default = NULL
		       const char *name)	// default = NULL
    : QDialog(parent, name)
{
    // Configuration choices:
    //
    //		0	==>  Toolbox option:  "Next" and "Quit" boxes.
    //		1	==>  Function option: "Done" box only.
    //
    // "Done" and "Quit" have the same action.

    node = NULL;
    slider_initialized = false;
    timeout = 1;	// progression: 1 (just in case) --> INIT --> STD
    seconds = 0;

    // Font:

    QFont font, font2;
    font.setPixelSize(FONT_SIZE);
    font2.setPixelSize(FONT2_SIZE);
    qApp->setFont(font, true);
    
    // Widgets:

    label = new QLabel("ID ", this);
    idline = new QLineEdit(this);
    slider = new QSlider(Qt::Horizontal, this);
    tlabel = new QLabel(this);
    if (config == 0) {
	next = new QPushButton("&Next snapshot", this);
	quit = new QPushButton("&Quit", this);
    } else
	quit = new QPushButton("&Done", this);
    text = new QTextEdit(this);
    timer = new QTimer(this);
    counter = new QTimer(this);

    setCaption("Qt_pp3");

    // Connections:

    // connect(idline, SIGNAL(textChanged(const QString &)),
    //	       this, SLOT(set_node_from_string(const QString &)));

    connect(idline, SIGNAL(returnPressed()),
	    this, SLOT(set_node_from_line()));
    connect(slider, SIGNAL(valueChanged(int)),
	    this, SLOT(set_node_from_index(int)));
    if (config == 0)
	connect(next, SIGNAL(clicked()), qApp, SLOT(quit()));
    connect(quit, SIGNAL(clicked()), this, SLOT(quit_pp3()));
    connect(timer, SIGNAL(timeout()), this, SLOT(timeout_quit_pp3()));
    connect(counter, SIGNAL(timeout()), this, SLOT(update_time()));

    // Attributes:

    idline->setFixedWidth(6*FONT_SIZE);		// width ~ height/2
    tlabel->setFont(font2);
    tlabel->setAlignment(AlignHCenter | AlignVCenter | ExpandTabs);
    if (config == 0)
	next->setAutoDefault(false);
    quit->setAutoDefault(false);
    text->setReadOnly(true);
    text->setMinimumSize(45*FONT_SIZE, 15*FONT_SIZE);

    start_timer();

    // Layout:

    QHBoxLayout *topleft = new QHBoxLayout;
    topleft->addWidget(label);
    topleft->addWidget(idline);

    QVBoxLayout *left = new QVBoxLayout;
    left->addLayout(topleft);
    left->addSpacing(10);
    left->addWidget(slider);
    left->addSpacing(20);
    left->addWidget(tlabel);
    left->addStretch(1);
    if (config == 0) left->addWidget(next);
    left->addWidget(quit);

    QHBoxLayout *main = new QHBoxLayout(this);
    main->addLayout(left);
    main->addWidget(text);
    main->setStretchFactor(left, 0);
    main->setStretchFactor(text, 1);
}

void pp3_widget::update_time()
{
    char tmp[64];
    seconds++;
    sprintf(tmp, "( timeout = %d s )", timeout-seconds);
    tlabel->setText(tmp);
}

void pp3_widget::start_timer()
{
    seconds = -1;
    update_time();
    timer->start(1000*timeout, true);
    counter->start(1000);		// 1 second increments
}

void pp3_widget::initialize_slider(hdyn *b)
{
    int i = b->get_index();
    if (i > 0) slider->setValue(i);	// won't generate valueChanged
					// if node not yet set
    b = b->get_root();
    int imin = 1000000000, imax = -1;
    for_all_leaves(hdyn, b, bb) {
	int index = bb->get_index();
	if (index < imin) imin = index;
	if (index > imax) imax = index;
    }

    slider->setMinValue(imin);
    slider->setMaxValue(imax);
    slider->setTickmarks(QSlider::Both);

    if (i <= 0) slider->setValue(imin);
    slider_initialized = true;
}

void pp3_widget::set_node(hdyn *b)			    // set all from b
{
    if (b && b->is_valid()) {

	if (b->is_root()) b = b->get_oldest_daughter();	    // avoid full dump

	if (!slider_initialized) initialize_slider(b);

	if (node != b) {

	    // Node pointer:

	    node = b;

	    // ID line (no additional signal generated):

	    idline->setText(b->format_label());

	    // Slider value (note that this will generate a sliderChanged
	    // signal):

	    if (b->get_index() >= 0) slider->setValue(b->get_index());

	    // Text window:

	    ostringstream s;
	    pp3(b, s);
	    text->setText(s.str());
	}

	// Reset the timer:

	if (timeout < INIT_TIMEOUT) timeout = INIT_TIMEOUT;
	start_timer();
    }
}

void pp3_widget::set_node_from_string(const QString &line)
{
    if (node) {
	char *name = new char[line.length()+1];
	strcpy(name, line.latin1());
	hdyn *new_node = (hdyn*)node_with_name(name, node->get_root());

	if (new_node) {			// "Enter" will reset the timer
	    timeout = STD_TIMEOUT;	// even if no node change occurs

	    if (new_node != node)
		set_node(new_node);
	    else
		start_timer();
	}
    }
}

void pp3_widget::set_node_from_line()
{
    set_node_from_string(idline->text());
}

void pp3_widget::set_node_from_index(int i)
{
    if (node) {
	hdyn *new_node = (hdyn*)node_with_index(i, node->get_root());
	if (new_node && new_node != node) {
	    timeout = STD_TIMEOUT;
	    set_node(new_node);
	}
    }
}

void pp3_widget::quit_pp3()
{
    qApp->exit(1);
}

void pp3_widget::timeout_quit_pp3()
{
    qApp->exit(2);
}

