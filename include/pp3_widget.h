#ifndef QTPP3
#define QTPP3

#include <qdialog.h>

class QLabel;
class QLineEdit;
class QSlider;
class QPushButton;
class QTextEdit;
class QTimer;
class hdyn;

class pp3_widget : public QDialog {

    Q_OBJECT

    public:
        pp3_widget(int config = 0, QWidget *parent = 0, const char *name = 0);
	void start_timer();
 	void initialize_slider(hdyn *b);
 	void set_node(hdyn *b);

    private slots:
	void set_node_from_string(const QString &line);
	void set_node_from_line();
	void set_node_from_index(int);
	void quit_pp3();
	void timeout_quit_pp3();
	void update_time();

    private:

	hdyn *node;
	bool slider_initialized;
	int timeout;
 	int seconds;

        QLabel *label;
	QLineEdit *idline;
	QSlider *slider;
	QLabel *tlabel;
	QPushButton *next;
	QPushButton *quit;
	QTextEdit *text;
	QTimer *timer;
	QTimer *counter;
};

#endif
