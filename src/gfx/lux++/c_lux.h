
/* lux.h: Definitions needed for the lux interface. */

#include "c_stdinc.h"

typedef unsigned long Window;
typedef unsigned long Color;

#define GRAPH_WINDOW    1
#define POPUP_WINDOW    2
#define DIALOG_WINDOW   3
#define BUTTON_WINDOW   4
#define INPUT_WINDOW    5
#define TEXT_WINDOW     6    

#define NO_TYPE         0
#define OK_BUTTON       1
#define CANCEL_BUTTON   2
#define CHECK_BUTTON    3
#define OK_KEEP_BUTTON  4

enum   {ok=1, ok_keep, cancel};
