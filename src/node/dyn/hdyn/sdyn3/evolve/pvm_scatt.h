#include "sigma3.h"

// PVM-related stuff:

// Message IDs:

#define SEND_DATA_MSG	1
#define RETURN_DATA_MSG	2
#define HANDSHAKE_MSG	9
#define TERMINATE_MSG	99
#define KILL_CONFIRM	998
#define FAILURE_MSG	999

void initialize_processors(int nproc, int debug);
void terminate_processors(int debug);

void pack(initial_state3 &init);
void pack(intermediate_state3 &inter);
void pack(final_state3 &final);
void unpack(initial_state3 &init);
void unpack(intermediate_state3 &inter);
void unpack(final_state3 &final);
