
#ifdef HAS_PVM

#include "scatter3.h"
#include "pvm3.h"

// Pack and unpack scattering structures (overloaded!).

// This version just sends the data as byte arrays.  It remains to
// be seen what happens if the data alignment required by the system
// leaves gaps in the structures...

void pack(initial_state3 &init)
{
    pvm_pkbyte((char*)&init, sizeof(init), 1);
}

void pack(intermediate_state3 &inter)
{
    pvm_pkbyte((char*)&inter, sizeof(inter), 1);
}

void pack(final_state3 &final)
{
    pvm_pkbyte((char*)&final, sizeof(final), 1);
}

void unpack(initial_state3 &init)
{
    pvm_upkbyte((char*)&init, sizeof(init), 1);
}

void unpack(intermediate_state3 &inter)
{
    pvm_upkbyte((char*)&inter, sizeof(inter), 1);
}

void unpack(final_state3 &final)
{
    pvm_upkbyte((char*)&final, sizeof(final), 1);
}

#endif
