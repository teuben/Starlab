//
// tdyn_tt.C: pdyn- and tdyn-specific node-handling functions.
//

#include "tdyn.h"

real _pdyn_::get_radius()
{
    if (oldest_daughter) return 0;
    if (temperature <= 0 || luminosity <= 0) return 0;

    // Use the radius-temperature-luminosity relation and return
    // radius in solar units.

    real radius = sqrt(luminosity) / pow(temperature/5800, 2);
    return radius;
}
