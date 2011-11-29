#include "lsoda.h"
#include "lsoda_internal.h"
#include "common.h"

/* newly added static variables */

int      mord[3] = {0, 12, 5};
double   sm1[13] = {0., 0.5, 0.575, 0.55, 0.45, 0.35, 0.25, 0.2, 0.15, 0.1, 0.075, 0.05, 0.025};

struct common_t common = {0};


