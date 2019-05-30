#include "main.h"
#include "poisson.h"
#include "fields.h"
#include "boundary.h"

typedef struct particle{
    REAL x;
    REAL y;
    REAL u;
    REAL v;
    int onScreen : 2;
} particle;
