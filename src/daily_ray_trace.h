#include <Windows.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <float.h>

#include "types.h"
#include "spectrum.h"
#include "rng.h"
#include "geometry.h"

#include "spectrum.c"
#include "rng.c"
#include "geometry.c"

#ifndef NAN
#error "NAN not supported, dingus!"
#endif
