#pragma once

/*
* This section defines the q1-compressed-compressed color key that is stored in the identifier / color index
*/

#define SECOND_QUARTILE_BOUNDARY (2 * maxN / 4 + 1) 

#define Q0DET_COLOR
#define CLASS_NAME   colorQ2_t
#define STORAGE_TYPE uintQ2_t
#define INDEX_TYPE   sizeQ2_t
#define BIT_LENGTH   (1 * maxN - SECOND_QUARTILE_BOUNDARY + 1)
#define LEX_INTEGER_COMPARATORS
#include "../byte.h"
