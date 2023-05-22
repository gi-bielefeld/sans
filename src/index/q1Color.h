#pragma once
/*
* This section defines the q1-compressed-compressed color key that is stored in the identifier / color index
*/

#define FIRST_QUARTILE_BOUNDARY (1*maxN/4) 

#define Q0DET_COLOR
#define CLASS_NAME   colorQ1_t
#define STORAGE_TYPE uintQ1_t
#define INDEX_TYPE   sizeQ1_t
#define BIT_LENGTH   (1 * maxN - FIRST_QUARTILE_BOUNDARY)
#define LEX_INTEGER_COMPARATORS
#include "../byte.h"