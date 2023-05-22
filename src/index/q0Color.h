#pragma once

/*
* This section defines the q0-uncompressed-compressed color key that is stored in the color index
*/
#define Q0DET_COLOR
#define CLASS_NAME   colorQ0_t
#define STORAGE_TYPE uintQ0_t
#define INDEX_TYPE   sizeQ0_t
#define BIT_LENGTH   (1*maxN)
#define LEX_INTEGER_COMPARATORS
#include "../byte.h"
