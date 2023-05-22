#pragma once
#include <iostream>
using namespace std;


/*
* The mod power defines the power of two to use as a basis for the 
* Kmer redistribution *
*/
#ifndef MOD_POWER
#define MOD_POWER 16
#endif



/*
* This section defines the compressed kmer key that is stored in the kmer / color index
* For a kmer using 2*k bits,
* The kmer_key uses 2*k - MOD_POWER bits
*/
#define KMER_KEY_BITS ((2*maxK - MOD_POWER) > 0 ? (2*maxK - MOD_POWER) : 0)

#define CLASS_NAME   kmerKey_t
#define STORAGE_TYPE uint2KM_t
#define INDEX_TYPE   size2KM_t
#define BIT_LENGTH   K_ENTRY_BITS
#define LEX_INTEGER_COMPARATORS
#include "../byte.h"


/*
* This section defines the identifier of a determining-color-tree and a color within that tree
*/
#define TREE_IDENTIFICATION_BITS 32
#define COLOR_IDENTIFICATION_BITS 32
