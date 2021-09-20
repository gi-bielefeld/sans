#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>

#include "util.h"
#include "translator.h"
#include "cleanliness.h"

#ifdef useBF
    #ifndef MAX_KMER_SIZE
        #define MAX_KMER_SIZE (((maxK-1)/32) + 1) * 32
    #endif
    #include <bifrost/CompactedDBG.hpp>
    #include <bifrost/ColoredCDBG.hpp>
#endif

using namespace std;

// Symmetric Alignment-free phylogeNomic Splits
// simple efficient re-implementation + filters
#define SANS_VERSION "2.1_09A"    // SANS serif

/**
 * This is the entry point of the program.
 *
 * @param argc number of cmd args
 * @param argv cmd args
 * @return exit status
 */ 
int main(int argc, char* argv[]);
