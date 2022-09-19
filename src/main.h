#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <thread>

#include "graph.h"
#include "util.h"

#ifndef MAX_KMER_SIZE
    #define MAX_KMER_SIZE maxK+1
#endif
#ifdef useBF
    #include <bifrost/CompactedDBG.hpp>
    #include <bifrost/ColoredCDBG.hpp>
#endif

using namespace std;

// Symmetric Alignment-free phylogeNomic Splits
// simple efficient re-implementation + filters
#define SANS_VERSION "2.0_10E"    // SANS serif

/**
 * This is the entry point of the program.
 *
 * @param argc number of cmd args
 * @param argv cmd args
 * @return exit status
 */
int main(int argc, char* argv[]);
