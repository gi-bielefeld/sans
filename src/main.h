#include <cstring>
#include <fstream>
#include <thread>

#include "graph.h"
#include "util.h"

#ifdef useBF
    #include <bifrost/CompactedDBG.hpp>
    #include <bifrost/ColoredCDBG.hpp>
    #ifndef MAX_KMER_SIZE
        #define MAX_KMER_SIZE (maxK+1)
    #endif
#endif

using namespace std;

// Symmetric Alignment-free phylogeNomic Splits
// simple efficient re-implementation + filters
#define SANS_VERSION "2.0_10E-KC-2.2_11D"    // SANS-KC

/**
 * This is the entry point of the program.
 *
 * @param argc number of cmd args
 * @param argv cmd args
 * @return exit status
 */
int main(int argc, char* argv[]);
