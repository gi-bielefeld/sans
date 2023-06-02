#include <vector>
#include <thread>
#include <mutex>
#include <array>

#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include "../tsl/sparse_map.h"
#include "../tsl/sparse_set.h"

template <typename K, typename V>
    // using hash_map = unordered_map<K,V>;
    using hash_map = tsl::sparse_pg_map<K,V>;
template <typename T>
    // using hash_set = unordered_set<T>;
    using hash_set = tsl::sparse_pg_set<T>;



// Index includes
#include "../kmer.h"
// Color includes
#include "idQueue.h"
#include "../color.h"

#ifndef MOD_POWER
#define MOD_POWER (16)
#endif

using namespace std;

class Index
{    
    public:

        // The number of hash-maps in the kmerMatrix (equal to the module 2**MOD_POWER + 1)
        static uint64_t bins;

        static vector<IDQueue> idQueue;
        
        // Kmer
        static vector<mutex> kmer_lock;
        static vector<hash_map<kmer_t, uint64_t>> kmerMatrix;
        // ID queue
        static vector<mutex> color_lock;

        static vector<hash_map<color_t, uint64_t>> id_by_color;
        static vector<hash_map<uint64_t, color_t>> color_by_id;
        static vector<hash_map<uint64_t, array<uint32_t, 2>>> support;

        static uint64_t relevant_bits; // number of targets
        static vector<uint64_t> color_period;
        static uint64_t color_period_sum;

        static void init(uint64_t num);

        static void add_colored_kmer(const kmer_t& kmer, uint64_t& kmer_bin, const uint64_t& color);

        /**
        *This method compresses a kmer_t into a k_entry vie increase of order after binning
        *@param kmer The kmer to compress
        *@param kmer_bin The target hash_table_id of this kmer
        *@return The compressed kmer_entry
        */
        static kmer_t compress_kmer(const kmer_t& kmer, uint64_t& kmer_bin);

        // Not implemented yet
        static kmer_t decompress_kmer(kmer_t& kmer_entry, uint64_t& kmer_bin);

};