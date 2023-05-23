#include <vector>
#include <thread>
#include <mutex>

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
#include "kEntry.h"

#include "../color.h"
#include "q0Color.h"
#include "q2Color.h"

#include "subtree.h"

using namespace std;

class Index
{    
    public:

        // The number of hash-maps in the kmerMatrix (equal to the module 2**MOD_POWER + 1)
        static uint64_t bins;

        static vector<IDQueue> idQueue;
        
        // Kmer
        static vector<mutex> kmer_lock;
        static vector<hash_map<kmer_t, uint32_t>> kmerMatrix;
        // ID queue
        static vector<mutex> color_lock;

        static vector<hash_map<color_t, uint32_t>> id_by_color;
        static vector<hash_map<uint32_t, color_t>> color_by_id;
        static vector<hash_map<uint32_t, uint32_t>> support;

        static vector<uint32_t> color_period;

        static void init();

        static void add_colored_kmer(const kmer_t& kmer, uint64_t& bin, const uint64_t& color);

        /**
        *This method compresses a kmer_t into a k_entry vie increase of order after binning
        *@param kmer The kmer to compress
        *@param bin The target hash_table_id of this kmer
        *@return The compressed kmer_entry
        */
        static kmer_t compress_kmer(const kmer_t& kmer, uint64_t& bin);

        // Not implemented yet
        static kmer_t decompress_kmer(kmer_t& kmer_entry, uint64_t& bin);

};