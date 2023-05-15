#include <vector>
#include <thread>
#include <mutex>

#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include "tsl/sparse_map.h"
#include "tsl/sparse_set.h"

template <typename K, typename V>
    // using hash_map = unordered_map<K,V>;
    using hash_map = tsl::sparse_pg_map<K,V>;
template <typename T>
    // using hash_set = unordered_set<T>;
    using hash_set = tsl::sparse_pg_set<T>;

#include "kmer.h"
#include "color.h"

using namespace std;

class Index
{    
    public:
        static hash_map<uint64_t, color_t> color_by_id;
        static hash_map<color_t, uint64_t> id_by_color;
        static hash_map<kmer_t, uint64_t> kmer_color;
        static uint64_t next_address;
        
        static void add_colored_kmer(const kmer_t& kmer, const uint64_t& color);
};