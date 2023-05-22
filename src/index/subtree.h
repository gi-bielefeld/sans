#pragma once 

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

#include "../kmer.h"
#include "../color.h"

// Index includes
#include "idQueue.h"
#include "q0Color.h"
#include "q1Color.h"
#include "q2Color.h"


/*
* The Subtree class implements the leave collection of a binary new_color tree routed,
* such that all colors stored here have the same amount of leading zeros
*/

template<typename T>
class subtree
{
    public:
        mutex subtreeLock;     // A mutex to lock the tree when changing verticess
        IDQueue idQueue;        // A queue for the creation of per tree-unique identifiers of the stored colors 

        hash_map<uint32_t, T> color_by_id;             // Maps ColorID->Color
        hash_map<T, uint32_t> id_by_color;             // Maps Color->ColorID        

        // This method creates a new entry for a
        uint32_t incrementColorSupport(const uint64_t& new_color);

        uint32_t incrementColorSupport(const uint64_t& new_color, uint32_t& current_color_id);
};

#include "subtree.cpp"