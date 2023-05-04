#include "mcl/concurrent_queue.h"
#include "tsl/sparse_set.h"
#include "tsl/sparse_map.h"

#include "tree.h"
#include "kmer.h"
using namespace std;

template <typename K, typename V>
  using queue = mcl::concurrent_queue<pair<K,V>>;
template <typename T>
  using hash_set = tsl::sparse_pg_set<T>;
template <typename K, typename V>
  using hash_map = tsl::sparse_pg_map<K,V>;

/**
 * This class manages the k-mer/color hash tables and splits calculation.
 */
class graph {

 private:

    /**
     * This is the min. coverage threshold for k-mers.
     */
    static uint64_t quality;

    /**
     * This is [Q] hash tables mapping k-mers to their colors/splits.
     */
    static vector<hash_map<kmer_t, color_t>> kmer_table;

    /**
     * This is [1] hash table mapping colors/splits to their weights.
     */
    static hash_map<color_t, uint32_t[2]> color_table;

    /**
     * This is [P] hash sets used to filter k-mers for coverage (q > 1).
     */
    static vector<hash_set<kmer_t>> quality_set;

    /**
     * This is [P] hash maps used to filter k-mers for coverage (q > 2).
     */
    static vector<hash_map<kmer_t, uint64_t>> quality_map;

    /**
     * This is [Q] queues used to synchronize k-mers from multiple threads.
     */
    static vector<queue<kmer_t, size1N_t>> thread_queue;

    /**
     * This is the number of k-mers to retrieve from the queue.
     */
    static uint64_t buffer;

 public:

    /**
     * This function initializes the thread queues and coverage threshold.
     *
     * @param quality coverage threshold
     * @param reverse merge complements
     * @param P number of file reading threads
     * @param Q number of queue hashing threads
     */
    static void init(const uint64_t& quality, const bool& reverse, const uint64_t& P, const uint64_t& Q);

    /**
     * This function extracts k-mers from a sequence and adds them to the hash table.
     *
     * @param T thread id [P]
     * @param str dna sequence
     * @param color color flag
     */
    static void add_kmers(const uint64_t& T, const string& str, const size1N_t& color);

    /**
     * This function extracts k-mer minimizers from a sequence and adds them to the hash table.
     *
     * @param T thread id [P]
     * @param str dna sequence
     * @param color color flag
     * @param window number of k-mers to minimize
     */
    static void add_minimizers(const uint64_t& T, const string& str, const size1N_t& color, const uint64_t& window);

    /**
     * This function extracts k-mers from a sequence and adds them to the hash table.
     *
     * @param T thread id [P]
     * @param str dna sequence
     * @param color color flag
     * @param max_iupac allowed number of ambiguous k-mers per position
     */
    static void add_kmers(const uint64_t& T, const string& str, const size1N_t& color, const uint64_t& max_iupac);

    /**
     * This function extracts k-mer minimizers from a sequence and adds them to the hash table.
     *
     * @param T thread id [P]
     * @param str dna sequence
     * @param color color flag
     * @param window number of k-mers to minimize
     * @param max_iupac allowed number of ambiguous k-mers per position
     */
    static void add_minimizers(const uint64_t& T, const string& str, const size1N_t& color, const uint64_t& window, const uint64_t& max_iupac);

    /**
     * This function iterates over the hash tables and outputs all the k-mer/color pairs.
     *
     * @param kmer string to store the k-mer
     * @param color string to store the color
     */
    static void lookup_kmer(const function<void(string&, string&)>& iterator);

    /**
     * This function iterates over the hash tables and outputs matching k-mer/color pairs.
     *
     * @param kmer string to store the k-mer
     * @param color string to store the color
     * @param query query sequence
     */
    static void lookup_kmer(const function<void(string&, string&)>& iterator, const string& query);

    /**
     * This function iterates over the hash tables and calculates the split weights.
     *
     * @param mean weight function
     * @param verbose print progress
     */
    static void calc_weights(const function<double(const uint32_t&, const uint32_t&)>& mean, const bool& verbose);

    /**
     * This function transfers k-mers & colors from the queue to the k-mer table.
     *
     * @param T thread id [Q]
     */
    static void merge_thread(const uint64_t& T);

    /**
     * This function clears a quality filter after full procession of a color.
     *
     * @param T thread id [P]
     */
    static void clear_filter(const uint64_t& T);

    /**
     * This function erases the quality filters and updates the distributor.
     */
    static void erase_filters();

 protected:

    /**
     * This function reverse complements a k-mer and applies a gap pattern.
     *
     * @param kmer bit sequence
     */
    static function<void(kmer_t&)> process_kmer;

    /**
     * This function restores a gap pattern for a right-compressed k-mer.
     *
     * @param kmer bit sequence
     */
    static function<void(kmer_t&)> restore_kmer;

    /**
     * This function qualifies a k-mer and places it into the hash table.
     *
     * @param T thread id [P]
     * @param kmer bit sequence
     * @param color color flag
     */
    static function<void(const uint64_t&, const kmer_t&, const size1N_t&)> emplace_kmer;

    /**
     * This function calculates the multiplicity of IUPAC k-mers.
     *
     * @param product overall multiplicity
     * @param factors per base multiplicity
     * @param chr IUPAC character
     */
    static void iupac_multiply(long double& product, vector<uint8_t>& factors, const char& chr);

    /**
     * This function shifts a base into a set of ambiguous IUPAC k-mers.
     *
     * @param prev set of k-mers
     * @param next set of k-mers
     * @param chr IUPAC character
     */
    static void iupac_shift(hash_set<kmer_t>& prev, hash_set<kmer_t>& next, const char& chr);

};
