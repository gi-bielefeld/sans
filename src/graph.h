#include "tsl/sparse_map.h"
#include "tsl/sparse_set.h"
#include <set>

#include "tree.h"
#include "kmer.h"

using namespace std;

template <typename K, typename V>
  using hash_map = tsl::sparse_pg_map<K,V>;
template <typename T>
  using hash_set = tsl::sparse_pg_set<T>;

/**
 * This class manages the k-mer/color hash tables and split list.
 */
class graph {

 private:

    /**
     * This is a hash table mapping k-mers to their colors/splits.
     */
    static vector<hash_map<kmer_t, color_t>> kmer_table;

    /**
     * This is a hash table mapping colors/splits to their weights.
     */
    static hash_map<color_t, array<uint32_t,2>> color_table;

 public:

    /**
     * This is the size of the top list.
     */
    static uint64_t t;

    /**
     * This function initializes the top list size.
     *
     * @param top_size top list size
     */
    static void init(const uint64_t& T, const uint64_t& top_size);

    /**
     * This function extracts k-mers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     * @param reverse merge complements
     */
    static void add_kmers(const uint64_t& T, const string& str, const size1N_t& color, const bool& reverse);

    /**
     * This function extracts k-mer minimizers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     * @param reverse merge complements
     * @param window number of k-mers to minimize
     */
    static void add_minimizers(const uint64_t& T, const string& str, const size1N_t& color, const bool& reverse, const uint64_t& window);

    /**
     * This function extracts k-mers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     * @param reverse merge complements
     * @param max_iupac allowed number of ambiguous k-mers per position
     */
    static void add_kmers(const uint64_t& T, const string& str, const size1N_t& color, const bool& reverse, const uint64_t& max_iupac);

    /**
     * This function extracts k-mer minimizers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     * @param reverse merge complements
     * @param window number of k-mers to minimize
     * @param max_iupac allowed number of ambiguous k-mers per position
     */
    static void add_minimizers(const uint64_t& T, const string& str, const size1N_t& color, const bool& reverse, const uint64_t& window, const uint64_t& max_iupac);

    /**
     * This function iterates over the hash table and outputs all the k-mer/color pairs.
     *
     * @param kmer string to store the k-mer
     * @param color string to store the color
     * @return iterator function
     */
    static function<bool(string&, string&)> lookup_kmer();

    /**
     * This function iterates over the hash table and outputs matching k-mer/color pairs.
     *
     * @param query query sequence
     * @param reverse merge complements
     * @param kmer string to store the k-mer
     * @param color string to store the color
     * @return iterator function
     */
    static function<bool(string&, string&)> lookup_kmer(const string& query, const bool& reverse);

    /**
     * This function iterates over the hash table and calculates the split weights.
     *
     * @param mean weight function
     * @param verbose print progress
     */
    static void calc_weights(const function<double(const uint32_t&, const uint32_t&)>& mean, const bool& verbose);

    /**
     * This function adds a single split (weight and colors) to the output list.
     *
     * @param weight split weight
     * @param color split colors
     */
    static void insert_split(const double& weight, const color_t& color);

    /**
     * This function merges two thread-separate hash tables.
     *
     * @param T1 first thread
     * @param T2 second thread
     */
    static void merge_threads(const uint64_t& T1, const uint64_t& T2);

    /**
     * This function destructs a thread-separate hash table.
     *
     * @param T thread index
     */
    static void erase_thread(const uint64_t& T);

 protected:

    /**
     * This function calculates the multiplicity of iupac k-mers.
     *
     * @param product overall multiplicity
     * @param factors per base multiplicity
     * @param chr iupac character
     */
    static void iupac_multiply(long double& product, vector<uint8_t>& factors, const char& chr);

    /**
     * This function shifts a base into a set of ambiguous iupac k-mers.
     *
     * @param prev set of k-mers
     * @param next set of k-mers
     * @param chr iupac character
     */
    static void iupac_shift(hash_set<kmer_t>& prev, hash_set<kmer_t>& next, const char& chr);

};
