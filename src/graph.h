#include <iostream>
#include <limits>
#include <algorithm>
#include <functional>
#include <utility>
#include <vector>

#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include "tsl/sparse_map.h"
#include "tsl/sparse_set.h"

using namespace std;

template <typename K, typename V>
    // using hash_map = unordered_map<K,V>;
    using hash_map = tsl::sparse_pg_map<K,V>;
template <typename T>
    // using hash_set = unordered_set<T>;
    using hash_set = tsl::sparse_pg_set<T>;

#include "kmer32.h"
#include "kmerXX.h"
#include "kmerAmino12.h"
#include "kmerAminoXX.h"

#if maxK > 32 // store k-mers in a bitset, allows larger k-mers
    typedef kmerXX kmer;
    typedef bitset<2*maxK> kmer_t;
#else // store k-mer bits in an integer, optimizes performance
    typedef kmer32 kmer;
    typedef uint64_t kmer_t;
#endif

#if maxK > 12 // store k-mers in a bitset, allows larger k-mers
    typedef kmerAminoXX kmerAmino;
    typedef bitset<5*maxK> kmerAmino_t;
#else // store k-mer bits in an integer, optimizes performance
    typedef kmerAmino12 kmerAmino;
    typedef uint64_t kmerAmino_t;
#endif



#include "color64.h"
#include "colorXX.h"

#if maxN > 64 // store colors in a bitset, allows more input files
    typedef colorXX color;
    typedef bitset<maxN> color_t;
#else // store color bits in an integer, optimizes performance
    typedef color64 color;
    typedef uint64_t color_t;
#endif

/**
 * A tree structure that is needed for generating a NEWICK string.
 */
struct node {
    color_t taxa;
    double weight;
    vector<node*> subsets;
};

/**
 * This class manages the k-mer/color hash tables and split list.
 */
class graph {

private:

    static bool isAmino;
    static bool considerOccurrences;

    /**
     * This is a hash table mapping k-mers to colors [O(1)].
     */
    static hash_map<kmer_t, color_t> kmer_table;
    /**
     * This is a hash table mapping k-mers to colors [O(1)].
     */
    static hash_map<kmerAmino_t, color_t> kmer_tableAmino;
    /**
     * This is a hash table mapping k-mers to occurrences-vectors.
     */
    static hash_map<kmer_t, vector<int>> copyNumber;
    /**
     * This is a hash table mapping k-mers to occurrences-vectors.
     */
    static hash_map<kmerAmino_t, vector<int>> copyNumberAmino;
    /**
     * This is a hash table mapping colors to weights [O(1)].
     */
    static hash_map<color_t, array<uint32_t,2>> color_table;
    /**
     * This is an int vector saving the occurrence of weights.
     */
    static vector<int> weight_Occ;

public:

    /**
     * This is an ordered tree collecting the splits [O(log n)].
     */
    static multimap<double, color_t, greater<>> split_list;

    /**
     * This is the size of the top list.
     */
    static uint64_t t;

    /**
    * These are the allowed chars.
    */
    static vector<char> allowedChars;


    /**
     * This function initializes the top list size and the allowed chars.
     *
     * @param t top list size
     * @param copyNumber count k-mer occurrences
     */
    static void init(uint64_t& top_size, bool isAmino, bool copyNumber);

    /**
     * This function extracts k-mers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     * @param reverse merge complements
     */
    static void add_kmers(string& str, uint64_t& color, bool& reverse);

    /**
     * This function extracts k-mer minimizers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     * @param reverse merge complements
     * @param m number of k-mers to minimize
     */
    static void add_minimizers(string& str, uint64_t& color, bool& reverse, uint64_t& m);

    /**
     * This function extracts k-mers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     * @param reverse merge complements
     * @param max_iupac allowed number of ambiguous k-mers per position
     */
    static void add_kmers(string& str, uint64_t& color, bool& reverse, uint64_t& max_iupac);

    /**
     * This function extracts k-mer minimizers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     * @param reverse merge complements
     * @param m number of k-mers to minimize
     * @param max_iupac allowed number of ambiguous k-mers per position
     */
    static void add_minimizers(string& str, uint64_t& color, bool& reverse, uint64_t& m, uint64_t& max_iupac);

    /**
     * This function calculates the split weight for a single entry of the hash table
     * 
     * @param mean weight function
     * @param verbose print progress
     * @param min_value the minimal weight represented in the top list
     * @return the new minimal weight represented in the top list
     */
    static double add_weight(color_t& color, double mean(uint32_t&, uint32_t&), double min_value, bool pos);

    /**
     * This function iterates over the hash table and calculates the split weights.
     *
     * @param mean weight function
     * @param verbose print progress
     * @param min_value the minimal weight currently represented in the top list
     */
    static void add_weights(double mean(uint32_t&, uint32_t&), double min_value, bool& verbose);

    /**
     * This function calculates the greatest common occurrence in all genomes from the occurrences-vector
     *
     * @param occurrences occurrences-vector for current k-mere
     * @return the weight that occurs in all genomes
     */
    static int get_greatestCommonWeight(vector<int>& occurrences);

    /**
     * This function iterates over the hash table of Copy Numbers and calculates the split weights.
     *
     * @param color split color
     * @param mean weight function
     * @param verbose print progress
     * @param min_value the minimal weight currently represented in the top list
     * @param occurrences occurrences-vector for current k-mere
     * @return the new minimal weight represented in the top list
     */
    static double add_weightsCopyNumber(color_t& color, double mean(uint32_t&, uint32_t&),
                                      double min_value, bool pos, vector<int>& occurrences);

    /**
     * This function adds a single split (weight and colors) to the output list.
     *
     * @param weight split weight
     * @param color split colors
     */
    static void add_split(double& weight, color_t& color);

    /**
     * This funtion adds a sigle split from a cdbg to the output list.
     * 
     * @param mean mean function
     * @param kmer_seq kmer
     * @param kmer_color split colors
     * @param min_value the minimal weight currently represented in the top list
     * @return The minimal weight currently represented in the top list
     */
     static double add_cdbg_colored_kmer(double mean(uint32_t&, uint32_t&), string kmer_seq, color_t& kmer_color, double min_value);       

    /**
     * This function filters a greedy maximum weight tree compatible subset.
     *
     * @param verbose print progress
     * @return the new minimal weight represented in the top list
     */
    static void filter_strict(bool& verbose);

    /**
     * This function filters a greedy maximum weight tree compatible subset and returns a newick string.
     *
     * @param map function that maps an integer to the original id, or null
     * @param verbose print progress
     */
    static string filter_strict(std::function<string(const uint64_t&)> map, bool& verbose);

    /**
     * This function filters a greedy maximum weight weakly compatible subset.
     *
     * @param verbose print progress
     */
    static void filter_weakly(bool& verbose);

    /**
     * This function filters a greedy maximum weight n-tree compatible subset.
     *
     * @param n number of trees
     * @param verbose print progress
     */
    static void filter_n_tree(uint64_t n, bool& verbose);

    /**
     * This function filters a greedy maximum weight n-tree compatible subset and returns a string with all trees in newick format.
     *
     * @param n number of trees
     * @param map function that maps an integer to the original id, or null
     * @param verbose print progress
     */
    static string filter_n_tree(uint64_t n, std::function<string(const uint64_t&)> map, bool& verbose);

protected:

    /**
     * This function tests if a split is compatible with an existing set of splits.
     *
     * @param color new split
     * @param color_set set of splits
     * @return true, if compatible
     */
    static bool test_strict(color_t& color, vector<color_t>& color_set);

    /**
     * This function tests if a split is weakly compatible with an existing set of splits.
     *
     * @param color new split
     * @param color_set set of splits
     * @return true, if weakly compatible
     */
    static bool test_weakly(color_t& color, vector<color_t>& color_set);

    /**
     * This function calculates the multiplicity of iupac k-mers.
     *
     * @param product overall multiplicity
     * @param factors per base multiplicity
     * @param input iupac character
     */
    static void iupac_calc(long double& product, vector<uint8_t>& factors, char& input);

    /**
     * This function shifts a base into a set of ambiguous iupac k-mers.
     *
     * @param prev set of k-mers
     * @param next set of k-mers
     * @param input iupac character
     */
    static void iupac_shift(hash_set<kmer_t>& prev, hash_set<kmer_t>& next, char& input);

    /**
   * This function shifts a base into a set of ambiguous iupac k-mers.
   *
   * @param prev set of k-mers
   * @param next set of k-mers
   * @param input iupac character
   */
    static void iupac_shift_amino(hash_set<kmerAmino_t>& prev, hash_set<kmerAmino_t>& next, char& input);

    /**
     * This function returns a tree structure (struct node) generated from the given list of color sets.
     *
     * @param color_set list of color sets
     * @return tree structure (struct node)
     */
    static node* build_tree(vector<color_t>& color_set);

    /**
     * This function recursively refines a given set/tree structure by a given split.
     *
     * @param current_set node of currently considered (sub-)set/tree structure
     * @param split color set to refine by
     * @return whether or not the given split is compatible with the set/tree structure
     */
    static bool refine_tree(node* current_set, color_t& split, color_t& allTaxa);

    /**
     * This function returns a newick string generated from the given tree structure (set).
     *
     * @param root root of the tree/set structure
     * @return newick string
     */
    static string print_tree(node* root, std::function<string(const uint64_t&)> map);

    /**
     * This function checks if the character at the given position is allowed.
     * @param pos position in str
     * @param str the current part of the sequence
     * @return true if allowed, false otherwise
     */
    static bool isAllowedChar(uint64_t pos, string &str);
};
