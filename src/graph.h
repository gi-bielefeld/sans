#include "tsl/sparse_map.h"
#include "tsl/sparse_set.h"
#include <set>

#include "tree.h"
#include "kmer.h"

#if 2*maxK > 64
    #define bmi2_pext(x,y) (x &= y)
    #define bmi2_pdep(x,y)

#elif defined(__BMI2__)
    #include <immintrin.h>
    #if 2*maxK > 32
        #define bmi2_pext(x,y) (x = _pext_u64(x,y))
        #define bmi2_pdep(x,y) (x = _pdep_u64(x,y))
    #else
        #define bmi2_pext(x,y) (x = _pext_u32(x,y))
        #define bmi2_pdep(x,y) (x = _pdep_u32(x,y))
    #endif
#else
    #define bmi2_pext(kmer, pattern) {\
        uint2K_t lmer = 0, mask = pattern;\
        for (uint2K_t pos = 1; mask; pos <<= 1) {\
            if (kmer & mask & -mask) lmer |= pos;\
            mask &= mask-1;\
        }   kmer = lmer;\
    }
    #define bmi2_pdep(kmer, pattern) {\
        uint2K_t lmer = 0, mask = pattern;\
        for (uint2K_t pos = 1; mask; pos <<= 1) {\
            if (kmer & pos) lmer |= mask & -mask;\
            mask &= mask-1;\
        }   kmer = lmer;\
    }
#endif

using namespace std;

template <typename K, typename V>
  using hash_map = tsl::sparse_pg_map<K,V>;
template <typename T>
  using hash_set = tsl::sparse_pg_set<T>;

/**
 * This class manages the k-mer/color hash tables and splits calculation.
 */
class graph {

 private:

    /**
     * This is the pattern of gapped k-mers.
     */
    static kmer_t pattern;

    /**
     * This is the min. coverage threshold for k-mers.
     */
    static uint64_t quality;

    /**
     * This is a hash table mapping k-mers to their colors/splits.
     */
    static vector<hash_map<kmer_t, color_t>> kmer_table;

    /**
     * This is a hash table mapping colors/splits to their weights.
     */
    static hash_map<color_t, array<uint32_t,2>> color_table;

    /**
     * This is a hash set used to filter k-mers for coverage (q > 1).
     */
    static vector<hash_set<kmer_t>> quality_set;

    /**
     * This is a hash map used to filter k-mers for coverage (q > 2).
     */
    static vector<hash_map<kmer_t, uint64_t>> quality_map;

 public:

    /**
     * This function initializes the k-mer pattern and coverage threshold.
     *
     * @param pattern gapped k-mer pattern
     * @param quality coverage threshold
     * @param reverse merge complements
     */
    static void init(const uint64_t& T, const string& pattern, const uint64_t& quality, const bool& reverse);

    /**
     * This function extracts k-mers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     */
    static void add_kmers(const uint64_t& T, const string& str, const size1N_t& color);

    /**
     * This function extracts k-mer minimizers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     * @param window number of k-mers to minimize
     */
    static void add_minimizers(const uint64_t& T, const string& str, const size1N_t& color, const uint64_t& window);

    /**
     * This function extracts k-mers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     * @param max_iupac allowed number of ambiguous k-mers per position
     */
    static void add_kmers(const uint64_t& T, const string& str, const size1N_t& color, const uint64_t& max_iupac);

    /**
     * This function extracts k-mer minimizers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     * @param window number of k-mers to minimize
     * @param max_iupac allowed number of ambiguous k-mers per position
     */
    static void add_minimizers(const uint64_t& T, const string& str, const size1N_t& color, const uint64_t& window, const uint64_t& max_iupac);

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
     * @param kmer string to store the k-mer
     * @param color string to store the color
     * @return iterator function
     */
    static function<bool(string&, string&)> lookup_kmer(const string& query);

    /**
     * This function iterates over the hash table and calculates the split weights.
     *
     * @param mean weight function
     * @param verbose print progress
     */
    static void calc_weights(const function<double(const uint32_t&, const uint32_t&)>& mean, const bool& verbose);

    /**
     * This function merges two thread-separate hash tables.
     *
     * @param T1 first thread
     * @param T2 second thread
     */
    static void merge_threads(const uint64_t& T1, const uint64_t& T2);

    /**
     * This function clears thread-separate temporary files.
     *
     * @param T thread index
     */
    static void clear_thread(const uint64_t& T);

    /**
     * This function destructs a thread-separate hash table.
     *
     * @param T thread index
     */
    static void erase_thread(const uint64_t& T);

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
