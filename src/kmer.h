#include <iostream>
#include <bitset>
using namespace std;

#ifndef maxK   // max. k-mer length defined as a preprocessor directive
    #define maxK 32
#endif

#if   2*maxK <= 8
    typedef uint_least8_t  uint2K_t;
#elif 2*maxK <= 16
    typedef uint_least16_t uint2K_t;
#elif 2*maxK <= 32
    typedef uint_least32_t uint2K_t;
#else
    typedef uint_least64_t uint2K_t;
#endif

#if   2*maxK <= 255
    typedef uint_fast8_t  size2K_t;
#elif 2*maxK <= 65535
    typedef uint_fast16_t size2K_t;
#elif 2*maxK <= 4294967295
    typedef uint_fast32_t size2K_t;
#else
    typedef uint_fast64_t size2K_t;
#endif

#if   2*maxK <= 64   // store k-mers as integers, optimizes performance
    typedef uint2K_t kmer_t;
    #define x_0b11u_(x) (x & 0b11u)
#else   // store k-mers in a bitset, allows arbitrarily large k-mers
    typedef bitset<2*maxK> kmer_t;
    #define x_0b11u_(x) (2*x[1]+x[0])

    namespace std {   // comparison function extending std::bitset
        template <size_t N>
        static inline bool operator<(const bitset<N>& x, const bitset<N>& y) {
            for (size_t i = N-1; i != -1; --i)
                if (x[i] ^ y[i]) return y[i];
            return false;
        }
    }
#endif

/**
 * This class contains functions for working with k-mer types.
 */
class kmer {

 private:

    /**
     * This is a bit-mask to erase all bits that exceed the k-mer length.
     */
    static kmer_t mask;

 public:

    /**
     * This is the length of a k-mer.
     */
    static size2K_t k;

    /**
     * This function initializes the k-mer length and bit-mask.
     *
     * @param kmer_length k-mer length
     */
    static void init(const size2K_t& kmer_length);

    /**
     * This function shifts a k-mer appending a new character to the right.
     *
     * @param kmer bit sequence
     * @param chr right character
     */
    static void shift(kmer_t& kmer, const char& chr);

    /**
     * This function unshifts a k-mer returning the character on the right.
     *
     * @param kmer bit sequence
     * @param chr right character
     */
    static void unshift(kmer_t& kmer, char& chr);

    /**
     * This function constructs the reverse complement of a given k-mer.
     *
     * @param kmer bit sequence
     */
    static void reverse_complement(kmer_t& kmer);

    /**
     * This function constructs the r.c. representative of a given k-mer.
     *
     * @param kmer bit sequence
     * @return 1 if inverted, 0 otherwise
     */
    static bool reverse_represent(kmer_t& kmer);

 protected:

    /**
     * This function encodes a single character to two bits.
     *
     * @param chr character
     * @return bit sequence
     */
    static uint2K_t char_to_bits(const char& chr);

    /**
     * This function decodes two bits to a single character.
     *
     * @param b bit sequence
     * @param chr character
     */
    static void bits_to_char(const uint2K_t& b, char& chr);

};
