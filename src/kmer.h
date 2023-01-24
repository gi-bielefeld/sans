#include <iostream>
using namespace std;

#ifndef maxK     // max. k-mer length defined
#define maxK 32  // as preprocessor directive
#endif

#define CLASS_NAME   kmer_t
#define STORAGE_TYPE uint2K_t
#define INDEX_TYPE   size2K_t
#define BIT_LENGTH   (2*maxK)
#define LEX_INTEGER_COMPARATORS
#include "byte.h"

/**
 * This class contains functions for working with k-mer types.
 */
class kmer {

 private:

    /**
     * This is a bit-mask to erase all bits that exceed the k-mer length.
     */
    static kmer_t mask;

    /**
     * This is a bit-mask to erase certain bits in a gapped k-mer pattern.
     */
    static kmer_t _mask;

 public:

    /**
     * This is the length of a k-mer (including gap positions).
     */
    static size2K_t k;

    /**
     * This is the length of a k-mer (not counting gap positions).
     */
    static size2K_t _k;

    /**
     * This function initializes the k-mer length and bit-mask.
     *
     * @param length k-mer length
     * @param pattern gapped k-mer pattern
     */
    static void init(const size2K_t& length, const string& pattern);

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

    /**
     * This function applies a gap pattern and right-compresses the k-mer.
     *
     * @param kmer bit sequence
     */
    static inline void extract_pattern(kmer_t& kmer) noexcept {
        kmer.pext(_mask);
    }

    /**
     * This function applies a gap pattern and left-decompresses the k-mer.
     *
     * @param kmer bit sequence
     */
    static inline void deposit_pattern(kmer_t& kmer) noexcept {
        kmer.pdep(_mask);
    }

 protected:

    /**
     * This function encodes a single character as two bits.
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
