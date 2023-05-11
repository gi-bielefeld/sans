#include "kmer.h"

/*
 * This class contains functions for working with k-mer types.
 */
size2K_t kmer::k;      // length of a k-mer (including gap positions)
kmer_t   kmer::mask;   // bit-mask to erase all bits that exceed the k-mer length

/**
 * This function initializes the k-mer length and bit-mask.
 *
 * @param length k-mer length
 */
void kmer::init(const size2K_t& length) {
    k = length; mask = 0b0u;
    for (size2K_t i = 0; i < k; ++i)  // fill all bits within the k-mer length with ones
        (mask <<= 02u) |= 0b11u;     // the remaining zero bits can be used to mask bits
}

/**
 * This function shifts a k-mer appending a new character to the right.
 *
 * @param kmer bit sequence
 * @param chr right character
 */
void kmer::shift(kmer_t& kmer, const char& chr) {
    kmer <<= 02u;    // shift all current bits to the left by two positions
    kmer |= char_to_bits(chr);    // encode the new rightmost character
    kmer &= mask;    // set all bits to zero that exceed the k-mer length
}

char kmer::shift_right(kmer_t& kmer, char& c) {
    uint64_t left = 2*kmer.test(2*k-1)+kmer.test(2*k-2);    // old leftmost character
    char c_left;

    uint64_t right = char_to_bits(c);    // new rightmost character

    kmer <<= 02u;    // shift all current bits to the left by two positions
    kmer |= right;    // encode the new character within the rightmost two bits
    kmer &= mask;    // set all bits to zero that exceed the k-mer length

    bits_to_char(left, c_left);    // return the dropped leftmost character
    return c_left;
}

/**
 * This function unshifts a k-mer returning the character on the right.
 *
 * @param kmer bit sequence
 * @param chr right character
 */
void kmer::unshift(kmer_t& kmer, char& chr) {
    bits_to_char(kmer & 0b11u, chr);    // return the rightmost character
    kmer >>= 02u;    // shift all current bits to the right by two positions
//  kmer &= mask;    // set all bits to zero that exceed the k-mer length
}

/**
 * This function constructs the reverse complement of a given k-mer.
 *
 * @param kmer bit sequence
 */
void kmer::reverse_complement(kmer_t& kmer) {
    kmer_t bits = ~kmer;    // flip the original k-mer
    kmer_t rcmp;
    for (size2K_t i = 0; i < k; ++i) {
        rcmp <<= 02u;    // shift in the first base
        rcmp |= bits & 0b11u;    // transfer the character
        bits >>= 02u;    // shift out the last base
    }
    kmer = rcmp;
}

/**
 * This function constructs the r.c. representative of a given k-mer.
 *
 * @param kmer bit sequence
 * @return 1 if inverted, 0 otherwise
 */
bool kmer::reverse_represent(kmer_t& kmer) {
    kmer_t bits = ~kmer;    // flip the original k-mer
    kmer_t rcmp;
    for (size2K_t i = 0; i < k; ++i) {
        rcmp <<= 02u;    // shift in the first base
        rcmp |= bits & 0b11u;    // transfer the character
        bits >>= 02u;    // shift out the last base
    }
    // return the lexicographically smaller
    if  (kmer < rcmp) return false;    // not reversed
    else kmer = rcmp; return true;    // reversed
}

/**
 * This function encodes a single character as two bits.
 *
 * @param chr character
 * @return bit encoding
 */
uint2K_t kmer::char_to_bits(const char& chr) {
    switch (chr) {
        case 'A': return 0b00u;
        case 'C': return 0b01u;
        case 'G': return 0b10u;
        case 'T': return 0b11u;
        default:
            cerr << "Error: invalid character: " << chr << endl;
            return -1;
    }
}

/**
 * This function decodes two bits to a single character.
 *
 * @param b bit encoding
 * @param chr character
 */
void kmer::bits_to_char(const uint2K_t& b, char& chr) {
    switch (b) {
        case 0b00u: chr = 'A'; break;
        case 0b01u: chr = 'C'; break;
        case 0b10u: chr = 'G'; break;
        case 0b11u: chr = 'T'; break;
        default:
            cerr << "Error: invalid bit encoding: " << b << endl;
    }
}
