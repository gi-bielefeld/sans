#include "kmer.h"
#include "ansi.h"

/*
 * This class contains functions for working with k-mer types.
 */
size2K_t kmer::k;      // length of a k-mer
kmer_t   kmer::mask;   // bit-mask to erase all bits that exceed the k-mer length

/**
 * This function initializes the k-mer length and bit-mask.
 *
 * @param kmer_length k-mer length
 */
void kmer::init(const size2K_t& kmer_length) {
    k = kmer_length; mask = 0b0u;
    for (size2K_t i = 0; i < k; ++i) {
        mask <<= 02u;    // fill all bits within the k-mer length with ones
        mask |= 0b11u;    // the remaining zero bits can be used to mask bits
    }
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
    kmer_t rcmp = 0b0u;    // empty reverse complement

    for (size2K_t i = 0; i < k; ++i) {
        rcmp <<= 02u;    // shift in the first base
        rcmp |= bits & 0b11u;
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
    kmer_t rcmp = 0b0u;    // empty reverse complement

    for (size2K_t i = 0; i < k; ++i) {
        rcmp <<= 02u;    // shift in the first base
        rcmp |= bits & 0b11u;
        bits >>= 02u;    // shift out the last base
    }
    // return the lexicographically smaller
    if  (kmer < rcmp) return false;    // not reversed
    else kmer = rcmp; return true;    // reversed
}

/**
 * This function applies a gap pattern and right-compresses the k-mer.
 *
 * @param kmer bit sequence
 * @param pattern bit mask
 */
void kmer::bmi2_pext(kmer_t& kmer, const kmer_t& pattern) {
    kmer.pext(pattern);
}

/**
 * This function applies a gap pattern and left-decompresses the k-mer.
 *
 * @param kmer bit sequence
 * @param pattern bit mask
 */
void kmer::bmi2_pdep(kmer_t& kmer, const kmer_t& pattern) {
    kmer.pdep(pattern);
}

/**
 * This function encodes a single character to two bits.
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
            $err << "Error: invalid character: " << chr << _end$$;
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
            $err << "Error: invalid bit encoding: " << b << _end$$;
    }
}
