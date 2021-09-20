#include "kmerAmino12.h"

/**
 * This is the length of a k-mer.
 */
uint64_t kmerAmino12::k;

/**
 * This is a bit-mask to erase all bits that exceed the k-mer length.
 */
uint64_t kmerAmino12::mask = -1;

/**
 * This function initializes the k-mer length and bit-mask.
 *
 * @param kmer_length k-mer length
 */
void kmerAmino12::init(uint64_t& kmer_length) {
    k = kmer_length; mask = 0;
    for (uint64_t i = 0; i < 5*k; ++i) {
        mask <<= 01u;    // fill all bits within the k-mer length with ones
        mask |= 01u;    // the remaining zero bits can be used to mask bits
    }
}

/**
 * This function shifts a k-mer adding a new character to the left.
 *
 * @param kmer bit sequence
 * @param c left character
 * @return right character
 */
char kmerAmino12::shift_left(uint64_t& kmer, char& c) {
    uint64_t left = util::amino_char_to_bits(c);    // new leftmost character
    uint64_t right = kmer & 0b11111u;    // old rightmost character

    kmer >>= 05u;    // shift all current bits to the right by five positions
    kmer |= left << (5*k-05u);    // encode the new character within the leftmost five bits
    kmer &= mask;    // set all bits to zero that exceed the k-mer length

    return util::amino_bits_to_char(right);    // return the dropped rightmost character
}

/**
 * This function shifts a k-mer adding a new character to the right.
 *
 * @param kmer bit sequence
 * @param c right character
 * @return left character
 */
char kmerAmino12::shift_right(uint64_t& kmer, char& c) {
    uint64_t left = (kmer >> (5*k-05u)) & 0b11111u;    // old leftmost character
    uint64_t right = util::amino_char_to_bits(c);    // new rightmost character

    kmer <<= 05u;    // shift all current bits to the left by five positions
    kmer |= right;    // encode the new character within the rightmost five bits
    kmer &= mask;    // set all bits to zero that exceed the k-mer length

    return util::amino_bits_to_char(left);    // return the dropped leftmost character
}

/**
 * This function constructs the reverse complement of a given k-mer.
 *
 * @param kmer bit sequence
 * @param minimize only invert, if smaller
 * @return 1 if inverted, 0 otherwise
 */
bool kmerAmino12::reverse_complement(uint64_t& kmer, bool minimize) {
   return false;
}

