#include "kmerAminoXX.h"

/**
 * This is the length of a k-mer.
 */
uint64_t kmerAminoXX::k;

/**
 * This is a bit-mask to erase all bits that exceed the k-mer length.
 */
bitset<5*maxK> kmerAminoXX::mask = mask.set();

/**
 * This function initializes the k-mer length and bit-mask.
 *
 * @param kmer_length k-mer length
 */
void kmerAminoXX::init(uint64_t& kmer_length) {
    k = kmer_length; mask.reset();
    for (uint64_t i = 0; i < 5*k; ++i) {
        mask <<= 01u;    // fill all bits within the k-mer length with ones
        mask |= 01u;    // the remaining zero bits can be used to mask bits
    }
}

/**
 * This function shifts a k-mer adding a new character to the left.
 * * http://www.cplusplus.com/forum/general/97378/
 * @param kmer bit sequence
 * @param c left character
 * @return right character
 */
char kmerAminoXX::shift_left(bitset<5 * maxK>& kmer, char& c) {
    uint64_t left = util::amino_char_to_bits(c);    // new leftmost character
    uint64_t right = 16*kmer[4]+8*kmer[3]+4*kmer[2]+2*kmer[1]+1*kmer[0];    // old rightmost character

    kmer >>= 05u;    // shift all current bits to the right by five positions

    for(int i = 4; i>=0; i--){   // encode the new character within the leftmost five bits
        kmer[5*k-(5-i)] = (left >> i) & 0x1;
    }

    kmer &= mask;    // set all bits to zero that exceed the k-mer length

    return util::amino_bits_to_char(right);    // return the dropped rightmost character
}


/**
 * This function shifts a k-mer adding a new character to the right.
 * http://www.cplusplus.com/forum/general/97378/
 * @param kmer bit sequence
 * @param c right character
 * @return left character
 */
char kmerAminoXX::shift_right(bitset<5 * maxK>& kmer, char& c) {
    uint64_t left = 16*kmer[5*k-1]+8*kmer[5*k-2]+4*kmer[5*k-3]+2*kmer[5*k-4]+1*kmer[5*k-5];    // old leftmost character
    uint64_t right = util::amino_char_to_bits(c);    // new rightmost character

    kmer <<= 05u;    // shift all current bits to the left by five positions
    for(int i = 4; i>=0; i--){// encode the new character within the rightmost five bits
        kmer[i] = (right >> i) & 0x1;
    }

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
bool kmerAminoXX::reverse_complement(bitset<5 * maxK>& kmer, bool minimize) {
    return false;
}


