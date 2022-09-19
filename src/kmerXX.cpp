#include "kmerXX.h"

/**
 * This is the length of a k-mer.
 */
uint64_t kmerXX::k;

/**
 * This is a bit-mask to erase all bits that exceed the k-mer length.
 */
bitset<2*maxK> kmerXX::mask = mask.set();

/**
 * This function initializes the k-mer length and bit-mask.
 *
 * @param kmer_length k-mer length
 */
void kmerXX::init(uint64_t& kmer_length) {
    k = kmer_length; mask.reset();
    for (uint64_t i = 0; i < 2*k; ++i) {
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
char kmerXX::shift_left(bitset<2*maxK>& kmer, char& c) {
    uint64_t left = char_to_bits(c);    // new leftmost character
    uint64_t right = 2*kmer[1]+kmer[0];    // old rightmost character

    kmer >>= 02u;    // shift all current bits to the right by two positions
    kmer[2*k-1] = left / 2;    // encode the new character within the leftmost two bits
    kmer[2*k-2] = left % 2;
    kmer &= mask;    // set all bits to zero that exceed the k-mer length

    return bits_to_char(right);    // return the dropped rightmost character
}

/**
 * This function shifts a k-mer adding a new character to the right.
 *
 * @param kmer bit sequence
 * @param c right character
 * @return left character
 */
char kmerXX::shift_right(bitset<2*maxK>& kmer, char& c) {
    uint64_t left = 2*kmer[2*k-1]+kmer[2*k-2];    // old leftmost character
    uint64_t right = char_to_bits(c);    // new rightmost character

    kmer <<= 02u;    // shift all current bits to the left by two positions
    kmer[1] = right / 2;    // encode the new character within the rightmost two bits
    kmer[0] = right % 2;
    kmer &= mask;    // set all bits to zero that exceed the k-mer length

    return bits_to_char(left);    // return the dropped leftmost character
}

/**
 * This function constructs the reverse complement of a given k-mer.
 *
 * @param kmer bit sequence
 * @param minimize only invert, if smaller
 * @return 1 if inverted, 0 otherwise
 */
bool kmerXX::reverse_complement(bitset<2*maxK>& kmer, bool minimize) {
    bitset<2*maxK> bits = kmer;    // copy the original k-mer
    bitset<2*maxK> rcmp;    // empty reverse complement

    for (uint64_t i = 0; i < k; ++i) {
        uint64_t base = 2*bits[1]+bits[0];
        bits >>= 02u;    // shift out the last base
        rcmp <<= 02u;    // shift in as the first base
        rcmp |= (~base & 0b11u);    // flip the bits
    }

    // if minimize == true, return the lexicographically smaller
    if (minimize) {
        for (uint64_t i = 2*k-1; i != -1; --i) {
            if (kmer[i] > rcmp[i]) {
                break;
            }
            if (kmer[i] < rcmp[i]) {
                return false;    // not reversed
            }
        }
    }
    kmer = rcmp;
    return true;    // reversed
}

/**
 * This function encodes a single character to two bits.
 *
 * @param c character
 * @return bit sequence
 */
uint64_t kmerXX::char_to_bits(char& c) {
    switch (c) {
        case 'A':
            return 0b00u;
        case 'C':
            return 0b01u;
        case 'G':
            return 0b10u;
        case 'T':
            return 0b11u;
        default:
            cerr << "Error: Invalid character " << c << "." << endl;
            return -1;
    }
}

/**
 * This function decodes two bits to a single character.
 *
 * @param b bit sequence
 * @return character
 */
char kmerXX::bits_to_char(uint64_t& b) {
    switch (b) {
        case 0b00u:
            return 'A';
        case 0b01u:
            return 'C';
        case 0b10u:
            return 'G';
        case 0b11u:
            return 'T';
        default:
            cerr << "Error: Invalid bit sequence " << b << "." << endl;
            return -1;
    }
}
