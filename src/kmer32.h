#include <iostream>

using namespace std;

/**
 * This class contains functions for working with k-mers (maxK <= 32).
 */
class kmer32 {

private:

    /**
     * This is a bit-mask to erase all bits that exceed the k-mer length.
     */
    static uint64_t mask;

public:

    /**
     * This is the length of a k-mer.
     */
    static uint64_t k;

    /**
     * This function initializes the k-mer length and bit-mask.
     *
     * @param kmer_length k-mer length
     */
    static void init(uint64_t& kmer_length);

    /**
     * This function shifts a k-mer adding a new character to the left.
     *
     * @param kmer bit sequence
     * @param c left character
     * @return right character
     */
    static char shift_left(uint64_t& kmer, char& c);

    /**
     * This function shifts a k-mer adding a new character to the right.
     *
     * @param kmer bit sequence
     * @param c right character
     * @return left character
     */
    static char shift_right(uint64_t& kmer, char& c);

    /**
     * This function constructs the reverse complement of a given k-mer.
     *
     * @param kmer bit sequence
     * @param minimize only invert, if smaller
     * @return 1 if inverted, 0 otherwise
     */
    static bool reverse_complement(uint64_t& kmer, bool minimize);

protected:

    /**
     * This function encodes a single character to two bits.
     *
     * @param c character
     * @return bit sequence
     */
    static uint64_t char_to_bits(char& c);

    /**
     * This function decodes two bits to a single character.
     *
     * @param b bit sequence
     * @return character
     */
    static char bits_to_char(uint64_t& b);

};
