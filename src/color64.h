#ifndef COLOR64_H
#define COLOR64_H
#include <iostream>
#include <cmath>

using namespace std;

/**
 * This class contains functions for working with color sets (maxN <= 64).
 */
class color64 {

private:

    /**
     * This is a bit-mask to erase all bits that exceed the color number.
     */
    static uint64_t mask;

public:

    /**
     * This is the number of colors.
     */
    static uint64_t n;

    /**
     * This function initializes the color number and bit-mask.
     *
     * @param color_number color number
     */
    static void init(uint64_t& color_number);

    /**
     * This function sets the bit at the given position to true.
     *
     * @param color bit sequence
     * @param pos position
     */
    static void set(uint64_t& color, uint64_t& pos);

    /**
     * This function returns the position of a single color.
     *
     * @param color bit sequence
     * @return position
     */
    static uint64_t pos(uint64_t& color);

    /**
     * This function sets the bit at the given position to false.
     *
     * @param color bit sequence
     * @param pos position
     */
    static void erase(uint64_t& color, uint64_t& pos);

    /**
     * This function tests if the bit at the given position is set.
     *
     * @param color bit sequence
     * @param pos position
     * @return 1 if bit is set, 0 otherwise
     */
    static bool test(uint64_t& color, uint64_t& pos);

    /**
     * This function inverts a color set to reduce the number of ones.
     *
     * @param color bit sequence
     * @param minimize only invert, if smaller
     * @return 1 if inverted, 0 otherwise
     */
    static bool complement(uint64_t& color, bool minimize);

    /**
     * This function returns the number of ones, or - if both is true - the number of zeros if it is larger.
     *
     * @param color bit sequence
     * @param both consider both zeros and ones and return max of both
     * @return number of ones (or zeros)
     */
    static int size(uint64_t& color, bool both);

    /**
     * This function tests if two splits of colors are compatible.
     *
     * @param c1 bit sequence
     * @param c2 bit sequence
     * @return true, if compatible
     */
    static bool is_compatible(uint64_t& c1, uint64_t& c2);

    /**
     * This function tests if three splits of colors are weakly compatible.
     *
     * @param c1 bit sequence
     * @param c2 bit sequence
     * @param c3 bit sequence
     * @return true, if weakly compatible
     */
    static bool is_weakly_compatible(uint64_t& c1, uint64_t& c2, uint64_t& c3);

    /**
     * This function returns the number of ones
     * @param color bit sequence
     * @return number of ones
     */
    static int numberOnes(uint64_t& color);

protected:

};

#endif