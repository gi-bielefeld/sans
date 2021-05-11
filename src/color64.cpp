#include "color64.h"

/**
 * This is the number of colors.
 */
uint64_t color64::n;

/**
 * This is a bit-mask to erase all bits that exceed the color number.
 */
uint64_t color64::mask = -1;

/**
 * This function initializes the color number and bit-mask.
 *
 * @param color_number color number
 */
void color64::init(uint64_t& color_number) {
    n = color_number; mask = 0;
    for (uint64_t i = 0; i < n; ++i) {
        mask <<= 01u;    // fill all bits within the color number with ones
        mask |= 01u;    // the remaining zero bits can be used to mask bits
    }
}

/**
 * This function sets the bit at the given position to true.
 *
 * @param color bit sequence
 * @param pos position
 */
void color64::set(uint64_t& color, uint64_t& pos) {
    color |= (0b1ull << pos) & mask;
}

/**
 * This function returns the position of a single color.
 *
 * @param color bit sequence
 * @return position
 */
uint64_t color64::pos(uint64_t& color) {
    return log2(color);
}

/**
 * This function sets the bit at the given position to false.
 *
 * @param color bit sequence
 * @param pos position
 */
void color64::erase(uint64_t& color, uint64_t& pos) {
    color &= ~(0b1ull << pos) & mask;
}

/**
 * This function tests if the bit at the given position is set.
 *
 * @param color bit sequence
 * @param pos position
 * @return 1 if bit is set, 0 otherwise
 */
bool color64::test(uint64_t& color, uint64_t& pos) {
    return (color >> pos) & 0b1u;
}

/**
 * This function inverts a color set to reduce the number of ones.
 *
 * @param color bit sequence
 * @param minimize only invert, if smaller
 * @return 1 if inverted, 0 otherwise
 */
bool color64::complement(uint64_t& color, bool minimize) {
    uint64_t bits = color;    // copy the original color
    uint64_t ones = 0;    // counter for the number of ones

    for (uint64_t i = 0; i < n; ++i) {
        ones += bits & 0b1u;    // count the last bit
        bits >>= 01u;    // shift to the next bit pos.
    }

    // if minimize == true, return the color set with fewer ones
    if (minimize && (2*ones < n || (2*ones == n && (color & 0b1u)))) {
        return false;    // not inverted
    } else {
        color = ~color & mask;    // flip the bits
        return true;    // inverted
    }
}

/**
 * This function returns the number of ones, or - if both is true - the number of zeros if it is larger.
 *
 * @param color bit sequence
 * @param both consider both zeros and ones and return max of both
 * @return number of ones (or zeros)
 */
int color64::size(uint64_t& color, bool both) {
    uint64_t bits = color;    // copy the original color
    uint64_t ones = 0;    // counter for the number of ones

    for (uint64_t i = 0; i < n; ++i) {
        ones += bits & 0b1u;    // count the last bit
        bits >>= 01u;    // shift to the next bit pos.
    }
    if (!both) { return ones; }

    uint64_t zeros = n - ones;
    if (ones > zeros) { return ones; }
    else { return zeros; }
}

/**
 * This function tests if two splits of colors are compatible.
 *
 * @param c1 bit sequence
 * @param c2 bit sequence
 * @return true, if compatible
 */
bool color64::is_compatible(uint64_t& c1, uint64_t& c2) {
    uint64_t n1 = ~c1 & mask, n2 = ~c2 & mask;
    return ((c1 & c2) == 0b0u || (c1 & n2) == 0b0u || (n1 & c2) == 0b0u || (n1 & n2) == 0b0u);
}

/**
 * This function tests if three splits of colors are weakly compatible.
 *
 * @param c1 bit sequence
 * @param c2 bit sequence
 * @param c3 bit sequence
 * @return true, if weakly compatible
 */
bool color64::is_weakly_compatible(uint64_t& c1, uint64_t& c2, uint64_t& c3) {
    uint64_t n1 = ~c1 & mask, n2 = ~c2 & mask, n3 = ~c3 & mask;
    return ((c1 & c2 & c3) == 0b0u || (c1 & n2 & n3) == 0b0u || (n1 & c2 & n3) == 0b0u || (n1 & n2 & c3) == 0b0u)
        && ((n1 & n2 & n3) == 0b0u || (n1 & c2 & c3) == 0b0u || (c1 & n2 & c3) == 0b0u || (c1 & c2 & n3) == 0b0u);
}

/**
 * This function returns the number of ones of the color.
 * @param color bit sequence
 * @return number of ones
 */
int color64::numberOnes(uint64_t &color) {
    uint64_t bits = color;    // copy the original color
    uint64_t ones = 0;    // counter for the number of ones

    for (uint64_t i = 0; i < n; ++i) {
        ones += bits & 0b1u;    // count the last bit
        bits >>= 01u;    // shift to the next bit pos.
    }
    return ones;
}

