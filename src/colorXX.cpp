#include "colorXX.h"

/**
 * This is the number of colors.
 */
uint64_t colorXX::n;

/**
 * This is a bit-mask to erase all bits that exceed the color number.
 */
bitset<maxN> colorXX::mask = mask.set();

/**
 * This function initializes the color number and bit-mask.
 *
 * @param color_number color number
 */
void colorXX::init(uint64_t& color_number) {
    n = color_number; mask.reset();
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
void colorXX::set(bitset<maxN>& color, uint64_t& pos) {
    color.set(pos) &= mask;
}

/**
 * This function returns the position of a single color.
 *
 * @param color bit sequence
 * @return position, or -1 if all zero
 */
uint64_t colorXX::pos(bitset<maxN>& color) {
    for (uint64_t i = 0; i < color.size(); i++) {
        if (color.test(i)) {
            return i;
        }
    }
    return -1;
}

/**
 * This function sets the bit at the given position to false.
 *
 * @param color bit sequence
 * @param pos position
 */
void colorXX::erase(bitset<maxN>& color, uint64_t& pos) {
    color.reset(pos) &= mask;
}

/**
 * This function tests if the bit at the given position is set.
 *
 * @param color bit sequence
 * @param pos position
 * @return 1 if bit is set, 0 otherwise
 */
bool colorXX::test(bitset<maxN>& color, uint64_t& pos) {
    return color.test(pos);
}

/**
 * This function inverts a color set to reduce the number of ones.
 *
 * @param color bit sequence
 * @param minimize only invert, if smaller
 * @return 1 if inverted, 0 otherwise
 */
bool colorXX::complement(bitset<maxN>& color, bool minimize) {
    uint64_t ones = color.count();    // count the number of ones
    // if minimize == true, return the color set with fewer ones
    if (minimize && (2*ones < n || (2*ones == n && color[0]))) {
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
int colorXX::size(bitset<maxN>& color, bool both) {
    uint64_t ones = color.count();    // count the number of ones
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
bool colorXX::is_compatible(bitset<maxN>& c1, bitset<maxN>& c2) {
    bitset<maxN> n1 = ~c1 & mask, n2 = ~c2 & mask;
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
bool colorXX::is_weakly_compatible(bitset<maxN>& c1, bitset<maxN>& c2, bitset<maxN>& c3) {
    bitset<maxN> n1 = ~c1 & mask, n2 = ~c2 & mask, n3 = ~c3 & mask;
    return ((c1 & c2 & c3) == 0b0u || (c1 & n2 & n3) == 0b0u || (n1 & c2 & n3) == 0b0u || (n1 & n2 & c3) == 0b0u)
        && ((n1 & n2 & n3) == 0b0u || (n1 & c2 & c3) == 0b0u || (c1 & n2 & c3) == 0b0u || (c1 & c2 & n3) == 0b0u);
}
