#include "color.h"

/*
 * This class contains functions for working with color types.
 */
size1N_t color::n;      // number of colors
color_t  color::mask;   // bit-mask to erase all bits that exceed the color number

/**
 * This function initializes the color number and bit-mask.
 *
 * @param color_number color number
 */
void color::init(const size1N_t& color_number) {
    n = color_number; mask = 0b0u;
    for (size1N_t i = 0; i < n; ++i) {
        mask <<= 01u;    // fill all bits within the color number with ones
        mask |= 0b1u;    // the remaining zero bits can be used to mask bits
    }
}

/**
 * This function sets the bit at the given position to one.
 *
 * @param color bit sequence
 * @param pos position
 */
void color::set(color_t& color, const size1N_t& pos) {
   #if maxN <= 64
     color |= ((uint1N_t) 0b1u << pos);
   #else
     color.set(pos);
   #endif
}

/**
 * This function sets the bit at the given position to zero.
 *
 * @param color bit sequence
 * @param pos position
 */
void color::unset(color_t& color, const size1N_t& pos) {
   #if maxN <= 64
     color &= ~((uint1N_t) 0b1u << pos);
   #else
     color.reset(pos);
   #endif
}

/**
 * This function tests if the bit at the given position is set.
 *
 * @param color bit sequence
 * @param pos position
 * @return 1 if bit is set, 0 otherwise
 */
bool color::test(const color_t& color, const size1N_t& pos) {
   #if maxN <= 64
     return ((color >> pos) & 0b1u);
   #else
     return color.test(pos);
   #endif
}

/**
 * This function returns the index of the first bit set to one.
 *
 * @param color bit sequence
 * @return position (or -1 if all zero)
 */
size1N_t color::index(const color_t& color) {
    if (color == 0b0u) return -1;
   #if maxN <= 64
     size1N_t index;    // counter for the position
     color_t _color_ = color^(color-1);
     for (index = -1; _color_; ++index)
         _color_ >>= 1;    // count last position and shift to next bit
     return index;
   #else
     for (size1N_t i = 0; i < n; ++i)
         if (color.test(i)) return i;    // iterate over each position
     return -1;
   #endif
}

/**
 * This function returns the number of bits that are set to one.
 *
 * @param color bit sequence
 * @return number of ones
 */
size1N_t color::count(const color_t& color) {
   #if maxN <= 64
     size1N_t count;    // counter for the number of ones
     color_t _color_ = color;
     for (count = 0; _color_; ++count)
         _color_ &= _color_-1;    // count last bit and shift to next pos
     return count;
   #else
     return color.count();    // count the number of ones directly
   #endif
}

/**
 * This function constructs the bit complement of a given color set.
 *
 * @param color bit sequence
 */
void color::complement(color_t& color) {
    color = ~color & mask;    // flip the bits
}

/**
 * This function constructs the representative of a given color set.
 *
 * @param color bit sequence
 * @return 1 if inverted, 0 otherwise
 */
bool color::represent(color_t& color) {
   #if maxN <= 64
     size1N_t count;    // counter for the number of ones
     color_t _color_ = color;
     for (count = 0; _color_; ++count)
         _color_ &= _color_-1;    // count last bit and shift to next pos
   #else
     size1N_t count = color.count();    // count the number directly
   #endif
    // return the color set with fewer ones to represent the split
    if (2*count < n || 2*count == n && x_0b1u_(color))
        return false;    // not inverted
    else color = ~color & mask;    // flip the bits
        return true;    // inverted
}

/**
 * This function tests if two splits of colors are compatible.
 *
 * @param c1 bit sequence
 * @param c2 bit sequence
 * @return true, if compatible
 */
bool color::is_compatible(const color_t& c1, const color_t& c2) {
    color_t n1 = ~c1 & mask, n2 = ~c2 & mask;
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
bool color::is_weakly_compatible(const color_t& c1, const color_t& c2, const color_t& c3) {
    color_t n1 = ~c1 & mask, n2 = ~c2 & mask, n3 = ~c3 & mask;
    return ((c1 & c2 & c3) == 0b0u || (c1 & n2 & n3) == 0b0u || (n1 & c2 & n3) == 0b0u || (n1 & n2 & c3) == 0b0u)
        && ((n1 & n2 & n3) == 0b0u || (n1 & c2 & c3) == 0b0u || (c1 & n2 & c3) == 0b0u || (c1 & c2 & n3) == 0b0u);
}
