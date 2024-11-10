#include "color.h"

/*
 * This class contains functions for working with color types.
 */
size1N_t color::n;      // number of colors
color_t  color::mask;   // bit-mask to erase all bits that exceed the color number

/**
 * This function initializes the color number and bit-mask.
 *
 * @param number color number
 */
void color::init(const size1N_t& number) {
    n = number; mask = 0b0u;
    for (size1N_t i = 0; i < n; ++i)  // fill all bits within the color number with ones
        (mask <<= 01u) |= 0b1u;      // the remaining zero bits can be used to mask bits
}

/**
 * This function shifts a color appending a new bit char to the right.
 *
 * @param color bit sequence
 * @param chr right color bit
 */
void color::shift(color_t& color, const char& chr) {
    color <<= 01u;    // shift all current bits to the left by one position
    color |= chr-48;    // encode the new rightmost color bit char
    color &= mask;    // set all bits to zero that exceed the color number
}

/**
 * This function unshifts a color returning the bit char on the right.
 *
 * @param color bit sequence
 * @param chr right color bit
 */
void color::unshift(color_t& color, char& chr) {
    chr = (color & 0b1u)+48;    // return the rightmost color bit char
    color >>= 01u;    // shift all current bits to the right by one position
//  color &= mask;    // set all bits to zero that exceed the color number
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
 * This function constructs the bit complement of a given color set with respect to a custom mask.
 *
 * @param color the bit sequence to invert
 * @param mask the mask to be used to erase all other '1' bits
 */
void color::complement(color_t& color, const color_t& mask) {
    color = ~color & mask;    // flip the bits
}

/**
 * This function constructs the representative of a given color set.
 *
 * @param color bit sequence
 * @return 1 if inverted, 0 otherwise
 */
bool color::represent(color_t& color) {
    return represent(color, mask, n);
}

/*bool color::old_represent(color_t& color) {
    size1N_t count = color.popcnt();
    // return the color set with fewer ones to represent the split
    if (2*count < n || 2*count == n && (color & 0b1u))
        return false;    // not inverted
    else color = ~color & mask;    // flip the bits
    return true;    // inverted
}*/

/**
 * This function constructs the representative of a given color set respective to the supplied mask.
 *
 * @param color bit sequence
 * @param mask the custom mask to be used to erase all other '1' bits
 * @param size the number of '1' bits inside the given mask. Will be calculated if not present
 * @return 1 if inverted, 0 otherwise
 */
bool color::represent(color_t& color, const color_t& mask, const size1N_t& size) {
    size1N_t count = color.popcnt();
    if (2*count > size) {
        // flip the bits -> inverted
        color = ~color & mask;
        return true;
    } else if (2*count == size && ((~color & mask) < color)) {
        // flip the bits if both have same no. of 1 bits but the complement represents a smaller
        // (decimal) number
        // this has to be done, as the partial mask would otherwise lead to flips in both directions!
        // example: Mask: 1111 0000
        //          color:0110 0000
        //        flipped:1001 0000
        // old criterion (AND with bin. 1) would flag true for both as the mask does not cover that
        // region ...
        color = ~color & mask; // inverted
        return true;
    }
    // the complement would have fewer 1 bits -> do not flip the color
    // not inverted
    return false;

}
bool color::represent(color_t& color, const color_t& mask) {
    size1N_t size = mask.popcnt();
    return represent(color, mask, size);
}

/**
 * This function tests if two splits of colors are compatible.
 *
 * @param c1 bit sequence
 * @param c2 bit sequence
 * @return true, if compatible
 */
bool color::is_compatible(const color_t& c1, const color_t& c2) {
    color_t n1 = ~c1 & mask, n2 = ~c2 & mask; using _ = color_t;
    return (_::disjoint(c1, c2) || _::disjoint(c1, n2) || _::disjoint(n1, c2) || _::disjoint(n1, n2));
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
    color_t n1 = ~c1 & mask, n2 = ~c2 & mask, n3 = ~c3 & mask; using _ = color_t;
    return (_::disjoint(c1, c2, c3) || _::disjoint(c1, n2, n3) || _::disjoint(n1, c2, n3) || _::disjoint(n1, n2, c3))
        && (_::disjoint(n1, n2, n3) || _::disjoint(n1, c2, c3) || _::disjoint(c1, n2, c3) || _::disjoint(c1, c2, n3));
}

/**
* This function tests whether a given color set is the complete set of colors.
* 
* @param c color set to test
* @return true, if color set equals all colors
*/
bool color::is_complete(const color_t& c){
	return (~c & mask)==0b0u;
}

/**
* This function tests whether a given color set is the complete set of colors given a custom mask.
*
* @param c color set to test
* @param mask the custom mask to be used to erase all other '1' bits
* @return true, if color set equals all colors
*/
bool color::is_complete(const color_t& c, const color_t& mask){
    return (~c & mask)==0b0u;
}


/**
	* This function tests whether a given color set is only a single color.
	* 
	* @param c color set to test
	* @return true, if color set contains exactly one color
	*/
bool color::is_singleton(const color_t& c){
	return c.popcnt()==1;
}

