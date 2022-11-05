#include <iostream>
#include <bitset>
using namespace std;

#ifndef maxN   // max. color number defined as a preprocessor directive
    #define maxN 64
#endif

#if   1*maxN <= 8
    typedef uint_least8_t  uint1N_t;
#elif 1*maxN <= 16
    typedef uint_least16_t uint1N_t;
#elif 1*maxN <= 32
    typedef uint_least32_t uint1N_t;
#else
    typedef uint_least64_t uint1N_t;
#endif

#if   1*maxN <= 255
    typedef uint_fast8_t  size1N_t;
#elif 1*maxN <= 65535
    typedef uint_fast16_t size1N_t;
#elif 1*maxN <= 4294967295
    typedef uint_fast32_t size1N_t;
#else
    typedef uint_fast64_t size1N_t;
#endif

#if   1*maxN <= 64   // store colors as integers, optimizes performance
    typedef uint1N_t color_t;
    #define x_0b1u_(x) (x & 0b1u)
#else   // store colors in a bitset, allows arbitrarily many colors
    typedef bitset<maxN> color_t;
    #define x_0b1u_(x) (x[0])
#endif

/**
 * This class contains functions for working with color types.
 */
class color {

 private:

    /**
     * This is a bit-mask to erase all bits that exceed the color number.
     */
    static color_t mask;

 public:

    /**
     * This is the number of colors.
     */
    static size1N_t n;

    /**
     * This function initializes the color number and bit-mask.
     *
     * @param color_number color number
     */
    static void init(const size1N_t& color_number);

    /**
     * This function sets the bit at the given position to one.
     *
     * @param color bit sequence
     * @param pos position
     */
    static void set(color_t& color, const size1N_t& pos);

    /**
     * This function sets the bit at the given position to zero.
     *
     * @param color bit sequence
     * @param pos position
     */
    static void unset(color_t& color, const size1N_t& pos);

    /**
     * This function tests if the bit at the given position is set.
     *
     * @param color bit sequence
     * @param pos position
     * @return 1 if bit is set, 0 otherwise
     */
    static bool test(const color_t& color, const size1N_t& pos);

    /**
     * This function shifts a color appending a new bit char to the right.
     *
     * @param color bit sequence
     * @param chr right color bit
     */
    static void shift(color_t& color, const char& chr);

    /**
     * This function unshifts a color returning the bit char on the right.
     *
     * @param color bit sequence
     * @param chr right color bit
     */
    static void unshift(color_t& color, char& chr);

    /**
     * This function returns the index of the first bit set to one.
     *
     * @param color bit sequence
     * @return position (or -1 if all zero)
     */
    static size1N_t index(const color_t& color);

    /**
     * This function returns the number of bits that are set to one.
     *
     * @param color bit sequence
     * @return number of ones
     */
    static size1N_t count(const color_t& color);

    /**
     * This function constructs the bit complement of a given color set.
     *
     * @param color bit sequence
     */
    static void complement(color_t& color);

    /**
     * This function constructs the representative of a given color set.
     *
     * @param color bit sequence
     * @return 1 if inverted, 0 otherwise
     */
    static bool represent(color_t& color);

    /**
     * This function tests if two splits of colors are compatible.
     *
     * @param c1 bit sequence
     * @param c2 bit sequence
     * @return true, if compatible
     */
    static bool is_compatible(const color_t& c1, const color_t& c2);

    /**
     * This function tests if three splits of colors are weakly compatible.
     *
     * @param c1 bit sequence
     * @param c2 bit sequence
     * @param c3 bit sequence
     * @return true, if weakly compatible
     */
    static bool is_weakly_compatible(const color_t& c1, const color_t& c2, const color_t& c3);

 protected:

};
