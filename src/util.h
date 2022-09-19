#include <iostream>
#include <chrono>
#include <cmath>

using namespace std;

/**
 * This class contains some helpful utility functions.
 */
class util {

private:

public:

    /**
     * This function calculates the arithmetic mean of two values.
     *
     * @param x first value
     * @param y second value
     * @return arithmetic mean
     */
    static double arithmetic_mean(uint32_t& x, uint32_t& y);

    /**
     * This function calculates the geometric mean of two values.
     *
     * @param x first value
     * @param y second value
     * @return geometric mean
     */
    static double geometric_mean(uint32_t& x, uint32_t& y);

    /**
     * This function calculates the geometric mean with pseudo counts.
     *
     * @param x first value
     * @param y second value
     * @return geometric mean
     */
    static double geometric_mean2(uint32_t& x, uint32_t& y);

    /**
     * This function displays a duration in a human readable format.
     *
     * @param time duration
     * @return formatted string
     */
    static string format_time(chrono::high_resolution_clock::duration time);

protected:

};
