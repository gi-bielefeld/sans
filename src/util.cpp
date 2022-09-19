#include "util.h"

/**
 * This function calculates the arithmetic mean of two values.
 *
 * @param x first value
 * @param y second value
 * @return arithmetic mean
 */
double util::arithmetic_mean(uint32_t& x, uint32_t& y) {
    return x / 2.0 + y / 2.0;
}

/**
 * This function calculates the geometric mean of two values.
 *
 * @param x first value
 * @param y second value
 * @return geometric mean
 */
double util::geometric_mean(uint32_t& x, uint32_t& y) {
    return sqrt(x) * sqrt(y);
}

/**
 * This function calculates the geometric mean with pseudo counts.
 *
 * @param x first value
 * @param y second value
 * @return geometric mean
 */
double util::geometric_mean2(uint32_t& x, uint32_t& y) {
    return sqrt(x+1) * sqrt(y+1) - 1;
}

/**
 * This function displays a duration in a human readable format.
 *
 * @param time duration
 * @return formatted string
 */
string util::format_time(chrono::high_resolution_clock::duration time) {
    double value = chrono::duration_cast<chrono::milliseconds>(time).count();
    string unit = "ms";
    if (value >= 1000) {
        value /= 1000;
        unit = "sec";
        if (value >= 200) {
            value /= 60;
            unit = "min";
            if (value >= 200) {
                value /= 60;
                unit = "h";
            }
        }
    }
    string number = to_string(value);
    return number.substr(0, number.find('.')+2) + ' ' + unit;
}
