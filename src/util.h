#ifndef SANS_UTIL_H
#define SANS_UTIL_H


#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <vector>
#include <regex>
#include <sys/stat.h>
#include <mutex>

using namespace std;


/** 
Makro definitions for advanced output options
(Transcribed from kc)
*/
/*
#define $log  std::cerr
#define $err  std::cerr << "\33[38;5;131m"
#define $warn std::cerr << "\33[38;5;137m"
#define $note std::cerr << "\33[38;5;132m"

#define $link std::cerr << "\33[38;5;67m"
#define $lite std::cerr << "\33[38;5;244m"

#define  end$              std::endl
#define  end$$             std::endl, exit(0)
#define _end$  "\33[0m" << std::endl
#define _end$$ "\33[0m" << std::endl, exit(1)

*/


/**
Minimal output synchronization routine
**/
#define $out  std::cout
#define  $                 std::flush
#define $_              << "\33[2K\r"
#define _$     "\33[0m" << std::flush

// This declaration allows to initialize the sync_mutex outside of 
// this header file. The instantiation is to be found in util.cpp 
extern std::mutex sync_mutex;
#define $SYNC(STREAM)  sync_mutex.lock(), STREAM, sync_mutex.unlock()

/**
 * This class contains some helpful utility functions.
 */
class util {

private:

public:

    /**
	* This function compares the number of input genomes (n) to the compile parameter DmaxN (maxN).
	* Produces new makefile and exits (code 3) if they disagree (too much) and a re-compilation is necessary (n>maxN) or recommended (n much smaller maxN).
	*
	* @param n number of input genomes
	* @return nothing
	*/ 
	static void check_n(uint64_t& n, string &path, const uint64_t& max_N);


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


    /**
     * This function encodes a single character to two bits.
     *
     * @param c character
     * @return bit sequence
     */
    static uint64_t char_to_bits(const char& c);

    /**
     * This function decodes two bits to a single character.
     *
     * @param b bit sequence
     * @return character
     */
    static char bits_to_char(const uint64_t& b);

    /**
     * This function encodes a single character to five bits.
     *
     * @param c character
     * @return bit sequence
     */
    static uint64_t amino_char_to_bits(const char& c);

    /**
     * This function decodes five bits to a single character.
     *
     * @param b bit sequence
     * @return character
     */
    static char amino_bits_to_char(const uint64_t& b);

    /**
     * Trims leading whitespaces from the string.
     * @param s
     */
    static void ltrim(string &s);
    /**
    * Trims trailing whitespaces from the string.
    * @param s
    */
    static void rtrim(string &s);

    /**
    * Trims leading and trailing whitespaces from the string.
    * @param s
    */
    static void trim(string &s);

    /**
    * Replaces all occurrences of replaceString with the new string.
    * @param s the String which should be changed
    * @param replaceString the string, which should be replaced
    * @param newString the new string instead of the replaced one
    */
    static void replaceAll(string &s, string replaceString, string newString);

    /**
    * Splits the given string s into parts by a given delimiter.
    * @param s the string which should be splitted
    * @param delimiter the delimiter
    */
    static vector<string> split(string &s, string delimiter);

    /**
     * Checks if the string entered is a number.
     * @param s the string
     * @return bool
     */
    static bool is_number(string &s);
	
	
	  /**
     * Checks if the path in the file string is an existing path.
     * @param filename the string conting a path/to/file
     * @return bool
     */
	static bool path_exist(const string &filename);

	
protected:
};

#endif
