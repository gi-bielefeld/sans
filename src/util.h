#include <iostream>
#include <chrono>
using namespace std;

/**
 * This class can display the help page of the program.
 */
class help {

 public:

    /**
     * This function prints a minimal version of the help page.
     */
    static void print_usage();

    /**
     * This function prints an extended version of the help page.
     */
    static void print_extended_usage();

};

/**
 * This class contains some helpful utility functions.
 */
class util {

 public:

    /**
     * This function converts a command line argument to a string.
     *
     * @param argc argument count
     * @param argv argument vector
     * @param i argument index
     * @return string
     */
    static string atos(int& argc, char* argv[], int& i);

    /**
     * This function converts a command line argument to a number.
     *
     * @param argc argument count
     * @param argv argument vector
     * @param i argument index
     * @return number
     */
    static uint64_t aton(int& argc, char* argv[], int& i);

    /**
     * This function converts a string argument to a number.
     *
     * @param param parameter
     * @param args argument
     * @return number
     */
    static uint64_t ston(const string& param, const string& args);

    /**
     * This function displays a duration in a human readable format.
     *
     * @param time duration
     * @return formatted string
     */
    static string format_time(const chrono::high_resolution_clock::duration& time);

};
