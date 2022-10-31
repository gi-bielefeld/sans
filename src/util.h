#include <iostream>
#include <chrono>
using namespace std;

/**
 * This class contains some helpful utility functions.
 */
class util {

 private:

 public:

    /**
     * This function prints a minimal version of the help page.
     */
    static void print_help();

    /**
     * This function prints an extended version of the help page.
     */
    static void print_extended_help();

    /**
     * This function displays a duration in a human readable format.
     *
     * @param time duration
     * @return formatted string
     */
    static string format_time(const chrono::high_resolution_clock::duration& time);

 protected:

};
