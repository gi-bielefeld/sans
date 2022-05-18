#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstring>

using namespace std;


class fileHandler {

public:
    /**
    * Try to extract the file format from the given file name 
    */
    static string get_format_by_file_name(string file_name);
    
};
