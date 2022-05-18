#include "file_handler.h"
#include "data.h"

#include <cstring>
#include <regex>

/**
* The file_hanlder class allows to import and export splits
* from and to tsv, nexus and newick files
*/

/**
 * Try to match the given file name to one of the supported formats
 * @param string; file_name; The file_name to process
 * @return string; format; The found format extension (empty string if not matched)
 */
string fileHandler::get_format_by_file_name(string file_name){
	string format;
	if (regex_match(file_name, regex("(.*).tsv$"))){format="tsv"; cout << "FormatIsTSV" << endl;}
	return format;
}

/**
 *
 *
 */

