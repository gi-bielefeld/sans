#include "converter_main.h"


/**
 * Converter entry point
 * @author Fabian Kolesch
 *
 *@param argc number of cmd args
 *@param argv cmd args
 *@return exit status
 */

int main(int argc, char*argv[]){
	
    /**
    * [Info]
    * --- Help page ---
    * - Print the help page to consile
    * - Describes the cmd arguments
    */
    if (argc <= 1 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0){
        cout << endl;
	// Header
	cout << "------------" << endl;
	cout << "SANS-convert" << endl;
	cout << "A standalone converter for TSV, Newick and Nexus files" << endl; 
	cout << "------------" << endl;
	// Usage
	cout << "Usage: SANS-convert [PARAMTERS]" << endl << endl;
	cout << "[Required arguments]" << endl;
	cout << "-i, --input <format> <file>	\tInput file" << endl;
	cout << "-o, --output <format> <file>	\tOutput file" << endl;
	cout << endl;
	cout << "[Optional arguments]" << endl;
	cout << "-t --taxa            		\tFile of files; List of taxa (required for .tsv files)" << endl;
   	cout << endl << endl;
	exit(0);
    }

    /**
    * [Meta]
    * --- Parameter parsing ---
    */
    cout << "Creating default" << endl; 
    string input_format;
    string input;

    string output;
    string output_format;
	
    string taxa;

    cout << "Parsing Params" << endl; 
    for (int i = 1; i < argc; ++i){
	// [Required arguments]
        // Input file
        if(strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input") == 0){
            input_format = argv[++i];
	    input = argv[++i];
        }
	// Output file
	else if(strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0){
	    output_format = argv[++i];
	    output = argv[++i];
	}
	// [Optional arguments]
	// Taxa file
	else if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--taxa") == 0){
	    taxa = argv[++i];
	}
    }
    cout << "Verifying Params" << endl;
    /**
     * --- Parameter verification ---
     */
    // Catch unknown input format
    if (!(input_format.compare("tsv")||input_format.compare("nexus")||input_format.compare("newick"))){
	cerr << "[Error] Unknown input format: " << input_format << endl;
	exit(1);
    }
    // Catch empty input
    else if(input.empty()){
	cerr << "[Error] No input specified" << endl;
	exit(1);
    }
    
    // Catch unknown output format
    else if (!((output_format.compare("tsv")==0)||output_format.compare("nexus")||output_format.compare("newick"))){
        cerr << "[Error] Unknown output format: " << output_format << endl;
	exit(1);
    }
    // Catch missing outout
    else if(output.empty()){
	cerr << "[Error] No output specified" << endl;
	exit(1);
    }
    
    /**
     * --- Dataset instanciation
     */
    
    // If there is no guide for taxa, build from input
    if (taxa.empty()){ data input_data(input, input_format);}

    // If there is a guide of taxa, build from input and guid
    else {data input_data(input, input_format, taxa);}
}
