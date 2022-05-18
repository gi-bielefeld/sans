#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;
/**
 *@author: Fabian Kolesch
 */

/**
 * The data class implements a generalized data interface
 * for tsv, newick and nexus files.
 * The Data object itself holds:
 * 	- Taxa
 * 	- Splits
 *
 */

class data{

private:
	string data_path; // The file to load from
	string data_format; // Format of the data file
	string taxa_path; // The file containing taxa (optional)

	vector<string> taxa; // The list of taxa
	vector<vector<uint64_t>> splits; // The set of splits;

	/* This method resolves the meta information of
	*  the input file and sets the corresponding instance fields.
	*  - format
	*  (This method not part of the constructor to avoid overloading)
	*/
	void test();

	// This method loads taxa from a file of files
	void load_taxa();
	// This method loads splits from a tsv file
	void load_tsv();
	// This method loads splits from a nexus file
	void load_nexus();
	// This method loads splits from a newick file
	void load_newick();

	// This method writs splits to a tsv file
	void save_tsv();
	// This method writes splits to a nexus file
	void save_nexus();
	// This method writes splits to a newick file
	void save_newick();

public:
	/* Taxa guided constructor of the data class
	*@Param: input_file: The file containing splits
	*@Param: string taxa_file: list of taxa
	*/
       	data(string input_file, string input_format, string taxa_path);
	
	/*
	 * Unguided constructor of the data class
	 * @Param: input_file: The file containing the splits
	 * NOTICE: If the taxa can not be extracted from the input file,
	 * they are named numerically!
	 */
	data(string input_file, string format);
	
	/*
	 * This method loads the splits and 
	 * taxa of the target dataset 
	 * into the data object.
	 */
	void load();
};

