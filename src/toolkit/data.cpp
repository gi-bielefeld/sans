#include "data.h"

#include <fstream>
#include <iostream>
#include <regex>
#include <unordered_map>

/*
 * @author: Fabian Kolesch
 */

/*
 * ---- The following code block implements the construction of a data object ----
 *
 */


/*
 * Constructor for taxa-less parsing
 */
data::data(string data_path, string data_format){
	// Throw an error if the taxa can not be inferred from the file
	if (data_format ==  "tsv"){
	cerr << "[Data Error] Taxa can not be inferred from a tsv file. Pleasy specify a list of taxa." << endl;
	exit(1);
	}

	this->data_path=data_path; // The data file
	this->data_format=data_format; // The format of the target file
	data::test(); // Ensure the data exists
}

/*
 * Constructor overload for taxa guided parsing
 *
 * This constructor allows for conversion of file formats 
 * that to not provide a list of taxa inherently.
 */
data::data(string data_path, string data_format, string taxa_path){
	this->data_path = data_path; // The path of the data file
	this->data_format = data_format; // The format of the data file
	this->taxa_path = taxa_path; // The path of the taxa file
	data::test(); // Ensure the data exists
}

/*
 * This method resolves the meta information of the 
 * input file. First it checks for existence and
 * readabilty. It extracts the file format by 
 * checking for a meaningful extension and file format.
 * It sets the fields:
 * - format
 *
 */
void data::test(){
	// --- Data file check
	// Check the data file exists
	ifstream file(this->data_path);
	if (!file.good()){
		cerr << "[Data Error] Could not read input file:\t" << this->data_path << endl;
		exit(1);
	}
	
	// Throw a warning if the data container file is empty
	string line;
	if (!getline(file, line)){
		cerr << "[Data Warning] Empty input file:\t" << data_path << endl;
	        file.close();
	}
	file.close();

	// --- Taxa file check
	// Ensure the taxa file exists
	if (!this->taxa_path.empty()){
		ifstream taxa_file(this->taxa_path);
		if(!taxa_file.good()){
			cerr << "[Data Error] Could not read taxa file:\t" << this->taxa_path << endl;
			exit(1);
		}
	}

}


/*
 * This method loads the 
 * splits and taxa into 
 * the data object.
 */
void data::load(){
}


/**
 * --- This codeblock implements the split loading 
 *  methods for the supported file formats
 */

/**
 * This method loads taxa from a file of files
 */
void data::load_taxa(){
}

/**
 * This method loads a set of splits from a .tsv file
 */
void data::load_tsv(){
}


