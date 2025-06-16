#include "main.h"
#include <algorithm>
#include <regex>
// gzstream imports
#include "gz/gzstream.h"

/**
 * This is the entry point of the program.
 *
 * @param argc number of cmd args
 * @param argv cmd args
 * @return exit status
 */
int main(int argc, char* argv[]) {

    /**
    * [help page definition]
    * - show the help page
    * - describes the software and argument usage
    */

    auto show_help_page = [](){
        cout << endl;
        cout << "SANS ambages | version " << SANS_VERSION << endl;
        cout << "Usage: SANS [PARAMETERS]" << endl;
        cout << endl;
        cout << "  Input arguments:" << endl;
        cout << endl;
        cout << "    -i, --input   \t Input file: file of files format" << endl;
        cout << "                  \t Either: one genome per line (space-separated for multifile genomes)" << endl;
        cout << "                  \t Or: kmtricks input format (see https://github.com/tlemane/kmtricks)" << endl;
        cout << endl;
        cout << "    -g, --graph   \t Graph file: load a Bifrost graph, file name prefix" << endl;
        cout << "                  \t optional: provide additional input file (format as for -i) to filter the graph" << endl;
        cout << "                  \t           (multiple colors in graph can be assigned to one genome ID)" << endl;
        #ifndef useBF
        cout << "                  \t (requires compiler flag -DuseBF, please edit makefile)" << endl;
        #endif
        cout << endl;
        cout << "    -s, --splits  \t Splits file: load an existing list of splits file" << endl;
        cout << "                  \t (allows to filter -t/-f, other arguments are ignored)" << endl;
        cout << endl;
        cout << "    -B, --blacklist\t File (Fasta, Fastq) of k-mers to be ignored" << endl;
        cout << endl;
        cout << "    (either --input and/or --graph, or --splits must be provided)" << endl;
        cout << endl;
        cout << "  Output arguments:" << endl;
        cout << endl;
        cout << "    -o, --output  \t Output TSV file: list of splits, sorted by weight desc." << endl;
        cout << endl;
        cout << "    -N, --newick  \t Output Newick file" << endl;
        cout << "                  \t (only applicable in combination with -f strict or n-tree)" << endl;
        cout << endl;
        cout << "    -X, --nexus  \t Output Nexus file" << endl;
        cout << "                 \t (Warning: Already existing files will be overwritten)" << endl;
        cout << "                 \t Attention: For a reliable visualization using SplitsTree,\n"
                "                 \t split weights in the Nexus file are scaled to the range 0 to 1." << endl;
        cout << endl;
        cout << "    -p, --pdf  \t\t Output network as PDF file" << endl;
        cout << "                 \t Requires SplitsTree in the PATH" << endl;
        cout << "                 \t Warning: Already existing files will be overwritten" << endl;
        cout << endl;
        cout << "    -S, --svg  \t\t Output network as SVG file" << endl;
        cout << "                 \t Requires SplitsTree in the PATH" << endl;
        cout << "                 \t Warning: Already existing files will be overwritten" << endl;
        cout << endl;
        cout << "    -r, --core  \t Output core k-mers in fasta file" << endl;
        cout << endl;
        cout << "    -R, --raw  \t Output both counts per split in TSV file" << endl;
        cout << endl;
        cout << "    (at least --output, --newick, --nexus, --pdf, --svg, --core, or --raw must be provided)" << endl;
        cout << endl;
        cout << "  Optional arguments:" << endl;
        cout << endl;
        cout << "    -k, --kmer    \t Length of k-mers (default: 31, or 10 for --amino and --code)" << endl;
        cout << endl;
        // cout << "    -w, --window  \t Number of k-mers per minimizer window (default: 1)" << endl;
        // cout << endl;
        cout << "    -t, --top     \t Number of splits in the output list (default: all)." << endl;
        cout << "                  \t Use -t <integer>n to limit relative to number of input files, or" << endl;
        cout << "                  \t use -t <integer> to limit by absolute value." << endl;
        cout << endl;
        cout << "    -m, --mean    \t Mean weight function to handle asymmetric splits" << endl;
        cout << "                  \t options: arith: arithmetic mean" << endl;
        cout << "                  \t          geom:  geometric mean" << endl;
        cout << "                  \t          geom2: geometric mean with pseudo-counts (default)" << endl;
        cout << endl;
        cout << "    -f, --filter  \t Output a greedy maximum weight subset" << endl;
        // cout << "                  \t additional output: (weighted) cleanliness, i.e., ratio of" << endl;
        // cout << "                  \t filtered splits w.r.t. original splits" << endl;
        cout << "                  \t options: strict: compatible to a tree" << endl;
        cout << "                  \t          weakly: weakly compatible network" << endl;
        cout << "                  \t          planar: compatible to a planar graph" << endl;
        cout << "                  \t                  (a.k.a. circular compatible, outer labeled planar)" << endl;
        cout << "                  \t          n-tree: compatible to a union of n trees" << endl;
        cout << "                  \t                  (where n is an arbitrary number, e.g. 2-tree)" << endl;
        cout << endl;
        cout << "    -x, --iupac   \t Extended IUPAC alphabet, resolve ambiguous bases or amino acids" << endl;
        cout << "                  \t Specify a number to limit the k-mers per position between" << endl;
        cout << "                  \t 1 (no ambiguity) and 4^k respectively 22^k (allows NNN...N)" << endl;
        cout << "                  \t Without --iupac respective k-mers are ignored" << endl;
        cout << endl;
        cout << "    -q, --qualify \t Discard k-mers with lower coverage than a threshold" << endl;
        cout << endl;
        cout << "    -n, --norev   \t Do not consider reverse complement k-mers" << endl;
        cout << endl;
        cout << "    -a, --amino   \t Consider amino acids: --input provides amino acid sequences" << endl;
        cout << "                  \t Implies --norev and a default k of 10" << endl;
        cout << endl;
        cout << "    -c, --code    \t Translate DNA: --input provides coding sequences" << endl;
        cout << "                  \t Implies --norev and a default k of 10" << endl;
        cout << "                  \t optional: ID of the genetic code to be used" << endl;
        cout << "                  \t Default: 1" << endl;
        cout << "                  \t Use 11 for Bacterial, Archaeal, and Plant Plastid Code" << endl;
        cout << "                  \t (See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for details.)" << endl;
        cout << endl;
        cout << "    -M, --maxN    \t Compare number of input genomes to compile paramter DmaxN" << endl;
        cout << "                  \t Add path/to/makefile (default is makefile in current working directory)." << endl;
        cout << endl;
        cout << "    -b, --bootstrap \t Perform bootstrapping with the specified number of replicates" << endl;
        cout << "                  \t optional: provide threshold to filter low support splits (e.g. 0.75)" << endl;
        cout << endl;
        cout << "    -C, --consensus\t Apply final filter w.r.t. support values" << endl;
        cout << "                  \t else: final filter w.r.t. split weights" << endl;
        cout << "                  \t optional: specify separate filter (see --filter for available filters.)" << endl;
        cout << endl;
        cout << "    -l, --label\t\t Color taxa according to given groups" << endl;
        cout << "                  \t (Requires SplitsTree in the PATH)" << endl;
        cout << "                  \t required file: file with name of taxon and group (tab separated)" << endl;
        cout << "                  \t optional: additional file with group and " << endl;
        cout << "                  \t color (rgb values, e.g. 90 0 255) (tab separated)" << endl;
        cout << "                  \t Only applicable together with -X or -p " << endl;
        cout << endl;
        cout << "    -v, --verbose \t Print information messages during execution" << endl;
        cout << endl;
        cout << "    -T, --threads \t The number of threads to spawn (default is all)" << endl;
        cout << endl;
        cout << "    -h, --help    \t Display this help page and quit" << endl;
        cout << endl;
        cout << "  Contact: pangenomics-service@cebitec.uni-bielefeld.de" << endl;
        cout << "  Evaluation: https://www.surveymonkey.de/r/denbi-service?sc=bigi&tool=sans" << endl;
        cout << endl;
    };

    // show the help page if no args are given or the arg is --help 
    if (argc <= 1 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0 ) { show_help_page(); return(0);}
    // show the help page if SANS is run from autoN-script without further arguments or --help
    if ((strcmp(argv[1],"-M") == 0 ) && (argc <= 3 || strcmp(argv[3], "-h") == 0 || strcmp(argv[3], "--help") == 0 )) { show_help_page(); return(0);}


    /**
    * [argument and meta variable initialization]
    * defaults are set here
    */

    // file definitions
    string input;    // name of input file
    string graph;    // name of graph file
    string graph_filter;    // name of graph filter file
    string splits;    // name of splits file
    string blacklistfile; // name of blacklist file
    string output;    // name of output file
    string newick;    // name of newick output file // Todo
    string nexus;   // name of nexus output file
    string pdf;     // name of PDF output file
    string svg;     // name of SVG output file
    string core;     // name of file for core k-mers
    string raw;  // name of file for raw count output
    string groups; // name of input file giving groups
    string coloring; // name of input file for using specified color
    string translate; // name of translate file

    // input
    uint64_t num = 0;    // number of input files

    // automatic recompilation
    string path = "./makefile"; // path to makefile for automatic recompilation
    bool check_n = false; // compare num (number of input genomes) to maxN (compile parameter DmaxN)

    // kmer args --
    bool userKmer = false; // is k-mer default or custom
    uint64_t kmer = 31;    // length of k-mers
    uint64_t window = 1;    // number of k-mers in a minimizer window

    // kmer preprocessing and filtering
    bool reverse = true;    // consider reverse complement k-mers
    uint64_t iupac = 1;    // allow extended iupac characters
    int quality = 1;    // min. coverage threshold for k-mers (if individual q values per file are given, this is the maximum among all)

    // amino processing
    bool amino = false;      // input files are amino acid sequences
    bool shouldTranslate = false;   // translate input files
    uint64_t code = 1;

    // split processing
    auto mean = util::geometric_mean2;    // weight function
    string filter;    // filter function

    // output limiter
    uint64_t top = -1;    // maximal number of splits to compute
    bool dyn_top = false; // use multiplicity of inputs as maximal number of splits

    // parallel hashing
    uint64_t threads = thread::hardware_concurrency(); // The number of threads to run on (default is #cores including smt / ht)

    // bootsrapping
    string consensus_filter; // filter function for filtering after bootstrapping
	uint32_t bootstrap_no=0; // = no bootstrapping
	float bootstrap_threshold=0; // threshold to filter low support splits
	hash_map<color_t, uint32_t> support_values; // hash_map for each original split with zero counts

    // qol
    bool verbose = false;    // print messages during execution
	chrono::high_resolution_clock::time_point end;

    // simple nexus, colored nexus, pdf, svg
    bool nexus_wanted = false;
    bool c_nexus_wanted = false;
    bool pdf_wanted = false;
    bool svg_wanted = false;
	bool raw_wanted = false;
	
	/**
	* Look-up set for k-mers that are ignored, i.e., not stored, counted etc.
	*/
	hash_set<kmer_t> blacklist;
	hash_set<kmerAmino_t> blacklist_amino;


    /**
     * [argument parser]
     * - parse the command line arguments and update the meta variables accordingly
     */

    // this lambda function throws an error, if a dependent argument is NULL
    auto catch_missing_dependent_args = [] (auto arg, string parent)
    {
        if(arg == NULL){cerr << "Error: Missing dependent argument after " << parent << endl;exit(1);}
    };

    // this lambda function throws an error if the target arg-string could not be parsed as int 
    auto catch_failed_stoi_cast = [] (auto arg, string parent)
    {
    try {stoi(arg);}
    catch( invalid_argument &excp )    {cerr << "Error: Invalid dependent argument: '" << arg << "' at: " << parent << " " << arg << endl; cerr << endl; exit(1);}
    };

    // this lambda function asks for user confirmation
    auto ask_for_user_confirmation = [] ()
    {
        string user_confirmation;
        int requests = 0;
        cout << "Please enter 'yes' or 'no'" << endl;
        while (requests <= 5)
        {
            getline(cin, user_confirmation);
            if (user_confirmation.compare("yes") == 0){return true;}
            if (user_confirmation.compare("no") == 0){return false;}
            requests++;
        }
        return false;
    };   

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input") == 0) {
            catch_missing_dependent_args(argv[i + 1], argv[i]);
            input = argv[++i];    // Input file: list of sequence files, one per line
        }
        else if (strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--graph") == 0) {
            catch_missing_dependent_args(argv[i + 1], argv[i]);
            graph = argv[++i];    // Graph file: load a Bifrost graph, file name prefix
            #ifndef useBF
                cerr << "Error: requires compiler flag -DuseBF" << endl;
                return 1;
            #endif
            if (i+1 < argc && argv[i+1][0]!='-') { // optional filter file
                graph_filter = argv[++i];
            }
        }
        else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--splits") == 0) {
            catch_missing_dependent_args(argv[i + 1], argv[i]);
            splits = argv[++i];    // Splits file: load an existing list of splits file
        }
        else if (strcmp(argv[i], "-B") == 0 || strcmp(argv[i], "--blacklist") == 0) {
            catch_missing_dependent_args(argv[i + 1], argv[i]);
            blacklistfile = argv[++i];    // Blacklist file: load kmers to be ignored
        }
        else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
            catch_missing_dependent_args(argv[i + 1], argv[i]);
            output = argv[++i];    // Output file: list of splits, sorted by weight desc.
            if (!util::path_exist(output)){
				cerr << "Error: output folder does not exist: "<< output << endl;
                return 1;
			}
        }
        else if (strcmp(argv[i], "-N") == 0 || strcmp(argv[i], "--newick") == 0) {
            catch_missing_dependent_args(argv[i + 1], argv[i]);
            newick = argv[++i];    // Output newick file
            if (!util::path_exist(newick)){
				cerr << "Error: output folder does not exist: "<< newick << endl;
                return 1;
			}
        }
        else if (strcmp(argv[i], "-X") == 0 || strcmp(argv[i], "--nexus") == 0) {
            catch_missing_dependent_args(argv[i + 1], argv[i]);
            nexus = argv[++i];    // Nexus output file
            nexus_wanted = true;
            if (!util::path_exist(nexus)){
                cerr << "Error: output folder does not exist: "<< nexus << endl;
                return 1;
            }
        }
        else if (strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--core") == 0) {
            catch_missing_dependent_args(argv[i + 1], argv[i]);
            core = argv[++i];    // fasta output file for core k-mers
            if (!util::path_exist(core)){
                cerr << "Error: output folder does not exist: "<< core << endl;
                return 1;
            }
        }
        else if (strcmp(argv[i], "-R") == 0 || strcmp(argv[i], "--raw") == 0) {
            catch_missing_dependent_args(argv[i + 1], argv[i]);
            raw = argv[++i];    // tsv output file for raw counts
            raw_wanted = true;
            if (!util::path_exist(raw)){
                cerr << "Error: output folder does not exist: "<< raw << endl;
                return 1;
            }
        }
        else if (strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "--label") == 0) {
            catch_missing_dependent_args(argv[i + 1], argv[i]);
            c_nexus_wanted = true;
            groups = argv[++i];

            ifstream file_stream(groups);
            if (!file_stream.good()) { // catch unreadable file
                cout << "\33[2K\r" << "\u001b[31m" << "(ERR)" << " Could not read file " <<  "<" << groups << ">" << "\u001b[0m" << endl;
                file_stream.close();
                return 1;
            } else { file_stream.close();}

            // optional scnd argument for specified coloring
            if (argv[i+1] != NULL && string(argv[i+1]).rfind("-", 0) != 0){
                coloring = argv[++i];

                ifstream file_stream(coloring);
                if (!file_stream.good()) { // catch unreadable file
                    cout << "\33[2K\r" << "\u001b[31m" << "(ERR)" << " Could not read file " <<  "<" << coloring << ">" << "\u001b[0m" << endl;
                    file_stream.close();
                    return 1;
                } else { file_stream.close();}

            }
        }
        else if (strcmp(argv[i], "-k") == 0 || strcmp(argv[i], "--kmer") == 0) {            
            catch_missing_dependent_args(argv[i + 1], argv[i]);

            kmer = stoi(argv[++i]);    // Length of k-mers (default: 31, 10 for amino acids)
            userKmer = true;
        }
        else if (strcmp(argv[i], "-w") == 0 || strcmp(argv[i], "--window") == 0) {
            window = stoi(argv[++i]);    // Number of k-mers (default: 1)
            if (window > 1) {
                cerr << "Warning: using experimental feature --window" << endl;
            }
        }
        else if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--top") == 0) {
            catch_missing_dependent_args(argv[i + 1], argv[i]);
			if (strcmp(argv[i+1],"all") != 0){ // if user selects "-t all", do nothing
				catch_failed_stoi_cast(argv[i + 1], argv[i]);
				i++;
				string top_str = argv[i];
				top = stoi(top_str); // Number of splits (default: all)
				if (top_str[top_str.size() - 1] == 'n'){ // Dynamic split num (default: false)
					dyn_top = true;
				}
			}else{
				i++;				
			}
        }
        else if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--mean") == 0) {
            catch_missing_dependent_args(argv[i + 1], argv[i]);
            string arg = argv[++i];    // Mean weight function to handle asymmetric splits
            if (arg == "arith") {
                mean = util::arithmetic_mean;
            }
            else if (arg == "geom") {
                mean = util::geometric_mean;
            }
            else if (arg == "geom2") {
                mean = util::geometric_mean2;
            }
            else {
                cerr << "Error: unknown argument: --mean " << arg << endl;
                return 1;
            }
        }
        else if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--filter") == 0) {
            catch_missing_dependent_args(argv[i + 1], argv[i]);
            filter = argv[++i];    // Filter a greedy maximum weight subset
            if (filter == "strict" || filter == "tree" | filter == "planar") {
                // compatible to a tree
            }
            else if (filter == "weakly") {
                // weakly compatible network
            }
            else if (filter.find("-tree") != -1 && filter.substr(filter.find("-tree")) == "-tree") {
                for (const char &c: filter.substr(0, filter.find("-tree"))){
                    if (!isdigit(c)){
                        cerr << "Error: unexpected argument: --filter " << filter << ". Please specify n (Example usage: --filter 2-tree)" << endl;
                        return 1;
                    }
                }
                stoi(filter.substr(0, filter.find("-tree")));
            }
            else if (filter.find("tree") != -1 && filter.substr(filter.find("tree")) == "tree") {
                for (const char &c: filter.substr(0, filter.find("-tree"))){
                    if (!isdigit(c)){
                        cerr << "Error: unexpected argument: --filter " << filter << ". Please specify n (Example usage: --filter 2-tree)" << endl;
                        return 1;
                    }
                }
                stoi(filter.substr(0, filter.find("tree")));
            }
            else {
                cerr << "Error: unknown argument: --filter " << filter << endl;
                return 1;
            }
        }
        else if (strcmp(argv[i], "-x") == 0 || strcmp(argv[i], "--iupac") == 0) {
            catch_missing_dependent_args(argv[i + 1], argv[i]);
            catch_failed_stoi_cast(argv[i + 1], argv[i]);
            iupac = stoi(argv[++i]);    // Extended IUPAC alphabet, resolve ambiguous bases
        }
        else if (strcmp(argv[i], "-q") == 0 || strcmp(argv[i], "--qualify") == 0) {
            catch_missing_dependent_args(argv[i + 1], argv[i]);
            catch_failed_stoi_cast(argv[i + 1], argv[i]);
            quality = stoi(argv[++i]);    // Discard k-mers below a min. coverage threshold
        }
        else if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--pdf") == 0) {
            catch_missing_dependent_args(argv[i + 1], argv[i]);
            pdf = argv[++i];    // PDF output file
            pdf_wanted = true;    // Output of tree as pdf
            if (!util::path_exist(pdf)){
                cerr << "Error: output folder does not exist: "<< pdf << endl;
                return 1;
            }
        }
        else if (strcmp(argv[i], "-S") == 0 || strcmp(argv[i], "--svg") == 0) {
            catch_missing_dependent_args(argv[i + 1], argv[i]);
            svg = argv[++i];    // SVG output file
            svg_wanted = true;    // Output of tree as svg
            if (!util::path_exist(svg)){
                cerr << "Error: output folder does not exist: "<< svg << endl;
                return 1;
            }
        }
        else if (strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--norev") == 0) {
            reverse = false;    // Do not consider reverse complement k-mers
        }
        else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose = true;    // Print messages during execution
        }
        else if (strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "--amino") == 0) {
            amino = true;   // Input provides amino acid sequences
        }
        else if (strcmp(argv[i], "-M") == 0 || strcmp(argv[i], "--maxN") == 0) {
            check_n = true; // compare num (number of input genomes) to maxN (compile parameter DmaxN)
            if (i+1 < argc && argv[i+1][0] != '-'){
				path = argv[++i]; // path to makefile
			}
		}
        else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--code") == 0) {
            if (i+1 < argc) {
                string param = argv[++i];
                if (!util::is_number(param)){
                    i--;
                } else {
                    code = stoi(param);     // Number of the genetic code to be used
                }

            }
            shouldTranslate = true;
        }
        // parallelization 
        else if (strcmp(argv[i], "-T") == 0 || strcmp(argv[i], "--threads") == 0){ 
            catch_missing_dependent_args(argv[i + 1], argv[i]);
            catch_failed_stoi_cast(argv[i + 1], argv[i]);
            threads = stoi(argv[++i]); // Number of threads to create
            // Catch negative thread count
            if (threads <= 0) 
            {
                cerr << "Error: processing requires at least one thread" << threads << endl;
            }
            // Ask for user confirmation if the number of threads exceeds 64
            else if(threads > 64)
            {
                cout << "Are you sure you want to create " << threads << " threads?" << endl;
                if (!ask_for_user_confirmation()){return 0;}
            }
        }
        // bootsrapping
        else if (strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--bootstrapping") == 0 || strcmp(argv[i], "--bootstrap") == 0) {
            catch_missing_dependent_args(argv[i + 1], argv[i]);
            catch_failed_stoi_cast(argv[i + 1], argv[i]);
            bootstrap_no = stoi(argv[++i]);
			if (i+1 < argc && argv[i+1][0]!='-') { // optional filter threshold
				try {
					bootstrap_threshold = std::stof(argv[++i]);
					if (bootstrap_threshold<0 || bootstrap_threshold>1){
					 cerr << "Error: bootstrap filter threshold must be between 0 and 1." << endl;
					 exit(1);
					}
				} catch (const std::exception& e) {
					 cerr << "Error: Could not read bootstrap filter threshold: " << argv[i] << endl;
					 exit(1);
				}
                
            }

        }
        else if (strcmp(argv[i], "-C") == 0 || strcmp(argv[i], "--consensus") == 0) {
			if (i+1 < argc && argv[i+1][0]!='-') {
				consensus_filter = argv[++i];    // Filter a greedy maximum weight subset
				if (consensus_filter == "strict" || consensus_filter == "tree") {
					// compatible to a tree
				}
				else if (consensus_filter == "weakly") {
					// weakly compatible network
				}
				else if (consensus_filter.find("-tree") != -1 && consensus_filter.substr(consensus_filter.find("-tree")) == "-tree") {
					for (const char &c: consensus_filter.substr(0, consensus_filter.find("-tree"))){
						if (!isdigit(c)){
							cerr << "Error: unexpected argument: --consensus_filter " << consensus_filter << ". Please specify n (Example usage: --consensus 2-tree)" << endl;
							return 1;
						}
					}
					stoi(consensus_filter.substr(0, consensus_filter.find("-tree")));
				}
				else if (consensus_filter.find("tree") != -1 && consensus_filter.substr(consensus_filter.find("tree")) == "tree") {
					for (const char &c: consensus_filter.substr(0, consensus_filter.find("-tree"))){
						if (!isdigit(c)){
							cerr << "Error: unexpected argument: --consensus " << consensus_filter << ". Please specify n (Example usage: --consensus 2-tree)" << endl;
							return 1;
						}
					}
					stoi(consensus_filter.substr(0, consensus_filter.find("tree")));
				}
				else {
					cerr << "Error: unknown argument: --consensus " << consensus_filter << endl;
					return 1;
				}
			} else {
				consensus_filter="SameAsFilter"; // set to the same as --filter later because --filter might not be set yet.
			}

        }
        else {
            cerr << "Error: unknown argument: " << argv[i] <<  "\t type --help" << endl;
            return 1;
        }
    }


    
    /**
     *   [version check]
     *   request the current SANS version from gitlab and check if this version is up to date (requires wget)
     */

    /* Turn off for CloWM
      if (verbose){cout << "Checking for updates" << endl;}
    bool version_checked = false;
    if (!system("wget --timeout=1 --tries=1 -qO- https://gitlab.ub.uni-bielefeld.de/gi/sans/raw/master/src/main.h | grep -q SANS_VERSION")){
        version_checked = true;
        if (system("wget --timeout=1 --tries=1 -qO- https://gitlab.ub.uni-bielefeld.de/gi/sans/raw/master/src/main.h | grep -q " SANS_VERSION)) {
        cout << "NEW VERSION AVAILABLE: https://gitlab.ub.uni-bielefeld.de/gi/sans" << endl;
        }
        else if(verbose){cout << "Version up to date" << endl;}
    }
    if (!version_checked && verbose) {cout << "Could not fetch version information" << endl;}
   */
    
    // set consensus filter to default (same as --filter) if necessary
	if (consensus_filter=="SameAsFilter"){
		consensus_filter=filter; 
	}


    /**
     *  [restriction check]
     * - Check if the given argument configuration does violate any run restrictions
     */ 
    if (input.empty() && !splits.empty()) // Input file is required for correct split loading
    {
        cerr << "Error: missing argument --input <file_name> is required for split loading" << endl;
        return 1;
    }
    
    if (input.empty() && graph.empty()) {
        cerr << "Error: missing argument: --input <file_name> or --graph <file_name>" << endl;
        return 1;
    }

    if (!input.empty() && !graph.empty() && !splits.empty()) {
        cerr << "Error: too many input arguments: --input, --graph, and --splits" << endl;
        return 1;
    }
    if (!graph.empty() && !splits.empty()) { // ---- Why not?
        cerr << "Error: too many input arguments: --graph and --splits" << endl;
        return 1;
    }

    if (input.empty() && amino) {
        cerr << "Error: missing argument: --input <file_name> for option --amino" << endl;
        return 1;
    }
    if (!splits.empty() && amino) {
        cerr << "Error: too many input arguments: --splits and --amino" << endl;
        return 1;
    }
    if (!graph.empty() && amino) {
        cerr << "Error: too many input arguments: --graph and --amino" << endl;
        return 1;
    }

    if (output.empty() && newick.empty() && nexus.empty() && pdf.empty() && svg.empty() && core.empty() && !raw_wanted) {
        cerr << "Error: missing argument: --output <file_name> or --newick <file_name> or --nexus <file_name> or --pdf <file_name> or --svg <file_name> or --core <file_name> or --raw <file_name>" << endl;
        return 1;
    }
	if (output.empty() && newick.empty() && nexus.empty() && pdf.empty() && svg.empty() && !core.empty()) {
		if(!filter.empty() || !consensus_filter.empty() || bootstrap_no>0 || mean != util::geometric_mean2 || top!=-1 ){
			cerr << "Warning: No output option for a phylogeny given. Only core k-mers are computed. Some given arguments only make sense for phylogeny construction and are redundant." << endl;
		}
    }
	if (!core.empty() && !splits.empty()) {
		cerr << "Error: From splits as input, no core k-mers can be determined." << endl;
		return 1;
    }
	if (!blacklist.empty() && input.empty() && graph.empty()) {
		cerr << "Error: Blacklist can only be applied when reading sequences as input, i.e. -i or -g." << endl;
		return 1;
    }
    if (kmer > maxK && splits.empty()) {
        cerr << "Error: k-mer length exceeds -DmaxK=" << maxK << endl;
        cerr << "Solution: Modify -DmaxK in makefile, run make, run SANS." << endl;
        return 1;
    }
    if (!newick.empty() && filter != "strict" && filter.find("tree") == -1 && consensus_filter.empty()) {
        cerr << "Error: Newick output only applicable in combination with -f strict or n-tree." << endl;
        return 1;
    }
    if (!newick.empty() && !consensus_filter.empty() && consensus_filter != "strict" && consensus_filter.find("tree") == -1) {
        cerr << "Error: Newick output only applicable in combination with -C strict or -C n-tree." << endl;
        return 1;
    }
    if (c_nexus_wanted && !nexus_wanted && !pdf_wanted && !svg_wanted){
        cerr << "Error: Labeled (colored) nexus output only applicable in combination with --nexus <filename> or --pdf <filename> or --svg <filename>." << endl;
        return 1;
    }
    if (c_nexus_wanted || pdf_wanted){
        const string splitstree_path = "SplitsTree";
        if (!nexus_color::program_in_path(splitstree_path)) {
            cerr << splitstree_path << " is not in the PATH." << endl;
        }
    }
    if (amino && quality > 1) {
        cerr << "Error: using --qualify with --amino is (currently) not supported" << endl;
        cerr << "       Please send us a message if you have this special use case" << endl;
        return 1;
    }
    if (!graph.empty()) {
        if (quality > 1) {
            cerr << "Warning: input from graph with --qualify can produce unexpected results" << endl;
        }
        else if (window > 1) {
            cerr << "Warning: input from graph with --window can produce unexpected results" << endl;
        }
    } else if (quality > 1 && window > 1) {
        cerr << "Warning: using --qualify with --window could produce unexpected results" << endl;
    }

    if (!input.empty() && !splits.empty()) {
        cerr << "Note: two input arguments --input and --splits were provided" << endl;
        cerr << "      --input is used for lookup only, no additional splits are inferred" << endl;
    }
    if (input.empty() && !splits.empty() && !newick.empty()) {
        cerr << "Note: Newick output from a list of splits, some taxa could be missing" << endl;
        cerr << "      --input can be used to provide the original list of taxa" << endl;
    }
    if (input.empty() && graph.empty() && bootstrap_no>0){
        cerr << "Error: Bootstrapping can only be applied with given sequence data (--input or --graph)" << endl;
		return 1;
	}
    if (bootstrap_no>0 && filter.empty()){
        cerr << "Error: Bootstrapping can only be applied when a filter is selected (--filter)" << endl;
		return 1;
	}
	if (bootstrap_no==0 && !consensus_filter.empty()){
        cerr << "Error: Filter on bootstrap values (--consensus) can only be chosen in combination with bootstrapping (--boostrapping)" << endl;
		return 1;
	}


    /*[processing setup]
    *
    * prepare the processing based on the input
    */

    // enable translation
    if (shouldTranslate) {
        amino = true;
        if (!translator::init(code)) {
            cerr << "Error: No translation data found" << translate << endl;
        }
    }
    // deduct default kmer size if not user defined
    if (!userKmer) {kmer = amino == true ? 10 : 31;}


    /**
     *  [indexing]
     * - Collect all target file names and locations
     * - Check whether the given files exist or not
     */ 

    // determine the folder the list is contained in
    string folder="";
	uint64_t found=input.find_last_of("/\\");
	if (found!=string::npos)
    { 
        folder=input.substr(0,found+1);
    }

    // parse the list of input sequence files
    hash_map<string, uint64_t> name_table; // the name to color map
    unordered_set<string> graph_filter_names; // storing the colors to filter a graph
    vector<string> denom_names; // storing the representative name per color
    vector<vector<string>> gen_files; // genome file collection
    vector<int> q_table; // q value (k-mer occurrence threshold) per genome/color

    
    if (!input.empty()) {
        // check the input file 
        ifstream file(input);
        if (!file.good()) {
            cerr << "Error: could not read input file: " << input << endl;
            return 1;
        }


        // parse the list of input sequence files
        string line; // the iterated input line
        string file_name; // the current file name
        bool is_first; // indicating the first filename of a line (For file list)
        bool has_files; // indicating if a line contains filenames
		int min_q=quality;
		int max_q=quality;
        
        getline(file, line);
        // check the file format
        std::smatch matches;

        // parse file of files
        while(true){
			std::regex_search(line, matches, std::regex("(:)"));
			int q=quality;
            vector<string> target_files; // container of the current target files
            if (!matches.empty()){ // parse kmt format
                // ensure the terminal signs " !" exists.
                if (line.find_first_of('!') != line.npos){
					int p=line.find_first_of('!');
					q = stoi(line.substr(p+1,line.length()));
					line = line.substr(0, line.find_first_of('!') + 1); // cut off tail
				}
                else if (line.back() == ' '){line += '!';} // append terminal sign if missing
                else {line += " !";} // append both terminal signs if missing

                string denom = line.substr(0, line.find_first_of(" ")); // get the dataset-id
                denom_names.push_back(denom); // add id to denominators
				name_table[denom] = num;
				
                line = line.substr(line.find_first_of(":") + 2, line.npos); // cut off the dataset-id

                std::smatch matches; // Match files
                while (std::regex_search(line, matches, std::regex("[ ; ]|[ !]"))){
                    file_name = matches.prefix().str(); // get filename from match
                    line=matches.suffix().str(); // update the line
                    if (file_name.length() == 0){continue;} // skip empty file name
                    else {target_files.push_back(file_name); name_table[file_name] = num;} // add the file name to target files and name table
				}
                num ++;
            }

            else{ // parse file list format
                is_first = true;
                has_files = false;
                string file_name = "";
                size_t it = 0;
                size_t line_length = line.length();
                for (auto x: line){ // iterate the line
                    it ++;
                    if (x == ' ' | it == line_length){ // checkout the file name if a space occurs or the line ends
                        if (it == line_length){file_name += x;} // add the last character to the last file name
                        if (file_name.length() == 0){file_name = ""; continue;} // skip continuous spaces
                        if (is_first){ // use first file name as denom name
                            has_files = true;
                            denom_names.push_back(file_name); // set denom name
                            is_first = false;
                        }
                        target_files.push_back(file_name); // add the file_name to the genome file vector
                        name_table[file_name] = num; // add the file tp the name_table
                        file_name = "";
                    }
                    else{file_name += x;}
                }
                if (has_files) {num++;}
            }

            q_table.push_back(q);
			if (q>max_q){max_q=q;}
			if (q<min_q){min_q=q;}

            // check files
			if (splits.empty()){
					for(string file_name: target_files){
						if(file_name[0]!='/'){ //no absolute path?
							file_name=folder+file_name;
						}
						ifstream file_stream = ifstream(file_name);
						if (!file_stream.good()) { // catch unreadable file
							cout << "\33[2K\r" << "\u001b[31m" << "(ERR)" << " Could not read file " <<  "<" << file_name << ">" << "\u001b[0m" << endl;
							file_stream.close();
							return 1;
						}
						else{ file_stream.close();}	
					}
			}
			
            gen_files.push_back(target_files); // add the files of the current genome to the genome collection
            if (!getline(file, line)) {break;}
        }
        
        quality=max_q;
        if(max_q==min_q){q_table.clear();} // all q_values the same (=quality)
    }
    int denom_file_count = denom_names.size();

	


    /**
     * --- [indexing CDBG input] --- 
     * - Collect all target sequence names from the Bifrost-CDBG
     */ 

#ifdef useBF
    // load an existing Bifrost graph
    ColoredCDBG<> cdbg;
    if (!graph.empty()) {
        if (cdbg.read(graph + ".gfa", graph + ".color.bfg", threads, verbose)) { // Allow parallel reading with new t parameter.

        } else {
            cerr << "Error: could not load Bifrost graph" << endl;
            return 1;
        }

        if (kmer != cdbg.getK()) {
		kmer = cdbg.getK();
	 	if (kmer > maxK) {
                cerr << "Error: k-mer length exceeds -DmaxK=" << maxK << endl;
		cerr << "Solution: modify -DmaxK in makefile, run make, run SANS." << maxK << endl;
                return 1;
            }
            cerr << "Warning: setting k-mer length to match the given Bifrost graph. New lenght: " << kmer << endl;
		}

        vector<string> cdbg_names = cdbg.getColorNames(); // color names of the cdbg compacted genomes.

		// read filter cdbg
		if (!graph_filter.empty()) {
			// check the input file 
			ifstream file(graph_filter);
			if (!file.good()) {
				cerr << "Error: could not read graph filter file: " << graph_filter << endl;
				return 1;
			}

			// parse the list of input sequence files
			string line; // the iterated input line
			string file_name; // the current file name
			bool is_first; // indicating the first filename of a line (For file list)
			bool has_files; // indicating if a line contains filenames
			
			getline(file, line);
			// check the file format
			std::smatch matches;

			// parse file of files
			while(true){
				std::regex_search(line, matches, std::regex("(:)"));
				vector<string> target_files; // container of the current target files
				if (!matches.empty()){ // parse kmt format
					// ensure the terminal signs " !" exists.
					if (line.find_first_of('!') != line.npos){
						int p=line.find_first_of('!');
						line = line.substr(0, line.find_first_of('!') + 1); // cut off tail
					}
					else if (line.back() == ' '){line += '!';} // append terminal sign if missing
					else {line += " !";} // append both terminal signs if missing

					string denom = line.substr(0, line.find_first_of(" ")); // get the dataset-id
					if (find(denom_names.begin(),denom_names.end(),denom) == denom_names.end()){
						denom_names.push_back(denom); // add id to denominators
						name_table[denom] = num;
					}
					else
					{
						cout << "Warning: " << denom << " exists in input and graph. It is treated as one sequence" << endl;
					}

					line = line.substr(line.find_first_of(":") + 2, line.npos); // cut off the dataset-id

					std::smatch matches; // Match files
					while (std::regex_search(line, matches, std::regex("[ ; ]|[ !]"))){
						file_name = matches.prefix().str(); // get filename from match
						line=matches.suffix().str(); // update the line
						if (file_name.length() == 0){continue;} // skip empty file name
						else {
							if (name_table.find(file_name) == name_table.end()){
								name_table[file_name] = num;
								target_files.push_back(file_name);
							}
							else
							{
								// compatible mapping?
 								if(denom_names[name_table[file_name]] != denom){
									cerr << "Error: Incompatible mapping of " << file_name << " in graph and sequence kmtricks file." << endl;
									exit(1);
								}
							}
						} // add the file name to target files and name table
						if (find(cdbg_names.begin(),cdbg_names.end(),file_name)==cdbg_names.end()){
								cerr << "Error: " << file_name << " in graph filter file but not in graph." << endl;
								exit(1);
						}
						graph_filter_names.insert(file_name);
					}
					num ++;
				}

				else{ // parse file list format
					is_first = true;
					has_files = false;
					string file_name = "";
					string denom="";
					size_t it = 0;
					size_t line_length = line.length();
					for (auto x: line){ // iterate the line
						it ++;
						if (x == ' ' | it == line_length){ // checkout the file name if a space occurs or the line ends
							if (it == line_length){file_name += x;} // add the last character to the last file name
							if (file_name.length() == 0){file_name = ""; continue;} // skip continuous spaces
							if (is_first){ // use first file name as denom name
								has_files = true;
								denom = file_name;
								if (find(denom_names.begin(),denom_names.end(),denom)==denom_names.end()){
									denom_names.push_back(denom); // set denom name
								} else
								{
									cout << "Warning: " << denom << " exists in input and graph. It is treated as one sequence" << endl;
								}
								is_first = false;
							}
							if (name_table.find(file_name) == name_table.end()){
								name_table[file_name] = num;
								target_files.push_back(file_name);
							}
							else
							{
								// compatible mapping?
 								if(denom_names[name_table[file_name]] != denom){
									cerr << "Error: Incompatible mapping of " << file_name << " in graph and sequence kmtricks file." << endl;
									exit(1);
								}
							}
							if (find(cdbg_names.begin(),cdbg_names.end(),file_name)==cdbg_names.end()){
									cerr << "Error: " << file_name << " in graph filter file but not in graph." << endl;
									exit(1);
							}
							graph_filter_names.insert(file_name);
							file_name = "";
						}
						else{file_name += x;}
					}
					if (has_files) {num++;}
				}

				gen_files.push_back(target_files); // add the files of the current genome to the genome collection
				if (!getline(file, line)) {break;}
			}// end reading filter
			
		} else { //no filter -> read all colors separately
			
		
			for (auto &col_name: cdbg_names){ // iterate the cdbg names and transcribe them to the name table
				if (name_table.find(col_name) == name_table.end()){
					name_table[col_name] = num++;
					denom_names.push_back(col_name);
					vector<string> dummy;	
					gen_files.push_back(dummy);
				}
				else
				{
					cout << "Warning: " << col_name << " exists in input and graph. It is treated as one sequence" << endl;
				}
			}
		}
		
        if (verbose) {
            cout << endl;
		}
	


    }

#endif


    /**
     *   [post indexing check]
     * - Update and check validity of input dependent meta variables
     */ 

    // check if the number of genomes is reasonably close the maximal storable color set
    if (check_n) {
       util::check_n(num,path,maxN);
	}

    // check if the number of genomes exceeds the maximal storable color set
    if (num > maxN) {
        cerr << "Error: number of input genomes ("<<num<<") exceeds -DmaxN=" << maxN << endl;
        cerr << "Solution: modify -DmaxN in makefile, run make, run SANS; or use SANS-autoN.sh." << endl;
        return 1;
    }
    if (maxN-num>=100) {
		cout << "Warning: number of input genomes ("<<num<<") much lower than -DmaxN=" << maxN << endl;
		cout << "Recommendation: modify -DmaxN in makefile, run make, run SANS; or use SANS-autoN.sh." << endl;
	}


    // Set dynamic top by filenum
    if (dyn_top){
        top = top * num;
    }
    if(verbose && top>-1){
	cout<<"Restricting output to "<<top<<" splits."<< endl;
    }

    /**
     * [input processing]
     * - transcribe each given input to the graph
     */ 

    chrono::high_resolution_clock::time_point begin = chrono::high_resolution_clock::now();    // time measurement

    kmer::init(kmer);      // initialize the k-mer length
    kmerAmino::init(kmer); // initialize the k-mer length
    color::init(num);    // initialize the color number
    graph::init(top, amino, q_table, quality, blacklist, blacklist_amino, threads); // initialize the toplist size and the allowed characters

	
	/**
	 * Read blacklist
	 */
	if(!blacklistfile.empty()){

        if (verbose) {
            cout << "Reading blacklist file... " << flush;
        }
 
        string sequence;    // read in the sequence files and extract the k-mers
		char c_name[(blacklistfile).length()+1]; // Create char array for c compatibilty
		strcpy(c_name, (blacklistfile).c_str()); // Transcire to char array

		igzstream file(c_name, ios::in);    // input file stream
				count::deleteCount();

				string appendixChars; 
				string line;    // read the file line by line
				while (getline(file, line)) {
					if (line.length() > 0) {
						if (line[0] == '>' || line[0] == '@') {    // FASTA & FASTQ header -> process
							graph::fill_blacklist(sequence, reverse);
							sequence.clear();
						}
						else if (line[0] == '+') {    // FASTQ quality values -> ignore
							getline(file, line);
						}
						else {
							transform(line.begin(), line.end(), line.begin(), ::toupper);
							string newLine = line;
							if (shouldTranslate) {
								if (appendixChars.length() >0 ) {
									newLine= appendixChars + newLine;
									appendixChars = "";
								}
								auto toManyChars = line.length() % 3;
								if (toManyChars > 0) {
									appendixChars = newLine.substr(line.length() - toManyChars, toManyChars);
									newLine = newLine.substr(0, line.length() - toManyChars);
								}
								newLine = translator::translate(newLine);
							}
							sequence += newLine;    // FASTA & FASTQ sequence -> read
						}
					}
				}
				if (verbose && count::getCount() > 0) {
					cerr << count::getCount()<< " triplets could not be translated while reading blacklist."<< endl;
				}
				graph::fill_blacklist(sequence, reverse);
				sequence.clear();

				file.close();
       if (verbose) {
            cout << graph::size_blacklist() << " k-mers read." << endl << flush;
        }
        if (graph::size_blacklist()==0){
			cerr << "Warning: Blacklist provided, but no k-mers read." << endl << flush;
		}
		graph::activate_blacklist();
	}


	
    /**
     *  ---> Split processing
     *  - transcibe all splits from the input split file
     */

    // iterate splits and add them to the toplist
    if (!splits.empty()) {
    ifstream file(splits);
    if (!file.good()) { // check if the target splits file exists
        cerr << "Error: could not read splits file: " << splits << endl;
        return 1;
    }
    string line;
    while (getline(file, line)) { // Iterate each split
        uint64_t curr = line.find('\t');
        double weight = stod(line.substr(0, curr));
        uint64_t next = curr + 1;

        color_t color = 0;
        do {
            curr = line.find('\t', next);
            string name = line.substr(next, curr-next);
            if (name_table.find(name) == name_table.end()) { // check if the splits genome names are already indexed
                cerr << "Error: unlisted file " << name << " in split file" << endl;
                return 1; 
            }
            color.set(name_table[name]);
            next = curr + 1;
        } while (curr != string::npos);

        graph::add_split(weight, color); // add the split to the graph
    }
    file.close();
    }

    
    
    
    /**
     * ---> sequence processing ---
     * - translate all given sequence k-mers
     * - transcribe all given sequence k-mers to the graph
     */ 
    

    if (!input.empty() && splits.empty()) {
        if (verbose) {
            cout << "Reading input files..." << endl << flush;
        }

        // Thread safe implementation of getting the index of the next input to preocess
        uint64_t index = 0;
        std::mutex index_mutex;
        auto index_lambda = [&] () { std::lock_guard<mutex> lg(index_mutex); return index++;};

        auto lambda = [&] (uint64_t T, vector<uint16_t> genome_ids, vector<uint16_t> file_ids){ // This lambda expression wraps the sequence-kmer hashing
            string sequence;    // read in the sequence files and extract the k-mers
            uint64_t i = index_lambda();
            while (i < genome_ids.size()) {
				string file_name = gen_files[genome_ids[i]][file_ids[i]]; // the filenames corresponding to the target  
				if(file_name[0]!='/'){ //no absolute path?
					file_name=folder+file_name;
				}

				char c_name[(file_name).length()+1]; // Create char array for c compatibilty
				strcpy(c_name, (file_name).c_str()); // Transcire to char array

				igzstream file(c_name, ios::in);    // input file stream
				if (verbose) {     // print progress
// 					cout << "\33[2K\r" << file_name;
					if (q_table.size()>0) {
						cout <<" q="<<q_table[genome_ids[i]];
					}
					cout << " (genome " << genome_ids[i]+1 << "/" << denom_file_count;
					if(genome_ids.size()>gen_files.size()){
						cout << "; file " << i+1 << "/" << genome_ids.size();
					}
					cout << ")" << endl;
				}
				count::deleteCount();

				string appendixChars; 
				string line;    // read the file line by line
				while (getline(file, line)) {
					if (line.length() > 0) {
						if (line[0] == '>' || line[0] == '@') {    // FASTA & FASTQ header -> process
							if (window > 1) {
								iupac > 1 ? graph::add_minimizers(T, sequence, genome_ids[i], reverse, window, iupac)
										: graph::add_minimizers(T, sequence, genome_ids[i], reverse, window);
							} else {
								iupac > 1 ? graph::add_kmers(T, sequence, genome_ids[i], reverse, iupac)
										: graph::add_kmers(T, sequence, genome_ids[i], reverse);
							}

							sequence.clear();

//                            if (verbose) {
//                                cout << "\33[2K\r" << line << flush << endl;    // print progress
//                            }
						}
						else if (line[0] == '+') {    // FASTQ quality values -> ignore
							getline(file, line);
						}
						else {
							transform(line.begin(), line.end(), line.begin(), ::toupper);
							string newLine = line;
							if (shouldTranslate) {
								if (appendixChars.length() >0 ) {
									newLine= appendixChars + newLine;
									appendixChars = "";
								}
								auto toManyChars = line.length() % 3;
								if (toManyChars > 0) {
									appendixChars = newLine.substr(line.length() - toManyChars, toManyChars);
									newLine = newLine.substr(0, line.length() - toManyChars);
								}

								newLine = translator::translate(newLine);
							}
							sequence += newLine;    // FASTA & FASTQ sequence -> read
						}
					}
				}
				if (verbose && count::getCount() > 0) {
					cerr << count::getCount()<< " triplets could not be translated."<< endl;
				}
				if (window > 1) {
					iupac > 1 ? graph::add_minimizers(T, sequence, genome_ids[i], reverse, window, iupac)
							: graph::add_minimizers(T, sequence, genome_ids[i], reverse, window);
				} else {
					iupac > 1 ? graph::add_kmers(T, sequence, genome_ids[i], reverse, iupac)
							: graph::add_kmers(T, sequence, genome_ids[i], reverse);
				}
				sequence.clear();

				
// 				if (verbose) {
// 					cout << "\33[2K\r" << flush;
// 				}
				file.close();
                graph::clear_thread(T);
                i = index_lambda();
            }
        }; // End of lambda expression

        // Driver code for multithreaded kmer hashing
		vector<uint16_t> genome_ids; //unfold multiple files per genome to two flat lists, one listing the genome ids and one listing the file ids.
		vector<uint16_t> file_ids;
		for (int g=0;g<gen_files.size();g++){
			for (int f=0;f<gen_files[g].size();f++){
				genome_ids.push_back(g);
				file_ids.push_back(f);
			}
		}
		vector<thread> thread_holder(threads);
        for (uint64_t thread_id = 0; thread_id < threads; ++thread_id){thread_holder[thread_id] = thread(lambda, thread_id, genome_ids, file_ids);}
        for (uint64_t thread_id = 0; thread_id < threads; ++thread_id){thread_holder[thread_id].join();}
        

    }

    /**
     * ---> bifrost CDBG processing ---
     * - iterate all colored k-mers from a CDBG
     * - compute the splits created by the CDBG k-mers given the graphs colore k-mer collection
     * (has to be executed after sequence processing)
     * // Todo: Bug
     */ 

double min_value = numeric_limits<double>::min(); // current minimal weight represented in the top list
#ifdef useBF
     if (!graph.empty()) {
		if (verbose){
            cout << "SANS::main(): Processing unitigs..." << endl;
		}
		
		uint64_t cur = 0, prog = 0, next;
		uint64_t max = cdbg.size();

		for (auto& unitig : cdbg) {
			if (verbose) {
				next = 100 * cur / max;
				if (prog < next)  cout << "\33[2K\r" << "Processed " << cur << " unitigs (" << next << "%) " << flush;
				prog = next; cur++;
			}
			auto sequence = unitig.mappedSequenceToString();
			auto* matrix = unitig.getData()->getUnitigColors(unitig);

			for (auto it = matrix->begin(unitig); it != matrix->end(); ++it) {
				auto kmer_sequence = sequence.substr(it.getKmerPosition(), kmer);
				auto c = cdbg.getColorName(it.getColorID());
				if(graph_filter.empty() || graph_filter_names.find(c)!=graph_filter_names.end()){ //
					const uint16_t& color = name_table[c]; // set the k-mer color
					graph::add_cdbg_colored_kmer(kmer_sequence, color);
				}
			}
		}
		if (verbose) {
			end = chrono::high_resolution_clock::now(); 
			cout<< "\33[2K\r" << "Processed " << max << " unitigs (100%)" << " (" << util::format_time(end - begin) << ")" << endl << flush;
		}
		denom_file_count = denom_names.size();
	}


#endif
       if(splits.empty() && (graph::number_singleton_kmers()+graph::number_kmers()==0)){
		cout << "no k-mers found." << endl;
	       exit(0);
       }
        
	if(verbose & ((!input.empty() && splits.empty()) || !graph.empty())){
		uint64_t s=graph::number_singleton_kmers();
		uint64_t all=s+graph::number_kmers();
		end = chrono::high_resolution_clock::now(); 
		cout << all << " k-mers read." << flush;
		cout << " (" << s << " / "<< (100*s/all) <<"% singleton k-mers)" << " (" << util::format_time(end - begin) << ")" << endl << flush;
	}

	///DEBUGING////
	// 	cout << "\nname_table:" << endl;
	// 	 for (const auto& pair : name_table) {
	//         cout << pair.first << " : " << pair.second << endl;
	//     }
	// 	
	// 	cout << "\ndenom_names:" << endl;
	//     for (const std::string& item : denom_names) {
	//         std::cout << item << std::endl;
	//     }
	// 	
	// 	cout << "\ngen_files:" << endl;
	//     for (const auto& row : gen_files) {
	//         for (const auto& item : row) {
	//             std::cout << item << " ";
	//         }
	//         std::cout << std::endl; // Newline after each row
	//     }



	/*
	 * [core k-mers]
	 */
	if(!core.empty()){
		// output file stream
		ofstream core_file(core);
		ostream core_stream(core_file.rdbuf());
		graph::output_core(core_stream,verbose);
		if(verbose){
			cout << " (" << util::format_time(end - begin) << ")" << endl << flush;
		}
	}
	

	// if only core-kmers are asked for, no further processing necessary
	if (!output.empty() || !newick.empty() || !nexus.empty() || !pdf.empty() || !svg.empty() || raw_wanted){ 
	
		/*
		* [graph processing]
		* - collect all colors and kmers into the color table
		* - weight the colors based on the weight function and the kmers that support them
		*/

		// function to map color position to file name
		std::function<string(const uint64_t&)> map=[=](uint64_t i) {
			if (i < denom_names.size()) return denom_names[i];
			cerr << "Error: color bit does not correspond to color name" << endl;
			exit(EXIT_FAILURE);
		};

		if (verbose) {
			cout << "Accumulating splits from non-singleton k-mers..."  << flush;
		}
		graph::add_weights(mean, min_value, verbose);  // accumulate split weights
		if (verbose) {
			end = chrono::high_resolution_clock::now();
			cout << "\33[2K\r" << "Accumulating splits from non-singleton k-mers... (" << util::format_time(end - begin) << ")" << endl;
			cout << "Accumulating splits from singleton k-mers..."  << flush;
		}
		graph::add_singleton_weights(mean, min_value, verbose);  // accumulate split weights
		if (verbose) {
			end = chrono::high_resolution_clock::now();
			cout << "\33[2K\r"  << "Accumulating splits from singleton k-mers... (" << util::format_time(end - begin) << ")" << endl;
		}


		
		

		/* [split processing]
		* - compute the splits
		*/ 

	
		if (verbose) {
			cout << "Compile split list..."  << flush;
		}
		graph::compile_split_list(mean, min_value);
		if (verbose) {
			end = chrono::high_resolution_clock::now();
			cout << " (" << util::format_time(end - begin) << ")" << endl;
		}


		/*
		* [bootstrap handling]
		*/


		if(bootstrap_no==0){ // if bootstrapping -> no initial filtering
			
			// NO BOOTSTRAPPING
				
			if(verbose){
				cout << "Filtering splits..." << flush;
			}
			apply_filter(filter,newick, map, graph::split_list,verbose, num);
			if (verbose) {
				end = chrono::high_resolution_clock::now();
				cout << "\33[2K\r" << "Filtering splits... (" << util::format_time(end - begin) << ")" << endl;
			}

		}else{
			
		// BOOTSTRAPPING
			

			bool verbose_orig=verbose;
			if(verbose){
// 				cout << "\n" << flush;
				verbose=false; // switch off output of filtering
			}
			
			// init bootstrap support value counting
			// remember original split weights for later
			hash_map<color_t, double> orig_weights;
			for (auto& it : graph::split_list){
				support_values.insert({it.second,0});
				orig_weights.insert({it.second,it.first});
			}
			


		// Thread safe implementation of getting the index of the next run
			uint64_t index = 0;
			std::mutex index_mutex;
			std::mutex count_mutex;
			auto index_lambda_bootstrap = [&] () {std::lock_guard<mutex> lg(index_mutex); return index++;};
			
			auto lambda_bootstrap_count = [&] (multimap_<double, color_t> split_list_bs, hash_map<color_t, uint32_t>& support_values) { std::lock_guard<mutex> lg(count_mutex); for (auto& it : split_list_bs){color_t colors = it.second;support_values[colors]++;}};

			auto lambda_bootstrap = [&] (uint64_t T, uint64_t max){ // This lambda expression wraps the sequence-kmer hashing
				uint64_t i = index_lambda_bootstrap();
				while (i<max){

					if (verbose_orig) {
						cout << "\33[2K\r" << "Bootstrapping... ("<<(i+1)<<"/"<<max<<")" << flush;
					}
					
					// create bootstrap replicate
					multimap_<double, color_t>  split_list_bs = graph::bootstrap(mean);
					apply_filter(filter,"", map, split_list_bs,verbose, num);
					
					// count conserved splits
					lambda_bootstrap_count(split_list_bs, support_values);
					
					// more to do for this thread?
					i = index_lambda_bootstrap();
				}
			};
			
			// Driver code for multithreaded bootstrapping
			vector<thread> thread_holder(threads);
			for (uint64_t thread_id = 0; thread_id < threads; ++thread_id){thread_holder[thread_id] = thread(lambda_bootstrap, thread_id, bootstrap_no);}
			for (uint64_t thread_id = 0; thread_id < threads; ++thread_id){thread_holder[thread_id].join();}

			
			verbose=verbose_orig; //switch back to verbose if originally set
			
			if (verbose) {
				end = chrono::high_resolution_clock::now();
				cout << "\33[2K\r" << "Bootstrapping... (" << util::format_time(end - begin) << ")" << endl << flush;
 				cout << "Filtering splits... "<< flush;
			}
			
			
			if(bootstrap_threshold>0){
				//erase low support splits
				auto it = graph::split_list.begin();
				while (it != graph::split_list.end()) {
					double conf=(1.0*support_values[it->second])/bootstrap_no;
					if (conf<bootstrap_threshold) {
						it = graph::split_list.erase(it);
					} else {
						++it;
					}
				}
			}

			if(consensus_filter.empty()) {
				// filter original splits by weight
				apply_filter(filter,newick, map, graph::split_list,&support_values,bootstrap_no,verbose, num);
			}else{
				// filter original splits by bootstrap value
				// compose a corresponding split list
				multimap_<double, color_t> split_list_conf;
				int max_weight = -1;
				for (auto& it : graph::split_list){
					int weight = it.first;
					color_t colors = it.second;
					if (max_weight == -1){ max_weight = weight; }
					double conf=(1.0*support_values[it.second]*max_weight)+weight;
					graph::add_split(conf,colors,split_list_conf);
				}

				//filter
	// 			apply_filter(consensus_filter,newick, map, split_list_conf,verbose);
				apply_filter(consensus_filter,newick, map, split_list_conf,&support_values,bootstrap_no,verbose, num);

				//apply result to original split list
				graph::split_list.clear();
				for (auto& it : split_list_conf){
					color_t colors = it.second;
					graph::add_split(orig_weights[colors],colors);
				}
			}
			if (verbose) {
				end = chrono::high_resolution_clock::now();
				cout << "\33[2K\r" << "Filtering splits... (" << util::format_time(end - begin) << ")" << endl << flush;
			}

		}


		/**
		* [write to output]
		*/ 

		if (verbose) {
			cout << "Writing output..." << endl << flush;
		}

		ofstream file;    // output file stream
		ostream stream(file.rdbuf());
		if (!output.empty()){
			file.open(output);
		}

		ofstream file_bootstrap;
		ostream stream_bootstrap(file_bootstrap.rdbuf());
		if (bootstrap_no>0){
			file_bootstrap.open(output+".bootstrap");    // output file stream
		}

		ofstream file_nexus; // output for nexus file
		ostream stream_nexus(file_nexus.rdbuf());

		ofstream file_raw;    // output for raw counts
		ostream stream_raw(file_raw.rdbuf());
		if(!raw.empty()){
			file_raw.open(raw);
		}
		
		if(nexus_wanted || pdf_wanted || svg_wanted){
			if(nexus.empty()){ // temporarily name nexus file to create pdf with it
				if (!pdf.empty()){
					nexus = pdf + ".nex";
				}else {
					nexus = svg + ".nex";
				}
				
			}
			file_nexus.open(nexus);
			// nexus format stuff
			stream_nexus << "#nexus\n\nBEGIN Taxa;\nDIMENSIONS ntax=" << denom_file_count << ";\nTAXLABELS" ;
			for(int i = 0; i < denom_file_count; ++i){
				string taxa = nexus_color::remove_extensions(denom_names[i]); // cutting off file extension
				stream_nexus << "\n[" << i+1 << "] '" << taxa << "'";
			}
			stream_nexus << "\n;\nEND; [TAXA]\n";
			if(bootstrap_no>0){ // confidence values
				stream_nexus << "\nBEGIN Splits;\nDIMENSIONS ntax=" << denom_file_count << " nsplits=" << graph::split_list.size() << ";\nFORMAT CONFIDENCES=YES;\nMATRIX";
			} else {
				stream_nexus << "\nBEGIN Splits;\nDIMENSIONS ntax=" << denom_file_count << " nsplits=" << graph::split_list.size() << ";\nMATRIX";
			}
		}

		uint64_t pos = 0;
		//cleanliness.setFilteredCount(graph::split_list.size());
		color_t split_color;

		int split_num = 0; // number of split
		int split_size; // #taxa in split
		string split_comp = ""; // save split components/taxa

		for (auto& split : graph::split_list) {

			if(nexus_wanted || pdf_wanted || svg_wanted){ // nexus
				++split_num;
				split_size = 0; // reset split size
				split_comp = ""; // reset split components
			}

			double weight = split.first;
			split_color = split.second;
		// cleanliness.setSmallestWeight(weight, split.second);
			if (!output.empty()){
				stream << weight;    // weight of the split
			}
			if (bootstrap_no>0) {
				stream_bootstrap << ((1.0 * support_values[split.second]) / bootstrap_no);
				// stream_bootstrap << support_values[split.second];
			}
			if (raw_wanted){
				array<uint32_t,2> weights = graph::color_table[split_color];
				stream_raw << weights[0] << '\t' << weights[1];
			}
			for (uint64_t i = 0; i < num; ++i) {

				if (split_color.test(pos)) {
					if (i < denom_names.size()) {
						if (!output.empty()){
							stream << '\t' << denom_names[i]; // name of the file
						}
						if (bootstrap_no > 0) {
							stream_bootstrap << '\t' << denom_names[i]; // name of the file
						}
						if (raw_wanted){
							stream_raw << '\t' << denom_names[i];
						}

						if(nexus_wanted || pdf_wanted || svg_wanted){ // nexus
							++split_size;
							if(split_size > 1) split_comp += " "; // " " only if not first value
							split_comp += to_string(i+1);
						}

					}
				}
				split_color >>= 01u;
			}

			if(nexus_wanted || pdf_wanted || svg_wanted){
				if(bootstrap_no>0){ // Adding bootstrap values
					stream_nexus << "\n[" << split_num << ", size=" << split_size << "]\t";
					stream_nexus << weight << "\t" << ((1.0 * support_values[split.second]) / bootstrap_no) << "\t" << split_comp << ",";
				} else {
					stream_nexus << "\n[" << split_num << ", size=" << split_size << "]\t" << weight << "\t" << split_comp << ",";
				}
			}
			if (!output.empty()){
				stream << endl;
			}
			if(bootstrap_no>0){
				stream_bootstrap<<endl;
			}
			if(raw_wanted){
				stream_raw<<endl;
			}
		}

		if (nexus_wanted || pdf_wanted || svg_wanted){ // nexus
			stream_nexus << "\n;\nEND; [Splits]\n";
			// filter = strict (=greedy),  weakly (=greedyWC), tree (=?greedy)
			string fltr = "none"; // filter used for SplitsTree
			stream_nexus << "\nBEGIN st_Assumptions;\nuptodate;\nsplitstransform=EqualAngle UseWeights = true RunConvexHull = true DaylightIterations = 0\nOptimizeBoxesIterations = 0 SpringEmbedderIterations = 0;\nSplitsPostProcess filter=";
			stream_nexus << fltr << ";\nexclude no missing;\nautolayoutnodelabels;\nEND; [st_Assumptions]";
		}

		//cleanliness.calculateWeightBeforeCounter();

		if (!output.empty()){
			file.close();
		}
		if(bootstrap_no>0){
			file_bootstrap.close();
		}
		if(raw_wanted){
			file_raw.close();
		}
		if(nexus_wanted || pdf_wanted || svg_wanted){
			file_nexus.close();
			// naming modified nexus output file
			//string modded_file = nexus_color::modify_filename(nexus, "labeled_");
			// scale weights to 0-1
			nexus_color::scale_nexus(nexus, verbose, nexus_wanted);

			if(c_nexus_wanted){
				// use scaled file to open, mod and save in SplitsTree
				nexus_color::open_in_splitstree(nexus, pdf, svg, verbose, true, nexus);

				if(verbose) cout << "Adding color..." << endl << flush;
				nexus_color::color_nexus(nexus, groups, coloring);
				if(pdf_wanted || svg_wanted){
					nexus_color::open_in_splitstree(nexus, pdf, svg, verbose, false);
				}
			} else if(pdf_wanted || svg_wanted){
				nexus_color::open_in_splitstree(nexus, pdf, svg, verbose, true, nexus);
			}

			// Delete nexus file if only pdf wanted
			if(!nexus_wanted){
				remove(nexus.c_str());
				//remove(modded_file.c_str());
			}
		}
	}


// time measurement
    if (verbose) {
        if (!filter.empty()) {
          // cleanliness.reportCleanliness();
        }
        cout << " Done!" << flush;    // print progress and time
		end = chrono::high_resolution_clock::now();
		cout << " (" << util::format_time(end - begin) << ")" << endl;
    }
    return 0;
}

/** This function applies the specified filter to the given split list.
 * 
 * @param filter string specifying the type of filter
 * @param newick string with a file name to write a newick output to (or empty string)
 * @param map function that maps an integer to the original id, or null
 * @param split_list the list of splits to be filtered, e.g. graph::split_list
 * (@param support_values a hash map storing the absolut support values for each color set)
 * (@param bootstrap_no the number of bootstrap replicates for computing the per centage support)
 * @param verbose print progress
 * 
 */
void apply_filter(string filter, string newick, std::function<string(const uint64_t&)> map, multimap_<double, color_t>& split_list, bool verbose, uint64_t& num){
	apply_filter(filter, newick, map, split_list, nullptr, 0, verbose, num);
}

void apply_filter(string filter, string newick, std::function<string(const uint64_t&)> map, multimap_<double, color_t>& split_list, hash_map<color_t, uint32_t>* support_values, const uint32_t& bootstrap_no, bool verbose, uint64_t& num){

		if (!filter.empty()) {    // apply filter

			if (filter == "strict" || filter == "tree") {
				if (!newick.empty()) {
					ofstream file(newick);    // output file stream
					ostream stream(file.rdbuf());
					stream << graph::filter_strict(map, split_list, support_values, bootstrap_no, verbose); // filter and output
					file.close();
				} else {
					graph::filter_strict(split_list, verbose);
				}
			}
			else if (filter == "weakly") {
				graph::filter_weakly(split_list, verbose);
			}
            else if (filter == "planar") {
                graph::filter_planar(split_list, verbose, num);
            }
			else if (filter.find("tree") != -1 && filter.substr(filter.find("tree")) == "tree") {
				if (!newick.empty()) {
					ofstream file(newick);    // output file stream
					ostream stream(file.rdbuf());
					auto ot = graph::filter_n_tree(stoi(filter.substr(0, filter.find("tree"))), map, split_list, support_values, bootstrap_no, verbose);
					stream <<  ot;
					file.close();
				} else {
					graph::filter_n_tree(stoi(filter.substr(0, filter.find("tree"))), split_list, verbose);
				}
			}
		}	

}
