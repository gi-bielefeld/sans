#include "main.h"



/**
 * This is the entry point of the program.
 *
 * @param argc number of cmd args
 * @param argv cmd args
 * @return exit status
 */
int main(int argc, char* argv[]) {

    // print a help message describing the program arguments
    if (argc <= 1 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        cout << endl;
        cout << "SANS serif | version " << SANS_VERSION << endl;
        cout << "Usage: SANS [PARAMETERS]" << endl;
        cout << endl;
        cout << "  Input arguments:" << endl;
        cout << endl;
        cout << "    -i, --input   \t Input file: list of sequence files, one per line" << endl;
        cout << endl;
        cout << "    -g, --graph   \t Graph file: load a Bifrost graph, file name prefix" << endl;
#ifndef useBF
        cout << "                  \t (requires compiler flag -DuseBF, please edit makefile)" << endl;
#endif
        cout << endl;
        cout << "    -s, --splits  \t Splits file: load an existing list of splits file" << endl;
        cout << "                  \t (allows to filter -t/-f, other arguments are ignored)" << endl;
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
        cout << "    (at least --output or --newick must be provided, or both)" << endl;
        cout << endl;
        cout << "  Optional arguments:" << endl;
        cout << endl;
        cout << "    -k, --kmer    \t Length of k-mers (default: 31, or 10 for --amino and --code)" << endl;
        cout << endl;
//        cout << "    -w, --window  \t Number of k-mers per minimizer window (default: 1)" << endl;
//        cout << endl;
        cout << "    -t, --top     \t Number of splits in the output list (default: all)." << endl;
        cout << "                  \t Use -t <integer>n to limit relative to number of input files, or" << endl;
        cout << "                  \t use -t <integer> to limit by absolute value." << endl;
        cout << endl;
        cout << "    -m, --mean    \t Mean weight function to handle asymmetric splits" << endl;
        cout << "                  \t options: arith: arithmetic mean" << endl;
        cout << "                  \t          geom:  geometric mean (default)" << endl;
        cout << "                  \t          geom2: geometric mean with pseudo-counts" << endl;
        cout << endl;
        cout << "    -f, --filter  \t Output a greedy maximum weight subset" << endl;
        // cout << "                  \t additional output: (weighted) cleanliness, i.e., ratio of" << endl;
        // cout << "                  \t filtered splits w.r.t. original splits" << endl;
        cout << "                  \t options: strict: compatible to a tree" << endl;
        cout << "                  \t          weakly: weakly compatible network" << endl;
        cout << "                  \t          n-tree: compatible to a union of n trees" << endl;
        cout << "                  \t                  (where n is an arbitrary number)" << endl;
        cout << endl;
        cout << "    -C, --cluster \t Clustering. Output clusters of genomes that are not split into given file." << endl;
        cout << endl;
        cout << "    -x, --iupac   \t Extended IUPAC alphabet, resolve ambiguous bases or amino acids" << endl;
        cout << "                  \t Specify a number to limit the k-mers per position between" << endl;
        cout << "                  \t 1 (no ambiguity) and 4^k respectively 22^k (allows NNN...N)" << endl;
        cout << "                  \t Without --iupac respective k-mers are ignored" << endl;
        cout << endl;
        cout << "    -n, --norev   \t Do not consider reverse complement k-mers" << endl;
        cout << endl;
        cout << "    -a, --amino   \t Consider amino acids: --input provides amino acid sequences" << endl;
        cout << "                  \t Implies --norev and a default k of 10" << endl;
        cout << endl;
        cout << "    -c, --code   \t Translate DNA: --input provides coding sequences" << endl;
        cout << "                 \t Implies --norev and a default k of 10" << endl;
        cout << "                 \t optional: ID of the genetic code to be used" << endl;
        cout << "                 \t Default: 1" << endl;
        cout << "                 \t Use 11 for Bacterial, Archaeal, and Plant Plastid Code" << endl;
        cout << "                 \t (See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for details.)" << endl;
        cout << endl;
        cout << "    -v, --verbose \t Print information messages during execution" << endl;
        cout << endl;
        cout << "    -h, --help    \t Display this help page and quit" << endl;
        cout << endl;
        cout << "  Contact: sans-service@cebitec.uni-bielefeld.de" << endl;
        cout << "  Evaluation: https://www.surveymonkey.de/r/denbi-service?sc=bigi&tool=sans" << endl;
        cout << endl;
        return 0;
    }

    string input;    // name of input file
    string graph;    // name of graph file
    string splits;    // name of splits file
    string output;    // name of output file
    string newick;    // name of newick output file
    string translate; // name of translate file

    uint64_t kmer = 31;    // length of k-mers
    uint64_t window = 1;    // number of k-mers
    uint64_t num = 0;    // number of input files
    uint64_t top = -1;    // number of splits
    bool dyn_top = false; // bind number of splits to num

    auto mean = util::geometric_mean;    // weight function
    string filter;    // filter function
    uint64_t iupac = 1;    // allow extended iupac characters
    string clusterfile;    // cluster genomes and write to this file
    bool reverse = true;    // consider reverse complement k-mers
    bool verbose = false;    // print messages during execution
    bool amino = false;      // input files are amino acid sequences
    bool shouldTranslate = false;   // translate input files
    bool userKmer = false; // is k-mer default or custom
    uint64_t code = 1;

    // parse the command line arguments and update the variables above
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input") == 0) {
            input = argv[++i];    // Input file: list of sequence files, one per line
        }
        else if (strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--graph") == 0) {
            graph = argv[++i];    // Graph file: load a Bifrost graph, file name prefix
            #ifndef useBF
                cerr << "Error: requires compiler flag -DuseBF" << endl;
                return 1;
            #endif
        }
        else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--splits") == 0) {
            splits = argv[++i];    // Splits file: load an existing list of splits file
        }
        else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
            output = argv[++i];    // Output file: list of splits, sorted by weight desc.
        }
        else if (strcmp(argv[i], "-N") == 0 || strcmp(argv[i], "--newick") == 0) {
            newick = argv[++i];    // Output newick file
        }
        else if (strcmp(argv[i], "-k") == 0 || strcmp(argv[i], "--kmer") == 0) {
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
            i++;
            string top_str = argv[i];
            top = stoi(top_str); // Number of splits (default: all)

            if (top_str[top_str.size() - 1] == 'n'){ // Dynamic split num (default: false)
                dyn_top = true;
            }
        }
        else if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--mean") == 0) {
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
            filter = argv[++i];    // Filter a greedy maximum weight subset
            if (filter == "strict" || filter == "tree") {
                // compatible to a tree
            }
            else if (filter == "weakly") {
                // weakly compatible network
            }
            else if (filter.find("tree") != -1 && filter.substr(filter.find("tree")) == "tree") {
                stoi(filter.substr(0, filter.find("tree")));
            }
            else {
                cerr << "Error: unknown argument: --filter " << filter << endl;
                return 1;
            }
        }
        else if (strcmp(argv[i], "-C") == 0 || strcmp(argv[i], "--clustter") == 0) {
            clusterfile = argv[++i]; // cluster genomes and write to this file
        }
        else if (strcmp(argv[i], "-x") == 0 || strcmp(argv[i], "--iupac") == 0) {
            iupac = stoi(argv[++i]);    // Extended IUPAC alphabet, resolve ambiguous bases
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
        else {
            cerr << "Error: unknown argument: " << argv[i] <<  "\t type --help" << endl;
            return 1;
        }
    }

    // check for a new version of SANS at program start.
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


    if (!userKmer) {
        if (!amino) {
            kmer = 31;
        } else {
            kmer = 10;
        }
    }

    if (input.empty() && graph.empty() && splits.empty()) {
        cerr << "Error: missing argument: --input <file_name> or --graph <file_name>" << endl;
        return 1;
    }
    if (!input.empty() && !graph.empty() && !splits.empty()) {
        cerr << "Error: too many input arguments: --input, --graph, and --splits" << endl;
        return 1;
    }
    if (!graph.empty() && !splits.empty()) {
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

    if (output.empty() && newick.empty() && clusterfile.empty()) {
        cerr << "Error: missing argument: --output <file_name> or --newick <file_name> or --cluster <file_name>" << endl;
        return 1;
    }
    if (kmer > maxK && splits.empty()) {
        cerr << "Error: k-mer length exceeds -DmaxK=" << maxK << endl;
        return 1;
    }
    if (!newick.empty() && filter != "strict" && filter.find("tree") == -1) {
        cerr << "Error: Newick output only applicable in combination with -f strict or n-tree" << endl;
        return 1;
    }
    if (!input.empty() && !splits.empty()) {
        cerr << "Note: two input arguments --input and --splits were provided" << endl;
        cerr << "      --input is used for lookup only, no additional splits are inferred" << endl;
    }
    if (input.empty() && !splits.empty() && !newick.empty()) {
        cerr << "Note: Newick output from a list of splits, some taxa could be missing" << endl;
        cerr << "      --input can be used to provide the original list of taxa" << endl;
    }

    // check if we need to init translation
    if (shouldTranslate) {
        amino = true;
        if (!translator::init(code)) {
            cerr << "Error: No translation data found" << translate << endl;
        }
    }

    // parse the list of input sequence files
    vector<string> files;
    hash_map<string, uint64_t> name_table;

    if (!input.empty()) {
        ifstream file(input);
        if (!file.good()) {
            cerr << "Error: could not read input file: " << input << endl;
            return 1;
        }
        string line;
        while (getline(file, line)) {
            files.emplace_back(line);
            name_table[line] = num++;
        }
        file.close();

        if (num > maxN) {
            cerr << "Error: number of files exceeds -DmaxN=" << maxN << endl;
            return 1;
        }
    }

    // Set dynamic top
    if (dyn_top){
        top = top * num;
    }

#ifdef useBF
    // load an existing Bifrost graph
    ColoredCDBG<> cdbg(kmer);
    if (!graph.empty()) {
        if (cdbg.read(graph + ".gfa", graph + ".bfg_colors", 1, verbose)) {
            num += cdbg.getNbColors();
        } else {
            cerr << "Error: could not load Bifrost graph" << endl;
            return 1;
        }

        if (kmer != cdbg.getK()) {
            kmer = cdbg.getK();
            if (kmer > maxK) {
                cerr << "Error: k-mer length exceeds -DmaxK=" << maxK << endl;
                return 1;
            }
            cerr << "Warning: setting k-mer length to " << kmer << endl;
        }
        if (num > maxN) {
            cerr << "Error: number of colors exceeds -DmaxN=" << maxN << endl;
            return 1;
        }
        if (verbose) {
            cout << endl;
        }
    }
#endif

    chrono::high_resolution_clock::time_point begin = chrono::high_resolution_clock::now();    // time measurement
    graph::init(top, amino); // initialize the toplist size and the allowed characters
    if (!splits.empty()) {
        ifstream file(splits);
        if (!file.good()) {
            cerr << "Error: could not read splits file: " << splits << endl;
            return 1;
        }
        string line;
        while (getline(file, line)) {
            uint64_t curr = line.find('\t');
            double weight = stod(line.substr(0, curr));
            uint64_t next = curr + 1;

            color_t color = 0;
            do {
                curr = line.find('\t', next);
                string name = line.substr(next, curr-next);
                if (name_table.find(name) == name_table.end()) {
                    files.emplace_back(name);
                    name_table[name] = num++;
                    if (num > maxN) {
                        cerr << "Error: number of files exceeds -DmaxN=" << maxN << endl;
                        return 1;
                    }
                }
                color::set(color, name_table[name]);
                next = curr + 1;
            } while (curr != string::npos);

            graph::add_split(weight, color);
        }
        file.close();
    }
    
    kmer::init(kmer);      // initialize the k-mer length
    kmerAmino::init(kmer); // initialize the k-mer length
    color::init(num);    // initialize the color number

    if (!input.empty() && splits.empty()) {
        if (verbose) {
            cout << "SANS::main(): Reading input files..." << flush;
        }
        string sequence;    // read in the sequence files and extract the k-mers

	// determine folder the list is contained in
        string folder="";
	uint64_t found=input.find_last_of("/\\");
	if (found!=string::npos){
	  folder=input.substr(0,found+1);
	}

        for (uint64_t i = 0; i < files.size(); ++i) {
            ifstream file(folder+files[i]);    // input file stream
            if (!file.good()) {
                cout << "\33[2K\r" << "\u001b[31m" << folder+files[i] << " (ERR)" << "\u001b[0m" << endl;    // could not read file
            }
            else if (verbose) {
                cout << "\33[2K\r" << folder+files[i] << " (" << i+1 << "/" << files.size() << ")" << endl;    // print progress
            }
            count::deleteCount();

            string appendixChars; 
            string line;    // read the file line by line
            while (getline(file, line)) {
                if (line.length() > 0) {
                    if (line[0] == '>' || line[0] == '@') {    // FASTA & FASTQ header -> process
                        if (window > 1) {
                            iupac > 1 ? graph::add_minimizers(sequence, i, reverse, window, iupac)
                                      : graph::add_minimizers(sequence, i, reverse, window);
                        } else {
                            iupac > 1 ? graph::add_kmers(sequence, i, reverse, iupac)
                                      : graph::add_kmers(sequence, i, reverse);
                        }

                        sequence.clear();

                        if (verbose) {
                            cout << "\33[2K\r" << line << flush << endl;    // print progress
                        }
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
                iupac > 1 ? graph::add_minimizers(sequence, i, reverse, window, iupac)
                          : graph::add_minimizers(sequence, i, reverse, window);
            } else {
                iupac > 1 ? graph::add_kmers(sequence, i, reverse, iupac)
                          : graph::add_kmers(sequence, i, reverse);
            }
            sequence.clear();

            if (verbose) {
                cout << "\33[2K\r" << flush;
            }
            file.close();
        }
    }

#ifdef useBF
    if (!graph.empty()) {
        if (verbose) {
            cout << "SANS::main(): Processing unitigs..." << flush;
        }
        uint64_t cur = 0, progress;
        uint64_t max = cdbg.size();

        for (auto& unitig : cdbg) {
            if (verbose) {
                cout << "\33[2K\r" << "Processed " << cur << " unitigs (" << 100*cur/max << "%) " << flush;
            }   cur++;

            auto sequence = unitig.mappedSequenceToString();
            auto *colors = unitig.getData()->getUnitigColors(unitig);
            auto it = colors->begin(unitig);
            auto end = colors->end();

            for (; it != end; ++it) {
                auto seq = sequence.substr(it.getKmerPosition(), kmer);
                auto col = files.size() + it.getColorID();
                graph::add_kmers(seq, col, reverse);
            }
        }
        if (verbose) {
            cout << "\33[2K\r" << "Processed " << max << " unitigs (100%)" << endl;
        }
    }
#endif

    // function to map color position to file name
    std::function<string(const uint64_t&)> map=[=](uint64_t i) {
        if (i < files.size()) return files[i];
        #ifdef useBF
        else return cdbg.getColorName(i-files.size());
        #endif
        cerr << "Error: color bit does not correspond to color name" << endl;
        exit(EXIT_FAILURE);
    };

    if (verbose) {
        cout << "Processing splits..." << flush;
    }
    graph::add_weights(mean, verbose);    // accumulate split weights

    if (verbose) {
        cout << "\33[2K\r" << "Filtering splits..." << flush;
    }

    //cleanliness cleanliness;
    if (!filter.empty()) {    // apply filter
        for (auto& split : graph::split_list) {
            //cleanliness.addWeightStateBefore(split.first, split.second);
        }

        if (filter == "strict" || filter == "tree") {
            if (!newick.empty()) {
                ofstream file(newick);    // output file stream
                ostream stream(file.rdbuf());
                stream << graph::filter_strict(map, verbose);    // filter and output
                file.close();
            } else {
               graph::filter_strict(verbose);
            }
        }
        else if (filter == "weakly") {
           graph::filter_weakly(verbose);
        }
        else if (filter.find("tree") != -1 && filter.substr(filter.find("tree")) == "tree") {
            if (!newick.empty()) {
                ofstream file(newick);    // output file stream
                ostream stream(file.rdbuf());
                auto ot = graph::filter_n_tree(stoi(filter.substr(0, filter.find("tree"))), map, verbose);
                stream <<  ot;
                file.close();
            } else {
               graph::filter_n_tree(stoi(filter.substr(0, filter.find("tree"))), verbose);
            }
        }
    }
    
    // CLUSTERING CLUSTERING CLUSTERING CLUSTERING CLUSTERING CLUSTERING CLUSTERING CLUSTERING CLUSTERING CLUSTERING
    if (!clusterfile.empty()) {
		// create one big cluster with all IDs
		std::set<uint64_t> initial_cluster;
		for (uint64_t i=0; i<color::n; i++) initial_cluster.insert(i);
		
		//create bag of clusters
		vector<set<uint64_t>> clusters;
		clusters.push_back(initial_cluster);
		uint64_t cur = 0;
		uint64_t max = graph::split_list.size();
		
		//iterate over all splits
		for (auto& split : graph::split_list) {
			//iterate over IDs in split
			uint64_t zero=0;
			color_t my_split = split.second;
			for (uint64_t i = 0; i < num; ++i) {
				if (color::test(my_split, zero)) { // we shift over the split, so always test at position 0
					//find cluster containing i
					set<uint64_t> new_c;
					for (std::vector<set<uint64_t>>::iterator it_c = clusters.begin(); it_c != clusters.end(); ++it_c) {
						set<uint64_t> & c = *it_c;

						std::set<uint64_t>::iterator it;
						uint64_t offset=0; // keep i at current pos, use offset instead
						it=c.find(i+offset);
						if (it!=c.end()) { //hit!
							// identify overlap
							while (it!=c.end()) {
								new_c.insert(*it);
								color::erase(my_split,offset);
								while(!color::test(my_split, offset)) {
									offset++;
									if (i+offset>=num) break;
								}
								it=c.find(i+offset);
							}
							if (new_c.size() == c.size()){
								break; //nothing to be done
							}
							//remove overlap from cluster
							std::set<uint64_t>::iterator it1 = c.begin();
							std::set<uint64_t>::iterator it2 = new_c.begin();
							while ( (it1 != c.end()) && (it2 != new_c.end()) ) {
								if (*it1 == *it2) {
									it1=c.erase(it1);
									++it2;
								} else if (*it2 < *it1) {
									++it2;
								} else { // *it1 < *it2
									++it1;
								}
							}
							// add new cluster to bag of clusters
							clusters.push_back(new_c);
							// continue with next element in split
							break;
						}
					}
				}
				// done with last splitting of cluster,
				// shift, i.e. check next position in split
				my_split >>= 01u;
			}
			if (verbose) {
				cout << "\33[2K\r" << "Processed " << cur++ << " splits (" << 100*cur/max << "%) " << flush;
			}
		}
		if (verbose) {
			cout << "\33[2K\r" << "Processed " << max << " splits (100%)" << endl;
		}
		// output
		ofstream file(clusterfile);    // output file stream
		ostream stream(file.rdbuf());
		for (std::set<uint64_t> const &c : clusters) {
			stream << c.size() << endl;
// 			    for(const uint64_t i : c){std::cout << i << flush << " ";}
// 			    std::cout << flush << endl;
			    
		}
		file.close();
	}
    // END_OF  CLUSTERING CLUSTERING CLUSTERING CLUSTERING CLUSTERING CLUSTERING CLUSTERING

    if (verbose) {
        cout << "\33[2K\r" << "Please wait..." << flush << endl;
    }
    ofstream file(output);    // output file stream
    ostream stream(file.rdbuf());

    uint64_t pos = 0;
    //cleanliness.setFilteredCount(graph::split_list.size());
    for (auto& split : graph::split_list) {
        double weight = split.first;
       // cleanliness.setSmallestWeight(weight, split.second);
        stream << weight;    // weight of the split
        for (uint64_t i = 0; i < num; ++i) {
            if (color::test(split.second, pos)) {
                if (i < files.size())
                    stream << '\t' << files[i];    // name of the file
                #ifdef useBF
                else
                    stream << '\t' << cdbg.getColorName(i-files.size());
                #endif
            }
            split.second >>= 01u;
        }
        stream << endl;
    }

    //cleanliness.calculateWeightBeforeCounter();

    file.close();

    chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();    // time measurement

    if (verbose) {
        if (!filter.empty()) {
          // cleanliness.reportCleanliness();
        }
        cout << " Done!" << flush << endl;    // print progress and time
        cout << " (" << util::format_time(end - begin) << ")" << endl;
    }
    return 0;
}



