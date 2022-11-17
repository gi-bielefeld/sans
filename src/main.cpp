#include "main.h"

/**
 * This is the entry point of the program.
 *
 * @param argc number of cmd args
 * @param argv cmd args
 * @return exit status
 */
int main(int argc, char* argv[]) {
    ios_base::sync_with_stdio(false);

    // check for a new version of SANS-KC at program start
    if (!system("wget --timeout=1 --tries=1 -qO- https://gitlab.ub.uni-bielefeld.de/gi/sans/raw/kc/src/main.h | grep -q SANS_VERSION")
      && system("wget --timeout=1 --tries=1 -qO- https://gitlab.ub.uni-bielefeld.de/gi/sans/raw/kc/src/main.h | grep -q " SANS_VERSION)) {
        cerr << "NEW VERSION AVAILABLE: https://gitlab.ub.uni-bielefeld.de/gi/sans/tree/kc" << endl;
    }
    // print a help message describing the program arguments
    if (argc <= 1 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        cout << endl;
        cout << "SANS-KC | version " << SANS_VERSION << endl;
        cout << "Usage: SANS [PARAMETERS]" << endl;
        if (argc <= 1) util::print_help();
        else           util::print_extended_help();
        cout << "  Contact: sans-service@cebitec.uni-bielefeld.de" << endl;
        cout << "  Evaluation: https://www.surveymonkey.de/r/denbi-service?sc=bigi&tool=sans" << endl;
        cout << endl;
        return 0;
    }

    string input;    // name of input file
    string graph;    // name of graph file
    string splits;    // name of splits file
    string output;    // name of output file
    string newick;    // name of newick file
    string counts;    // name of counts file
    bool interactive = false;    // lookup k-mer

    uint64_t kmer = 0;    // length of k-mers
    string pattern;    // pattern of gapped k-mers
    uint64_t window = 1;    // number of k-mers per window
    uint64_t iupac = 1;    // allow extended IUPAC characters
    uint64_t quality = 1;    // min. coverage threshold for k-mers
    bool reverse = false;    // consider reverse complement k-mers

    uint64_t num = 0;    // number of input files
    uint64_t top = -1;    // number of top splits
    string filter;    // filter function of splits
    uint64_t T = 1;    // number of parallel threads
    bool verbose = false;    // print messages during execution

    auto arithmetic_mean = [] (const uint32_t& x, const uint32_t& y) { return x / 2.0 + y / 2.0; };
    auto geometric_mean  = [] (const uint32_t& x, const uint32_t& y) { return sqrt(x) * sqrt(y); };
    auto geometric_mean2 = [] (const uint32_t& x, const uint32_t& y) { return sqrt(x+1) * sqrt(y+1) - 1; };
    function<double(const uint32_t&, const uint32_t&)> default_mean  = geometric_mean;   // weight function

    // parse the command line arguments and update the variables above
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input") == 0) {
            input = argv[++i];    // Input FASTA files: list of sequence files, one per line
        }
        else if (strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--graph") == 0) {
            graph = argv[++i];    // Input Graph file: load a Bifrost graph, filename prefix
            #ifndef useBF
                cerr << "Error: requires compiler flag -DuseBF, please see makefile" << endl;
                return 1;
            #endif
        }
        else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--splits") == 0) {
            splits = argv[++i];    // Input Splits file: load an existing list of splits file
        }
        else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
            output = argv[++i];    // Output TSV file: list of splits, sorted by weight desc.
        }
        else if (strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--newick") == 0) {
            newick = argv[++i];    // Output Newick file: convert splits to a tree topology
        }
        else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--counts") == 0) {
            counts = argv[++i];    // Output K-mer file: list k-mer occurrence per input file
        }
        else if (strcmp(argv[i], "-C") == 0) {    /* hidden option */
            interactive = true;    // Interactive mode: list k-mer occurrence per input file
        }

        else if (strcmp(argv[i], "-k") == 0 || strcmp(argv[i], "--kmer") == 0) {
            kmer = stoi(argv[++i]);    // Length of k-mers (default: 31)
        }
        else if (strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "--gapped") == 0) {
            pattern = argv[++i];    // Pattern of gapped k-mers (default: no gaps)
        }
        else if (strcmp(argv[i], "-w") == 0 || strcmp(argv[i], "--window") == 0) {
            window = stoi(argv[++i]);    // Number of k-mers per window (default: 1)
        }
        else if (strcmp(argv[i], "-x") == 0 || strcmp(argv[i], "--iupac") == 0) {
            iupac = stoi(argv[++i]);    // Extended IUPAC alphabet, resolve ambiguous bases
        }
        else if (strcmp(argv[i], "-q") == 0 || strcmp(argv[i], "--qualify") == 0) {
            quality = stoi(argv[++i]);    // Discard k-mers below a min. coverage threshold
        }
        else if (strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--reverse") == 0) {
            reverse = true;    // Keep one repr. for reverse complement k-mers
        }

        else if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--top") == 0) {
            top = stoi(argv[++i]);    // Number of top splits to output (default: all)
        }
        else if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--mean") == 0) {
            string arg = argv[++i];    // Mean weight function to handle asymmetric splits
            if      (arg == "arith") default_mean = arithmetic_mean;
            else if (arg == "geom")  default_mean = geometric_mean;
            else if (arg == "geom2") default_mean = geometric_mean2;
            else {
                cerr << "Error: unknown argument: --mean " << arg << endl;
                cerr << "       type --help to see a list of supported options" << endl;
                return 1;
            }
        }
        else if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--filter") == 0) {
            filter = argv[++i];    // Output a greedy maximum weight subset of splits
            if      (filter == "strict" || filter == "tree"); // compatible to a tree
            else if (filter == "weakly");                     // weakly compatible network
            else if (filter.find("tree") != -1 && filter.substr(filter.find("tree")) == "tree")
                stoi(filter.substr(0, filter.find("tree")));
            else {
                cerr << "Error: unknown argument: --filter " << filter << endl;
                cerr << "       type --help to see a list of supported options" << endl;
                return 1;
            }
        }

        else if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--threads") == 0) {
            string newT = argv[++i];
            if (newT == "max" || newT == "MAX") {
                uint64_t maxT = thread::hardware_concurrency();
                if (maxT == 0) {
                    cerr << "Error: could not determine number of threads" << endl;
                    return 1;
                } else {
                    T = maxT;
                    cerr << "Note: running SANS-KC using " << T << " threads" << endl;
                }
            } else {
                T = stoi(newT);    // Number of parallel threads (default: 1)
            }
        }
        else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose = true;    // Print information messages during execution
        }
        else if (strcmp(argv[i], "-rv") == 0 || strcmp(argv[i], "-vr") == 0) {
            reverse = true;    // Keep one repr. for reverse complement k-mers
            verbose = true;    // Print information messages during execution
        }
        else if (! (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) ) {
            cerr << "Error: unknown argument: " << argv[i] << endl;
            cerr << "       type --help to see the full list of parameters" << endl;
            return 1;
        }
    }

    if (input.empty() && graph.empty() && splits.empty()) {
        cerr << "Error: missing argument: --input, --graph, or --splits <file_name>" << endl;
        return 1;
    }
    if (!graph.empty() && !splits.empty()) {
        if (!input.empty())
             cerr << "Error: too many arguments: --input, --graph, and --splits <file_name>" << endl;
        else cerr << "Error: too many input arguments: --graph and --splits <file_name>" << endl;
        return 1;
    }
    if (output.empty() && newick.empty() && (counts.empty() && !interactive)) {
        cerr << "Error: missing argument: --output, --newick, or --counts <file_name>" << endl;
        return 1;
    }
    if (!splits.empty() && !(counts.empty() && !interactive)) {
        cerr << "Error: k-mer counts cannot be calculated if the input is a list of splits" << endl;
        return 1;
    }
    if (!newick.empty() && filter != "strict" && filter.find("tree") == -1) {
        cerr << "Error: Newick output only applicable in combination with -f strict/n-tree" << endl;
        return 1;
    }

    if (!pattern.empty()) {
        uint64_t lmer = 0;    // length without gaps
        for (size_t i = 0; i < pattern.length(); ++i) {
            if (pattern[i] == '1') {
                ++lmer; continue;
            }
            if (pattern[i] != '0') {
                cerr << "Error: pattern must be a sequence of 0s and 1s, where 0 means a gap" << endl;
                return 1;
            }
            if (reverse && pattern[i] != pattern[pattern.length()-1-i]) {
                cerr << "Error: pattern should be symmetric to work with reverse complements" << endl;
                return 1;
            }
        }
        if (kmer != 0 && kmer != lmer && kmer != pattern.length()) {
            cerr << "Error: pattern length does not match the given k-mer length" << endl;
            return 1;
        }
        kmer = pattern.length();

        if (lmer == 0) {
            cerr << "Error: pattern must be a sequence of 0s and 1s, with at least one 1" << endl;
            return 1;
        }
        if (lmer == pattern.length()) {
            pattern.clear();
        }
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
    if (!splits.empty()) {
        if (!input.empty()) {
            cerr << "Note: two input arguments --input and --splits were provided" << endl;
            cerr << "      --input is used for lookup only, no additional splits are inferred" << endl;
        }
        else if (!newick.empty()) {
            cerr << "Note: Newick output from a list of splits, some taxa could be missing" << endl;
            cerr << "      --input can be used to provide the original list of taxa" << endl;
        }
    }

    // parse the list of input sequence files
    vector<string> files;
    hash_map<string, size1N_t> name_table;
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
    }

#ifdef useBF
    // load an existing Bifrost graph file
    ColoredCDBG<> cdbg(kmer);
    if (!graph.empty()) {
        if (cdbg.read(graph + ".gfa", graph + ".bfg_colors", 1, verbose)) {
            num += cdbg.getNbColors();
            if (verbose) cerr << endl;
        } else {
            cerr << "Error: could not load Bifrost graph: " << graph << endl;
            return 1;
        }
        if (kmer != 0 && kmer != cdbg.getK()) {
            if (pattern.empty()) {
                cerr << "Warning: graph file does not match the given k-mer length" << endl;
            } else {
                cerr << "Error: graph file does not match the given pattern length" << endl;
                return 1;
            }
        }
        kmer = cdbg.getK();
    }
#endif

    // load & parse an existing list of splits
    if (!splits.empty()) {
        ifstream file(splits);
        if (!file.good()) {
            cerr << "Error: could not read splits file: " << splits << endl;
            return 1;
        }
        string line;
        while (getline(file, line)) {
            size_t curr = line.find('\t');
            size_t next = curr + 1;
            do {
                curr = line.find('\t', next);
                string name = line.substr(next, curr-next);
                if (name_table.find(name) == name_table.end()) {
                    files.emplace_back(name);
                    name_table[name] = num++;
                }
                next = curr + 1;
            } while (curr != string::npos);
        }
        file.close();
    }

    if (kmer == 0 && splits.empty()) {
        kmer = min(maxK, 31);    // k-mer length not specified & not estimated by graph/pattern
        //->default
    }
    if (kmer > maxK) {
        cerr << "Error: k-mer length exceeds -DmaxK=" << maxK << ", please see makefile" << endl;
        return 1;
    }
    if (num > maxN) {
        cerr << "Error: number of files exceeds -DmaxN=" << maxN << ", please see makefile" << endl;
        return 1;
    }

    auto begin = chrono::high_resolution_clock::now();    // time measurement
    graph::init(T, pattern, quality, reverse);    // initialize the graph
    tree::init(top);    // initialize the splits list size
    kmer::init(kmer);    // initialize the k-mer length
    color::init(num);    // initialize the color number

    if (!input.empty() && splits.empty()) {
        if (verbose)
            cerr << "\33[2K\r" << "Reading input files..." << flush;

       auto lambda = [&] (const uint64_t& T, const size1N_t& lower_bound, const size1N_t& upper_bound) {
        for (size1N_t i = lower_bound; i < upper_bound; ++i) {
            ifstream file(files[i]);    // input file stream
            if (!file.good())
                cerr << "\33[2K\r" << "\u001b[31m" << files[i] << " (ERR)" << "\u001b[0m" << endl;    // could not read file
            else if (verbose)
                cerr << "\33[2K\r" << files[i] << " (" << i+1 << "/" << files.size() << ")" << endl;    // print progress

            string line;    // read the file line by line
            string sequence;    // read in the sequence files and extract the k-mers
            while (getline(file, line)) {
                if (line.length() > 0) {
                    if (line[0] == '>' || line[0] == '@') {    // FASTA & FASTQ header -> process
                        if (window <= 1)
                             iupac <= 1 ? graph::add_kmers(T, sequence, i)
                                        : graph::add_kmers(T, sequence, i, iupac);
                        else iupac <= 1 ? graph::add_minimizers(T, sequence, i, window)
                                        : graph::add_minimizers(T, sequence, i, window, iupac);
                        sequence.clear();

                        if (verbose)
                            cerr << "\33[2K\r" << line << flush;    // print progress
                    }
                    else if (line[0] == '+') {    // FASTQ quality values -> ignore
                        getline(file, line);
                    } else {
                        transform(line.begin(), line.end(), line.begin(), ::toupper);
                        sequence += line;    // FASTA & FASTQ sequence -> read
                    }
                }
            }
            if (window <= 1)
                 iupac <= 1 ? graph::add_kmers(T, sequence, i)
                            : graph::add_kmers(T, sequence, i, iupac);
            else iupac <= 1 ? graph::add_minimizers(T, sequence, i, window)
                            : graph::add_minimizers(T, sequence, i, window, iupac);
            sequence.clear();

            if (verbose)
                cerr << "\33[2K\r" << flush;
            graph::clear_thread(T);
            file.close();
        }};

        const size1N_t MAX = files.size();
        const size1N_t SIZE = MAX / T;
        const size1N_t CARRY = MAX % T;

        vector<size1N_t> size(T);
        for (uint64_t x = 0; x < T; ++x) {
            size[x] = SIZE;
        }
        for (uint64_t x = 0; x < CARRY; ++x) {
            size[x] += 1;
        }
        size1N_t curr_size = 0;

        vector<thread> thr(T);
        for (uint64_t x = 0; x < T; ++x) {
            thr[x] = thread(lambda, x, curr_size, curr_size + size[x]);
            curr_size += size[x];
        }
        for (uint64_t x = 0; x < T; ++x) {
            thr[x].join();
        }
        size.clear();

        while (thr.size() > 1) {
            for (uint64_t x = 0; x < thr.size()-1; x+=2) {
                thr[x] = thread(graph::merge_threads, x, x+1);
                if (verbose) cerr << ":" << flush;
            }   if (verbose) cerr << "\r" << flush;

            for (uint64_t x = 0; x < thr.size()-1; x+=2) {
                thr[x].join();
                if (verbose) cerr << "." << flush;
            }   if (verbose) cerr << "\r" << flush;

            uint64_t next_size = thr.size() / 2;
            for (uint64_t x = 1; x <= next_size; ++x) {
                thr.erase(thr.begin()+x);
                graph::erase_thread(x);
                if (verbose) cerr << "." << flush;
            }   if (verbose) cerr << "\33[2K\r" << flush;
        }
    }

#ifdef useBF
    if (!graph.empty()) {
        if (verbose)
            cerr << "\33[2K\r" << "Processing unitigs..." << flush;
        uint64_t cur = 0, prog = 0, next;
        uint64_t max = cdbg.size();

        for (auto& unitig : cdbg) {
            if (verbose) {
                next = 100 * cur / max;
                 if (prog < next)  cerr << "\33[2K\r" << "Processed " << cur << " unitigs (" << next << "%) " << flush;
                prog = next; cur++;
            }
            auto sequence = unitig.mappedSequenceToString();
            auto* colors = unitig.getData()->getUnitigColors(unitig);

            for (auto it = colors->begin(unitig); it != colors->end(); ++it) {
                auto substr = sequence.substr(it.getKmerPosition(), kmer);
                auto color = files.size() + it.getColorID();
                graph::add_kmers(0, substr, color);
            }
        }
        if (verbose)
            cerr << "\33[2K\r" << "Processed " << max << " unitigs (100%)" << endl;
        graph::clear_thread(0);
    }
#endif

    if (!splits.empty()) {
        if (verbose)
            cerr << "\33[2K\r" << "Reading splits file..." << flush;
        ifstream file(splits);

        string line;
        while (getline(file, line)) {
            size_t curr = line.find('\t');
            double weight = stod(line.substr(0, curr));
            color_t color = 0b0u;
            size_t next = curr + 1;
            do {
                curr = line.find('\t', next);
                string name = line.substr(next, curr-next);
                color::set(color, name_table[name]);
                next = curr + 1;
            } while (curr != string::npos);

            tree::insert_split(weight, color);
        }
        if (verbose)
            cerr << "\33[2K\r" << flush;
        file.close();
    }

    // function to map a color position to a file name
    function<string(const size1N_t&)> color_name = [&] (const size1N_t& i) {
        if (i < files.size()) return files[i];
       #ifdef useBF
         else return cdbg.getColorName(i-files.size());
       #endif
        cerr << "Error: color bit does not correspond to a color name" << endl;
        exit(1);
    };

    if (interactive) {    // lookup k-mer, interactive mode
        string query;
        cerr << ">>> " << flush;    // display a command line prompt
        while (getline(cin, query)) {    // wait for user to enter a k-mer
            cerr << "\33[2K\r" << flush;
                                            // remove leading whitespaces
            auto l1 = find_if(query.begin(), query.end(), ::isspace);
            auto l2 = find_if_not(query.begin(), query.end(), ::isspace);
             if (l1 < l2) query.erase(l1, l2);    // remove trailing whitespaces
            auto r1 = find_if_not(query.rbegin(), query.rend(), ::isspace);
            auto r2 = find_if(query.rbegin(), query.rend(), ::isspace);
             if (r1.base() < r2.base()) query.erase(r1.base(), r2.base());

            if (query.empty()) { cerr << ">>> " << flush; continue; }
            transform(query.begin(), query.end(), query.begin(), ::toupper);
            replace(query.begin(), query.end(), '.', 'N');    // allow . gaps
            replace(query.begin(), query.end(), '-', 'N');    // allow - gaps
            replace(query.begin(), query.end(), '*', 'N');    // allow * gaps

            cerr << "Searching..." << flush;
            auto iterator = graph::lookup_kmer(query);
            cerr << "\33[2K\r" << flush;
            string kmer, color;    // iterate over the hash table
            while (iterator(kmer, color)) {
                for (size_t i = 0; i < pattern.length(); ++i)
                    if (pattern[i] == '0') kmer[i] = '-';
                cerr << "... " << flush;
                cout << kmer << ": " << color << endl;
                cerr << "\33[2K\r" << flush;
            }
            cerr << ">>> " << flush;
        }
        cerr << "\33[2K\r" << flush;
    }
    if (!counts.empty()) {    // lookup k-mer, output to file
        if (verbose)
            cerr << "\33[2K\r" << "Please wait... " << flush;
        ofstream file(counts);    // output file stream
        ostream stream(file.rdbuf());
        auto iterator = graph::lookup_kmer();
        string kmer, color;    // iterate over the hash table
        while (iterator(kmer, color)) {
            for (size_t i = 0; i < pattern.length(); ++i)
                if (pattern[i] == '0') kmer[i] = '-';
            stream << kmer << ": " << color << endl;
        }
        file.close();
    }

    if (!output.empty() || !newick.empty()) {
        if (verbose)
            cerr << "\33[2K\r" << "Processing splits..." << flush;
        graph::calc_weights(default_mean, verbose);    // accumulate split weights

        if (!filter.empty()) {    // apply filter
            if (verbose)
                cerr << "\33[2K\r" << "Filtering splits..." << flush;
            if (filter == "strict" || filter == "tree")
                tree::filter_strict(verbose);
            else if (filter == "weakly")
                tree::filter_weakly(verbose);
            else if (filter.find("tree") != -1 && filter.substr(filter.find("tree")) == "tree")
                tree::filter_n_tree(stoi(filter.substr(0, filter.find("tree"))), verbose);
        }

        if (verbose)
            cerr << "\33[2K\r" << "Please wait... " << flush;
        if (!output.empty()) {
            ofstream file(output);    // output file stream
            ostream stream(file.rdbuf());
            for (auto& split : tree::splits) {
                stream << split.first;    // weight of the split
                for (size1N_t i = 0; i < num; ++i) {
                    if (color::test(split.second, i)) {
                        stream << '\t' << color_name(i);    // name of the file
                    }
                }
                stream << endl;
            }
            file.close();
        }
        if (!newick.empty()) {
            ofstream file(newick);    // output file stream
            ostream stream(file.rdbuf());
            stream << tree::build_string(color_name);    // filter and print tree
            file.close();
        }
    }

    auto end = chrono::high_resolution_clock::now();    // time measurement
    if (verbose)    // print progress and time
        cerr << "Done! (" << util::format_time(end - begin) << ")" << endl;
    return 0;
}
