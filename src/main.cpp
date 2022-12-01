#include "main.h"
#include "ansi.h"

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
        $link << "NEW VERSION AVAILABLE: https://gitlab.ub.uni-bielefeld.de/gi/sans/tree/kc" << _end$;
    }
    // print a help message describing the program arguments
    if (argc <= 1 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        $out << end$;
        $out << "SANS-KC | version " << SANS_VERSION << end$;
        $out << "Usage: SANS [PARAMETERS]" << end$;
        if (argc <= 1) util::print_help();
        else           util::print_extended_help();
        $out << "  Contact: sans-service@cebitec.uni-bielefeld.de" << end$;
        $out << "  Evaluation: https://www.surveymonkey.de/r/denbi-service?sc=bigi&tool=sans" << end$;
        $out << end$;
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
                $err << "Error: requires compiler flag -DuseBF, please see makefile" << _end$;
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
                $err << "Error: unknown argument: --mean " << arg << _end$;
                $err << "       type --help to see a list of supported options" << _end$;
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
                $err << "Error: unknown argument: --filter " << filter << _end$;
                $err << "       type --help to see a list of supported options" << _end$;
                return 1;
            }
        }

        else if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--threads") == 0) {
            string newT = argv[++i];
            if (newT == "max" || newT == "MAX") {
                uint64_t maxT = thread::hardware_concurrency();
                if (maxT == 0) {
                    $err << "Error: could not determine number of threads" << _end$;
                    return 1;
                } else {
                    T = maxT;
                    $note << "Note: running SANS-KC using " << T << " threads" << _end$;
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
            $err << "Error: unknown argument: " << argv[i] << _end$;
            $err << "       type --help to see the full list of parameters" << _end$;
            return 1;
        }
    }

    if (input.empty() && graph.empty() && splits.empty()) {
        $err << "Error: missing argument: --input, --graph, or --splits <file_name>" << _end$;
        return 1;
    }
    if (!graph.empty() && !splits.empty()) {
        if (!input.empty())
             $err << "Error: too many arguments: --input, --graph, and --splits <file_name>" << _end$;
        else $err << "Error: too many input arguments: --graph and --splits <file_name>" << _end$;
        return 1;
    }
    if (output.empty() && newick.empty() && (counts.empty() && !interactive)) {
        $err << "Error: missing argument: --output, --newick, or --counts <file_name>" << _end$;
        return 1;
    }
    if (!splits.empty() && !(counts.empty() && !interactive)) {
        $err << "Error: k-mer counts cannot be calculated if the input is a list of splits" << _end$;
        return 1;
    }
    if (!newick.empty() && filter != "strict" && filter.find("tree") == -1) {
        $err << "Error: Newick output only applicable in combination with -f strict/n-tree" << _end$;
        return 1;
    }

    if (!pattern.empty()) {
        uint64_t lmer = 0;    // length without gaps
        for (size_t i = 0; i < pattern.length(); ++i) {
            if (pattern[i] == '1') {
                ++lmer; continue;
            }
            if (pattern[i] != '0') {
                $err << "Error: pattern must be a sequence of 0s and 1s, where 0 means a gap" << _end$;
                return 1;
            }
            if (reverse && pattern[i] != pattern[pattern.length()-1-i]) {
                $err << "Error: pattern should be symmetric to work with reverse complements" << _end$;
                return 1;
            }
        }
        if (kmer != 0 && kmer != lmer && kmer != pattern.length()) {
            $err << "Error: pattern length does not match the given k-mer length" << _end$;
            return 1;
        }
        kmer = pattern.length();

        if (lmer == 0) {
            $err << "Error: pattern must be a sequence of 0s and 1s, with at least one 1" << _end$;
            return 1;
        }
        if (lmer == pattern.length()) {
            pattern.clear();
        }
    }

    if (!graph.empty()) {
        if (quality > 1) {
            $warn << "Warning: input from graph with --qualify can produce unexpected results" << _end$;
        }
        else if (window > 1) {
            $warn << "Warning: input from graph with --window can produce unexpected results" << _end$;
        }
    } else if (quality > 1 && window > 1) {
        $warn << "Warning: using --qualify with --window could produce unexpected results" << _end$;
    }
    if (!splits.empty()) {
        if (!input.empty()) {
            $note << "Note: two input arguments --input and --splits were provided" << _end$;
            $note << "      --input is used for lookup only, no additional splits are inferred" << _end$;
        }
        else if (!newick.empty()) {
            $note << "Note: Newick output from a list of splits, some taxa could be missing" << _end$;
            $note << "      --input can be used to provide the original list of taxa" << _end$;
        }
    }

    // parse the list of input sequence files
    vector<string> files;
    hash_map<string, size1N_t> name_table;
    if (!input.empty()) {
        ifstream file(input);
        if (!file.good()) {
            $err << "Error: could not read input file: " << input << _end$;
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
            if (verbose) $log << end$;
        } else {
            $err << "Error: could not load Bifrost graph: " << graph << _end$;
            return 1;
        }
        if (kmer != 0 && kmer != cdbg.getK()) {
            if (pattern.empty()) {
                $warn << "Warning: graph file does not match the given k-mer length" << _end$;
            } else {
                $err << "Error: graph file does not match the given pattern length" << _end$;
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
            $err << "Error: could not read splits file: " << splits << _end$;
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
        $err << "Error: k-mer length exceeds -DmaxK=" << maxK << ", please see makefile" << _end$;
        return 1;
    }
    if (num > maxN) {
        $err << "Error: number of files exceeds -DmaxN=" << maxN << ", please see makefile" << _end$;
        return 1;
    }

    auto begin = chrono::high_resolution_clock::now();    // time measurement
    graph::init(T, pattern, quality, reverse);    // initialize the graph
    tree::init(top);    // initialize the splits list size
    kmer::init(kmer);    // initialize the k-mer length
    color::init(num);    // initialize the color number

    if (!input.empty() && splits.empty()) {
        if (verbose)
            $log $_ << "Reading input files..." << $;

       auto lambda = [&] (const uint64_t& T, const uint64_t& maxT) {
        for (size1N_t i = T; i < files.size(); i += maxT) {
            ifstream file(files[i]);    // input file stream
            if (!file.good())
                $err $_ << files[i] << " (ERR)" << _end$;    // could not read file
            else if (verbose)
                $log $_ << files[i] << " (" << i+1 << "/" << files.size() << ")" << end$;    // print progress

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
                            $lite $_ << line << _$;    // print progress
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
                $log $_ << $;
            graph::clear_thread(T);
            file.close();
        }};

        vector<thread> thr(T);
        for (uint64_t x = 0; x < T; ++x)
            thr[x] = thread(lambda, x, T);
        for (uint64_t x = 0; x < T; ++x)
            thr[x].join();

        while (thr.size() > 1) {
            for (uint64_t x = 0; x < thr.size()-1; x+=2) {
                thr[x] = thread(graph::merge_threads, x, x+1);
                if (verbose) $log << ":" << $;
            }   if (verbose) $log << "\r" << $;

            for (uint64_t x = 0; x < thr.size()-1; x+=2) {
                thr[x].join();
                if (verbose) $log << "." << $;
            }   if (verbose) $log << "\r" << $;

            uint64_t next_size = thr.size() / 2;
            for (uint64_t x = 1; x <= next_size; ++x) {
                graph::erase_thread(x);
                thr.erase(thr.begin()+x);
                if (verbose) $lite << "." << _$;
            }   if (verbose) $lite $_ << _$;
        }
    }

#ifdef useBF
    if (!graph.empty()) {
        if (verbose)
            $log $_ << "Processing unitigs..." << $;
        uint64_t cur = 0, prog = 0, next;
        uint64_t max = cdbg.size();

        for (auto& unitig : cdbg) {
            if (verbose) {
                next = 100 * cur / max;
                 if (prog < next)  $log $_ << "Processed " << cur << " unitigs (" << next << "%) " << $;
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
            $log $_ << "Processed " << max << " unitigs (100%)" << end$;
        graph::clear_thread(0);
    }
#endif

    if (!splits.empty()) {
        if (verbose)
            $log $_ << "Reading splits file..." << $;
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
            $log $_ << $;
        file.close();
    }

    // function to map a color position to a file name
    function<string(const size1N_t&)> color_name = [&] (const size1N_t& i) {
        if (i < files.size()) return files[i];
       #ifdef useBF
         else return cdbg.getColorName(i-files.size());
       #endif
        $err << "Error: color bit does not correspond to a color name" << _end$;
        exit(1);
    };

    if (interactive) {    // lookup k-mer, interactive mode
        string query;
        $link << ">>> " << _$;    // display a command line prompt
        while (getline(cin, query)) {    // wait for user to enter a k-mer
            $link $_ << _$;
                                            // remove leading whitespaces
            auto l1 = find_if(query.begin(), query.end(), ::isspace);
            auto l2 = find_if_not(query.begin(), query.end(), ::isspace);
             if (l1 < l2) query.erase(l1, l2);    // remove trailing whitespaces
            auto r1 = find_if_not(query.rbegin(), query.rend(), ::isspace);
            auto r2 = find_if(query.rbegin(), query.rend(), ::isspace);
             if (r1.base() < r2.base()) query.erase(r1.base(), r2.base());

            if (query.empty()) { $link << ">>> " << _$; continue; }
            transform(query.begin(), query.end(), query.begin(), ::toupper);
            replace(query.begin(), query.end(), '.', 'N');    // allow . gaps
            replace(query.begin(), query.end(), '-', 'N');    // allow - gaps
            replace(query.begin(), query.end(), '*', 'N');    // allow * gaps

            $lite << "Searching..." << _$;
            auto iterator = graph::lookup_kmer(query);
            $lite $_ << _$;
            string kmer, color;    // iterate over the hash table
            while (iterator(kmer, color)) {
                for (size_t i = 0; i < pattern.length(); ++i)
                    if (pattern[i] == '0') kmer[i] = '-';
                $lite << "... " << _$;
                $out << kmer << ": " << color << end$;
                $lite $_ << _$;
            }
            $link << ">>> " << _$;
        }
        $link $_ << _$;
    }
    if (!counts.empty()) {    // lookup k-mer, output to file
        if (verbose)
            $log $_ << "Please wait... " << $;
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
            $log $_ << "Processing splits..." << $;
        graph::calc_weights(default_mean, verbose);    // accumulate split weights

        if (!filter.empty()) {    // apply filter
            if (verbose)
                $log $_ << "Filtering splits..." << $;
            if (filter == "strict" || filter == "tree")
                tree::filter_strict(verbose);
            else if (filter == "weakly")
                tree::filter_weakly(verbose);
            else if (filter.find("tree") != -1 && filter.substr(filter.find("tree")) == "tree")
                tree::filter_n_tree(stoi(filter.substr(0, filter.find("tree"))), verbose);
        }

        if (verbose)
            $log $_ << "Please wait... " << $;
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
        $log << "Done! (" << util::format_time(end - begin) << ")" << end$;
    return 0;
}
