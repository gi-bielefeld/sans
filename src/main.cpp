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
        $out << end$$;
    }

    file input;    // Input FASTA files: list of sequence files, one per line
    file index;    // Input Index file: load a k-mer index, e.g. counts table
    file graph;    // Input Graph file: load a Bifrost graph, filename prefix
    file splits;   // Input Splits file: load an existing list of splits file

    file output;     // Output TSV file: list of splits, sorted by weight desc.
    file newick;     // Output Newick file: convert splits to a tree topology
    file counts;     // Output K-mer file: list k-mer occurrence per input file
    file diff;       // Print the difference between two index or splits files
    bool query = 0;  // Interactive: list k-mer occurrence per input file

    uint64_t kmer = 0;     // Length of k-mers (default: 31)
    string   pattern;      // Pattern of gapped k-mers (default: no gaps)
    uint64_t window = 1;   // Number of k-mers per minimizer window (default: 1)
    uint64_t iupac = 1;    // Extended IUPAC alphabet, resolve ambiguous bases
    uint64_t quality = 1;  // Discard k-mers with lower coverage than a threshold
    bool     reverse = 0;  // Keep one repr. for reverse complement k-mers

    uint64_t num = 0;      // Number of colors, determined from input files
    uint64_t top = -1;     // Number of splits in the output list (default: all)
    string   filter;       // Output a greedy maximum weight subset of splits
    uint64_t T = -1;       // Number of parallel threads (default: auto)
    bool     verbose = 0;  // Print some information messages during execution
    bool     VERBOSE = 0;  // Print more information messages during execution

    auto arithmetic_mean = [] (const uint32_t& x, const uint32_t& y) { return x / 2.0 + y / 2.0; };
    auto geometric_mean  = [] (const uint32_t& x, const uint32_t& y) { return sqrt(x) * sqrt(y); };
    auto geometric_mean2 = [] (const uint32_t& x, const uint32_t& y) { return sqrt(x+1) * sqrt(y+1) - 1; };
    function<double(const uint32_t&, const uint32_t&)> default_mean  = geometric_mean;

    // parse the command line arguments and update the variables above
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input") == 0)
            input = file(util::atos(argc, argv, ++i), type::INPUT_FILE, mode::READ);
        else if (strcmp(argv[i], "-j") == 0 || strcmp(argv[i], "--index") == 0)
            index = file(util::atos(argc, argv, ++i), type::INDEX_FILE, mode::READ);
        else if (strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--graph") == 0)
            graph = file(util::atos(argc, argv, ++i), type::GRAPH_FILE, mode::READ);
        else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--splits") == 0)
            splits = file(util::atos(argc, argv, ++i), type::SPLITS_FILE, mode::READ);

        else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0)
            output = file(util::atos(argc, argv, ++i), type::SPLITS_FILE, mode::WRITE);
        else if (strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--newick") == 0)
            newick = file(util::atos(argc, argv, ++i), type::NEWICK_FILE, mode::WRITE);
        else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--counts") == 0)
            counts = file(util::atos(argc, argv, ++i), type::INDEX_FILE, mode::WRITE);
        else if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--diff") == 0)
            diff = file(util::atos(argc, argv, ++i), type::INDEX_OR_SPLITS_FILE, mode::READ);
        else if (strcmp(argv[i], "-C") == 0)    /* hidden option */
            query = true;    // Interactive: list k-mer occurrence per input file

        else if (strcmp(argv[i], "-k") == 0 || strcmp(argv[i], "--kmer") == 0)
            kmer = util::aton(argc, argv, ++i);    // Length of k-mers (default: 31)
        else if (strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "--gapped") == 0)
            pattern = util::atos(argc, argv, ++i);    // Pattern of gapped k-mers (default: no gaps)
        else if (strcmp(argv[i], "-w") == 0 || strcmp(argv[i], "--window") == 0)
            window = util::aton(argc, argv, ++i);    // Number of k-mers per window (default: 1)
        else if (strcmp(argv[i], "-x") == 0 || strcmp(argv[i], "--iupac") == 0)
            iupac = util::aton(argc, argv, ++i);    // Extended IUPAC alphabet, resolve ambiguous bases
        else if (strcmp(argv[i], "-q") == 0 || strcmp(argv[i], "--qualify") == 0)
            quality = util::aton(argc, argv, ++i);    // Discard k-mers below a min. coverage threshold
        else if (strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--reverse") == 0)
            reverse = true;    // Keep one repr. for reverse complement k-mers

        else if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--top") == 0)
            top = util::aton(argc, argv, ++i);    // Number of top splits to output (default: all)
        else if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--mean") == 0) {
            string mean = util::atos(argc, argv, ++i);    // Mean weight function to handle asymmetric splits
            if      (mean == "arith") default_mean = arithmetic_mean;
            else if (mean == "geom")  default_mean = geometric_mean;
            else if (mean == "geom2") default_mean = geometric_mean2;
            else {
                $err << "Error: unknown argument: " << argv[i-1] << ' ' << mean << _end$;
                $err << "       type --help to see a list of supported options" << _end$$;
           }}
        else if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--filter") == 0) {
            filter = util::atos(argc, argv, ++i);    // Output a greedy maximum weight subset of splits
            if      (filter == "strict" || filter == "tree"); // compatible to a tree
            else if (filter == "weakly");                     // weakly compatible network
            else if (filter.find("tree") != -1 && filter.substr(filter.find("tree")) == "tree")
                util::ston(argv[i-1], filter);
            else {
                $err << "Error: unknown argument: " << argv[i-1] << ' ' << filter << _end$;
                $err << "       type --help to see a list of supported options" << _end$$;
           }}

        else if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--threads") == 0)
            T = util::aton(argc, argv, ++i);    // Number of parallel threads (default: auto)
        else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0)
            verbose = true;    // Print information messages during execution
        else if (strcmp(argv[i], "-rv") == 0 || strcmp(argv[i], "-vr") == 0)
            reverse = true,    // Keep one repr. for reverse complement k-mers
            verbose = true;    // Print information messages during execution
        else if (strcmp(argv[i], "-V") == 0)    /* hidden option */
            verbose = true,    // Print some information messages during execution
            VERBOSE = true;    // Print more information messages during execution
        else if (!(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)) {
            $err << "Error: unknown argument: " << argv[i] << _end$;
            $err << "       type --help to see the full list of parameters" << _end$$;
        }
    }

    if (input.empty() && index.empty() && graph.empty() && splits.empty())
        $err << "Error: missing argument: --input, --index, --graph, or --splits <file_name>" << _end$$;
    if (!index.empty() && !graph.empty() && !splits.empty())
        $err << "Error: too many input arguments: --index, --graph, and --splits <file_name>" << _end$$;
    if (!index.empty() && !graph.empty())
        $err << "Error: too many input arguments: --index and --graph <file_name>" << _end$$;
    if (!index.empty() && !splits.empty())
        $err << "Error: too many input arguments: --index and --splits <file_name>" << _end$$;
    if (!graph.empty() && !splits.empty())
        $err << "Error: too many input arguments: --graph and --splits <file_name>" << _end$$;

    if (output.empty() && newick.empty() && (counts.empty() && !query) && diff.empty())
        $err << "Error: missing argument: --output, --newick, --counts, or --diff <file_name>" << _end$$;
    if (!output.empty() && !newick.empty() && !(diff.empty() || diff.is(type::SPLITS_FILE)))
        $err << "Error: too many output arguments: --output, --newick, and --diff <file_name>" << _end$$;
    if (!output.empty() && !(diff.empty() || diff.is(type::SPLITS_FILE)))
        $err << "Error: too many output arguments: --output and --diff <file_name>" << _end$$;
    if (!newick.empty() && !(diff.empty() || diff.is(type::SPLITS_FILE)))
        $err << "Error: too many output arguments: --newick and --diff <file_name>" << _end$$;

    if (!splits.empty() && !(counts.empty() && !query))
        $err << "Error: k-mer counts cannot be calculated if the input is a list of splits" << _end$$;
    if (!splits.empty() && !(diff.empty() || diff.is(type::SPLITS_FILE)))
        $err << "Error: k-mer counts cannot be calculated if the input is a list of splits" << _end$$;
    if (!newick.empty() && filter != "strict" && filter.find("tree") == -1)
        $err << "Error: Newick output only applicable in combination with -f strict/n-tree" << _end$$;

    if (!pattern.empty()) {
        uint64_t lmer = 0;    // length without gaps
        for (size_t i = 0; i < pattern.length(); ++i) {
            if (pattern[i] == '1')
               { ++lmer; continue; }
            if (pattern[i] != '0')
                $err << "Error: pattern must be a sequence of 0s and 1s, where 0 means a gap" << _end$$;
            if (reverse && pattern[i] != pattern[pattern.length()-1-i])
                $err << "Error: pattern should be symmetric to work with reverse complements" << _end$$;
        }
        if (kmer != 0 && kmer != lmer && kmer != pattern.length())
            $err << "Error: pattern length does not match the given k-mer length" << _end$$;
        kmer = pattern.length();

        if (lmer == 0)
            $err << "Error: pattern must be a sequence of 0s and 1s, with at least one 1" << _end$$;
        if (lmer == pattern.length())
            pattern.clear();
    }

    if (!index.empty()) {
        if (!input.empty() && verbose)
            $note << "Note: --input is used for lookup only, no additional k-mers are counted" << _end$;
        if (input.empty() && !output.empty())
            $warn << "Warning: splits output from a k-mer index, color names might be missing" << _end$,
            $warn << "         --input can be used to provide the original list of files" << _end$;
        if (input.empty() && output.empty() && !newick.empty())
            $warn << "Warning: Newick output from a k-mer index, color names might be missing" << _end$,
            $warn << "         --input can be used to provide the original list of files" << _end$;
    }
    if (!splits.empty()) {
        if (!input.empty() && verbose)
            $note << "Note: --input is used for lookup only, no additional splits are inferred" << _end$;
        if (input.empty() && !newick.empty())
            $warn << "Warning: Newick output from a list of splits, some taxa could be missing" << _end$,
            $warn << "         --input can be used to provide the original list of files" << _end$;
    }

    // parse the list of input sequence files
    vector<string> color_name;
    hash_map<string, size1N_t> color_index;
    if (!input.empty()) {
        ifstream file(input);
        if (!file.good()) {
            $err << "Error: could not read input file: " << input << _end$$;
        }
        string line;
        while (getline(file, line)) {
            if (index.empty() && splits.empty()) {
                ifstream fastx(line);
                if (!fastx.good()) {
                    $err << "Error: could not read sequence file: " << line << _end$$;
                }
                fastx.close();
            }
            color_name.emplace_back(line);
            color_index[line] = num++;
        }
        file.close();
    }
    uint64_t file_num = color_name.size();

#ifdef useBF
    // load an existing Bifrost graph file
    ColoredCDBG<> cdbg(kmer);
    if (!graph.empty()) {
        if (cdbg.read(graph + ".gfa", graph + ".bfg_colors", 1, verbose)) {
            if (verbose) $log << end$;
        } else {
            $err << "Error: could not load Bifrost graph: " << graph << _end$$;
        }
        if (kmer != 0 && kmer != cdbg.getK()) {
            if (pattern.empty()) {
                $warn << "Warning: graph file does not match the given k-mer length" << _end$;
            } else {
                $err << "Error: graph file does not match the given pattern length" << _end$$;
            }
        }
        for (uint64_t i = 0; i < cdbg.getNbColors(); ++i) {
            const string& name = cdbg.getColorName(i);
            if (color_index.find(name) == color_index.end()) {
                color_name.emplace_back(name);
                color_index[name] = num++;
            }
        }
        kmer = cdbg.getK();
    }
#endif

    // load & parse an existing list of splits
    if (!splits.empty()) {
        ifstream file(splits);
        if (!file.good()) {
            $err << "Error: could not read splits file: " << splits << _end$$;
        }
        string line; size_t curr, next;
        while (getline(file, line)) {
            curr = line.find('\t');
            next = curr + 1;
            do {
                curr = line.find('\t', next);
                const string& name = line.substr(next, curr-next);
                if (color_index.find(name) == color_index.end()) {
                    color_name.emplace_back(name);
                    color_index[name] = num++;
                }
                next = curr + 1;
            } while (curr != string::npos);
        }
        file.close();
    }

    if (kmer == 0)
        kmer = min<uint64_t>(31, maxK);    // k-mer length not specified & not defined by pattern
    if (kmer > maxK)
        $err << "Error: k-mer length exceeds -DmaxK=" << maxK << ", please see makefile" << _end$$;
    if (num > maxN)
        $err << "Error: number of files exceeds -DmaxN=" << maxN << ", please see makefile" << _end$$;

    uint64_t maxT = thread::hardware_concurrency()-1; T = min<uint64_t>({T, num, maxT});

    auto begin = chrono::high_resolution_clock::now();  // time measurement
    color::init(num);  // initialize the color number
    kmer::init(kmer, pattern);  // initialize the k-mer length
    tree::init(top);  // initialize the splits list size
    graph::init(T, quality, reverse);  // initialize the graph

    tree::color_name = [&] (const size1N_t& index)
     { return color_name[index]; };  // map color position to file name
    tree::color_index = [&] (const string& name)
     { return color_index[name]; };  // map file name to color position

    if (!input.empty() && index.empty() && splits.empty()) {
        if (T <= 0 && verbose) $note << "Note: running SANS-KC using 1(+0) threads" << _end$;
        if (T >= 1 && verbose) $note << "Note: running SANS-KC using " << T << "(+1) threads" << _end$;
        if (T >= 2 && quality > 1)
            $warn << "Warning: using --qualify on multiple threads could increase memory usage" << _end$,
            $warn << "         In case of problems, please specify a smaller number of threads" << _end$;

        vector<thread> list_thread(T); atomic<uint64_t> active(T);    // parallel file reading threads
             if (VERBOSE) $log $_ << "Reading input files..." << $;
        else if (verbose) $log $_ << "Reading input files..." << end$;

        auto lambda = [&] (const uint64_t& T, const uint64_t& maxT) {
            for (size1N_t i = T; i < file_num; i += maxT) {
                ifstream file(color_name[i]);    // input file stream & print progress
                     if (VERBOSE) $log $_ << color_name[i] << " (" << i+1 << "/" << file_num << ")" << end$;
                else if (verbose) $log $_ << color_name[i] << " (" << i+1 << "/" << file_num << ")" << $;

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

                                 if (VERBOSE) $lite $_ << line << _$;    // print more fine-grained progress
                            else if (verbose) $lite $_ << color_name[i] << " (" << i+1 << "/" << file_num << ")" << _$;
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
                file.close();

                graph::clear_thread(T);
            } --active;
        };

        if (T <= 0) lambda(0, 1);
        else {
            for (uint64_t x = 0; x < T; ++x)
                list_thread[x] = thread(lambda, x, T);
            while (active) {   // wait for threads
                graph::merge_threads();
                this_thread::sleep_for(1ns);
            }   graph::merge_threads();
            for (uint64_t x = 0; x < T; ++x)
                list_thread[x].join();
        }
    }
    graph::erase_threads();

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
                const string& name = cdbg.getColorName(it.getColorID());
                graph::add_kmers(0, substr, color_index[name]);
            }
        }
        if (verbose)
            $log $_ << "Processed " << max << " unitigs (100%)" << end$;
    }
#endif

    if (!splits.empty()) {
        if (verbose)
            $log $_ << "Reading splits file..." << $;
        ifstream file(splits);

        string line; size_t curr, next;
        double weight; color_t color;
        while (getline(file, line)) {
            curr = line.find('\t');
            weight = stod(line.substr(0, curr));
            color = 0b0u;
            next = curr + 1;
            do {
                curr = line.find('\t', next);
                const string& name = line.substr(next, curr-next);
                color.set(color_index[name]);
                next = curr + 1;
            } while (curr != string::npos);

            tree::insert_split(weight, color);
        }
        if (verbose)
            $log $_ << $;
        file.close();
    }

    if (query) {    // lookup k-mer, interactive query mode
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
                tree::filter_n_tree(stol(filter), verbose);
        }

        if (verbose)
            $log $_ << "Please wait... " << $;
        if (!output.empty()) {
            ofstream file(output);    // output file stream
            ostream stream(file.rdbuf());
            for (auto& split : tree::splits) {
                stream << split.first;    // weight of the split
                for (size1N_t i = 0; i < num; ++i) {
                    if (split.second.test(i)) {
                        stream << '\t' << color_name[i];    // name of the file
                    }
                }
                stream << endl;
            }
            file.close();
        }
        if (!newick.empty()) {
            ofstream file(newick);    // output file stream
            ostream stream(file.rdbuf());
            stream << tree::build_string();    // filter and print tree
            file.close();
        }
    }

    auto end = chrono::high_resolution_clock::now();  // time measurement
    if (verbose)  // print progress and time
        $log << "Done! (" << util::format_time(end - begin) << ")" << end$;
}
