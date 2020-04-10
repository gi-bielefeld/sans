#include <bifrost/CompactedDBG.hpp>
#include <bifrost/ColoredCDBG.hpp>
#include "SansOpt.h"
#include "TopSplits.cpp"

#define SANS_VERSION "1.1"

using namespace std;



// use:  PrintVersion();
// post: The version of the program has been printed to cout
void PrintVersion() {
    cout << SANS_VERSION << endl;
}

// use:  PrintCite();
// post: Information of how to cite this software has been printed to cerr
void PrintCite() {
    cout << "The paper describing this software has been published at WABI 2019: " << endl;
    cout << "Wittler, R.: Alignment- and reference-free phylogenomics with colored de Bruijn graphs. In: Huber, K. and Gusfield, D. (eds.) Proceedings of WABI 2019. LIPIcs. 143, Schloss Dagstuhl--Leibniz-Zentrum fuer Informatik, Dagstuhl, Germany (2019)." << endl;
}

void PrintUsage() {
    cout << endl;
    cout << "SANS " << SANS_VERSION << endl;
    cout << "Symmetric Alignment-free phylogeNomic Splits" << endl;
    cout << endl;
    cout << "Usage: SANS [PARAMETERS]" << endl;
    cout << endl;
    cout << " > Mandatory with required argument:" << endl;
    cout << endl;
    cout << "  -s, --input-seq-files   Input sequence files (FASTA/FASTQ possibly gzipped)" << endl;
    cout << "                          Input files can be provided as a list in a TXT file (one file per line)" << endl;
    cout << "                          K-mers with exactly 1 occurrence in the input files will be discarded" << endl;
    cout << "  -r, --input-ref-files   Input reference files (FASTA/FASTQ possibly gzipped and GFA)" << endl;
    cout << "                          Input files can be provided as a list in a TXT file (one file per line)" << endl;
    cout << "                          All k-mers of the input reference files are used" << endl;
    cout << "  -o, --output-file       name of output file" << endl;
    cout << endl;
    cout << " > Optional with required argument:" << endl;
    cout << endl;
    cout << "  -t, --threads           Number of threads (default: 1)" << endl;
    cout << "  -T, --top               Output the top T splits sorted by weight descending (default: all)" << endl;
    cout << "  -k, --kmer-length       Length of k-mers (default: 31)" << endl;
    cout << "  -m, --min-length        Length of minimizers (auto-adjusted by default: see verbose output)" << endl;
    cout << "  -b, --bloom-bits        Number of Bloom filter bits per k-mer with 1+ occurrences (default: 14)" << endl;
    cout << "  -B, --bloom-bits2       Number of Bloom filter bits per k-mer with 2+ occurrences (default: 14)" << endl;
    cout << "  -l, --load-mbbf         Input Blocked Bloom Filter file, skips filter step (default: no input)" << endl;
    cout << "  -w, --write-mbbf        Output Blocked Bloom Filter file (default: no output)" << endl;
    cout << "  -u, --chunk-size        Read chunk size per thread (default: 64)" << endl;
    cout << "  -f, --filter            Output a greedy maximum weight subset" << endl;
    cout << "                          options: 1-tree: compatible to a tree" << endl;
    cout << "                                   2-tree: compatible to union of two trees (network)" << endl;
    cout << endl;
    cout << " > Optional with no argument:" << endl;
    cout << endl;
    cout << "  -i, --clip-tips         Clip tips shorter than k k-mers in length" << endl;
    cout << "  -d, --del-isolated      Delete isolated contigs shorter than k k-mers in length" << endl;
    cout << "  -y, --keep-mercy        Keep low coverage k-mers connecting tips" << endl;
    cout << "  -a, --allow-asym        Do not discard asymmetric splits completely" << endl;
    cout << "  -S, --output-sequences  Output the conserved subsequences the splits are derived from" << endl;
    cout << "  -v, --verbose           Print information messages during execution" << endl;
    cout << endl;
}

void parse_ProgramOptions(int argc, char **argv, SANS_opt& opt) {

    int option_index = 0, c;
    const char *opt_string = "s:r:o:t:T:k:m:b:B:l:w:u:f:idyaSv";

    static struct option long_options[] = {
            {"input-seq-files",  required_argument, 0, 's'},
            {"input-ref-files",  required_argument, 0, 'r'},
            {"output-file",      required_argument, 0, 'o'},
            {"threads",          required_argument, 0, 't'},
            {"top",              required_argument, 0, 'T'},
            {"kmer-length",      required_argument, 0, 'k'},
            {"min-length",       required_argument, 0, 'm'},
            {"bloom-bits",       required_argument, 0, 'b'},
            {"bloom-bits2",      required_argument, 0, 'B'},
            {"load-mbbf",        required_argument, 0, 'l'},
            {"write-mbbf",       required_argument, 0, 'w'},
            {"chunk-size",       required_argument, 0, 'u'},
            {"filter",           required_argument, 0, 'f'},
            {"clip-tips",        no_argument,       0, 'i'},
            {"del-isolated",     no_argument,       0, 'd'},
            {"keep-mercy",       no_argument,       0, 'y'},
            {"allow-asym",       no_argument,       0, 'a'},
            {"output-sequences", no_argument,       0, 'S'},
            {"verbose",          no_argument,       0, 'v'},
            {0, 0,                                  0, 0}
    };

    opt.build = true;
    opt.outputColors = true;
    opt.g=0;

    while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {
        switch (c) {
            case 's':
                for (--optind; (optind < argc) && (*argv[optind] != '-'); ++optind) {
                    opt.filename_seq_in.push_back(argv[optind]);
                }
                break;
            case 'r':
                for (--optind; (optind < argc) && (*argv[optind] != '-'); ++optind) {
                    opt.filename_ref_in.push_back(argv[optind]);
                }
                break;
            case 'o':
                opt.prefixFilenameOut = optarg;
                break;
            case 't':
                opt.nb_threads = atoi(optarg);
                break;
            case 'T':
                opt.top_splits = atoi(optarg);
                break;
            case 'k':
                opt.k = atoi(optarg);
                break;
            case 'm':
                opt.g = atoi(optarg);
                break;
            case 'b':
                opt.nb_bits_unique_kmers_bf = atoi(optarg);
                break;
            case 'B':
                opt.nb_bits_non_unique_kmers_bf = atoi(optarg);
                break;
            case 'l':
                opt.inFilenameBBF = optarg;
                break;
            case 'w':
                opt.outFilenameBBF = optarg;
                break;
            case 'u':
                opt.read_chunksize = atoi(optarg);
                break;
            case 'f':
                opt.filter = optarg;
                break;
            case 'i':
                opt.clipTips = true;
                break;
            case 'd':
                opt.deleteIsolated = true;
                break;
            case 'y':
                opt.useMercyKmers = true;
                break;
            case 'a':
                opt.allow_asym = true;
                break;
            case 'S':
                opt.output_sequences = true;
                break;
            case 'v':
                opt.verbose = true;
                break;
            default: break;
        }
    }
}


bool check_ProgramOptions(SANS_opt& opt) {

    bool ret = true;

    size_t max_threads = std::thread::hardware_concurrency();

    auto check_files = [&](vector<string>& v_files) {

        vector<string> files_tmp;

        char* buffer = new char[4096]();

        for (const auto& file : v_files) {

            if (!check_file_exists(file)) {

                cerr << "Error: File " << file << " not found." << endl;
                ret = false;
            }
            else {

                const string s_ext = file.substr(file.find_last_of('.') + 1);

                if ((s_ext == "txt")){

                    FILE* fp = fopen(file.c_str(), "r");

                    if (fp != nullptr) {

                        fclose(fp);

                        ifstream ifs_file_txt(file);
                        istream i_file_txt(ifs_file_txt.rdbuf());

                        while (i_file_txt.getline(buffer, 4096)){

                            fp = fopen(buffer, "r");

                            if (fp == nullptr) {

                                cerr << "Error: Could not open file " << buffer << " for reading." << endl;
                                ret = false;
                            }
                            else {

                                fclose(fp);
                                files_tmp.push_back(string(buffer));
                            }
                        }

                        ifs_file_txt.close();
                    }
                    else {

                        cerr << "Error: Could not open file " << file << " for reading." << endl;
                        ret = false;
                    }
                }
                else files_tmp.push_back(file);
            }
        }

        v_files = move(files_tmp);

        delete[] buffer;
    };

    // Check general parameters

    if (opt.nb_threads <= 0){

        cerr << "Error: Number of threads cannot be less than or equal to 0." << endl;
        ret = false;
    }

    if (opt.nb_threads > max_threads){

        cerr << "Error: Number of threads cannot be greater than or equal to " << max_threads << "." << endl;
        ret = false;
    }

    if (opt.k <= 3){

        cerr << "Error: Length k of k-mers cannot be less than or equal to 3." << endl;
        ret = false;
    }

    if (opt.k >= MAX_KMER_SIZE){

        cerr << "Error: Length k of k-mers cannot exceed or be equal to " << MAX_KMER_SIZE << "." << endl;
        ret = false;
    }

    // set g appropriately
    if (opt.g == 0) {
        if (opt.k >= 29) opt.g = 23;
        else if (opt.k > 14) opt.g = opt.k - 8;
        else opt.g = opt.k/2;
        if (opt.verbose) {
            cout << "Length of minimizer auto-adjusted to " << opt.g << "." << endl;
        }

    }

    if (opt.g <= 0){

        cerr << "Error: Length m of minimizers cannot be less than or equal to 0." << endl;
        ret = false;
    }

    if (opt.g > opt.k - 2){

        cerr << "Error: Length m of minimizers cannot exceed k - 2 (" << (opt.k - 2) << ")." << endl;
        ret = false;
    }

    const string out = opt.prefixFilenameOut;

    FILE* fp = fopen(out.c_str(), "w");

    if (fp == nullptr) {

        cerr << "Error: Could not open file for writing output graph in GFA format: " << out << "." << endl;
        ret = false;
    }
    else {

        fclose(fp);
        if (remove(out.c_str()) != 0) cerr << "Error: Could not remove temporary file " << out << "." << endl;
    }

    if ((opt.filename_seq_in.size() + opt.filename_ref_in.size()) == 0) {

        cerr << "Error: Missing input files." << endl;
        ret = false;
    }
    else {

        check_files(opt.filename_seq_in);
        check_files(opt.filename_ref_in);
    }


    if (opt.read_chunksize <= 0) {

        cerr << "Error: Chunk size of reads to share among threads cannot be less than or equal to 0." << endl;
        ret = false;
    }

    if (opt.outFilenameBBF.length() != 0){

        FILE* fp = fopen(opt.outFilenameBBF.c_str(), "wb");

        if (fp == nullptr) {

            cerr << "Error: Could not open Blocked Bloom filter file " << opt.outFilenameBBF << " for writing." << endl;
            ret = false;
        }
        else {

            fclose(fp);

            if (remove(opt.outFilenameBBF.c_str()) != 0){

                cerr << "Error: Could not remove temporary file " << opt.outFilenameBBF << "." << endl;
            }
        }
    }

    if (opt.inFilenameBBF.length() != 0){

        if (check_file_exists(opt.inFilenameBBF)){

            FILE* fp = fopen(opt.inFilenameBBF.c_str(), "rb");

            if (fp == nullptr) {

                cerr << "Error: Could not read input Blocked Bloom filter file " << opt.inFilenameBBF << "." << endl;
                ret = false;
            }
            else fclose(fp);
        }
        else {

            cerr << "Error: Input Blocked Bloom filter " << opt.inFilenameBBF << " file does not exist." << endl;
            ret = false;
        }
    }

    return ret;
}


int main(int argc, char **argv){

    if (argc < 2) PrintUsage();
    else {

        SANS_opt opt;


        parse_ProgramOptions(argc, argv, opt); // Parse input parameters

        if (!check_ProgramOptions(opt)) return 0; // Check if input parameters are valid

        ColoredCDBG<> cdbg(opt.k, opt.g);

        cdbg.buildGraph(opt);
        cdbg.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);
        cdbg.buildColors(opt);

        if (opt.verbose) {
            cout << "SANS.cpp(): traversing graph and extracting..." << endl;
        }

        searchGraph(cdbg, opt);

    }
}
