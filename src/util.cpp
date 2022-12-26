#include "util.h"
#include "ansi.h"

/**
 * This function prints a minimal version of the help page.
 */
void util::print_help() {
    $out << end$;
    $out << "  Input arguments:" << end$;
    $out << end$;
    $out << "    -i, --input   \t Input FASTA files: list of sequence files, one per line" << end$;
    $out << "    -j, --index   \t Input Index file: load a k-mer index, e.g. counts table" << end$;
    $out << "    -g, --graph   \t Input Graph file: load a Bifrost graph, filename prefix" << end$;
    $out << "    -s, --splits  \t Input Splits file: load an existing list of splits file" << end$;
    $out << end$;
    $out << "  Output arguments:" << end$;
    $out << end$;
    $out << "    -o, --output  \t Output TSV file: list of splits, sorted by weight desc." << end$;
    $out << "    -n, --newick  \t Output Newick file: convert splits to a tree topology" << end$;
    $out << "    -c, --counts  \t Output K-mer file: list k-mer occurrence per input file" << end$;
    $out << "    -d, --diff    \t Print the difference between two index or splits files" << end$;
    $out << end$;
    $out << "  K-mer options:" << end$;
    $out << end$;
    $out << "    -k, --kmer    \t Length of k-mers (default: 31)" << end$;
    $out << "    -l, --gapped  \t Pattern of gapped k-mers (default: no gaps)" << end$;
    $out << "    -w, --window  \t Number of k-mers per minimizer window (default: 1)" << end$;
    $out << "    -x, --iupac   \t Extended IUPAC alphabet, resolve ambiguous bases" << end$;
    $out << "    -q, --qualify \t Discard k-mers with lower coverage than a threshold" << end$;
    $out << "    -r, --reverse \t Keep one repr. for reverse complement k-mers" << end$;
    $out << end$;
    $out << "  Filter options:" << end$;
    $out << end$;
    $out << "    -t, --top     \t Number of splits in the output list (default: all)" << end$;
    $out << "    -m, --mean    \t Mean weight function to handle asymmetric splits" << end$;
    $out << "    -f, --filter  \t Output a greedy maximum weight subset of splits" << end$;
    $out << end$;
    $out << "  Other settings:" << end$;
    $out << end$;
    $out << "    -p, --threads \t Number of parallel threads (default: auto)" << end$;
    $out << "    -v, --verbose \t Print information messages during execution" << end$;
    $out << "    -h, --help    \t Display an extended help page and quit" << end$;
    $out << end$;
}

/**
 * This function prints an extended version of the help page.
 */
void util::print_extended_help() {
    $out << " ___________________________________________________________________________" << end$;
    $out << end$;
    $out << "  Input arguments:" << end$;
    $out << end$;
    $out << "    -i, --input   \t Input FASTA files: list of sequence files, one per line" << end$;
    $out << end$;
    $out << "    -j, --index   \t Input Index file: load a k-mer index, e.g. counts table" << end$;
    $out << "                  \t (provide list -i to lookup names or extend the index)" << end$;
    $out << end$;
    $out << "    -g, --graph   \t Input Graph file: load a Bifrost graph, filename prefix" << end$;
    $out << "                  \t (requires compiler flag -DuseBF, please see makefile)" << end$;
    $out << end$;
    $out << "    -s, --splits  \t Input Splits file: load an existing list of splits file" << end$;
    $out << "                  \t (allows to filter -t/-f, other arguments are ignored)" << end$;
    $out << end$;
    $out << "    (either --input, --index, --graph, or --splits must be provided)" << end$;
    $out << " ___________________________________________________________________________" << end$;
    $out << end$;
    $out << "  Output arguments:" << end$;
    $out << end$;
    $out << "    -o, --output  \t Output TSV file: list of splits, sorted by weight desc." << end$;
    $out << end$;
    $out << "    -n, --newick  \t Output Newick file: convert splits to a tree topology" << end$;
    $out << "                  \t (only applicable in combination with -f strict/n-tree)" << end$;
    $out << end$;
    $out << "    -c, --counts  \t Output K-mer file: list k-mer occurrence per input file" << end$;
    $out << "                  \t (cannot be calculated if the input is a list of splits)" << end$;
    $out << end$;
    $out << "    -d, --diff    \t Print the difference between two index or splits files" << end$;
    $out << "                  \t (input a file with matching format, output to stdout)" << end$;
    $out << end$;
    $out << "    (either --output, --newick, --counts, or --diff must be provided)" << end$;
    $out << " ___________________________________________________________________________" << end$;
    $out << end$;
    $out << "  K-mer options:" << end$;
    $out << end$;
    $out << "    -k, --kmer    \t Length of k-mers (default: 31)" << end$;
    $out << end$;
    $out << "    -l, --gapped  \t Pattern of gapped k-mers (default: no gaps)" << end$;
    $out << "                  \t e.g. 1110111, where 0 means a gap/skip base" << end$;
    $out << end$;
    $out << "    -w, --window  \t Number of k-mers per minimizer window (default: 1)" << end$;
    $out << "                  \t Consider only the smallest k-mer among w neighbors" << end$;
    $out << "                  \t Reduces memory usage, but also lowers the accuracy" << end$;
    $out << end$;
    $out << "    -x, --iupac   \t Extended IUPAC alphabet, resolve ambiguous bases" << end$;
    $out << "                  \t Specify a number to limit the k-mers per position" << end$;
    $out << "                  \t between 1 (no ambiguity) and 4^k (allows NNN...N)" << end$;
    $out << end$;
    $out << "    -q, --qualify \t Discard k-mers with lower coverage than a threshold" << end$;
    $out << "                  \t Only keep k-mers with multiple occurrences per file" << end$;
    $out << end$;
    $out << "    -r, --reverse \t Keep one repr. for reverse complement k-mers" << end$;
    $out << "                  \t Reverse and consider the lex. smaller of both" << end$;
    $out << " ___________________________________________________________________________" << end$;
    $out << end$;
    $out << "  Filter options:" << end$;
    $out << end$;
    $out << "    -t, --top     \t Number of splits in the output list (default: all)" << end$;
    $out << end$;
    $out << "    -m, --mean    \t Mean weight function to handle asymmetric splits" << end$;
    $out << "                  \t options: arith: arithmetic mean" << end$;
    $out << "                  \t          geom:  geometric mean (default)" << end$;
    $out << "                  \t          geom2: geometric mean with pseudo-counts" << end$;
    $out << end$;
    $out << "    -f, --filter  \t Output a greedy maximum weight subset of splits" << end$;
    $out << "                  \t options: strict: compatible to a tree" << end$;
    $out << "                  \t          weakly: weakly compatible network" << end$;
    $out << "                  \t          n-tree: compatible to a union of n trees" << end$;
    $out << "                  \t                  (where n is an arbitrary number)" << end$;
    $out << " ___________________________________________________________________________" << end$;
    $out << end$;
    $out << "  Other settings:" << end$;
    $out << end$;
    $out << "    -p, --threads \t Number of parallel threads (default: auto)" << end$;
    $out << "                  \t Reduces runtime, but could increase memory usage" << end$;
    $out << "                  \t In case of problems, please try a smaller number" << end$;
    $out << end$;
    $out << "    -v, --verbose \t Print information messages during execution" << end$;
    $out << end$;
    $out << "    -h, --help    \t Display an extended help page and quit" << end$;
    $out << " ___________________________________________________________________________" << end$;
    $out << end$;
}

/**
 * This function converts a command line argument to a string.
 *
 * @param argc argument count
 * @param argv argument vector
 * @param i argument index
 * @return string
 */
string util::atos(int& argc, char* argv[], int& i) {
    if (i < argc) {
        return argv[i];
    } else {
        $err << "Error: incomplete argument: " << argv[i-1] << " (expecting a string)" << _end$$;
    }
}

/**
 * This function converts a command line argument to a number.
 *
 * @param argc argument count
 * @param argv argument vector
 * @param i argument index
 * @return number
 */
uint64_t util::aton(int& argc, char* argv[], int& i) {
    if (i < argc) {
        try {
            return stoull(argv[i]);
        } catch (...) {
            $err << "Error: malformed argument: " << argv[i-1] << ' ' << argv[i] << " (expected a number)" << _end$$;
        }
    } else {
        $err << "Error: incomplete argument: " << argv[i-1] << " (expecting a number)" << _end$$;
    }
}

/**
 * This function converts a string argument to a number.
 *
 * @param param parameter
 * @param args argument
 * @return number
 */
uint64_t util::ston(const string& param, const string& args) {
    try {
        return stoull(args);
    } catch (...) {
        $err << "Error: malformed argument: " << param << ' ' << args << " (expected a number)" << _end$$;
    }
}

/**
 * This function displays a duration in a human readable format.
 *
 * @param time duration
 * @return formatted string
 */
string util::format_time(const chrono::high_resolution_clock::duration& time) {
    double value = chrono::duration_cast<chrono::milliseconds>(time).count();
    string unit = "ms";
    if (value >= 1000) {
        value /= 1000; unit = "sec";
        if (value >= 200) {
            value /= 60; unit = "min";
            if (value >= 200) {
                value /= 60; unit = "h";
            }
        }
    }
    string number = to_string(value);
    return number.substr(0, number.find('.')+2) + ' ' + unit;
}
