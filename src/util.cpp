#include "util.h"

/**
 * This function prints a minimal version of the help page.
 */
void util::print_help() {
    cout << endl;
    cout << "  Input arguments:" << endl;
    cout << endl;
    cout << "    -i, --input   \t Input FASTA files: list of sequence files, one per line" << endl;
    cout << "    -g, --graph   \t Input Graph file: load a Bifrost graph, filename prefix" << endl;
    cout << "    -s, --splits  \t Input Splits file: load an existing list of splits file" << endl;
    cout << endl;
    cout << "  Output arguments:" << endl;
    cout << endl;
    cout << "    -o, --output  \t Output TSV file: list of splits, sorted by weight desc." << endl;
    cout << "    -n, --newick  \t Output Newick file: convert splits to a tree topology" << endl;
    cout << "    -c, --counts  \t Output K-mer file: list k-mer occurrence per input file" << endl;
    cout << endl;
    cout << "  K-mer options:" << endl;
    cout << endl;
    cout << "    -k, --kmer    \t Length of k-mers (default: 31)" << endl;
    cout << "    -l, --gapped  \t Pattern of gapped k-mers (default: no gaps)" << endl;
    cout << "    -w, --window  \t Number of k-mers per minimizer window (default: 1)" << endl;
    cout << "    -x, --iupac   \t Extended IUPAC alphabet, resolve ambiguous bases" << endl;
    cout << "    -q, --qualify \t Discard k-mers with lower coverage than a threshold" << endl;
    cout << "    -r, --reverse \t Keep one repr. for reverse complement k-mers" << endl;
    cout << endl;
    cout << "  Filter options:" << endl;
    cout << endl;
    cout << "    -t, --top     \t Number of splits in the output list (default: all)" << endl;
    cout << "    -m, --mean    \t Mean weight function to handle asymmetric splits" << endl;
    cout << "    -f, --filter  \t Output a greedy maximum weight subset of splits" << endl;
    cout << endl;
    cout << "  Other settings:" << endl;
    cout << endl;
    cout << "    -p, --threads \t Number of parallel threads (default: 1)" << endl;
    cout << "    -v, --verbose \t Print information messages during execution" << endl;
    cout << "    -h, --help    \t Display an extended help page and quit" << endl;
    cout << endl;
}

/**
 * This function prints an extended version of the help page.
 */
void util::print_extended_help() {
    cout << " ___________________________________________________________________________" << endl;
    cout << endl;
    cout << "  Input arguments:" << endl;
    cout << endl;
    cout << "    -i, --input   \t Input FASTA files: list of sequence files, one per line" << endl;
    cout << endl;
    cout << "    -g, --graph   \t Input Graph file: load a Bifrost graph, filename prefix" << endl;
    cout << "                  \t (requires compiler flag -DuseBF, please see makefile)" << endl;
    cout << endl;
    cout << "    -s, --splits  \t Input Splits file: load an existing list of splits file" << endl;
    cout << "                  \t (allows to filter -t/-f, other arguments are ignored)" << endl;
    cout << endl;
    cout << "    (either --input and/or --graph, or --splits must be provided)" << endl;
    cout << " ___________________________________________________________________________" << endl;
    cout << endl;
    cout << "  Output arguments:" << endl;
    cout << endl;
    cout << "    -o, --output  \t Output TSV file: list of splits, sorted by weight desc." << endl;
    cout << endl;
    cout << "    -n, --newick  \t Output Newick file: convert splits to a tree topology" << endl;
    cout << "                  \t (only applicable in combination with -f strict/n-tree)" << endl;
    cout << endl;
    cout << "    -c, --counts  \t Output K-mer file: list k-mer occurrence per input file" << endl;
    cout << "                  \t (cannot be calculated if the input is a list of splits)" << endl;
    cout << endl;
    cout << "    (either --output and/or --newick, or --counts must be provided)" << endl;
    cout << " ___________________________________________________________________________" << endl;
    cout << endl;
    cout << "  K-mer options:" << endl;
    cout << endl;
    cout << "    -k, --kmer    \t Length of k-mers (default: 31)" << endl;
    cout << endl;
    cout << "    -l, --gapped  \t Pattern of gapped k-mers (default: no gaps)" << endl;
    cout << "                  \t e.g. 1110111, where 0 means a gap/skip base" << endl;
    cout << endl;
    cout << "    -w, --window  \t Number of k-mers per minimizer window (default: 1)" << endl;
    cout << "                  \t Consider only the smallest k-mer among w neighbors" << endl;
    cout << "                  \t Reduces memory usage, but also lowers the accuracy" << endl;
    cout << endl;
    cout << "    -x, --iupac   \t Extended IUPAC alphabet, resolve ambiguous bases" << endl;
    cout << "                  \t Specify a number to limit the k-mers per position" << endl;
    cout << "                  \t between 1 (no ambiguity) and 4^k (allows NNN...N)" << endl;
    cout << endl;
    cout << "    -q, --qualify \t Discard k-mers with lower coverage than a threshold" << endl;
    cout << "                  \t Only keep k-mers with multiple occurrences per file" << endl;
    cout << endl;
    cout << "    -r, --reverse \t Keep one repr. for reverse complement k-mers" << endl;
    cout << "                  \t Reverse and consider the lex. smaller of both" << endl;
    cout << " ___________________________________________________________________________" << endl;
    cout << endl;
    cout << "  Filter options:" << endl;
    cout << endl;
    cout << "    -t, --top     \t Number of splits in the output list (default: all)" << endl;
    cout << endl;
    cout << "    -m, --mean    \t Mean weight function to handle asymmetric splits" << endl;
    cout << "                  \t options: arith: arithmetic mean" << endl;
    cout << "                  \t          geom:  geometric mean (default)" << endl;
    cout << "                  \t          geom2: geometric mean with pseudo-counts" << endl;
    cout << endl;
    cout << "    -f, --filter  \t Output a greedy maximum weight subset of splits" << endl;
    cout << "                  \t options: strict: compatible to a tree" << endl;
    cout << "                  \t          weakly: weakly compatible network" << endl;
    cout << "                  \t          n-tree: compatible to a union of n trees" << endl;
    cout << "                  \t                  (where n is an arbitrary number)" << endl;
    cout << " ___________________________________________________________________________" << endl;
    cout << endl;
    cout << "  Other settings:" << endl;
    cout << endl;
    cout << "    -p, --threads \t Number of parallel threads (default: 1)" << endl;
    cout << "                  \t Reduces runtime, but increases memory usage" << endl;
    cout << endl;
    cout << "    -v, --verbose \t Print information messages during execution" << endl;
    cout << endl;
    cout << "    -h, --help    \t Display an extended help page and quit" << endl;
    cout << " ___________________________________________________________________________" << endl;
    cout << endl;
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
