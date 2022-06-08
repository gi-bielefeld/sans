#include "graph.h"

/**
 * This is the size of the top list.
 */
uint64_t graph::t;

bool graph::isAmino;

/**
 * This is a hash table mapping k-mers to colors [O(1)].
 */
hash_map<kmer_t, color_t> graph::kmer_table;

hash_map<kmerAmino_t, color_t> graph::kmer_tableAmino;

/**
 * This is a hash table mapping colors to weights [O(1)].
 */
hash_map<color_t, array<uint32_t,2>> graph::color_table;

/**
 * This is an ordered tree collecting the splits [O(log n)].
 */
multimap<double, color_t, greater<>> graph::split_list;

/**
* These are the allowed chars.
*/
vector<char> graph::allowedChars;


/**
 * This is a comparison function extending std::bitset.
 */ 
#if (maxK > 12 || maxN > 64)
namespace std {
    template <size_t N>
    bool operator<(const bitset<N>& x, const bitset<N>& y) {
        for (uint64_t i = N-1; i != -1; --i) {
            if (x[i] ^ y[i]) return y[i];
        }
        return false;
    }
}
#endif


/**
 * Initializes a new node struct.
 *
 * @param taxa color_t coding all taxa beneath this node
 * @param subsets list of subsets
 */
struct node* newSet(color_t taxa, double weight, vector<node*> subsets) {
    // declare and allocate new node
    auto* node = new struct node();
    node->taxa = std::move(taxa);
    node->weight = std::move(weight);
    node->subsets = std::move(subsets);
    return(node);
}

/**
 * This function initializes the top list size and the allowed chars.
 *
 * @param t top list size
 */
void graph::init(uint64_t& top_size, bool amino) {
    t = top_size;
    isAmino = amino;

    if(!isAmino){
        graph::allowedChars.push_back('A');
        graph::allowedChars.push_back('C');
        graph::allowedChars.push_back('G');
        graph::allowedChars.push_back('T');
    }else{
        graph::allowedChars.push_back('A');
        //graph::allowedChars.push_back('B');
        graph::allowedChars.push_back('C');
        graph::allowedChars.push_back('D');
        graph::allowedChars.push_back('E');
        graph::allowedChars.push_back('F');
        graph::allowedChars.push_back('G');
        graph::allowedChars.push_back('H');
        graph::allowedChars.push_back('I');
        //graph::allowedChars.push_back('J');
        graph::allowedChars.push_back('K');
        graph::allowedChars.push_back('L');
        graph::allowedChars.push_back('M');
        graph::allowedChars.push_back('N');
        graph::allowedChars.push_back('O');
        graph::allowedChars.push_back('P');
        graph::allowedChars.push_back('Q');
        graph::allowedChars.push_back('R');
        graph::allowedChars.push_back('S');
        graph::allowedChars.push_back('T');
        graph::allowedChars.push_back('U');
        graph::allowedChars.push_back('V');
        graph::allowedChars.push_back('W');
        //graph::allowedChars.push_back('X');
        graph::allowedChars.push_back('Y');
        //graph::allowedChars.push_back('Z');
        graph::allowedChars.push_back('*');
    }


}

/**
 * This function extracts k-mers from a sequence and adds them to the hash table.
 *
 * @param str dna sequence
 * @param color color flag
 * @param reverse merge complements
 */
void graph::add_kmers(string& str, uint64_t& color, bool& reverse) {
    if (str.length() < kmer::k) return;    // not enough characters

    uint64_t pos;    // current position in the string, from 0 to length
    kmer_t kmer;    // create a new empty bit sequence for the k-mer
    kmer_t rcmer;    // create a bit sequence for the reverse complement

    kmerAmino_t kmerAmino;    // create a new empty bit sequence for the k-mer

    uint64_t begin = 0;
next_kmer:
    pos = begin;

    for (; pos < str.length(); ++pos) {    // collect the bases from the string
        if (!isAllowedChar(pos, str)) {
            begin = pos+1;    // str = str.substr(pos+1, string::npos);
            goto next_kmer;    // unknown base, start a new k-mer from the beginning
        }
        if (!isAmino) {
            kmer::shift_right(kmer, str[pos]);    // shift each base into the bit sequence

            if (pos+1 - begin >= kmer::k) {
                rcmer = kmer;
                if (reverse) kmer::reverse_complement(rcmer, true);    // invert the k-mer, if necessary
                color::set(kmer_table[rcmer], color);    // update the k-mer with the current color
            }
        } else {
            kmerAmino::shift_right(kmerAmino, str[pos]);    // shift each base into the bit sequence

            if (pos+1 - begin >= kmerAmino::k) {
                color::set(kmer_tableAmino[kmerAmino], color);    // update the k-mer with the current color
            }
        }


    }
}

/**
 * This function extracts k-mer minimizers from a sequence and adds them to the hash table.
 *
 * @param str dna sequence
 * @param color color flag
 * @param reverse merge complements
 * @param m number of k-mers to minimize
 */
void graph::add_minimizers(string& str, uint64_t& color, bool& reverse, uint64_t& m) {
    if (str.length() < (!isAmino ? kmer::k : kmerAmino::k)) return;    // not enough characters

    vector<kmer_t> sequence_order;    // k-mers ordered by their position in sequence
    multiset<kmer_t> value_order;    // k-mers ordered by their lexicographical value

    vector<kmerAmino_t> sequence_order_Amino;    // k-mers ordered by their position in sequence
    multiset<kmerAmino_t> value_order_Amino;    // k-mers ordered by their lexicographical value

    uint64_t pos;    // current position in the string, from 0 to length
    kmer_t kmer;    // create a new empty bit sequence for the k-mer
    kmer_t rcmer;    // create a bit sequence for the reverse complement

    kmerAmino_t kmerAmino;    // create a new empty bit sequence for the k-mer

    uint64_t begin = 0;
next_kmer:
    pos = begin;
    sequence_order.clear();
    sequence_order_Amino.clear();
    value_order.clear();
    sequence_order_Amino.clear();

    for (; pos < str.length(); ++pos) {    // collect the bases from the string
        if (!isAllowedChar(pos, str)) {
            begin = pos+1;    // str = str.substr(pos+1, string::npos);
            goto next_kmer;    // unknown base, start a new k-mer from the beginning
        }
        if (!isAmino) {
            kmer::shift_right(kmer, str[pos]);    // shift each base into the bit sequence

            if (pos+1 - begin >= kmer::k) {
                rcmer = kmer;
                if (reverse) kmer::reverse_complement(rcmer, true);    // invert the k-mer, if necessary

                if (sequence_order.size() == m) {
                    value_order.erase(*sequence_order.begin());    // remove k-mer outside the window
                    sequence_order.erase(sequence_order.begin());
                }
                value_order.emplace(rcmer);    // insert k-mer ordered by its lexicographical value
                sequence_order.emplace_back(rcmer);

                if (sequence_order.size() == m) {
                    color::set(kmer_table[*value_order.begin()], color);    // update the k-mer with the current color
                }
            }
        } else {
            kmerAmino::shift_right(kmerAmino, str[pos]);    // shift each base into the bit sequence

            if (pos+1 - begin >= kmerAmino::k) {
                if (sequence_order.size() == m) {
                    value_order_Amino.erase(*sequence_order_Amino.begin());    // remove k-mer outside the window
                    sequence_order_Amino.erase(sequence_order_Amino.begin());
                }
                value_order_Amino.emplace(kmerAmino);    // insert k-mer ordered by its lexicographical value
                sequence_order_Amino.emplace_back(kmerAmino);

                if (sequence_order_Amino.size() == m) {

                    color::set(kmer_tableAmino[*value_order_Amino.begin()], color);    // update the k-mer with the current color
                }
            }
        }

    }
}
/**
 * This function checks if the character at the given position is allowed.
 * @param pos position in str
 * @param str the current part of the sequence
 * @return true if allowed, false otherwise
 */
bool graph::isAllowedChar(uint64_t pos, string &str) {
    bool allowed = false;
    char &currentChar = str[pos];

    for (int i = 0; i < graph::allowedChars.size() && !allowed; i++){
        allowed =  graph::allowedChars.at(i) == currentChar;
    }

    return allowed;
}

/**
 * This function extracts k-mers from a sequence and adds them to the hash table.
 *
 * @param str dna sequence
 * @param color color flag
 * @param reverse merge complements
 * @param max_iupac allowed number of ambiguous k-mers per position
 */
void graph::add_kmers(string& str, uint64_t& color, bool& reverse, uint64_t& max_iupac) {
    if (str.length() < (!isAmino ? kmer::k : kmerAmino::k)) return;    // not enough characters

    if (!isAmino) {
        hash_set<kmer_t> ping;    // create a new empty set for the k-mers
        hash_set<kmer_t> pong;    // create another new set for the k-mers
        bool ball; bool wait;    // indicates which of the two sets should be used

        vector<uint8_t> factors;    // stores the multiplicity of iupac bases
        long double product;    // stores the overall multiplicity of the k-mers

        uint64_t pos;    // current position in the string, from 0 to length
        kmer_t kmer;    // create an empty bit sequence for the initial k-mer
        kmer_t rcmer;    // create a bit sequence for the reverse complement

        uint64_t begin = 0;
        next_kmer:
        pos = begin;

        ping.clear(); pong.clear(); factors.clear();
        ball = true; wait = false; product = 1;
        (ball ? ping : pong).emplace(kmer);

        for (; pos < str.length(); ++pos) {    // collect the bases from the string
            if (str[pos] == '.' || str[pos] == '-') {
                begin = pos+1;    // str = str.substr(pos+1, string::npos);
                goto next_kmer;    // gap character, start a new k-mer from the beginning
            }
            iupac_calc(product, factors, str[pos]);

            if (product <= max_iupac) {    // check if there are too many ambiguous k-mers
                if (wait) {
                    begin = pos-kmer::k+1;    // str = str.substr(pos-kmer::k+1, string::npos);
                    goto next_kmer;    // start a new k-mer from the beginning
                }
                iupac_shift(ball ? ping : pong, !ball ? ping : pong, str[pos]);
                ball = !ball;    // shift each base in, resolve iupac character
            } else { wait = true; continue; }

            if (pos+1 - begin >= kmer::k) {
                for (auto& kmer : (ball ? ping : pong)) {    // iterate over the current set of ambiguous k-mers
                    rcmer = kmer;
                    if (reverse) kmer::reverse_complement(rcmer, true);    // invert the k-mer, if necessary
                    color::set(kmer_table[rcmer], color);    // update the k-mer with the current color
                }
            }
        }
    } else {
        hash_set<kmerAmino_t> ping;    // create a new empty set for the k-mers
        hash_set<kmerAmino_t> pong;    // create another new set for the k-mers
        bool ball; bool wait;    // indicates which of the two sets should be used

        vector<uint8_t> factors;    // stores the multiplicity of iupac bases
        long double product;    // stores the overall multiplicity of the k-mers

        uint64_t pos;    // current position in the string, from 0 to length
        kmerAmino_t kmer;    // create an empty bit sequence for the initial k-mer

        uint64_t begin = 0;
        next_kmerAmino:
        pos = begin;

        ping.clear(); pong.clear(); factors.clear();
        ball = true; wait = false; product = 1;
        (ball ? ping : pong).emplace(kmer);

        for (; pos < str.length(); ++pos) {    // collect the bases from the string
            if (str[pos] == '.' || str[pos] == '-') {
                begin = pos+1;    // str = str.substr(pos+1, string::npos);
                goto next_kmerAmino;    // gap character, start a new k-mer from the beginning
            }
            iupac_calc(product, factors, str[pos]);

            if (product <= max_iupac) {    // check if there are too many ambiguous k-mers
                if (wait) {
                    begin = pos-kmerAmino::k+1;    // str = str.substr(pos-kmer::k+1, string::npos);
                    goto next_kmerAmino;    // start a new k-mer from the beginning
                }
                iupac_shift_amino(ball ? ping : pong, !ball ? ping : pong, str[pos]);
                ball = !ball;    // shift each base in, resolve iupac character
            } else { wait = true; continue; }

            if (pos+1 - begin >= kmerAmino::k) {
                for (auto& kmer : (ball ? ping : pong)) {    // iterate over the current set of ambiguous k-mers
                    color::set(kmer_tableAmino[kmer], color);    // update the k-mer with the current color
                }
            }
        }
    }


}

/**
 * This function extracts k-mer minimizers from a sequence and adds them to the hash table.
 *
 * @param str dna sequence
 * @param color color flag
 * @param reverse merge complements
 * @param m number of k-mers to minimize
 * @param max_iupac allowed number of ambiguous k-mers per position
 */
void graph::add_minimizers(string& str, uint64_t& color, bool& reverse, uint64_t& m, uint64_t& max_iupac) {
    if (str.length() < (!isAmino ? kmer::k : kmerAmino::k)) return;    // not enough characters

   if (!isAmino) {
       vector<kmer_t> sequence_order;    // k-mers ordered by their position in sequence
       multiset<kmer_t> value_order;    // k-mers ordered by their lexicographical value
       multiset<kmer_t> inner_value_order;

       hash_set<kmer_t> ping;    // create a new empty set for the k-mers
       hash_set<kmer_t> pong;    // create another new set for the k-mers
       bool ball; bool wait;    // indicates which of the two sets should be used

       vector<uint8_t> factors;    // stores the multiplicity of iupac bases
       long double product;    // stores the overall multiplicity of the k-mers

       uint64_t pos;    // current position in the string, from 0 to length
       kmer_t kmer;    // create an empty bit sequence for the initial k-mer
       kmer_t rcmer;    // create a bit sequence for the reverse complement

       uint64_t begin = 0;
       next_kmer:
       pos = begin;
       sequence_order.clear();
       value_order.clear();

       ping.clear(); pong.clear(); factors.clear();
       ball = true; wait = false; product = 1;
       (ball ? ping : pong).emplace(kmer);

       for (; pos < str.length(); ++pos) {    // collect the bases from the string
           if (str[pos] == '.' || str[pos] == '-') {
               begin = pos+1;    // str = str.substr(pos+1, string::npos);
               goto next_kmer;    // gap character, start a new k-mer from the beginning
           }
           iupac_calc(product, factors, str[pos]);

           if (product <= max_iupac) {    // check if there are too many ambiguous k-mers
               if (wait) {
                   begin = pos-kmer::k+1;    // str = str.substr(pos-kmer::k+1, string::npos);
                   goto next_kmer;    // start a new k-mer from the beginning
               }
               iupac_shift(ball ? ping : pong, !ball ? ping : pong, str[pos]);
               ball = !ball;    // shift each base in, resolve iupac character
           } else { wait = true; continue; }

           if (pos+1 - begin >= kmer::k) {
               for (auto& kmer : (ball ? ping : pong)) {    // iterate over the current set of ambiguous k-mers
                   rcmer = kmer;
                   if (reverse) kmer::reverse_complement(rcmer, true);    // invert the k-mer, if necessary
                   inner_value_order.emplace(rcmer);
               }

               if (sequence_order.size() == m) {
                   value_order.erase(*sequence_order.begin());    // remove k-mer outside the window
                   sequence_order.erase(sequence_order.begin());
               }
               value_order.emplace(*inner_value_order.begin());    // insert k-mer ordered by its lexicographical value
               sequence_order.emplace_back(*inner_value_order.begin());
               inner_value_order.clear();

               if (sequence_order.size() == m) {
                   color::set(kmer_table[*value_order.begin()], color);    // update the k-mer with the current color
               }
           }
       }
   } else {
       vector<kmerAmino_t> sequence_order;    // k-mers ordered by their position in sequence
       multiset<kmerAmino_t> value_order;    // k-mers ordered by their lexicographical value
       multiset<kmerAmino_t> inner_value_order;

       hash_set<kmerAmino_t> ping;    // create a new empty set for the k-mers
       hash_set<kmerAmino_t> pong;    // create another new set for the k-mers
       bool ball; bool wait;    // indicates which of the two sets should be used

       vector<uint8_t> factors;    // stores the multiplicity of iupac bases
       long double product;    // stores the overall multiplicity of the k-mers

       uint64_t pos;    // current position in the string, from 0 to length
       kmerAmino_t kmer;    // create an empty bit sequence for the initial k-mer

       uint64_t begin = 0;
       next_kmerAmino:
       pos = begin;
       sequence_order.clear();
       value_order.clear();

       ping.clear(); pong.clear(); factors.clear();
       ball = true; wait = false; product = 1;
       (ball ? ping : pong).emplace(kmer);

       for (; pos < str.length(); ++pos) {    // collect the bases from the string
           if (str[pos] == '.' || str[pos] == '-') {
               begin = pos+1;    // str = str.substr(pos+1, string::npos);
               goto next_kmerAmino;    // gap character, start a new k-mer from the beginning
           }
           iupac_calc(product, factors, str[pos]);

           if (product <= max_iupac) {    // check if there are too many ambiguous k-mers
               if (wait) {
                   begin = pos-kmerAmino::k+1;    // str = str.substr(pos-kmer::k+1, string::npos);
                   goto next_kmerAmino;    // start a new k-mer from the beginning
               }
               iupac_shift_amino(ball ? ping : pong, !ball ? ping : pong, str[pos]);
               ball = !ball;    // shift each base in, resolve iupac character
           } else { wait = true; continue; }

           if (pos+1 - begin >= kmerAmino::k) {
               for (auto& kmer : (ball ? ping : pong)) {    // iterate over the current set of ambiguous k-mers
                   inner_value_order.emplace(kmer);
               }

               if (sequence_order.size() == m) {
                   value_order.erase(*sequence_order.begin());    // remove k-mer outside the window
                   sequence_order.erase(sequence_order.begin());
               }
               value_order.emplace(*inner_value_order.begin());    // insert k-mer ordered by its lexicographical value
               sequence_order.emplace_back(*inner_value_order.begin());
               inner_value_order.clear();

               if (sequence_order.size() == m) {
                   color::set(kmer_tableAmino[*value_order.begin()], color);    // update the k-mer with the current color
               }
           }
       }
   }
}

/**
 * This function calculates the multiplicity of iupac k-mers.
 *
 * @param product overall multiplicity
 * @param factors per base multiplicity
 * @param input iupac character
 */
void graph::iupac_calc(long double& product, vector<uint8_t>& factors, char& input) {

    if(!isAmino){
        switch (input) {
            case 'A': case 'C': case 'G': case 'T':
                product *= 1;
                factors.emplace_back(1);
        }
        switch (input) {
            case 'R': case 'Y': case 'S': case 'W': case 'K': case 'M':
                product *= 2;
                factors.emplace_back(2);
        }
        switch (input) {
            case 'B': case 'D': case 'H': case 'V':
                product *= 3;
                factors.emplace_back(3);
        }
        switch (input) {
            case 'N':
                product *= 4;
                factors.emplace_back(4);
        }
        if (factors.size() > kmer::k) {
            product /= *factors.begin();
            factors.erase(factors.begin());
        }
    }else{
        switch (input) {
            case 'A': case 'C': case 'D': case 'E': case 'F': case 'G':
            case 'H': case 'I': case 'K': case 'L': case 'M': case 'N':
            case 'O': case 'P': case 'Q': case 'R': case 'S': case 'T':
            case 'U': case 'V': case 'W': case 'Y': case '*':
                product *= 1;
                factors.emplace_back(1);
        }
        switch (input) {
            case 'B': case 'Z': case 'J':
                product *= 2;
                factors.emplace_back(2);
        }
        switch (input) {
            case 'X':
                product *= 22;
                factors.emplace_back(22);
        }
        if (factors.size() > kmerAmino::k) {
            product /= *factors.begin();
            factors.erase(factors.begin());
        }
    }



}

/**
 * This function shifts a base into a set of ambiguous iupac k-mers.
 *
 * @param prev set of k-mers
 * @param next set of k-mers
 * @param input iupac character
 */
void graph::iupac_shift(hash_set<kmer_t>& prev, hash_set<kmer_t>& next, char& input) {
    kmer_t temp; char base;
    while (!prev.empty()) {    // extend each previous k-mer
        switch (input) {
            case 'A': case 'R': case 'W': case 'M':
            case 'D': case 'H': case 'V': case 'N':
                temp = *prev.begin(); base = 'A';
                kmer::shift_right(temp, base);
                next.emplace(temp);
        }
        switch (input) {
            case 'C': case 'Y': case 'S': case 'M':
            case 'B': case 'H': case 'V': case 'N':
                temp = *prev.begin(); base = 'C';
                kmer::shift_right(temp, base);
                next.emplace(temp);
        }
        switch (input) {
            case 'G': case 'R': case 'S': case 'K':
            case 'B': case 'D': case 'V': case 'N':
                temp = *prev.begin(); base = 'G';
                kmer::shift_right(temp, base);
                next.emplace(temp);
        }
        switch (input) {
            case 'T': case 'Y': case 'W': case 'K':
            case 'B': case 'D': case 'H': case 'N':
                temp = *prev.begin(); base = 'T';
                kmer::shift_right(temp, base);
                next.emplace(temp);
        }
        prev.erase(prev.begin());
    }
}

/**
 * This function shifts an amino acid into a set of ambiguous iupac k-mers.
 *
 * @param prev set of k-mers
 * @param next set of k-mers
 * @param input iupac character
 */
void graph::iupac_shift_amino(hash_set<kmerAmino_t>& prev, hash_set<kmerAmino_t>& next, char& input) {
    string acidsUnique = "ACDEFGHIKLMNOPQRSTUVWY*";
    string acidsB = "DN";
    string acidsZ = "EQ";
    string acidsJ = "LI";

    kmerAmino_t temp; char acid;
    while (!prev.empty()) {    // extend each previous k-mer

        //first: all unique acids OR X for all
        for (char i : acidsUnique) {
            acid = i;
            if (acid == input || input == 'X') {
                temp = *prev.begin();
                kmerAmino::shift_right(temp, acid);
                next.emplace(temp);
            }
        }

        //second: all ambigious acids

        if (input == 'B') {
            for (char i : acidsB) {
                temp = *prev.begin();
                acid = i;
                kmerAmino::shift_right(temp, acid);
                next.emplace(temp);
            }
        }

        if (input == 'Z') {
            for (char i : acidsZ) {
                temp = *prev.begin();
                acid = i;
                kmerAmino::shift_right(temp, acid);
                next.emplace(temp);
            }
        }

        if (input == 'J'){
            for (char i : acidsJ) {
                temp = *prev.begin();
                acid = i;
                kmerAmino::shift_right(temp, acid);
                next.emplace(temp);
            }
        }
        prev.erase(prev.begin());
    }
}

/**
 * This function iterates over the hash table and calculates the split weights.
 *
 * @param mean weight function
 * @param verbose print progress
 */
void graph::add_weights(double mean(uint32_t&, uint32_t&), bool& verbose) {
    double min_value = numeric_limits<double>::min();    // current min. weight in the top list (>0)
    uint64_t cur = 0, prog = 0, next;
    uint64_t max = !isAmino ? kmer_table.size() : kmer_tableAmino.size();

    if (!isAmino) {
        for (auto it = kmer_table.begin(); it != kmer_table.end(); ++it) {    // iterate over k-mer hash table
            if (verbose) {
                next = 100*cur/max;
                if (prog < next)  cout << "\33[2K\r" << "Processing splits... " << next << "%" << flush;
                prog = next; cur++;
            }

            color_t& color = it.value();    // get the color set for each k-mer
            bool pos = color::complement(color, true);    // invert the color set, if necessary
            if (color == 0) continue;    // ignore empty splits
            array<uint32_t,2>& weight = color_table[color];    // get the weight and inverse weight for the color set

            double old_value = mean(weight[0], weight[1]);    // calculate the old mean value
            if (old_value >= min_value) {    // if it is greater than the min. value, find it in the top list
                auto range = split_list.equal_range(old_value);    // get all color sets with the given weight
                for (auto it = range.first; it != range.second; ++it) {
                    if (it->second == color) {    // iterate over the color sets to find the correct one
                        split_list.erase(it);    // erase the entry with the old weight
                        break;
                    }
                }
            }
            weight[pos]++;    // update the weight or the inverse weight of the current color set

            double new_value = mean(weight[0], weight[1]);    // calculate the new mean value
            if (new_value >= min_value) {    // if it is greater than the min. value, add it to the top list
                split_list.emplace(new_value, color);    // insert it at the correct position ordered by weight
                if (split_list.size() > t) {
                    split_list.erase(--split_list.end());    // if the top list exceeds its limit, erase the last entry
                    min_value = split_list.rbegin()->first;    // update the min. value for the next iteration
                }
            }
        }
    } else {
        for (auto it = kmer_tableAmino.begin(); it != kmer_tableAmino.end(); ++it) {    // iterate over k-mer hash table
            if (verbose) {
                next = 100*cur/max;
                if (prog < next)  cout << "\33[2K\r" << "Processing splits... " << next << "%" << flush;
                prog = next; cur++;
            }
            color_t& color = it.value();    // get the color set for each k-mer
            bool pos = color::complement(color, true);    // invert the color set, if necessary
            if (color == 0) continue;    // ignore empty splits
            array<uint32_t,2>& weight = color_table[color];    // get the weight and inverse weight for the color set

            double old_value = mean(weight[0], weight[1]);    // calculate the old mean value
            if (old_value >= min_value) {    // if it is greater than the min. value, find it in the top list
                auto range = split_list.equal_range(old_value);    // get all color sets with the given weight
                for (auto it = range.first; it != range.second; ++it) {
                    if (it->second == color) {    // iterate over the color sets to find the correct one
                        split_list.erase(it);    // erase the entry with the old weight
                        break;
                    }
                }
            }
            weight[pos]++;    // update the weight or the inverse weight of the current color set

            double new_value = mean(weight[0], weight[1]);    // calculate the new mean value
            if (new_value >= min_value) {    // if it is greater than the min. value, add it to the top list
                split_list.emplace(new_value, color);    // insert it at the correct position ordered by weight
                if (split_list.size() > t) {
                    split_list.erase(--split_list.end());    // if the top list exceeds its limit, erase the last entry
                    min_value = split_list.rbegin()->first;    // update the min. value for the next iteration
                }
            }
        }
    }




}

/**
 * This function adds a single split (weight and colors) to the output list.
 *
 * @param weight split weight
 * @param color split colors
 */
void graph::add_split(double& weight, color_t& color) {
    split_list.emplace(weight, color);    // insert it at the correct position ordered by weight
    if (split_list.size() > t) {
        split_list.erase(--split_list.end());    // if the top list exceeds its limit, erase the last entry
    }
}

/**
 * This function filters a greedy maximum weight tree compatible subset.
 *
 * @param verbose print progress
 */
void graph::filter_strict(bool& verbose) {
    filter_strict(nullptr, verbose);
}

/**
 * This function filters a greedy maximum weight tree compatible subset and returns a newick string.
 *
 * @param map function that maps an integer to the original id, or null if no newick output wanted
 * @param verbose print progress
 */
string graph::filter_strict(std::function<string(const uint64_t&)> map, bool& verbose) {
    auto tree = vector<color_t>();    // create a set for compatible splits
    auto it = split_list.begin();
    uint64_t cur = 0, prog = 0, next;
    uint64_t max = split_list.size();
loop:
    while (it != split_list.end()) {
        if (verbose) {
            next = 100*cur/max;
             if (prog < next)  cout << "\33[2K\r" << "Filtering splits... " << next << "%" << flush;
            prog = next; cur++;
        }
        if (test_strict(it->second, tree)) {
            tree.emplace_back(it->second);
            ++it; goto loop;    // if compatible, add the new split to the set
        }
        it = split_list.erase(it);    // otherwise, remove split
    }
    if (map) {
        node* root = build_tree(tree);
        return print_tree(root, map) + ";\n";
    } else {
        return "";
    }
}

/**
 * This function filters a greedy maximum weight weakly compatible subset.
 *
 * @param verbose print progress
 */
void graph::filter_weakly(bool& verbose) {
    auto network = vector<color_t>();    // create a set for compatible splits
    auto it = split_list.begin();
    uint64_t cur = 0, prog = 0, next;
    uint64_t max = split_list.size();
loop:
    while (it != split_list.end()) {
        if (verbose) {
            next = 100*(cur*cur)/(max*max);
             if (prog < next)  cout << "\33[2K\r" << "Filtering splits... " << next << "%" << flush;
            prog = next; cur++;
        }
        if (test_weakly(it->second, network)) {
            network.emplace_back(it->second);
            ++it; goto loop;    // if compatible, add the new split to the set
        }
        it = split_list.erase(it);    // otherwise, remove split
    }
}

/**
 * This function filters a greedy maximum weight n-tree compatible subset.
 *
 * @param n number of trees
 * @param verbose print progress
 */
void graph::filter_n_tree(uint64_t n, bool& verbose) {
    filter_n_tree(n, nullptr, verbose);
}

/**
 * This function filters a greedy maximum weight n-tree compatible subset and returns a string with all trees in newick format.
 *
 * @param n number of trees
 * @param map function that maps an integer to the original id, or null
 * @param verbose print progress
 */
string graph::filter_n_tree(uint64_t n, std::function<string(const uint64_t&)> map, bool& verbose) {
    auto forest = vector<vector<color_t>>(n);    // create a set for compatible splits
    auto it = split_list.begin();
    uint64_t cur = 0, prog = 0, next;
    uint64_t max = split_list.size();
loop:
    while (it != split_list.end()) {
        if (verbose) {
            next = 100*cur/max;
             if (prog < next)  cout << "\33[2K\r" << "Filtering splits... " << next << "%" << flush;
            prog = next; cur++;
        }
       for (auto& tree : forest)
        if (test_strict(it->second, tree)) {
            tree.emplace_back(it->second);
            ++it; goto loop;    // if compatible, add the new split to the set
        }
        it = split_list.erase(it);    // otherwise, remove split
    }
    // output
    string s;
    if (map) {
        for (auto& tree : forest) {
            node* root = build_tree(tree);
            s += print_tree(root, map) + ";\n";
        }
    }
    return s;
}

/**
 * This function tests if a split is compatible with an existing set of splits.
 *
 * @param color new split
 * @param color_set set of splits
 * @return true, if compatible
 */
bool graph::test_strict(color_t& color, vector<color_t>& color_set) {
    for (auto& elem : color_set) {
        if (!color::is_compatible(elem, color)) {
            return false;    // compare to each split in the set
        }
    }
    return true;
}

/**
 * This function tests if a split is weakly compatible with an existing set of splits.
 *
 * @param color new split
 * @param color_set set of splits
 * @return true, if weakly compatible
 */
bool graph::test_weakly(color_t& color, vector<color_t>& color_set) {
    for (auto& elem1 : color_set) {
        for (auto& elem2 : color_set) {
            if (elem1 != elem2) {
                if (!color::is_weakly_compatible(elem1, elem2, color)) {
                    return false;    // compare to each split in the set
                }
            }
        }
    }
    return true;
}

/**
 * This function recursively refines a given set/tree structure by a given split.
 *
 * @param current_set node of currently considered (sub-)set/tree structure
 * @param split color set to refine by
 * @return whether or not the given split is compatible with the set/tree structure
 */
bool graph::refine_tree(node* current_set, color_t& split, color_t& allTaxa) {
    // possible cases:
    // splitsize <2: nothing has to be done
    // split equals one subset -> warning: split twice
    // split is fully contained in one subset -> recurse
    // inverse split ... (i.e. split covers one subset partially) -> recurse with inverse
    // split covers several subsets completely -> introduce new split
    if (color::size(split, false) < 2 || color::size(allTaxa, false) - color::size(split, false) < 2) { return true; }

    vector<node*> *subsets = &current_set->subsets;
    vector<node*> fullycoveredsubsets = {};
    node* partiallycoveredsubset = nullptr;

    for (node* subset : *subsets) {
        color_t subtaxa = subset->taxa;
        if (split == subtaxa) {
            return true;
        }
        // split.issubset(subtaxa)?
        if ((split & subtaxa) == split) { return refine_tree(subset, split, allTaxa); }
        // subtaxa.issubset(split):
        if ((subtaxa & split) == subtaxa) { fullycoveredsubsets.push_back(subset); }
        // elif not subtaxa.isdisjoint(split): # does intersect
        else if ((subtaxa & split) != 0b0u) {
            // if partiallycoveredsubset:
            if (partiallycoveredsubset != nullptr) { return false; } //there cannot be more than one
            else { partiallycoveredsubset = subset; }
        }
    }

    if (partiallycoveredsubset != nullptr) {
        if (fullycoveredsubsets.size() == subsets->size()-1) {
            // recurse into this subset with inverse split
            color_t inversesplit = color::complement(split, false);
            // if inversesplit.issubset(partiallycoveredsubset[1]):
            if ((inversesplit & partiallycoveredsubset->taxa) == inversesplit) {
                return refine_tree(partiallycoveredsubset, inversesplit, allTaxa);
            } else { return false; }
        } else { return false; }
    } else if (fullycoveredsubsets.size() > 1) {
        // introduce new split
        color_t newsubtaxa = 0b0u;
        for(node* subset : fullycoveredsubsets) { newsubtaxa |= subset->taxa; }
        // get weight of split
        double weight = 0;
        auto it = split_list.begin();
        while (it != split_list.end()) {
            if (it->second == split){
                weight = it->first;
                break;
            }
            it++;
        }
        node* newset = newSet(newsubtaxa, weight, fullycoveredsubsets);
        // remove old sets
        for(node* subset : fullycoveredsubsets) {
            // subsets.remove(subset)
            subsets->erase(std::remove(subsets->begin(), subsets->end(), subset), subsets->end());
        }
        // add new set
        subsets->push_back(newset);
        return true;
    } else {
        std::cerr << "ERROR: this cannot be: just one fully covered subset and nothing else!?" << endl;
        exit(EXIT_FAILURE);
    }
}

/**
 * This function returns a tree structure (struct node) generated from the given list of color sets.
 *
 * @param color_set list of color sets
 * @return tree structure (struct node)
 */
node* graph::build_tree(vector<color_t>& color_set) {
    //initialize set of trivial splits
    vector<node*> subsets = {};
    color_t allTaxa = 0b0u;

    for (uint64_t i = 0; i < color::n; i++) {
        color_t leaf = 0b0u;
        color::set(leaf, i);
        color::set(allTaxa, i);
        vector<node*> emptyset = {};
        // get weight
        double weight = 0;
        auto it = split_list.begin();
        while (it != split_list.end()) {
            if (it->second == leaf){
                weight = it->first;
                break;
            }
            it++;
        }
        node* newset = newSet(leaf, weight, emptyset);
        subsets.push_back(newset);
    }
    node* sets = newSet(allTaxa, 0, subsets);

    for (color_t split : color_set) {
        // split if possible
        if (!refine_tree(sets, split, allTaxa)) {
            std::cerr << "ERROR: splits are incompatible" << endl;
            exit(EXIT_FAILURE);
        }
    }
    return sets;
}

/**
 * This function returns a newick string generated from the given tree structure (set).
 *
 * @param root root of the tree/set structure
 * @return newick string
 */
string graph::print_tree(node* root, std::function<string(const uint64_t&)> map) {
    vector<node*> subsets = root->subsets;
    color_t taxa = root->taxa;

    if (subsets.empty()){    // leaf set
        if (color::size(taxa, false) == 0) {
            std::cerr << "ERROR: child with no taxon!?" << endl;
            exit(EXIT_FAILURE);
        } else if (color::size(taxa, false) == 1) {
            return map(color::pos(taxa)) + ":" + to_string(root->weight);
        } else {
            std::cerr << "ERROR: child with more than one taxon!?" << endl;
            exit(EXIT_FAILURE);
        }
    }
    else {
        string s = "(";
        for (node* subset : subsets) {
            s += print_tree(subset, map);
            if (subset != subsets.back()) { s += ","; }
        }
        s += "):";
        s += to_string(root->weight);
        return s;
    }
}

/**
 * This function calculates the k-mer entropy for the input genomes
 * @param n the number of genomes in the input
 */
void graph::entropy(int n){
    double ent = 0;
    //cout << kmer_table.size() << endl;
    /*for(auto test: kmer_table){
        color_t col = test.second;
        int freq = color64::size(col, false);
        cout << "Frequency : "<< freq << endl;
        ent = ent + freq * log2(freq * 1.0);
    }*/
    for(auto it = kmer_table.begin(); it != kmer_table.end(); ++it){
        color_t col = it.value();
        int freq = color::size(col, false);
        double rel = (freq *1.0)/n;
        //cout << "p_i" << rel << endl;
        ent = ent - rel * log2(rel);
    }

    //ent = - ent;

    cout << "entropy = " << ent << endl;
}

/**
 * This function calculates the Shannon diversity for the input genomes
 * @param n the number of genomes in the input
 */
void graph::diversity(int n){

    double div = 0.0;

    int ok = 0;

    for(int i = 1; i <= n; i++){
        int ci = 0;
        if(isAmino){
            ok = kmer_tableAmino.size();
            for(auto it = kmer_tableAmino.begin(); it != kmer_tableAmino.end(); ++it){
                color_t col = it.value();
                int freq = color::size(col, false);
                if(freq == i){
                    ci = ci+1;
                }
            }
        }
        else{
            ok = kmer_table.size();
            for(auto it = kmer_table.begin(); it != kmer_table.end(); ++it){
                color_t col = it.value();
                int freq = color::size(col, false);
                if(freq == i){
                    ci = ci+1;
                }
            }
        }
        //cout << ci << endl;
        if(ci != 0){
            double rel = (ci*1.0)/(ok*1.0);
            div = div + rel * log2(rel);
        }
    }

    div = -div;
    cout << "shannon diversity: " << div << endl;

}
