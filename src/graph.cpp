#include "graph.h"

/*
 * This class manages the k-mer/color hash tables and split list.
 */
uint64_t                             graph::t;             // size of the top list
vector<hash_map<kmer_t, color_t>>    graph::kmer_table;    // hash table mapping k-mers to their colors/splits
hash_map<color_t, array<uint32_t,2>> graph::color_table;   // hash table mapping colors/splits to their weights

/**
 * This function initializes the top list size.
 *
 * @param top_size top list size
 */
void graph::init(const uint64_t& T, const uint64_t& top_size) {
    t = top_size;
    kmer_table.resize(T);
}

/**
 * This function extracts k-mers from a sequence and adds them to the hash table.
 *
 * @param str dna sequence
 * @param color color flag
 * @param reverse merge complements
 */
void graph::add_kmers(const uint64_t& T, const string& str, const size1N_t& color, const bool& reverse) {
    if (str.length() < kmer::k) return;    // not enough characters

    size_t pos;    // current position in the string, from 0 to length
    kmer_t kmer;    // create a new empty bit sequence for the k-mer
    kmer_t rcmer;    // create a bit sequence for the reverse complement

    size_t begin = 0;
next_kmer:
    pos = begin;

    for (; pos < str.length(); ++pos) {    // collect the bases from the string
        if (str[pos] != 'A' && str[pos] != 'C' && str[pos] != 'G' && str[pos] != 'T') {
            begin = pos+1;    // str = str.substr(pos+1, string::npos);
            goto next_kmer;    // unknown base, start a new k-mer from the beginning
        }
        kmer::shift(kmer, str[pos]);    // shift each base into the bit sequence

        if (pos+1 - begin >= kmer::k) {
            rcmer = kmer;
            if (reverse) kmer::reverse_represent(rcmer);    // invert the k-mer, if necessary
            color::set(kmer_table[T][rcmer], color);    // update the k-mer with the current color
        }
    }
}

/**
 * This function extracts k-mer minimizers from a sequence and adds them to the hash table.
 *
 * @param str dna sequence
 * @param color color flag
 * @param reverse merge complements
 * @param window number of k-mers to minimize
 */
void graph::add_minimizers(const uint64_t& T, const string& str, const size1N_t& color, const bool& reverse, const uint64_t& window) {
    if (str.length() < kmer::k) return;    // not enough characters

    vector<kmer_t> sequence_order;    // k-mers ordered by their position in sequence
    multiset<kmer_t> value_order;    // k-mers ordered by their lexicographical value

    size_t pos;    // current position in the string, from 0 to length
    kmer_t kmer;    // create a new empty bit sequence for the k-mer
    kmer_t rcmer;    // create a bit sequence for the reverse complement

    size_t begin = 0;
next_kmer:
    pos = begin;

    sequence_order.clear();
    value_order.clear();

    for (; pos < str.length(); ++pos) {    // collect the bases from the string
        if (str[pos] != 'A' && str[pos] != 'C' && str[pos] != 'G' && str[pos] != 'T') {
            begin = pos+1;    // str = str.substr(pos+1, string::npos);
            goto next_kmer;    // unknown base, start a new k-mer from the beginning
        }
        kmer::shift(kmer, str[pos]);    // shift each base into the bit sequence

        if (pos+1 - begin >= kmer::k) {
            rcmer = kmer;
            if (reverse) kmer::reverse_represent(rcmer);    // invert the k-mer, if necessary

            if (sequence_order.size() == window) {
                value_order.erase(*sequence_order.begin());    // remove k-mer outside the window
                sequence_order.erase(sequence_order.begin());
            }
            value_order.emplace(rcmer);    // insert k-mer ordered by its lexicographical value
            sequence_order.emplace_back(rcmer);

            if (sequence_order.size() == window) {
                color::set(kmer_table[T][*value_order.begin()], color);    // update the k-mer with the current color
            }
        }
    }
}

/**
 * This function extracts k-mers from a sequence and adds them to the hash table.
 *
 * @param str dna sequence
 * @param color color flag
 * @param reverse merge complements
 * @param max_iupac allowed number of ambiguous k-mers per position
 */
void graph::add_kmers(const uint64_t& T, const string& str, const size1N_t& color, const bool& reverse, const uint64_t& max_iupac) {
    if (str.length() < kmer::k) return;    // not enough characters

    hash_set<kmer_t> ping;    // create a new empty set for the k-mers
    hash_set<kmer_t> pong;    // create another new set for the k-mers
    bool ball; bool wait;    // indicates which of the two sets should be used

    vector<uint8_t> factors;    // stores the multiplicity of iupac bases
    long double product;    // stores the overall multiplicity of the k-mers

    size_t pos;    // current position in the string, from 0 to length
    kmer_t kmer;    // create an empty bit sequence for the initial k-mer
    kmer_t rcmer;    // create a bit sequence for the reverse complement

    size_t begin = 0;
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
        iupac_multiply(product, factors, str[pos]);

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
                if (reverse) kmer::reverse_represent(rcmer);    // invert the k-mer, if necessary
                color::set(kmer_table[T][rcmer], color);    // update the k-mer with the current color
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
 * @param window number of k-mers to minimize
 * @param max_iupac allowed number of ambiguous k-mers per position
 */
void graph::add_minimizers(const uint64_t& T, const string& str, const size1N_t& color, const bool& reverse, const uint64_t& window, const uint64_t& max_iupac) {
    if (str.length() < kmer::k) return;    // not enough characters

    vector<kmer_t> sequence_order;    // k-mers ordered by their position in sequence
    multiset<kmer_t> value_order;    // k-mers ordered by their lexicographical value
    multiset<kmer_t> inner_value_order;

    hash_set<kmer_t> ping;    // create a new empty set for the k-mers
    hash_set<kmer_t> pong;    // create another new set for the k-mers
    bool ball; bool wait;    // indicates which of the two sets should be used

    vector<uint8_t> factors;    // stores the multiplicity of iupac bases
    long double product;    // stores the overall multiplicity of the k-mers

    size_t pos;    // current position in the string, from 0 to length
    kmer_t kmer;    // create an empty bit sequence for the initial k-mer
    kmer_t rcmer;    // create a bit sequence for the reverse complement

    size_t begin = 0;
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
        iupac_multiply(product, factors, str[pos]);

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
                if (reverse) kmer::reverse_represent(rcmer);    // invert the k-mer, if necessary
                inner_value_order.emplace(rcmer);
            }

            if (sequence_order.size() == window) {
                value_order.erase(*sequence_order.begin());    // remove k-mer outside the window
                sequence_order.erase(sequence_order.begin());
            }
            value_order.emplace(*inner_value_order.begin());    // insert k-mer ordered by its lexicographical value
            sequence_order.emplace_back(*inner_value_order.begin());
            inner_value_order.clear();

            if (sequence_order.size() == window) {
                color::set(kmer_table[T][*value_order.begin()], color);    // update the k-mer with the current color
            }
        }
    }
}

/**
 * This function calculates the multiplicity of iupac k-mers.
 *
 * @param product overall multiplicity
 * @param factors per base multiplicity
 * @param chr iupac character
 */
void graph::iupac_multiply(long double& product, vector<uint8_t>& factors, const char& chr) {
    switch (chr) {
        case 'A': case 'C': case 'G': case 'T':
            factors.emplace_back(1);
            product *= 1; break;
        case 'R': case 'Y': case 'S': case 'W': case 'K': case 'M':
            factors.emplace_back(2);
            product *= 2; break;
        case 'B': case 'D': case 'H': case 'V':
            factors.emplace_back(3);
            product *= 3; break;
        case 'N':
            factors.emplace_back(4);
            product *= 4; break;
        default:
            cerr << "Error: invalid character: " << chr << endl;
            exit(1);
    }
    if (factors.size() > kmer::k) {
        product /= *factors.begin();
        factors.erase(factors.begin());
    }
}

/**
 * This function shifts a base into a set of ambiguous iupac k-mers.
 *
 * @param prev set of k-mers
 * @param next set of k-mers
 * @param chr iupac character
 */
void graph::iupac_shift(hash_set<kmer_t>& prev, hash_set<kmer_t>& next, const char& chr) {
    kmer_t temp;
    while (!prev.empty()) {    // extend each previous k-mer
        switch (chr) {
            case 'A': case 'R': case 'W': case 'M':
            case 'D': case 'H': case 'V': case 'N':
                temp = *prev.begin();
                kmer::shift(temp, 'A');
                next.emplace(temp);
        }
        switch (chr) {
            case 'C': case 'Y': case 'S': case 'M':
            case 'B': case 'H': case 'V': case 'N':
                temp = *prev.begin();
                kmer::shift(temp, 'C');
                next.emplace(temp);
        }
        switch (chr) {
            case 'G': case 'R': case 'S': case 'K':
            case 'B': case 'D': case 'V': case 'N':
                temp = *prev.begin();
                kmer::shift(temp, 'G');
                next.emplace(temp);
        }
        switch (chr) {
            case 'T': case 'Y': case 'W': case 'K':
            case 'B': case 'D': case 'H': case 'N':
                temp = *prev.begin();
                kmer::shift(temp, 'T');
                next.emplace(temp);
        }
        prev.erase(prev.begin());
    }
}

/**
 * This function iterates over the hash table and outputs all the k-mer/color pairs.
 *
 * @param kmer string to store the k-mer
 * @param color string to store the color
 * @return iterator function
 */
function<bool(string&, string&)> graph::lookup_kmer() {
    auto it = kmer_table[0].begin();    // create iterator and a function to advance
    return [&, it] (string& kmer, string& color) mutable {
        if (it == kmer_table[0].begin()) {    // initialize placeholder strings
            kmer = string(kmer::k, 'N');    // reserve enough space for characters
            color = string(color::n, '0');    // reserve enough space for colors
        }
        if (it != kmer_table[0].end()) {
            kmer_t kmer_bits = it->first;    // shift the characters into the string
            for (size2K_t i = 0; i != kmer::k; ++i)
                kmer::unshift(kmer_bits, kmer[kmer::k-i-1]);
            color_t color_bits = it->second;    // shift the colors into the string
            for (size1N_t i = 0; i != color::n; ++i)
                color::unshift(color_bits, color[i]);
            ++it; return true;    // suspend and return the results to the caller
        }
        kmer.clear(); color.clear(); return false;    // clean up and finish
    };
}

/**
 * This function iterates over the hash table and outputs matching k-mer/color pairs.
 *
 * @param query query sequence
 * @param reverse merge complements
 * @param kmer string to store the k-mer
 * @param color string to store the color
 * @return iterator function
 */
function<bool(string&, string&)> graph::lookup_kmer(const string& query, const bool& reverse) {
    string str;    // fill shorter queries with N up to the k-mer length
    if (query.length() < kmer::k)
        for (size2K_t i = 0; i < kmer::k-query.size(); ++i)
            str += 'N';
    str += query;    // allows to search shorter strings within a k-mer
    if (query.length() < kmer::k)
        for (size2K_t i = 0; i < kmer::k-query.size(); ++i)
            str += 'N';

    kmer_table.resize(2);    // add a new table to store the results
    hash_set<kmer_t> ping, pong;    // create two sets for the k-mers
    bool ball;    // indicates which of the two sets should be used

    size_t pos;    // current position in the string, from 0 to length
    kmer_t kmer;    // create an empty bit sequence for the initial k-mer
    kmer_t rcmer;    // create a bit sequence for the reverse complement

    size_t begin = 0;
next_kmer:
    pos = begin;

    ping.clear(); pong.clear(); ball = true;
    (ball ? ping : pong).emplace(kmer);

    for (; pos < str.length(); ++pos) {    // collect the bases from the string
        iupac_shift(ball ? ping : pong, !ball ? ping : pong, str[pos]);
        ball = !ball;    // shift each base in, resolve iupac character
        if ((ball ? ping : pong).empty()) {
            begin = pos+1;    // str = str.substr(pos+1, string::npos);
            goto next_kmer;    // invalid character, start a new k-mer from next position
        }

        if (pos+1 - begin >= kmer::k) {
            for (auto& kmer : (ball ? ping : pong)) {    // iterate over the current set of ambiguous k-mers
                rcmer = kmer;
                if (reverse) kmer::reverse_represent(rcmer);    // invert the k-mer, if necessary
                if (kmer_table[0].find(rcmer) != kmer_table[0].end()) {    // lookup the k-mer/color
                    kmer_table[1][rcmer] = kmer_table[0][rcmer];    // transfer to the results table
                }
            }
        }
    }

    auto it = kmer_table[1].begin();    // create iterator and a function to advance
    return [&, it] (string& kmer, string& color) mutable {
        if (it == kmer_table[1].begin()) {    // initialize placeholder strings
            kmer = string(kmer::k, 'N');    // reserve enough space for characters
            color = string(color::n, '0');    // reserve enough space for colors
        }
        if (it != kmer_table[1].end()) {
            kmer_t kmer_bits = it->first;    // shift the characters into the string
            for (size2K_t i = 0; i != kmer::k; ++i)
                kmer::unshift(kmer_bits, kmer[kmer::k-i-1]);
            color_t color_bits = it->second;    // shift the colors into the string
            for (size1N_t i = 0; i != color::n; ++i)
                color::unshift(color_bits, color[i]);
            ++it; return true;    // suspend and return the results to the caller
        }
        kmer_table[1].clear(); kmer_table.resize(1);    // remove results table
        kmer.clear(); color.clear(); return false;    // clean up and finish
    };
}

/**
 * This function iterates over the hash table and calculates the split weights.
 *
 * @param mean weight function
 * @param verbose print progress
 */
void graph::calc_weights(const function<double(const uint32_t&, const uint32_t&)>& mean, const bool& verbose) {
    double min_value = numeric_limits<double>::min();    // current min. weight in the top list (>0)
    uint64_t cur = 0, prog = 0, next;
    uint64_t max = kmer_table[0].size();

    for (auto it = kmer_table[0].begin(); it != kmer_table[0].end(); ++it) {    // iterate over k-mer hash table
        if (verbose) {
            next = 100 * cur / max;
             if (prog < next)  cerr << "\33[2K\r" << "Processing splits... " << next << "%" << flush;
            prog = next; cur++;
        }
        color_t& color = it.value();    // get the color set for each k-mer
        bool pos = color::represent(color);    // invert the color set, if necessary
        if (color == 0b0u) continue;    // ignore empty splits
        array<uint32_t,2>& weight = color_table[color];    // get the weight and inverse weight for the color set

        double old_value = mean(weight[0], weight[1]);    // calculate the old mean value
        if (old_value >= min_value) {    // if it is greater than the min. value, find it in the top list
            auto range = tree::splits.equal_range(old_value);    // get all color sets with the given weight
            for (auto it = range.first; it != range.second; ++it) {
                if (it->second == color) {    // iterate over the color sets to find the correct one
                    tree::splits.erase(it);    // erase the entry with the old weight
                    break;
                }
            }
        }
        weight[pos]++;    // update the weight or the inverse weight of the current color set

        double new_value = mean(weight[0], weight[1]);    // calculate the new mean value
        if (new_value >= min_value) {    // if it is greater than the min. value, add it to the top list
            tree::splits.emplace(new_value, color);    // insert it at the correct position ordered by weight
            if (tree::splits.size() > t) {
                tree::splits.erase(--tree::splits.end());    // if the top list exceeds its limit, erase the last entry
                min_value = tree::splits.rbegin()->first;    // update the min. value for the next iteration
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
void graph::insert_split(const double& weight, const color_t& color) {
    tree::splits.emplace(weight, color);    // insert it at the correct position ordered by weight
    if (tree::splits.size() > t) {
        tree::splits.erase(--tree::splits.end());    // if the top list exceeds its limit, erase the last entry
    }
}

/**
 * This function merges two thread-separate hash tables.
 *
 * @param T1 first thread
 * @param T2 second thread
 */
void graph::merge_threads(const uint64_t& T1, const uint64_t& T2) {
    for (auto it = kmer_table[T2].begin(); it != kmer_table[T2].end();) {
        if (kmer_table[T1].find(it->first) != kmer_table[T1].end()) {
            kmer_table[T1][it->first] |= it->second;
        } else {
            kmer_table[T1][it->first] = it->second;
        }
        it = kmer_table[T2].erase(it);
    }
}

/**
 * This function destructs a thread-separate hash table.
 *
 * @param T thread index
 */
void graph::erase_thread(const uint64_t& T) {
    kmer_table[T].clear();
    kmer_table.erase(kmer_table.begin()+T);
}
