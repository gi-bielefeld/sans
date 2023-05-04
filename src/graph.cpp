#include "graph.h"
#include "ansi.h"

/*
 * This class manages the k-mer/color hash tables and splits calculation.
 */
uint64_t graph::quality;   // min. coverage threshold for k-mers
uint64_t graph::buffer;    // number of k-mers to retrieve from the queue

vector<hash_map<kmer_t, color_t>>  graph::kmer_table;    // [Q] hash tables mapping k-mers to their colors/splits
hash_map<color_t, uint32_t[2]>     graph::color_table;   // [1] hash table mapping colors/splits to their weights
vector<hash_set<kmer_t>>           graph::quality_set;   // [P] hash sets used to filter k-mers for coverage (q > 1)
vector<hash_map<kmer_t, uint64_t>> graph::quality_map;   // [P] hash maps used to filter k-mers for coverage (q > 2)
vector<queue<kmer_t, size1N_t>>    graph::thread_queue;  // [Q] queues used to synchronize k-mers from multiple threads

function<void(kmer_t&)>                                         graph::process_kmer;   // reverse complement a k-mer and apply a gap pattern
function<void(kmer_t&)>                                         graph::restore_kmer;   // restore a gap pattern for a right-compressed k-mer
function<void(const uint64_t&, const kmer_t&, const size1N_t&)> graph::emplace_kmer;   // qualify a k-mer and place it into the hash table

/**
 * This function initializes the thread queues and coverage threshold.
 *
 * @param quality coverage threshold
 * @param reverse merge complements
 * @param P number of file reading threads
 * @param Q number of queue hashing threads
 */
void graph::init(const uint64_t& quality, const bool& reverse, const uint64_t& P, const uint64_t& Q) {

    if (kmer::k == kmer::_k && !reverse) {
        process_kmer = [&] (kmer_t& kmer) {};
        restore_kmer = [&] (kmer_t& kmer) {};
    }
    if (kmer::k == kmer::_k && reverse) {
        process_kmer = [&] (kmer_t& kmer)
          { kmer::reverse_represent(kmer); };
        restore_kmer = [&] (kmer_t& kmer) {};
    }
    if (kmer::k != kmer::_k && !reverse) {
        process_kmer = [&] (kmer_t& kmer)
          { kmer::extract_pattern(kmer); };
        restore_kmer = [&] (kmer_t& kmer)
          { kmer::deposit_pattern(kmer); };
    }
    if (kmer::k != kmer::_k && reverse) {
        process_kmer = [&] (kmer_t& kmer)
          { kmer::extract_pattern(kmer);
            kmer::reverse_represent(kmer); };
        restore_kmer = [&] (kmer_t& kmer)
          { kmer::deposit_pattern(kmer); };
    }
    graph::quality = quality;

    if (Q == 0) { /* single-CPU: read one file at a time, don't use queues, put into one table directly */
        kmer_table.resize(1);
        switch (quality) {
            case 1:
                case 0: /* no quality check */
                emplace_kmer = [&] (const uint64_t& T, const kmer_t& kmer, const size1N_t& color) {
                    kmer_table[0][kmer].set(color);
                };  break;
            case 2:
                quality_set.resize(1);
                emplace_kmer = [&] (const uint64_t& T, const kmer_t& kmer, const size1N_t& color) {
                    if (quality_set[0].find(kmer) == quality_set[0].end()) {
                        quality_set[0].emplace(kmer);
                    } else {
                        quality_set[0].erase(kmer);
                        kmer_table[0][kmer].set(color);
                    }
                };  break;
            default:
                quality_map.resize(1);
                emplace_kmer = [&] (const uint64_t& T, const kmer_t& kmer, const size1N_t& color) {
                    if (quality_map[0][kmer] < quality-1) {
                        quality_map[0][kmer]++;
                    } else {
                        quality_map[0].erase(kmer);
                        kmer_table[0][kmer].set(color);
                    }
                };  break;
        }
    } else { /* multi-CPU: read multiple files in parallel, distribute to multiple queues & tables */
        graph::buffer = max<uint64_t>(P * (1048576 / sizeof(pair<kmer_t, size1N_t>)) / Q, 4);
        for (uint64_t i = 0; i < Q; ++i) thread_queue.emplace_back(graph::buffer);
        graph::buffer = max<uint64_t>(P * (1024 / sizeof(pair<kmer_t, size1N_t>)) / Q, 4);

        kmer_table.resize(Q);
        switch (quality) {
            case 1:
                case 0: /* no quality check */
                emplace_kmer = [&] (const uint64_t& T, const kmer_t& kmer, const size1N_t& color) {
                    const pair<const kmer_t&, const size1N_t&> item(kmer, color);
                    while (!thread_queue[kmer % Q].try_push(item));  // wait for enough space
                };  break;
            case 2:
                quality_set.resize(P);
                emplace_kmer = [&] (const uint64_t& T, const kmer_t& kmer, const size1N_t& color) {
                    if (quality_set[T].find(kmer) == quality_set[T].end()) {
                        quality_set[T].emplace(kmer);
                    } else {
                        quality_set[T].erase(kmer);
                        const pair<const kmer_t&, const size1N_t&> item(kmer, color);
                        while (!thread_queue[kmer % Q].try_push(item));  // wait for enough space
                    }
                };  break;
            default:
                quality_map.resize(P);
                emplace_kmer = [&] (const uint64_t& T, const kmer_t& kmer, const size1N_t& color) {
                    if (quality_map[T][kmer] < quality-1) {
                        quality_map[T][kmer]++;
                    } else {
                        quality_map[T].erase(kmer);
                        const pair<const kmer_t&, const size1N_t&> item(kmer, color);
                        while (!thread_queue[kmer % Q].try_push(item));  // wait for enough space
                    }
                };  break;
        }
    }
}

/**
 * This function extracts k-mers from a sequence and adds them to the hash table.
 *
 * @param T thread id [P]
 * @param str dna sequence
 * @param color color flag
 */
void graph::add_kmers(const uint64_t& T, const string& str, const size1N_t& color) {
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
            process_kmer(rcmer);    // invert the k-mer & apply gap pattern, if necessary
            emplace_kmer(T, rcmer, color);    // update the k-mer with the current color
        }
    }
}

/**
 * This function extracts k-mer minimizers from a sequence and adds them to the hash table.
 *
 * @param T thread id [P]
 * @param str dna sequence
 * @param color color flag
 * @param window number of k-mers to minimize
 */
void graph::add_minimizers(const uint64_t& T, const string& str, const size1N_t& color, const uint64_t& window) {
    if (str.length() < kmer::k) return;    // not enough characters

    vector<kmer_t> sequence_order;    // k-mers ordered by their position in sequence
    multiset<kmer_t> value_order;    // k-mers ordered by their lexicographical value
    multiset<kmer_t>::iterator temp;    // temporary pointer to erased/inserted k-mer

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
            process_kmer(rcmer);    // invert the k-mer & apply gap pattern, if necessary
            bool changed = false;    // indicates if the last min. value k-mer has changed

            if (sequence_order.size() == window) {
                temp = value_order.erase(value_order.find(*sequence_order.begin()));
                changed |= temp == value_order.begin();    // check if old min. value k-mer was removed
                sequence_order.erase(sequence_order.begin());    // remove left k-mer from the window
            } else changed = true;

            temp = value_order.emplace(rcmer);    // insert k-mer ordered by its lexicographical value
            changed |= temp == value_order.begin();    // check if inserted value is new min. value k-mer
            sequence_order.emplace_back(rcmer);    // add right k-mer to the window

            if (sequence_order.size() == window && changed) {    // check if the min. value k-mer has changed
                emplace_kmer(T, *value_order.begin(), color);    // update the k-mer with the current color
            }
        }
    }
}

/**
 * This function extracts k-mers from a sequence and adds them to the hash table.
 *
 * @param T thread id [P]
 * @param str dna sequence
 * @param color color flag
 * @param max_iupac allowed number of ambiguous k-mers per position
 */
void graph::add_kmers(const uint64_t& T, const string& str, const size1N_t& color, const uint64_t& max_iupac) {
    if (str.length() < kmer::k) return;    // not enough characters

    hash_set<kmer_t> ping;    // create a new empty set for the k-mers
    hash_set<kmer_t> pong;    // create another new set for the k-mers
    bool ball; bool wait;    // indicates which of the two sets should be used

    vector<uint8_t> factors;    // stores the multiplicity of IUPAC bases
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
            ball = !ball;    // shift each base in, resolve IUPAC character
        } else { wait = true; continue; }

        if (pos+1 - begin >= kmer::k) {
            for (auto& kmer : (ball ? ping : pong)) {    // iterate over the current set of ambiguous k-mers
                rcmer = kmer;
                process_kmer(rcmer);    // invert the k-mer & apply gap pattern, if necessary
                emplace_kmer(T, rcmer, color);    // update the k-mer with the current color
            }
        }
    }
}

/**
 * This function extracts k-mer minimizers from a sequence and adds them to the hash table.
 *
 * @param T thread id [P]
 * @param str dna sequence
 * @param color color flag
 * @param window number of k-mers to minimize
 * @param max_iupac allowed number of ambiguous k-mers per position
 */
void graph::add_minimizers(const uint64_t& T, const string& str, const size1N_t& color, const uint64_t& window, const uint64_t& max_iupac) {
    if (str.length() < kmer::k) return;    // not enough characters

    vector<kmer_t> sequence_order;    // k-mers ordered by their position in sequence
    multiset<kmer_t> value_order;    // k-mers ordered by their lexicographical value
    multiset<kmer_t> inner_value_order;    // current IUPAC k-mers ordered by value
    multiset<kmer_t>::iterator temp;    // temporary pointer to erased/inserted k-mer

    hash_set<kmer_t> ping;    // create a new empty set for the k-mers
    hash_set<kmer_t> pong;    // create another new set for the k-mers
    bool ball; bool wait;    // indicates which of the two sets should be used

    vector<uint8_t> factors;    // stores the multiplicity of IUPAC bases
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
            ball = !ball;    // shift each base in, resolve IUPAC character
        } else { wait = true; continue; }

        if (pos+1 - begin >= kmer::k) {
            for (auto& kmer : (ball ? ping : pong)) {    // iterate over the current set of ambiguous k-mers
                rcmer = kmer;
                process_kmer(rcmer);    // invert the k-mer & apply gap pattern, if necessary
                inner_value_order.emplace(rcmer);    // find smallest at the current position
            }
            bool changed = false;    // indicates if the last min. value k-mer has changed

            if (sequence_order.size() == window) {
                temp = value_order.erase(value_order.find(*sequence_order.begin()));
                changed |= temp == value_order.begin();    // check if old min. value k-mer was removed
                sequence_order.erase(sequence_order.begin());    // remove left k-mer from the window
            } else changed = true;

            temp = value_order.emplace(*inner_value_order.begin());    // insert k-mer ordered by value
            changed |= temp == value_order.begin();    // check if inserted value is new min. value k-mer
            sequence_order.emplace_back(*inner_value_order.begin());    // add right k-mer to the window
            inner_value_order.clear();

            if (sequence_order.size() == window && changed) {    // check if the min. value k-mer has changed
                emplace_kmer(T, *value_order.begin(), color);    // update the k-mer with the current color
            }
        }
    }
}

/**
 * This function calculates the multiplicity of IUPAC k-mers.
 *
 * @param product overall multiplicity
 * @param factors per base multiplicity
 * @param chr IUPAC character
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
            $err << "Error: invalid character: " << chr << _end$$;
    }
    if (factors.size() > kmer::k) {
        product /= *factors.begin();
        factors.erase(factors.begin());
    }
}

/**
 * This function shifts a base into a set of ambiguous IUPAC k-mers.
 *
 * @param prev set of k-mers
 * @param next set of k-mers
 * @param chr IUPAC character
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
 * This function iterates over the hash tables and outputs all the k-mer/color pairs.
 *
 * @param kmer string to store the k-mer
 * @param color string to store the color
 */
void graph::lookup_kmer(const function<void(string&, string&)>& iterator) {
    string kmer_string(kmer::k, 'N');    // reserve enough space for characters
    string color_string(color::n, '0');    // reserve enough space for color bits

    for (auto& _kmer_table : kmer_table) {
        for (auto it = _kmer_table.begin(); it != _kmer_table.end(); ++it) {
            kmer_t kmer_byte = it->first;    // shift the characters into the string
            restore_kmer(kmer_byte);    // restore the gap pattern, if necessary
            for (size2K_t i = 0; i != kmer::k; ++i)
                kmer::unshift(kmer_byte, kmer_string[kmer::k-i-1]);
            color_t color_byte = it->second;    // shift the color bits into the string
            for (size1N_t i = 0; i != color::n; ++i)
                color::unshift(color_byte, color_string[i]);
            iterator(kmer_string, color_string);    // return the results to the caller
        }
    }
}

/**
 * This function iterates over the hash tables and outputs matching k-mer/color pairs.
 *
 * @param kmer string to store the k-mer
 * @param color string to store the color
 * @param query query sequence
 */
void graph::lookup_kmer(const function<void(string&, string&)>& iterator, const string& query) {
    string str;    // fill shorter queries with N up to the k-mer length
 // if (query.length() < kmer::k)
 //     for (size2K_t i = 0; i < kmer::k-query.size(); ++i)
 //         str += 'N';
    str += query;    // allows to search shorter strings within a k-mer
 // if (query.length() < kmer::k)
 //     for (size2K_t i = 0; i < kmer::k-query.size(); ++i)
 //         str += 'N';
    if (str.length() < kmer::k) return;    // not enough characters

    hash_set<kmer_t> ping, pong;    // create two sets for all possible IUPAC k-mers
    hash_set<kmer_t> players;    // create another set for the unique r.c. k-mers
    bool ball;    // indicates which of the first two sets should be used

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
        ball = !ball;    // shift each base in, resolve IUPAC character
        if ((ball ? ping : pong).empty()) {
            begin = pos+1;    // str = str.substr(pos+1, string::npos);
            goto next_kmer;    // invalid character, start a new k-mer from next position
        }
        if (pos+1 - begin >= kmer::k) {
            for (auto& kmer : (ball ? ping : pong)) {    // iterate over the current set of ambiguous k-mers
                rcmer = kmer;
                process_kmer(rcmer);    // invert the k-mer & apply gap pattern, if necessary
                players.emplace(rcmer);    // add the r.c. k-mer to the set of unique candidates
            }
        }
    }

    string kmer_string(kmer::k, 'N');    // reserve enough space for characters
    string color_string(color::n, '0');    // reserve enough space for color bits
    const uint64_t& Q = kmer_table.size();    // number of thread-separate hash tables

    for (auto& rcmer : players) {
        const auto& query_table = kmer_table[rcmer % Q];
        const auto& it = query_table.find(rcmer);
        if (it != query_table.end()) {    // lookup the k-mer/color in the correct table
            kmer_t kmer_byte = it->first;    // shift the characters into the string
            restore_kmer(kmer_byte);    // restore the gap pattern, if necessary
            for (size2K_t i = 0; i != kmer::k; ++i)
                kmer::unshift(kmer_byte, kmer_string[kmer::k-i-1]);
            color_t color_byte = it->second;    // shift the color bits into the string
            for (size1N_t i = 0; i != color::n; ++i)
                color::unshift(color_byte, color_string[i]);
            iterator(kmer_string, color_string);    // return the results to the caller
        }
    }
}

/**
 * This function iterates over the hash tables and calculates the split weights.
 *
 * @param mean weight function
 * @param verbose print progress
 */
void graph::calc_weights(const function<double(const uint32_t&, const uint32_t&)>& mean, const bool& verbose) {
    double min_value = numeric_limits<double>::min();    // current min. weight in the splits list
    uint64_t cur = 0, max = 0, prog = 0, next;
    for (auto& _kmer_table : kmer_table)
        max += _kmer_table.size();

    for (auto& _kmer_table : kmer_table) {
        for (auto it = _kmer_table.begin(); it != _kmer_table.end(); ++it) {    // iterate over the k-mer table
            if (verbose) {
                next = 100 * cur / max;
                 if (prog < next)  $log $_ << "Processing splits... " << next << "%" << $;
                prog = next; cur++;
            }
            color_t& color = it.value();    // get the color set for each k-mer
            bool pos = color::represent(color);    // invert the color set, if necessary
            if (color != 0b0u) color_table[color][pos]++;    // update the weight or inverse weight
        }
    }
    for (auto it = color_table.begin(); it != color_table.end(); ++it) {    // iterate over the color table
        const auto& color = it.key();    // get the color set of the split
        const auto& weight = it.value();    // get the weights for each split
        double new_value = mean(weight[0], weight[1]);    // calculate the new mean value

        if (new_value >= min_value) {    // if it is greater than the min. value, add it to the top list
            tree::splits.emplace(new_value, color);    // insert it at the correct position ordered by weight
            if (tree::splits.size() > tree::t) {      // if the splits list exceeds its limit, erase the last entry
                tree::splits.erase(--tree::splits.end());
                min_value = tree::splits.rbegin()->first;    // update the min. value for the next iteration
            }
        }
    }
}

/**
 * This function transfers k-mers & colors from the queue to the k-mer table.
 *
 * @param T thread id [Q]
 */
void graph::merge_thread(const uint64_t& T) {
    pair<kmer_t, size1N_t> items[buffer]; size_t items_size;
    do { items_size = thread_queue[T].try_pop_bulk(items, buffer);
         for (uint64_t i = 0; i != items_size; ++i)
             kmer_table[T][items[i].first].set(items[i].second);
    }  while (items_size != 0);   // repeat until queue is empty
}

/**
 * This function clears a quality filter after full procession of a color.
 *
 * @param T thread id [P]
 */
void graph::clear_filter(const uint64_t& T) {
    switch (quality) {
        case 1:  case 0: break;
        case 2:  quality_set[T].clear(); break;
        default: quality_map[T].clear(); break;
    }
}

/**
 * This function erases the quality filters and updates the distributor.
 */
void graph::erase_filters() {
    switch (quality) {
        case 1:  case 0: break;
        case 2:  quality_set.clear(); break;
        default: quality_map.clear(); break;
    }
    const uint64_t& Q = thread_queue.size();    // number of thread-separate hash queues

    if (Q == 0) { /* single-CPU: read one file at a time, don't use queues, put into one table directly */
        emplace_kmer = [&] (const uint64_t& T, const kmer_t& kmer, const size1N_t& color) {
            kmer_table[0][kmer].set(color);
        };
    } else { /* multi-CPU: read multiple files in parallel, distribute to multiple queues & tables */
        emplace_kmer = [&] (const uint64_t& T, const kmer_t& kmer, const size1N_t& color) {
            const pair<const kmer_t&, const size1N_t&> item(kmer, color);
            while (!thread_queue[kmer % Q].try_push(item));  // wait for enough space
        };
    }
}
