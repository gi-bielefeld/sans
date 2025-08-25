#include "graph.h"
#include "util.h"
#include "pctree/PCTree.h"
#include <mutex>
#include <thread>
#include <algorithm>

/**
 * This is the size of the top list.
 */
uint64_t graph::t;

/**
 * This is the min. coverage threshold for k-mers.
 */
vector<int> graph::q_table;
int graph::quality;


bool graph::isAmino;

/*
*
* [Parallelization]
*
*/

uint64_t graph::table_count;

/**
 * This is a vecotr of spinlocks protecting the hash maps 
 */
vector<spinlock> graph::lock;

/**
 * This vector holds the carries of 2**i % table_count for fast distribution of binary represented kmers
 */
vector<uint_fast32_t> graph::period;

/**
 * This is vector of hash tables mapping k-mers to colors [O(1)].
 */
vector<hash_map<kmer_t, color_t>> graph::kmer_table;

/**
 * This is the amino equivalent.
 */ 
vector<hash_map<kmerAmino_t, color_t>> graph::kmer_tableAmino;

/**
 * This is a hash table mapping colors to weights [O(1)].
 */
hash_map<color_t, array<uint32_t,2>> graph::color_table;

/**
 * This is a hash set used to filter k-mers for coverage (q > 1).
 */
vector<hash_set<kmer_t>> graph::quality_set;

/**
 * This is a hash map used to filter k-mers for coverage (q > 2).
 */
vector<hash_map<kmer_t, uint16_t>> graph::quality_map;

/**
 * Look-up set for k-mers that are ignored, i.e., not stored, counted etc.
 */
hash_set<kmer_t> graph::blacklist;
hash_set<kmerAmino_t> graph::blacklist_amino;


/**
 * This is vector of hash tables mapping k-mers to genomes to buffer a k-mer before adding to the kmer_table. If it is seen a second time, it is added. Otherwise the singleton k-mer is ignored
 */
vector<hash_map<kmer_t, uint16_t>> graph::singleton_kmer_table;
vector<hash_map<kmerAmino_t, uint16_t>> graph::singleton_kmer_tableAmino;
uint64_t graph::singleton_counters[maxN];
spinlock graph::singleton_counters_locks[maxN];

/**
 * Count unique kmers per genome
 */
uint64_t graph::kmer_counters[maxN];



/**
 * This is a hash set used to filter k-mers for coverage (q > 1).
 */
vector<hash_set<kmerAmino_t>> graph::quality_setAmino;

/**
 * This is a hash map used to filter k-mers for coverage (q > 2).
 */
vector<hash_map<kmerAmino_t, uint16_t>> graph::quality_mapAmino;

/**
 * This is an ordered tree collecting the splits [O(log n)].
 */
// https://itecnote.com/tecnote/c-sorting-multimap-with-both-keys-and-values/
// This is necessairy to create a sorted output
multimap_<double, color_t> graph::split_list;

/**
* These are the allowed chars.
*/
vector<char> graph::allowedChars;

/**
 * This function qualifies a k-mer and places it into the hash table.
 */
function<void(const uint64_t& T, uint_fast32_t& bin, const kmer_t&, const uint16_t&)> graph::emplace_kmer;
function<void(const uint64_t& T, uint_fast32_t& bin, const kmer_t&, const uint16_t&)> graph::emplace_kmer_tmp;
function<void(const uint64_t& T, uint_fast32_t& bin, const kmerAmino_t&, const uint16_t&)> graph::emplace_kmer_amino;
function<void(const uint64_t& T, uint_fast32_t& bin, const kmerAmino_t&, const uint16_t&)> graph::emplace_kmer_amino_tmp;

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
 * This function activates the use of the blacklist when inserting k-mers. It has to be separated from the init function, because when the latter is called, the blacklist is still empty. 
 */
void graph::activate_blacklist(){
    // Black list for kmers given?
    if ((!isAmino && !blacklist.empty()) || (isAmino && !blacklist_amino.empty())) {
		emplace_kmer = [&] (const uint64_t& T, uint_fast32_t& bin, const kmer_t& kmer, const uint16_t& color) {
			// only add if kmer not in blacklist
			if (graph::blacklist.find(kmer)==graph::blacklist.end()) {
				emplace_kmer_tmp(T, bin, kmer, color);
			}
		};
		emplace_kmer_amino = [&] (const uint64_t& T, uint_fast32_t& bin, const kmerAmino_t& kmer, const uint16_t& color) {
			// only add if kmer not in blacklist
			if (graph::blacklist_amino.find(kmer)==graph::blacklist_amino.end()) {
				emplace_kmer_amino_tmp(T, bin, kmer, color);
			}
		};
	}
}


/**
* --- [Hash map access] ---
* The following methods are used to access the entries of the vectorized hash maps
*/ 

/**
 * This method shift-updates the bin of a kmer
 */

uint_fast32_t graph::shift_update_bin(uint_fast32_t& bin, uint_fast8_t& left, uint_fast8_t& right)
{
        return (8 * table_count // Bias
            + 4 * (bin - period[2*kmer::k-1] * (left / 2) - (left % 2) * period[2*kmer::k - 2]) // Shift
            + period[1] * (right / 2) + period[0] * (right % 2)) // Update   
            % table_count; // Mod
}

/**
 * This method shift updates the reverse complement bin of a kmer
 */
uint_fast32_t graph::shift_update_rc_bin(uint_fast32_t& rc_bin, uint_fast8_t& left, uint_fast8_t& right)
{   
    // Remove
    rc_bin += 8 * table_count - period[1] * (!(left / 2 )) - period[0] * (!(left % 2 ));
    // First shift
    if (rc_bin & 0b1u) {rc_bin += table_count;}
    rc_bin >>= 1;
    // Second shift
    if (rc_bin & 0b1u) {rc_bin += table_count;}
    rc_bin >>= 1;
    // Update
    rc_bin += (period[2*kmer::k-1] * (!(right / 2)) + period[2*kmer::k-2] * (!(right % 2)));
    rc_bin %= table_count;
    return rc_bin;
}


/**
 * This method shift-updates the bin of an amino kmer
 */
uint_fast32_t graph::shift_update_amino_bin(uint_fast32_t& bin, kmerAmino_t& kmer, uint_fast8_t& right)
{
    // update the binning carry (solution of the shift-update-carry equation)
    // shift
    bin = 160 * table_count + // Bias 
                32 * (bin // Shift
                - kmer.test(5*kmerAmino::k - 1) * period[5*kmerAmino::k - 1]
                - kmer.test(5*kmerAmino::k - 2) * period[5*kmerAmino::k - 2]
                - kmer.test(5*kmerAmino::k - 3) * period[5*kmerAmino::k - 3]
                - kmer.test(5*kmerAmino::k - 4) * period[5*kmerAmino::k - 4]
                - kmer.test(5*kmerAmino::k - 5) * period[5*kmerAmino::k - 5]);
    // update
    for(int i = 4; i>=0; i--){
        bin += period[i] * ((right >> i) & 0b1u);
    }
    // mod
    bin %= table_count;
    return bin;
}

/**
 * This method computes the bin of a given kmer(slower than shift update)
 * @param kmer The target kmer
 * @return uint64_t The bin
 */
#if (maxK <= 32)
uint_fast32_t graph::compute_bin(const kmer_t& kmer)
{
    return kmer % table_count;
}
#else
uint_fast32_t graph::compute_bin(const kmer_t& kmer)
{
	uint_fast32_t carry = 1;
	uint_fast32_t rest = 0;
	if (kmer.test(0)){rest++;} // Test the last bit
	for (uint_fast32_t it=1; it < 2 * kmer::k; it++){
	    carry = (2*carry) % table_count;
	    if (kmer.test(it)){rest += carry;}
	}
	return rest % table_count;
}
#endif

#if maxK <= 12
    uint_fast32_t graph::compute_amino_bin(const kmerAmino_t& kmer)
    {
        return kmer % table_count;
    }
#else
    uint_fast32_t graph::compute_amino_bin(const kmerAmino_t& kmer)
    {
	    uint_fast32_t carry = 1;
	    uint_fast32_t rest = 0;

	    if (kmer.test(0)){rest++;} // Test the last bit
	    for (uint_fast32_t it=1; it < 5* kmerAmino::k; it--){
	        carry = (2*carry) % table_count;
	        if (kmer.test(it)){rest += carry;}
	    }
	    return rest % table_count;
    }
#endif



/**
 * This function searches the corresponding hash table for the given kmer
 * @param kmer The kmer to search
 */
bool graph::search_kmer(const kmer_t& kmer)
{    
    return kmer_table[compute_bin(kmer)].contains(kmer);
}

/** 
 * This function searches the correstponding hash table for the given amino kmer
 * @param kmer The amino kmer to search
 * @return bool Kmer exists as key in the corresponding map
 */
bool graph::search_kmer_amino(const kmerAmino_t& kmer)
{
    return kmer_tableAmino[compute_amino_bin(kmer)].contains(kmer);
}


/**
* This function returns the stored colores of the given kmer
* @param kmer The target kmer
* @return color_t The stored colores
*/
color_t graph::get_color(const kmer_t& kmer, bool reversed){
    return kmer_table[compute_bin(kmer)][kmer];
}


/**
 * This function returns the stored color vector of the given amino kmer
 * @param kmer The target amino kmer
 * return color_t The stored color vector
 */
color_t graph::get_color_amino(const kmerAmino_t& kmer){
    return kmer_tableAmino[compute_amino_bin(kmer)][kmer];
}

/**
 * This function erases a stored kmer from its corresponding hash table
 * @param kmer The kmer to remove
 */
void graph::remove_kmer(const kmer_t& kmer, bool reversed){
    kmer_table[compute_bin(kmer)].erase(kmer);
}


/**
 * This function erases a stored amino kmer from its corresponding hash table
 * @param kmer The kmer to remove
 */
void graph::remove_kmer_amino(const kmerAmino_t& kmer){
    kmer_tableAmino[compute_amino_bin(kmer)].erase(kmer);
}

/*
*
* [Sequence processing]
*
*/


/**
 * This function extracts k-mers from a sequence and adds them to the black list.
 *
 * @param str sequence
 * @param reverse merge complements
 */
void graph::fill_blacklist(string& str, bool& reverse) {
    if (str.length() < kmer::k) return;    // not enough characters

    uint64_t pos;    // current position in the string, from 0 to length

    kmer_t kmer;    // create a new empty bit sequence for the k-mer
    kmer_t rcmer; // create a bit sequence for the reverse complement


    uint_fast8_t left;  // The character that is shifted out 
    uint_fast8_t right; // The binary code of the character that is shifted in

    kmerAmino_t kmerAmino=0;    // create a new empty bit sequence for the k-mer

    uint64_t begin = 0;
    
    next_kmer:

    pos = begin;
    for (; pos < str.length(); ++pos) {    // collect the bases from the string
        if (!isAllowedChar(pos, str)) {
            begin = pos+1;    // str = str.substr(pos+1, string::npos);
            goto next_kmer;    // unknown base, start a new k-mer from the beginning
        }
        // DNA processing 
        if (!isAmino) {

            right = util::char_to_bits(str[pos]);
            #if maxK <= 32
                kmer::shift(kmer, right); // shift each base into the bit sequence
                rcmer = kmer;
	            if (reverse){
                    kmer::reverse_complement(rcmer); // invert the k-mer
                }
            #else
                kmer::shift(kmer, right); // shift each base into the bit sequence
                rcmer = kmer;
                if (reverse){
                    kmer::reverse_complement(rcmer);
                }
            #endif
             // If the current word is a k-mer
            if (pos+1 - begin >= kmer::k) {
                rcmer < kmer ? blacklist.emplace(rcmer) : blacklist.emplace(kmer);
            }
        
        // Amino processing
        } else {
            #if maxK <= 12
                kmerAmino::shift_right(kmerAmino, str[pos]);    // shift each base into the bit sequence
            #else
                kmerAmino::shift_right(kmerAmino, str[pos]);
            #endif
            // The current word is a k-mer
            if (pos+1 - begin >= kmerAmino::k) {
                // Insert the k-mer
                blacklist_amino.emplace(kmerAmino);
            }
        }
    }
}


/**
* This function tells how many k-mers are in the black list.
*/
uint64_t graph::size_blacklist(){
	return isAmino ? blacklist_amino.size() : blacklist.size();
}


/**
 * This function extracts k-mers from a sequence and adds them to the hash table.
 *
 * @param str dna sequence
 * @param color color flag
 * @param reverse merge complements
 */
void graph::add_kmers(uint64_t& T, string& str, uint16_t& color, bool& reverse) {
    if (str.length() < kmer::k) return;    // not enough characters

    uint_fast32_t bin = 0; // current hash_map vector index
    uint_fast32_t rc_bin = 0; // current reverse hash_map vector index

    uint64_t pos;    // current position in the string, from 0 to length

    kmer_t kmer;    // create a new empty bit sequence for the k-mer
    kmer_t rcmer; // create a bit sequence for the reverse complement


    uint_fast8_t left;  // The character that is shifted out 
    uint_fast8_t right; // The binary code of the character that is shifted in

    #if maxK > 32
    if (!isAmino){
        for (int i =0; i < 2* kmer::k; i++){rc_bin += period[i];}
        rc_bin %= table_count;
    }
    #endif

    kmerAmino_t kmerAmino=0;    // create a new empty bit sequence for the k-mer

    uint64_t begin = 0;
    
    next_kmer:

    pos = begin;
    for (; pos < str.length(); ++pos) {    // collect the bases from the string
        if (!isAllowedChar(pos, str)) {
            begin = pos+1;    // str = str.substr(pos+1, string::npos);
            goto next_kmer;    // unknown base, start a new k-mer from the beginning
        }
        // DNA processing 
        if (!isAmino) {
            right = util::char_to_bits(str[pos]);
            #if maxK <= 32
                kmer::shift(kmer, right); // shift each base into the bit sequence
                rcmer = kmer;

                bin = kmer % table_count; // update the forward bin
	            if (reverse){
                    kmer::reverse_complement(rcmer); // invert the k-mer
                    rc_bin = rcmer % table_count;
                }
            #else
                left = 2*kmer.test(2*kmer::k-1)+kmer.test(2*kmer::k-2); // old leftmost character
                bin = shift_update_bin(bin, left, right); // Shift update the forward complement bin
                
                kmer::shift(kmer, right); // shift each base into the bit sequence
                rcmer = kmer;
                if (reverse){
                    kmer::reverse_complement(rcmer);
                    rc_bin = shift_update_rc_bin(rc_bin, left, right);  // Update the reverse complement table index
                }
            #endif
             // If the current word is a k-mer
            if (pos+1 - begin >= kmer::k) {
                rcmer < kmer ? emplace_kmer(T, rc_bin, rcmer, color) : emplace_kmer(T, bin, kmer, color);
            }
        
        // Amino processing
        } else {
            right = util::amino_char_to_bits(str[pos]);
            #if maxK <= 12
                kmerAmino::shift_right(kmerAmino, str[pos]);    // shift each base into the bit sequence
                bin = kmerAmino % table_count;
            #else
                bin = shift_update_amino_bin(bin, kmerAmino, right);
                kmerAmino::shift_right(kmerAmino, str[pos]);
            #endif
            // The current word is a k-mer
            if (pos+1 - begin >= kmerAmino::k) {
                // shift update the bin
                // Insert the k-mer into its table
                emplace_kmer_amino(T, bin, kmerAmino, color);  // update the k-mer with the current color
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
void graph::add_minimizers(uint64_t& T, string& str, uint16_t& color, bool& reverse, uint64_t& m) {
    if (str.length() < (!isAmino ? kmer::k : kmerAmino::k)) return;    // not enough characters

    vector<kmer_t> sequence_order;    // k-mers ordered by their position in sequence
    multiset<kmer_t> value_order;    // k-mers ordered by their lexicographical value

    vector<kmerAmino_t> sequence_order_Amino;    // k-mers ordered by their position in sequence
    multiset<kmerAmino_t> value_order_Amino;    // k-mers ordered by their lexicographical value

    uint64_t pos;    // current position in the string, from 0 to length
    kmer_t kmer;    // create a new empty bit sequence for the k-mer
    kmer_t rcmer;    // create a bit sequence for the reverse complement

    kmerAmino_t kmerAmino=0;    // create a new empty bit sequence for the k-mer

    uint64_t begin = 0;
next_kmer:
    pos = begin;
    sequence_order.clear();
    sequence_order_Amino.clear();
    value_order.clear();
    sequence_order_Amino.clear();
    
    uint_fast32_t bin = 0;
    uint_fast32_t amino_bin = 0;

    for (; pos < str.length(); ++pos) {    // collect the bases from the string
        if (!isAllowedChar(pos, str)) {
            begin = pos+1;    // str = str.substr(pos+1, string::npos);
            goto next_kmer;    // unknown base, start a new k-mer from the beginning
        }
        if (!isAmino) {
            kmer::shift(kmer, str[pos]);    // shift each base into the bit sequence

            if (pos+1 - begin >= kmer::k) {
                rcmer = kmer;
		        // Test for multitables
		        bool reversed = false;
                if (reverse) {reversed = kmer::reverse_represent(rcmer);}    // invert the k-mer, if necessary

                if (sequence_order.size() == m) {
                    value_order.erase(*sequence_order.begin());    // remove k-mer outside the window
                    sequence_order.erase(sequence_order.begin());
                }
                value_order.emplace(rcmer);    // insert k-mer ordered by its lexicographical value
                sequence_order.emplace_back(rcmer);

                if (sequence_order.size() == m) {
                    bin = compute_bin(*value_order.begin());
                    emplace_kmer(T, bin, *value_order.begin(), color);    // update the k-mer with the current color
                }
            }
        } else {
            kmerAmino::shift_right(kmerAmino, str[pos]);    // shift each base into the bit sequence
            bin = compute_amino_bin(kmerAmino);
            if (pos+1 - begin >= kmerAmino::k) {
                if (sequence_order.size() == m) {
                    value_order_Amino.erase(*sequence_order_Amino.begin());    // remove k-mer outside the window
                    sequence_order_Amino.erase(sequence_order_Amino.begin());
                }
                value_order_Amino.emplace(kmerAmino);    // insert k-mer ordered by its lexicographical value
                sequence_order_Amino.emplace_back(kmerAmino);

                if (sequence_order_Amino.size() == m) {
                    // Update the minimizer in the corresponding table
                    emplace_kmer_amino(T, bin, *value_order_Amino.begin(), color);    // update the k-mer with the current color
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
void graph::add_kmers(uint64_t& T, string& str, uint16_t& color, bool& reverse, uint64_t& max_iupac) {
    if (str.length() < (!isAmino ? kmer::k : kmerAmino::k)) return;    // not enough characters

    uint_fast32_t bin = 0;

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

                    if (reverse) kmer::reverse_represent(rcmer);    // invert the k-mer, if necessary
                    bin = compute_bin(rcmer);
                    emplace_kmer(T, bin, rcmer, color);    // update the k-mer with the current color
                }
            }
        }
    } 
    else {
        hash_set<kmerAmino_t> ping;    // create a new empty set for the k-mers
        hash_set<kmerAmino_t> pong;    // create another new set for the k-mers
        bool ball; bool wait;    // indicates which of the two sets should be used

        vector<uint8_t> factors;    // stores the multiplicity of iupac bases
        long double product;    // stores the overall multiplicity of the k-mers

        uint64_t pos;    // current position in the string, from 0 to length
        kmerAmino_t kmer=0;    // create an empty bit sequence for the initial k-mer

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
                    bin = compute_amino_bin(kmer);  
                    emplace_kmer_amino(T, bin, kmer, color);  // update the k-mer with the current color
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
void graph::add_minimizers(uint64_t& T, string& str, uint16_t& color, bool& reverse, uint64_t& m, uint64_t& max_iupac) {
    if (str.length() < (!isAmino ? kmer::k : kmerAmino::k)) return;    // not enough characters

    uint_fast32_t bin = 0;

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
		   if (reverse) kmer::reverse_represent(rcmer);    // invert the k-mer, if necessary
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
                    bin = compute_bin(*value_order.begin());
                    emplace_kmer(T, bin, *value_order.begin(), color);    // update the k-mer with the current color
               }
           }
       }
   } 
   else {
       vector<kmerAmino_t> sequence_order;    // k-mers ordered by their position in sequence
       multiset<kmerAmino_t> value_order;    // k-mers ordered by their lexicographical value
       multiset<kmerAmino_t> inner_value_order;

       hash_set<kmerAmino_t> ping;    // create a new empty set for the k-mers
       hash_set<kmerAmino_t> pong;    // create another new set for the k-mers
       bool ball; bool wait;    // indicates which of the two sets should be used

       vector<uint8_t> factors;    // stores the multiplicity of iupac bases
       long double product;    // stores the overall multiplicity of the k-mers

       uint64_t pos;    // current position in the string, from 0 to length
       kmerAmino_t kmer=0;    // create an empty bit sequence for the initial k-mer

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
                   // Todo: Get the target hash map index from the kmer bits
                   bin = compute_amino_bin(*value_order.begin());
                   emplace_kmer_amino(T, bin, *value_order.begin(), color);    // update the k-mer with the current color
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
                kmer::shift(temp, base);
                next.emplace(temp);
        }
        switch (input) {
            case 'C': case 'Y': case 'S': case 'M':
            case 'B': case 'H': case 'V': case 'N':
                temp = *prev.begin(); base = 'C';
                kmer::shift(temp, base);
                next.emplace(temp);
        }
        switch (input) {
            case 'G': case 'R': case 'S': case 'K':
            case 'B': case 'D': case 'V': case 'N':
                temp = *prev.begin(); base = 'G';
                kmer::shift(temp, base);
                next.emplace(temp);
        }
        switch (input) {
            case 'T': case 'Y': case 'W': case 'K':
            case 'B': case 'D': case 'H': case 'N':
                temp = *prev.begin(); base = 'T';
                kmer::shift(temp, base);
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
        } //second: all ambigious acids

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
 * This function clears color-related temporary files.
 */
void graph::clear_thread(uint64_t& T) {
    switch (quality) {
        case 1:  case 0: break;
        case 2:  quality_set[T].clear(); break;
        default: quality_map[T].clear(); break;
    }
}

/**
 * This function computes a color table entry from the current kmer map and a cdbg colored kmer. 
 * (To call befor add_weights)
 * @param seq kmer
 * @param kmer_color the split colors
 * @param min_value the minimal weight represented in the top list
 */
void graph::add_cdbg_colored_kmer(string kmer_seq, const uint16_t& kmer_color){
    
        kmer_t kmer; // create a kmer to search in the set of tables

        for (int pos=0; pos < kmer_seq.length(); ++pos) // collect the bases from the k-mer sequence.
        {
            kmer::shift(kmer, kmer_seq[pos]);
        }

	    bool reversed = kmer::reverse_represent(kmer);

		uint_fast32_t bin = compute_bin(kmer);
		hash_kmer<true>(bin, kmer, kmer_color);    // update the k-mer with the current color

}

/*
*
* [Color table processing]
*
*/

/**
 * This function iterates over the hash table and calculates the split weights.
 * 
 * @param mean weight function
 * @param min_value the minimal weight represented in the top list
 * @param verbose print progess
 */
void graph::add_weights(double mean(uint32_t&, uint32_t&), double min_value, bool& verbose) {
	
	
	
    //double min_value = numeric_limits<double>::min(); // current min. weight in the top list (>0)
    uint64_t cur=0, prog=0, next;

    // check table (Amino or base)
    uint64_t max = 0; // table size
    if (isAmino){for (auto table: kmer_tableAmino){max += table.size();}} // use the sum of amino table sizes
    else {for (auto table: kmer_table){max+=table.size();}} // use the sum of base table sizes

    // If the tables are empty, there is nothing to be done	    
    if (max==0){
        return;
    }
    // The iterators for the tables
    hash_map<kmer_t, color_t>::iterator base_it;
    hash_map<kmerAmino_t, color_t>::iterator amino_it;

    // Iterate the tables
    for (int i = 0; i < graph::table_count; i++) // Iterate all tables
    {
        if (!isAmino){base_it = kmer_table[i].begin();} // base table iterator
        else {amino_it = kmer_tableAmino[i].begin();} // amino table iterator

        while (true) { // process splits
            // show progress
            if (verbose) { 
                next = 100*cur/max;
                if (prog < next)  cout << "\33[2K\r" << "Accumulating splits from non-singleton k-mers... " << next << "%" << flush;
                prog = next; cur++;
            }
            // update the iterator
            color_t* color_ref; // reference of the current color
            if (isAmino) { // if the amino table is used, update the amino iterator
                
                if (amino_it == kmer_tableAmino[i].end()){break;} // stop iterating if done
                else{color_ref = &amino_it.value(); ++amino_it;} // iterate the amino table
                }
            else { // if the base tables is used update the base iterator
                // Todo: Get the target hash map index from the kmer bits
                if (base_it == kmer_table[i].end()){break;} // stop itearating if done
                else {color_ref = &base_it.value(); ++base_it;} // iterate the base table
                }
            // process
            color_t& color = *color_ref;
            bool pos = color::represent(color);    // invert the color set, if necessary
            if (color == 0) continue;    // ignore empty splits
            // add_weight(color, mean, min_value, pos);
			array<uint32_t,2>& weight = color_table[color];    // get the weight and inverse weight for the color set
			weight[pos]++; // update the weight or the inverse weight of the current color set
		}
    }
}




/**
 * This function iterates over the singleton tables and adds the split weights.
 * 
 * @param mean weight function
 * @param min_value the minimal weight represented in the top list
 * @param verbose print progess
 */
void graph::add_singleton_weights(double mean(uint32_t&, uint32_t&), double min_value, bool& verbose) {
	
	// not needed anymore
	singleton_kmer_table.clear();
	singleton_kmer_tableAmino.clear();	
	
    //double min_value = numeric_limits<double>::min(); // current min. weight in the top list (>0)
    uint64_t cur=0, prog=0, next;

    // check table (Amino or base)
    uint64_t max = number_singleton_kmers();

    // If the tables are empty, there is nothing to be done	    
    if (max==0){
        return;
    }
    
    // The iterators for the tables
	color_t color;
		
    // Iterate the counters
    for (int i = 0; i < maxN; i++) // Iterate all tables
    {
            // show progress
            if (verbose) { 
                next = 100*cur/max;
                if (prog < next)  cout << "\33[2K\r" << "Accumulating splits from singleton k-mers... " << next << "%" << flush;
                prog = next; cur+=singleton_counters[i];
            }
            if (singleton_counters[i]==0) continue;
// 			cerr << i << ": " << singleton_counters[i] << " " << endl << flush;
			color = 0b0u;
			color.set(i);
            // process
            // add_weight(color, mean, min_value, pos);
			array<uint32_t,2>& weight = color_table[color];    // get the weight and inverse weight for the color set
			weight[0]+=singleton_counters[i]; // update the weight or the inverse weight of the current color set
    }
}



/**
 * This function calculates the weight for all splits and puts them into the split_ölist
 * @param mean weight function
 * @param min_value the minimal weight represented in the top list
 */
void graph::compile_split_list(double mean(uint32_t&, uint32_t&), double min_value)
{
	// Iterating over the map using Iterator till map end.
	hash_map<color_t, array<uint32_t,2>>::iterator it = color_table.begin();
	while (it != color_table.end())	{

		// Accessing the key
		color_t colors = it->first;
		
		// Accessing the value
		array<uint32_t,2> weights = it->second;
		
		//insert into split list
		double new_mean = mean(weights[0], weights[1]);    // calculate the mean value
		if (new_mean >= min_value) {    // if it is greater than the min. value, add it to the top list
			split_list.emplace(new_mean, colors);    // insert it at the correct position ordered by weight
			if (split_list.size() > t) {
				split_list.erase(--split_list.end());    // if the top list exceeds its limit, erase the last entry
				min_value = split_list.rbegin()->first;    // update the min. value for the next iteration (only necessary of t is exceeded, otherwise min_value does not play a role.
			}
		}
		
		// iterator incremented to point next item
		it++;
	}
}

/**
 * This function determines the core k-mers, i.e., all k-mers present in all genomes.
 * Core k-mers are output to given file in fasta format, one k-mer per entry
 * @param file output file stream
 * @param verbose print progess
 */
void graph::output_core(ostream& file, bool& verbose)
{
    uint64_t cur=0, prog=0, next, core_count=0, all_count=0, singletons_count=0;

    // check table (Amino or base)
    uint64_t max = 0; // table size
    if (isAmino){for (auto table: kmer_tableAmino){max += table.size();}} // use the sum of amino table sizes
    else {for (auto table: kmer_table){max+=table.size();}} // use the sum of base table sizes

    // If the tables are empty, there is nothing to be done	    
    if (max==0){
        return;
    }
    // The iterators for the tables
    hash_map<kmer_t, color_t>::iterator base_it;
    hash_map<kmerAmino_t, color_t>::iterator amino_it;

    // Iterate the tables
    for (int i = 0; i < graph::table_count; i++) // Iterate all tables
    {
        if (!isAmino){base_it = kmer_table[i].begin();} // base table iterator
        else {amino_it = kmer_tableAmino[i].begin();} // amino table iterator

        while (true) { // process splits
            // show progress
            if (verbose) { 
                next = 100*cur/max;
                if (prog < next)  cout << "\33[2K\r" << "Collecting core k-mers... " << next << "%" << flush;
                prog = next; cur++;
            }
            // update the iterator
            color_t* color_ref; // reference of the current color
            kmer_t kmer;
			kmerAmino_t kmerAmino;
            if (isAmino) { // if the amino table is used, update the amino iterator
                if (amino_it == kmer_tableAmino[i].end()){break;} // stop iterating if done
                else{kmerAmino = amino_it.key(); color_ref = &amino_it.value(); ++amino_it;} // iterate the amino table
            }
            else { // if the base tables is used update the base iterator
                // Todo: Get the target hash map index from the kmer bits
                if (base_it == kmer_table[i].end()){break;} // stop itearating if done
                else {kmer = base_it.key(); color_ref = &base_it.value(); ++base_it;} // iterate the base table
            }
            // process
            color_t& color = *color_ref;
			all_count++;
			// is core?
			if(color::is_complete(color)){
				core_count++;
				//output
				file << ">" << endl;
 				file << (isAmino?(kmerAmino::kmer_to_string(kmerAmino)):(kmer::kmer_to_string(kmer))) << endl;
			}
			if(color::is_singleton(color)){
				singletons_count++;
			}
		}
    }
	if (verbose) { 
		cout  << "\33[2K\r" << "Collecting core k-mers... (" << core_count << " / "<< (100*core_count/all_count) << "%)"<< flush;
	}
}



/**
 * Get the number of k-mers in all tables.
 * @return number of k-mers in all tables.
 */
uint64_t graph::number_kmers(){
	uint64_t num=0;
	if (isAmino){ // use the sum of amino table sizes
		for (auto table: kmer_tableAmino){num += table.size();}
	} else { // use the sum of base table sizeskmer_table.size(); 
		for (auto table: kmer_table){num+=table.size();}
	}
	return num;
}


/**
 * Get the number of singleton k-mers in all tables.
 * @return number of k-mers in all singleton kmer tables.
 */
uint64_t graph::number_singleton_kmers(){
	uint64_t num=0;
	for (uint16_t g=0;g<maxN-1;g++){num += singleton_counters[g];}
	return num;
}


/**
 * Get the number of singleton k-mers in a given genome.
 * @param genome the id (rf. denom_names) of the genome
 * @return number of singleton k-mers read in that genome.
 */
uint64_t graph::number_singleton_kmers(uint16_t& genome){
	return singleton_counters[genome];
}


/**
 * Get the number of unique k-mers in a given genome.
 * @param genome the id (rf. denom_names) of the genome
 * @return number of unique k-mers read in that genome.
 */
uint64_t graph::number_kmers(uint16_t& genome){
	return kmer_counters[genome];
}


/**
 * @param n number of genomes
 * @return mean value of number of k-mers per genome
 */
double graph::mean_number_kmers(int n){
    return util::mean(kmer_counters,n);
}

/**
 * @param mu mean value
 * @param n number of genomes
 * @return standard deviation of number of k-mers per genome
 */
double graph::stdev_number_kmers(double mu, int n){
    return util::stdev(kmer_counters,mu,n);
}




/**
 * This function generates a bootstrap replicate. We mimic drawing n k-mers at random with replacement from all n observed k-mers. Say a k-mer would be drawn x times. Instead, we calculate x for each k-mer (in each split in color_table) from a binomial distribution (n repetitions, 1/n success rate) and calculate a new split weight according to the new number of k-mers.
 * @param mean weight function
 * @return the new list of splits of length at least t ordered by weight as usual
 */
multimap_<double, color_t> graph::bootstrap(double mean(uint32_t&, uint32_t&)) {

	uint64_t max = graph::number_kmers();

	std::random_device rd;
	std::mt19937 gen(rd());

	multimap_<double, color_t> sl;
	double min_value=0;

	// perform n time max trials, each succeeds 1/max
	std::binomial_distribution<> d(max, 1.0/max);
	
	// Iterating over the map using Iterator till map end.
	hash_map<color_t, array<uint32_t,2>>::iterator it = color_table.begin();
	while (it != color_table.end())	{

		// Accessing the key
		color_t colors = it->first;
		
		// Accessing the value
		array<uint32_t,2> weights = it->second;
		
		// bootstrap the number of kmer occurrences for split and inverse
		array<uint32_t,2> new_weights;
		new_weights[0]=0;
		new_weights[1]=0;
		for (int i=0;i<2;i++) {
			for (int r=0;r<weights[i];r++){
				new_weights[i] += d(gen);
			}
// 			uint64_t n = weights[i]*max;
// 			std::binomial_distribution<> dn(n, 1.0/max);
// 			cout << weights[i] << "\t" << dn(gen) << "\t" << new_weights[i] << "\n" << flush;
		}
		
		//insert into new split list
		double new_mean = mean(new_weights[0], new_weights[1]);    // calculate the new mean value
		if (new_mean >= min_value) {    // if it is greater than the min. value, add it to the top list
			sl.emplace(new_mean, colors);    // insert it at the correct position ordered by weight
			if (sl.size() > t) {
				sl.erase(--sl.end());    // if the top list exceeds its limit, erase the last entry
				min_value = sl.rbegin()->first;    // update the min. value for the next iteration (only necessary of t is exceeded, otherwise min_value does not play a role.
			}
		}
		
		// iterator incremented to point next item
		it++;
		
	}

	return sl;
}


/**
 * This function calls add_split onthe global split_list
 *
 * @param weight split weight
 * @param color split colors
 */
void graph::add_split(double& weight, color_t& color) {
	add_split(weight, color, split_list);
}

/**
 * This function adds a single split (weight and colors) to the output list.
 *
 * @param weight split weight
 * @param color split colors
 * @param split_list list of splits
 */
void graph::add_split(double& weight, color_t& color, multimap_<double, color_t>& split_list) {
    split_list.emplace(weight, color);    // insert it at the correct position ordered by weight
    if (split_list.size() > t) {
        split_list.erase(--split_list.end());    // if the top list exceeds its limit, erase the last entry
    }
}

/*
*
* [Filtering]
*
*/


/**
 * This function filters a greedy maximum weight tree compatible subset and returns a newick string.
 *
 * (@param map function that maps an integer to the original id, or null if no newick output wanted)
 * @param split_list list of splits to be filtered
 * (@param support_values a hash map storing the absolut support values for each color set)
 * (@param bootstrap_no the number of bootstrap replicates for computing the per centage support)
 * @param verbose print progress
 */

void graph::filter_strict(multimap_<double, color_t>& split_list, bool& verbose) {
    filter_strict(nullptr, split_list, nullptr, 0, verbose);
}

string graph::filter_strict(std::function<string(const uint16_t&)> map, multimap_<double, color_t>& split_list, hash_map<color_t, uint32_t>* support_values, const uint32_t& bootstrap_no, bool& verbose) {
    auto tree = vector<color_t>();    // create a set for compatible splits
    color_t col;
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
        col = it->second;
        if (test_strict(col, tree)) {
            tree.emplace_back(it->second);
            ++it; goto loop;    // if compatible, add the new split to the set
        }
        it = split_list.erase(it);  // otherwise, remove split
    }
    if (map) {
        node* root = build_tree(tree);
        return print_tree(root, map, support_values, bootstrap_no) + ";\n";
    } else {
        return "";
    }
}

/**
 * This function filters a greedy maximum weight weakly compatible subset.
 *
 * @param split_list list of splits to be filtered
 * @param verbose print progress
 */
void graph::filter_weakly(multimap_<double, color_t>& split_list, bool& verbose) {
    auto network = vector<color_t>();    // create a set for compatible splits
    color_t col;
    auto it = split_list.begin();
    uint64_t cur = 0, prog = 0, next;
    uint64_t max = split_list.size();

loop:
    while (it != split_list.end()) {
        if (verbose) {
            next = 100 * (cur * sqrt(cur)) / (max * sqrt(max));
             if (prog < next)  cout << "\33[2K\r" << "Filtering splits... " << next << "%" << flush;
            prog = next; cur++;
        }
        col = it -> second;
        if (test_weakly(col, network)) {
            network.emplace_back(it->second);
            ++it; goto loop;    // if compatible, add the new split to the set
        }
        it = split_list.erase(it);    // otherwise, remove split
    }
}

/**
 * This function filters a planar graph compatible subset.
 *
 * @param split_list list of splits to be filtered
 * @param verbose print progress
 */
void graph::filter_planar(multimap_<double, color_t>& split_list, bool& verbose, uint64_t& num) {

    auto it = split_list.begin();
    color_t split_copy;

    // single bit of split
    int state = 0;

    // construct PC-Tree with maxN nodes
    vector<pc_tree::PCNode*> leaves; 
	pc_tree::PCTree tree(num, &leaves);

   
loop:
    while(it != split_list.end()) {

        vector<pc_tree::PCNode*> consecutiveLeaves = {}; 
    
        // current split
        split_copy = it->second;

        // read all bits of the current split
        for(int j=0;j<maxN;j++) {
            state = (split_copy & 0b1u);
            split_copy >>= 01u;

            //collect all leaves which should be consecutive
            if (state == 1) {
                consecutiveLeaves.push_back(leaves.at(j));
            }
        }

        // if possible, insert the new split in pc-tree
        if (tree.makeConsecutive(consecutiveLeaves)) { 
            ++it; goto loop;
        }

        // delete the split, if it is not compatible to the pc-tree
        it = split_list.erase(it); 
                

    }
    
}


/**
 * This function filters a greedy maximum weight n-tree compatible subset and returns a string with all trees in newick format.
 *
 * @param n number of trees
 * (@param map function that maps an integer to the original id, or null)
 * @param split_list list of splits to be filtered
 * (@param support_values a hash map storing the absolut support values for each color set)
 * (@param bootstrap_no the number of bootstrap replicates for computing the per centage support)
 * @param verbose print progress
 */
void graph::filter_n_tree(uint64_t n, multimap_<double, color_t>& split_list, bool& verbose) {
    filter_n_tree(n, nullptr, split_list, nullptr, 0, verbose);
}

string graph::filter_n_tree(uint64_t n, std::function<string(const uint16_t&)> map, multimap_<double, color_t>& split_list, hash_map<color_t, uint32_t>* support_values, const uint32_t& bootstrap_no, bool& verbose) {
    auto forest = vector<vector<color_t>>(n);    // create a set for compatible splits
    color_t col;
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
        col = it-> second; 
        for (auto& tree : forest)
        if (test_strict(col, tree)) {
            tree.emplace_back(col);
            ++it; goto loop;    // if compatible, add the new split to the set
        }
        it = split_list.erase(it);    // otherwise, remove split
    }
    // output
    string s;
    if (map) {
        for (auto& tree : forest) {
            node* root = build_tree(tree);
            s += print_tree(root, map, support_values,bootstrap_no) + ";\n";
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
        if (!color::is_compatible(elem1, color)) {
            for (auto& elem2 : color_set) {
                if (!color::is_weakly_compatible(elem1, elem2, color)) {
                    return false;    // compare to each pair of splits in the set
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
    if (split.popcnt() < 2 || allTaxa.popcnt() - split.popcnt() < 2) { return true; }

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
			color_t inversesplit = split;
			color::complement(inversesplit);
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

    for (uint16_t i = 0; i < color::n; i++) {
        color_t leaf = 0b0u;
        leaf.set(i);
        allTaxa.set(i);
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
 * @param map function that maps an integer to the original id, or null
 * (@param support_values a hash map storing the absolut support values for each color set)
 * (@param bootstrap_no the number of bootstrap replicates for computing the per centage support)
 * @return newick string
 */
string graph::print_tree(node* root, std::function<string(const uint16_t&)> map) {
	return print_tree(root, map,nullptr,0);
}

string graph::print_tree(node* root, std::function<string(const uint16_t&)> map, hash_map<color_t, uint32_t>* support_values, const uint32_t& bootstrap_no) {
    vector<node*> subsets = root->subsets;
    color_t taxa = root->taxa;

    if (subsets.empty()){    // leaf set
        if (taxa.popcnt() == 0) {
            std::cerr << "ERROR: child with no taxon!?" << endl;
            exit(EXIT_FAILURE);
        } else if (taxa.popcnt() == 1) {
            return map(taxa.tzcnt()) + ":" + to_string(root->weight);
        } else {
            std::cerr << "ERROR: child with more than one taxon!?" << endl;
            exit(EXIT_FAILURE);
        }
    }
    else {
        string s = "(";
        for (node* subset : subsets) {
            s += print_tree(subset, map, support_values, bootstrap_no);
            if (subset != subsets.back()) { s += ","; }
        }
        s += ")";
		if(support_values!=nullptr){
			s+=to_string(((1.0*(*support_values)[taxa])/bootstrap_no));
		}
		s += ":";
        s += to_string(root->weight);
        return s;
    }
}


