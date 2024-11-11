#include "graph.h"
#include <algorithm>
#include <bits/fs_fwd.h>
#include <float.h>
#include <mutex>
#include <thread>

#include "util.h"

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


hash_map<color_t, double> graph::split_list_colors;

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
 * This function initializes the top list size, coverage threshold, and allowed characters.
 *
 * @param t top list size
 * @param q_table coverage thresholds
 * @param quality global q or maximum among all q values
 */

void graph::init(uint64_t& top_size, bool amino, vector<int>& q_table, int& quality, hash_set<kmer_t>& blacklist, hash_set<kmerAmino_t>& blacklist_amino, uint64_t& thread_count) {
    t = top_size;
    isAmino = amino;
    if(!isAmino){

        // Automatic table count
        //    table_count = 45 * thread_count - 33; // Estimated scaling
        //    table_count = table_count % 2 ? table_count : table_count + 1; // Ensure the table count is odd
        table_count = (0b1u << 14) + 1;
        

        // Init base tables
	    kmer_table = vector<hash_map<kmer_t, color_t>> (table_count);
	    singleton_kmer_table = vector<hash_map<kmer_t, uint16_t>> (table_count);

        // Init the lock vector
	    lock = vector<spinlock> (table_count);

        // Precompute the period for fast shift update kmer binning in bitset representation 
        #if (maxK > 32)     
        uint_fast32_t last = 1 % table_count;
        for (int i = 1; i <= 2*(kmer::k); i++)
        {
            // cout << last << endl;
	        period.push_back(last);
	        last = (2 * last) % table_count;
        }
        #endif

	    graph::allowedChars.push_back('A');
        graph::allowedChars.push_back('C');
        graph::allowedChars.push_back('G');
        graph::allowedChars.push_back('T');
    }else{
        // Automatic table count
        //    table_count = 33 * thread_count + 33; // Estimated scaling
        //    table_count = table_count % 2 ? table_count : table_count + 1; // Ensure the table count is odd
        table_count = (0b1u << 14) + 1;

        // Init amino tables
        kmer_tableAmino = vector<hash_map<kmerAmino_t, color_t>> (table_count);
        singleton_kmer_tableAmino = vector<hash_map<kmerAmino_t, uint16_t>> (table_count);
		
        // Init the mutex lock vector
        lock = vector<spinlock> (table_count);

        // Precompute the period for fast shift update kmer binning in bitset representation 
        #if (maxK > 12)     
        uint64_t last = 1 % table_count;
        for (int i = 1; i <= 5*(kmerAmino::k); i++)
        {
	        period.push_back(last);
	        last = (2 * last) % table_count;
        }
        #endif

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

    graph::quality = quality;
    graph::q_table = q_table;
	graph::blacklist = blacklist;
	graph::blacklist_amino = blacklist_amino;
    switch (quality) {
    case 1:
	case 0: /* no quality check */
        emplace_kmer_tmp = [&] (const uint64_t& T, uint_fast32_t& bin, const kmer_t& kmer, const uint16_t& color) {
            hash_kmer(bin, kmer, color);
        };
        emplace_kmer_amino_tmp = [&] (const uint64_t& T, uint_fast32_t& bin, const kmerAmino_t& kmer, const uint16_t& color) {
            hash_kmer_amino(bin, kmer, color);
        };
        break;

    case 2:
        isAmino ? quality_setAmino.resize(thread_count) : quality_set.resize(thread_count);
        if (q_table.size()>0){
            emplace_kmer_tmp = [&] (const uint64_t& T, uint_fast32_t& bin, const kmer_t& kmer, const uint16_t& color) {
                if (q_table[color]==1){
                    hash_kmer(bin, kmer, color);
                } else if (quality_set[T].find(kmer) == quality_set[T].end()) {
                    quality_set[T].emplace(kmer);
                } else {
                    quality_set[T].erase(kmer);
                    hash_kmer(bin, kmer, color);
                }
            };
            emplace_kmer_amino_tmp = [&] (const uint64_t& T, uint_fast32_t& bin, const kmerAmino_t& kmer, const uint16_t& color) {
                if (q_table[color]==1){
                    hash_kmer_amino(bin, kmer, color);
                } else if (quality_setAmino[T].find(kmer) == quality_setAmino[T].end()) {
                    quality_setAmino[T].emplace(kmer);
                } else {
                    quality_setAmino[T].erase(kmer);
                    hash_kmer_amino(bin, kmer, color);
                }
            };
        } else { // global quality value (one if-clause fewer)
            emplace_kmer_tmp = [&] (const uint64_t& T, uint_fast32_t& bin, const kmer_t& kmer, const uint16_t& color) {
                if (quality_set[T].find(kmer) == quality_set[T].end()) {
                    quality_set[T].emplace(kmer);
                } else {
                    quality_set[T].erase(kmer);
                    hash_kmer(bin, kmer, color);
                }
            };
            emplace_kmer_amino_tmp = [&] (const uint64_t& T, uint_fast32_t& bin, const kmerAmino_t& kmer, const uint16_t& color) {
                if (quality_setAmino[T].find(kmer) == quality_setAmino[T].end()) {
                    quality_setAmino[T].emplace(kmer);
                } else {
                    quality_setAmino[T].erase(kmer);
                    hash_kmer_amino(bin, kmer, color);
                }
            };
        }
        break;
    default:
        isAmino ? quality_mapAmino.resize(thread_count) : quality_map.resize(thread_count);
        if (q_table.size()>0){
            emplace_kmer_tmp = [&] (const uint64_t& T, uint_fast32_t& bin, const kmer_t& kmer, const uint16_t& color) {
                if (quality_map[T][kmer] < q_table[color]-1) {
                    quality_map[T][kmer]++;
                } else {
                    quality_map[T].erase(kmer);
                    hash_kmer(bin, kmer, color);
                }
            };
            emplace_kmer_amino_tmp = [&] (const uint64_t& T, uint_fast32_t& bin, const kmerAmino_t& kmer, const uint16_t& color) {
                if (quality_mapAmino[T][kmer] < q_table[color]-1) {
                    quality_mapAmino[T][kmer]++;
                } else {
                    quality_mapAmino[T].erase(kmer);
                    hash_kmer_amino(bin, kmer, color);
                }
            };
        }else { // global quality value
            emplace_kmer_tmp = [&] (const uint64_t& T, uint_fast32_t& bin, const kmer_t& kmer, const uint16_t& color) {
                if (quality_map[T][kmer] < quality-1) {
                    quality_map[T][kmer]++;
                } else {
                    quality_map[T].erase(kmer);
                    hash_kmer(bin, kmer, color);
                }
            };
            emplace_kmer_amino_tmp = [&] (const uint64_t& T, uint_fast32_t& bin, const kmerAmino_t& kmer, const uint16_t& color) {
                if (quality_mapAmino[T][kmer] < quality-1) {
                    quality_mapAmino[T][kmer]++;
                } else {
                    quality_mapAmino[T].erase(kmer);
                    hash_kmer_amino(bin, kmer, color);
                }
            };

        }
        break;
    }
	emplace_kmer = [&] (const uint64_t& T, uint_fast32_t& bin, const kmer_t& kmer, const uint16_t& color) {
		emplace_kmer_tmp(T, bin, kmer, color);
	};
	emplace_kmer_amino = [&] (const uint64_t& T, uint_fast32_t& bin, const kmerAmino_t& kmer, const uint16_t& color) {
		emplace_kmer_amino_tmp(T, bin, kmer, color);
	};
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
* This function hashes a k-mer and stores it in the correstponding hash table.
* The corresponding table is chosen by the carry of the encoded k-mer given the number of tables as module.
*  @param kmer The kmer to store
*  @param color The color to store 
*/
void graph::hash_kmer(uint_fast32_t& bin, const kmer_t& kmer, const uint16_t& color)
{
    lock[bin].lock();
	hash_map<kmer_t,color_t>::iterator entry=kmer_table[bin].find(kmer); 
	// already in the kmer table? -> add
	if(entry != kmer_table[bin].end()){
		entry.value().set(color);
	}
	// not yet in the kmer table?
	else{
		hash_map<kmer_t,uint16_t>::iterator s_entry = singleton_kmer_table[bin].find(kmer);
 		//seen once before? -> add to kmer table / remove from singleton table
		if(s_entry != singleton_kmer_table[bin].end()){
			if(s_entry.value() != color){
				kmer_table[bin][kmer].set(s_entry.value());
				kmer_table[bin][kmer].set(color);
				singleton_counters_locks[s_entry.value()].lock();
				singleton_counters[s_entry.value()]--;
				singleton_counters_locks[s_entry.value()].unlock();
				singleton_kmer_table[bin].erase(s_entry);
			}
		}
		// not seen before -> add to singleton_table
		else{
			singleton_kmer_table[bin][kmer]=color;
			singleton_counters_locks[color].lock();
			singleton_counters[color]++;
			singleton_counters_locks[color].unlock();
		}
	}
    lock[bin].unlock();
}


/**
 * This function hashes an amino k-mer and stores it in the corresponding hash table.
 * The correspontind table is chosen by the carry of the encoded k-mer bitset by the bit-module function.
 * @param kmer The kmer to store
 * @param color The color to store
 */
void graph::hash_kmer_amino(uint_fast32_t& bin, const kmerAmino_t& kmer, const uint16_t& color)
{
    lock[bin].lock();
	hash_map<kmerAmino_t,color_t>::iterator entry=kmer_tableAmino[bin].find(kmer); 
	// already in the kmer table? -> add
	if(entry != kmer_tableAmino[bin].end()){
		entry.value().set(color);
	}
	// not yet in the kmer table?
	else{
		hash_map<kmerAmino_t,uint16_t>::iterator s_entry = singleton_kmer_tableAmino[bin].find(kmer);
 		//seen once before? -> add to kmer table / remove from singleton table
		if(s_entry != singleton_kmer_tableAmino[bin].end()){
			if(s_entry.value() != color){
				kmer_tableAmino[bin][kmer].set(s_entry.value());
				kmer_tableAmino[bin][kmer].set(color);
				singleton_counters_locks[s_entry.value()].lock();
				singleton_counters[s_entry.value()]--;
				singleton_counters_locks[s_entry.value()].unlock();
				singleton_kmer_tableAmino[bin].erase(s_entry);
			}
		}
		// not seen before -> add to singleton_table
		else{
			singleton_kmer_tableAmino[bin][kmer]=color;
			singleton_counters_locks[color].lock();
			singleton_counters[color]++;
			singleton_counters_locks[color].unlock();
		}
	}	
    lock[bin].unlock();
}

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
void graph::add_cdbg_colored_kmer(double mean(uint32_t&, uint32_t&), string kmer_seq, color_t& kmer_color, double min_value){
    
    bool has_kmers = false; 
    for (auto table : kmer_table){if(!table.empty()){has_kmers = true; break;}} // Check if any entries exist

    if (has_kmers){ // check if the kmer is already stored
        kmer_t kmer; // create a kmer to search in the set of tables

        for (int pos=0; pos < kmer_seq.length(); ++pos) // collect the bases from the k-mer sequence.
        {
            kmer::shift(kmer, kmer_seq[pos]);
        }

	    bool reversed = kmer::reverse_represent(kmer);

        if (search_kmer(kmer)){ // Check if additional colors are stored for this kmer
            // Get the colors stored for this kmer
            color_t hashed_color = get_color(kmer, reversed); // the currently stored colores of the kmer
            for (uint64_t pos=0; pos < maxN; pos++){ // transcribe hashed colores to the cdbg color set
                    if(hashed_color.test(pos) && !kmer_color.test(pos)){ // test if the color is set in the stored color set
                        kmer_color.set(pos);
                    }
           }
           // Remove the kmer from the hash table
           remove_kmer(kmer, reversed); // remove the kmer from the table
	   }
    }            
    //Todo 
    bool pos = color::represent(kmer_color);  // invert the color set, if necessary
    if (kmer_color == 0) return; // ignore empty splits
    // min_value = add_weight(kmer_color, mean, min_value, pos); // compute weight
	array<uint32_t,2>& weight = color_table[kmer_color]; // get the weight and inverse weight for the color set
	weight[pos]++; // update the weight or the inverse weight of the current color set
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
 * This function calculates the weight for all splits and puts them into the split_list
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
 * This function calculates the weight for all splits based on the supplied color table and
 * puts them into the given split_map. All colors keys are saved as shared pointers to save
 * memory.
 * Note that this function does not allow for top list support !!
 *
* The split_map allows for efficient weight queries for given splits.
 * To get high scoring splits, use the ordered split_list set from compile_split_list instead!
 * @param mean weight function
 * @param color_table the custom color table to be used for the split calculation
 * @param mask the mask specifying all genomes that are processed in this split map
 * @param split_map the output parameter for the 'std::shared_ptr<color_t> -> double' split mapping
 * @param best_split the output parameter for the highest weighted split of all generated splits
 */
void graph::compile_split_map(double mean(uint32_t&, uint32_t&), hash_map<color_t, array<uint32_t, 2>>& color_table, const color_t& mask, hash_map<color_t, double>& split_map, pair<double, color_t>& best_split) {
    //best_split.first = -1; // make sure the best split is not set beforehand
    for (const auto& it : color_table) {
        // get the current color entry
        color_t color = it.first;
        array<uint32_t,2> occs = it.second;

        // check if the color needs to be flipped to be the represent (fewer 1 bits)
        // NOTE: This should never happen to begin with, but is included for the sake of defensive
        // programming, if this function is called with bad defined color_table inputs :)
        if (color::represent(color, mask)) {
            cerr << "WARNING: A representative color was actually not the representative while generating the split map. Self fixing error...";
            // make a new shared pointer to the color, as i cannot be inside of the color_table
            assert(color_table.find(color) == color_table.end());
            // this shared pointer has to be created before erasing the old pointer in order to
            // ensure that the reference count is at least 1. Otherwise, the actual value might be deleted
            //shared_ptr<color_t> new_shared_color = make_shared<color_t>(color);
            color_table.erase(color);
            color_table.emplace(color, occs);
            cerr << "done!" << endl;
        }

        // insert new split into the split list
        double new_mean = mean(occs[0], occs[1]);    // calculate the mean value
        // test if the complement or the color itself of the split is already in the list/map
        // NOTE: This can happen, because the gdac filter generates 2 new color tables in each
        // recursion which can include complementary color pairs if their mask has exactly twice
        // as many 1 bits as the colors themselves!
        auto map_change_it = split_map.find(color);
        if (map_change_it != split_map.end()) {
            // replace the old weight with the new mean (weight) if it is larger
            if (map_change_it->second < new_mean) {
                assert(map_change_it->first == color); // todo remove
                // update the weight inside the split map
                map_change_it.value() = new_mean;
                cerr << "ERROR: Collected already existing split in compile_split_map, its weight was smaller. Updated its weight with the higher value" << endl;
            } else if (map_change_it->second < new_mean) {
                cerr << "ERROR: Collected already existing split in compile_split_map, but its weight was equal. No updates needed." << endl;
            } else {
                cerr << "ERROR: Collected already existing split in compile_split_map, but its weight was greater. No updates needed." << endl;
            }

            // update the best split if the current weight is higher
            if (best_split.first < new_mean) {
                best_split.first = new_mean;
                best_split.second = map_change_it->first;
            }
            // skip as the color already exists. NOTE: THIS SHOULD NEVER BE THE CASE!!
            continue;
        } else { // check if complement split is already present
            color_t complement = color;
            color::complement(complement, mask);
            map_change_it = split_map.find(complement);
            if (map_change_it != split_map.end()) {
                if (map_change_it->second < new_mean) {
                    assert(map_change_it->first == color); // todo remove
                    // update the weight inside the split map
                    map_change_it.value() = new_mean;
                    cerr << "ERROR: Collected already existing complement split in compile_split_map, its weight was smaller. Updated its weight with the higher value" << endl;
                } else if (map_change_it->second < new_mean) {
                    cerr << "ERROR: Collected already existing complement split in compile_split_map, but its weight was equal. No updates needed." << endl;
                } else {
                    cerr << "ERROR: Collected already existing complement split in compile_split_map, but its weight was greater. No updates needed." << endl;
                }

                // update the best split if the current weight is higher
                if (best_split.first < new_mean) {
                    best_split.first = new_mean;
                    best_split.second = map_change_it->first;
                }
                // skip as the color already exists. NOTE: THIS SHOULD NEVER BE THE CASE!!
                continue;
            }
        }
        // split color (or its complement) is not in the list -> add the new split
        // CHECK IF COLOR IN TABLE OCCURS MULTIPLE TIMES
        auto pls_dont_exist = split_map.find(color);
        if (pls_dont_exist != split_map.end()) {
            cerr << "color was already present while iterating over color table in compile_split_list: " << endl;
        }
        // NOTE: the iterator is used again to ensure all shared pointers only point to ONE color (to safe memory)
        split_map.emplace(it.first, new_mean);

        // update the best split if the current weight is higher
        if (best_split.first < new_mean) {
            best_split.first = new_mean;
            best_split.second = it.first;
        }

        // remove the smallest weighted split if the top list limit was exceeded
        // NOT SUPPORTED WITH THE MAP VERSION FOR OBVIOUS REASONS...
    }
    if (color_table.size() > split_map.size()) {
        cerr << "ERROR: Color table size was greater than its split map after split compilation." << color_table.size() << " (color table) vs. " << split_map.size() << " (split map)" << endl;
    } else if (color_table.size() > split_map.size()) {
        cerr << "ERROR: Color table size was smaller than its split map after split compilation." << color_table.size() << " (color table) vs. " << split_map.size() << " (split map)" << endl;
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

// Entrypoint to GDAC implementation.

/**
 * This function filters the splitset using a greedy divide and conquer approach.
 * NOTE: This will currently force the filtered splits to be a compatible tree.
 *
 * (@param map function that maps an integer to the original id, or null if no newick output wanted)
 * @param split_list list of splits to be filtered
 * (@param support_values a hash map storing the absolut support values for each color set)
 * (@param bootstrap_no the number of bootstrap replicates for computing the per centage support)
 * @param verbose print progress
 */
void graph::filter_gdac(multimap_<double, color_t>& split_list, bool& verbose) {
    // TODO: actually make
    filter_gdac(nullptr, split_list, nullptr, 0, verbose);
}


string graph::filter_gdac(std::function<string(const uint16_t&)> map, multimap_<double, color_t>& split_list, hash_map<color_t, uint32_t>* support_values, const uint32_t& bootstrap_no, bool& verbose) {

    if (split_list.empty()) {
        // if the list is empty, there is nothing to do
        cerr << "ERROR: Tried to use the gdac filter on an initially empty split set." << endl;
        return ""; // TODO: replace with user-friendly error message to console
    }

    // convert the input split list set to use shared pointers in order to save memory for the colors
    // as most of the colors are not chaning through the recursion! Just their weights change ;)

    // compute the split map that maps from the shared colors to the entries of the split list set

    //shared_multimap_<double, color_t> shared_split_list;
    //shared_hash_map<color_t, double> split_shared_map;

    //multimap_shared<double, color_t> shared_split_list;

    // compute the color to weight split mapping that points with shared color pointers to safe memory on subsequent recursions
    hash_map<color_t, double> split_map;
    auto best_split = *split_list.begin();
    for (const auto& it : split_list) {
        //auto split_list_entry = make_pair(it.first, make_shared<color_t>(map_entry));
        //shared_split_list.emplace(split_list_entry);
        split_map.emplace(it.second, it.first);
    }

    // compute the color to weight split mapping that points to the split list shared pointers
    //hash_map<color_t, double> split_map;

    //shared_hash_map<color_t, double> split_map;
    //for (const auto& it : split_list) {
    //    split_map.emplace(it.second, it.first);
    //}

    // compute the input color mask
    color_t other_color = best_split.second;
    color::complement(other_color);
    color_t mask = mask | other_color;

    // apply the recursive filtering method
    auto filtered_splits_map = helper_filter_gdac(best_split, split_map, color_table, mask, verbose);

    // convert the hashmap back to an ordered set such that SANS can output it correctly :)
    split_list.clear();
    for (const auto& it : filtered_splits_map) {
        split_list.emplace(it.second, it.first);
    }

    //split_list = filtered_splits;
    // construct the newick string if necessary
    if (map) {
        // build the tree
        auto tree = vector<color_t>();
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
        node* root = build_tree(tree);
        return print_tree(root, map, support_values, bootstrap_no) + ";\n";
    }
    return "";


    // TODO DELETE
    //
    // // filtered result splits for the currently handled genomes
    // auto filtered_splits = vector<color_t>();
    // // BASE CASES
    // // TODO: Base Case genome set has 1 genome x .. -> filtered_splits = color corresponding to x | {}
    // // TODO: Base Case genome set has 2 genomes x,y .. -> filtered_splits = color corresponding to x | y
    //
    //
    // // split list contains splits of genome set T
    //
    // // get the highest weighted split of the current split set
    // auto whole_split_set_iter = split_list.begin();
    // // best split 'A | B'
    // color_t best_split_A_B = whole_split_set_iter->second;
    // // add the best split to the splits list
    // filtered_splits.push_back(best_split_A_B);
    //
    // // for easier reading purposes only. this is always equivalent to best_split_A_B
    // // TODO: maybe rename best_split_color to not waste unnecessary memory
    // color_t a_color = best_split_A_B;
    // // compute the complementary color representing set B of the 'A | B' split
    // color_t b_color = best_split_A_B;
    // color::complement(b_color);
    //
    // // collect all kmers corresponding to A and B respectively
    // auto a_kmers = vector<kmer_t>();
    // auto b_kmers = vector<kmer_t>();
    // for (const auto& kmer_map : kmer_table) {
    //     // iterate over the current bin map
    //     for (const auto& kmer_entry : kmer_map) {
    //         // sort the current kmer into the a list if the color matches
    //         if (kmer_entry.second == a_color) {
    //             a_kmers.emplace_back(kmer_entry.first);
    //         } else {
    //             // kmer does not correspond to the A set of the split, i.e. it must correspond to
    //             // the b set...
    //             b_kmers.emplace_back(kmer_entry.first);
    //         }
    //     }
    // }
    //
    // // TODO: construct new graphs from the kmer sets
    // // need q_table for graph::init? recounting kmer occurrence??
    // // main.cpp refrences q_table as the occurrence threshold per color
    // // -> kmer occ in color >= threshold -> kmer is in color? :(
    // // and more importantly, this would absolutely nuke the parallelization as the number of
    // // spinlocks would skyrocket...
    //
    // // call filter_gdac on these graphs recursively
    //
    //
    // cout << "Hallelujah" << "\n";
    //
    //
    //
    // // RECURSION
    // // somehow use compile_split_list with best_split_A (A)
    // multimap_<double, color_t> split_list_A {}; // TODO: DUMMY
    // // vice versa with split B (i.e. color::complement(best_split_A))
    // multimap_<double, color_t> split_list_B {}; // TODO: DUMMY
    //
    // // infer the A split (x^A, y^A) that maximizes the total weight of (x^A, B cup y^A) + (y^A, x^A cup B)
    // auto splits_iter_A = split_list_A.begin();
    // color_t cut_split_A {};
    // while (splits_iter_A != split_list_A.end()) {
    //     color_t x_A = splits_iter_A->second;
    //     color_t y_A = x_A; // deep copy, because color_t is primitive
    //     // TODO: how to combine two colors, i.e. splits that do not share the same genomes..
    //     // we need weight of color corresponding to split s_x_A = 'x^A | y^A cup B' and s_y_A = 'y^A | x^A cup B'
    //     color_t s_x_A {}; //dummy...
    //     color_t s_y_A {}; // dummy...
    //
    // }
    // // same inference for B split
    // auto splits_iter_B = split_list_B.begin();
    // color_t cut_split_B {};
    // while (splits_iter_B != split_list_A.end()) {
    //     color_t x_B = splits_iter_B->second;
    //     color_t y_B = x_B; // deep copy, because color_t is primitive
    //     // TODO: how to combine two colors, i.e. splits that do not share the same genomes..
    //     // we need weight of color corresponding to split s_x_B = 'x^B | y^B cup A' and s_y_B = 'y^B | x^B cup A'
    //     color_t s_x_B {}; //dummy...
    //     color_t s_y_B {}; // dummy...
    // }
    // // add the splits of A and B to the splits list
    // filtered_splits.push_back({}); // dummy...
    // // TODO: apply the filtering for the subset splits while calculating the best suited cut points for inclusion in here!
    //
    //
    // // TODO: add support for actually printing a newick string with bootstrap support
    // // auto tree = vector<color_t>();    // create a set for compatible splits
    // // color_t col;
    // // auto it = split_list.begin();
    // // uint64_t cur = 0, prog = 0, next;
    // // uint64_t max = split_list.size();
    // // loop:
    // //     while (it != split_list.end()) {
    // //         if (verbose) {
    // //             next = 100*cur/max;
    // //             if (prog < next)  cout << "\33[2K\r" << "Filtering splits... " << next << "%" << flush;
    // //             prog = next; cur++;
    // //         }
    // //         col = it->second;
    // //         if (test_strict(col, tree)) {
    // //             tree.emplace_back(it->second);
    // //             ++it; goto loop;    // if compatible, add the new split to the set
    // //         }
    // //         it = split_list.erase(it);  // otherwise, remove split
    // //     }
    // // if (map) {
    // //     node* root = build_tree(tree);
    // //     return print_tree(root, map, support_values, bootstrap_no) + ";\n";
    // // } else {
    // //     return "";
    // // }
}

hash_map<color_t, double> graph::helper_filter_gdac(pair<double, color_t>& best_split, hash_map<color_t, double>& split_map, const hash_map<color_t, array<uint32_t, 2>>& color_table, color_t& mask, bool verbose) {
    // split list contains splits of genome set T
    // get the highest weighted split of the current split set
    auto whole_split_set_iter = split_map.begin();

    // DEBUG BREAKPOINT ENTRY...
    if (split_map.size() == 637) {
        cerr << "reached breakpoint 1b2b3b" << endl;
    }

    // ==== SPLIT LIST SIZE BASE CASES ====
    if (split_map.size() <= 1) { // empty or lists with only one split cannot be filtered any further
        return split_map;
    } else if (split_map.size() == 2) {
        // check if the two splits of the list are each others complement, i.e. redundant
        //color_t best_color = whole_split_set_iter->first;
        color_t inverted_best = whole_split_set_iter->first;
        double best_weight = whole_split_set_iter->second;
        color::complement(inverted_best, mask);
        ++whole_split_set_iter;
        if (whole_split_set_iter->first == inverted_best) {
            // both splits can be reduced to the best split, as it holds the same information as both
            // i.e. 1. split = 'A | B' and 2. split 'B | A' => simplify to just 'A | B'
            if (best_weight < whole_split_set_iter->second) {
                // delete first split as its weight is smaller than the weight of its other identical color split
                split_map.erase(split_map.begin());
            } else {
                split_map.erase(whole_split_set_iter);
            }
            if (verbose) {
                cerr << "ERROR: simplified split set of size two, as both splits had complementary colors" << endl;
            }
            // no further refinement possible, as only one split is left (see first base case)
            return split_map;
        } else {
            // reset the split set iterator to point to its the start, i.e highest split
            whole_split_set_iter = split_map.begin();
        }
    }

    assert(whole_split_set_iter == split_map.begin()); // should always be the case, see above

    // ==== END SPLIT LIST SIZE BASE CASES ====

    // filtered result splits for the currently handled genomes
    auto filtered_splits = hash_map<color_t, double>();

    // ==== RETRIEVAL OF SUITABLE BEST SPLIT ====

    // test if highest split is empty or complete
    color_t best_color_A = best_split.second;
    double best_weight_AB = best_split.first;
    if (color::is_complete(best_color_A, mask) || best_color_A == 0b0u) {
        // NOTE: This is done just for the sake of defensive programming and should never be executed!

        // try next best split until it's not complete or empty (i.e. actual has additional info)
        ++whole_split_set_iter; // can never be the end as the list has at least 2 splits
        // (see split_list.size() <= 1 base case)
        best_weight_AB = -1; // ensure that the value is updated
        while (whole_split_set_iter != split_map.end() && (color::is_complete(best_color_A, mask) || best_color_A == 0b0u)) {
            if (whole_split_set_iter->second > best_weight_AB) {
                best_color_A = whole_split_set_iter->first;
                best_weight_AB = whole_split_set_iter->second;
            }
            ++whole_split_set_iter;
        }
        assert(best_weight_AB != -1);
        if (color::is_complete(best_color_A, mask) || best_color_A == 0b0u) {
            // the split set cannot be refined any further if no suitable best split is available
            if (verbose) {
                cerr << "ERROR: could not find a best split that is not empty/complete. This is not allowed to happen" << endl;
            }
            exit(EXIT_FAILURE);
            //return split_map;
        }
    }

    // ==== END RETRIEVAL OF SUITABLE BEST SPLIT IF ====

    // ==== DIVIDE STEP ====

    // best split 'A | B'
    //double best_weight_AB = whole_split_set_iter->first;
    //color_t best_color_A = whole_split_set_iter->second;
    // calculate the corresponding B part
    color_t best_color_B = best_color_A;
    color::complement(best_color_B, mask);
    // retrieve the used section of the bitmask
    color_t best_mask = best_color_B | best_color_A;

    size1N_t a_size = best_color_A.popcnt();
    size1N_t b_size = best_color_B.popcnt();
    size1N_t mask_size = mask.popcnt();

    // ==== CASE EXTINCTION FOR LEAVES (SINGLETON SPLITS) ====
    // test if the chosen best split connects two leaves (i.e. A and B both represent one taxon)
    if (a_size == 1 && mask_size == 2) { // BASE CASE
        // the chosen best split should be the only split in the refined split list
        if (verbose) {
            cout << "the size of A and B was 1, skipping further refinement as both are leaves." << endl;
        }
        filtered_splits.emplace(best_color_A,best_weight_AB);
        return filtered_splits;
    }

    // calculate new isolated split list for A and B from the given kmer color table

    auto color_table_size = color_table.size();
    hash_map<color_t, array<uint32_t,2>> color_table_A;
    hash_map<color_t, array<uint32_t,2>> color_table_B;

    if (a_size == 1) {
        // and A cup B (the mask) > 2... see leaf base case above
        // color table of a can be skipped, as the A color is a singleton, i.e. has no 'inner' splits...
        // but we still need the cut point in set B to connect A | B

        // collect new B splits from the A | B color table
        for (const auto& it : color_table) {
            // vice versa for the B part
            // erase all '1' bits that don't correspond to the B part
            color_t c_b_color = it.first & best_color_B;
            // check if new B color can be represented with fewer '1' bits
            array<uint32_t,2> c_b_occs = it.second;
            if (color::represent(c_b_color, best_color_B, b_size)) {
                // occurrence tuple has to be swapped if the color was inverted
                std::swap(c_b_occs[0], c_b_occs[1]);
            }

            if (c_b_color == 0b0u) {
                // skip the color if its empty
                continue;
            }
            // check if the representative color was already encountered
            auto it_b = color_table_B.find(c_b_color);
            if (it_b != color_table_B.end()) {
                // occurrences of the current color have to be added to the existing entry
                it_b.value()[0] += c_b_occs[0];
                it_b.value()[1] += c_b_occs[1];
            } else {
                // add the new color to the table
                // NOTE: the corresponding color has to be searched in the color table to get a new shared pointer.
                //auto color_it = color_table.find(c_b_color);
                //if (color_it != color_table.end()) {
                color_table_B.emplace(c_b_color, c_b_occs);
                //} else {
                // create shared pointer to the local variable if the color does not exist in the color_table
                // (can happen as the masking with one of the best split halves can introduce new colors)
                //color_table_B.emplace(make_shared<color_t>(c_b_color), c_b_occs);
                //}
            }
        }

    } else if (a_size > b_size) {
        // should never be reached, just for safety reasons / defensive programming.
        // assert would also be feasible here!
        if (verbose) {
            cerr << "ERROR: the number of 1 bits of the represent A of split (A, B) was greater than B, which is not allowed to happen" << endl;
        }
    } else {
        // both A and B do not represent leaves, i.e. we need to infer a cut point in both sets
        // to connect them with the chosen split 'A | B'!

        // collect the new isolated splits of sets A and B from the 'A cup B' color table
        for (const auto& it : color_table) {
            //if (a_size != 1) { // only process A if it is not a leaf
            // erase all '1' bits that don't correspond to the A part
            color_t c_a_color = it.first & best_color_A;
            // check if new A color can be represented with fewer '1' bits
            array<uint32_t,2> c_a_occs = it.second;
            if (color::represent(c_a_color, best_color_A, a_size)) {
                // occurrence tuple has to be swapped if the color was inverted
                std::swap(c_a_occs[0], c_a_occs[1]);
            }
            // ANOTHER TEST
            if (c_a_color != 0b0u) { // skip color if its empty
                // check if the representative color was already encountered
                auto it_a = color_table_A.find(c_a_color);
                if (it_a != color_table_A.end()) {
                    // occurrences of the current color have to be added to the existing entry
                    it_a.value()[0] += c_a_occs[0];
                    it_a.value()[1] += c_a_occs[1];
                } else {
                    // add the new color to the table
                    color_table_A.emplace(c_a_color, c_a_occs);
                }
            }
            // vice versa for the B part
            // erase all '1' bits that don't correspond to the B part
            color_t c_b_color = it.first & best_color_B;
            // check if new B color can be represented with fewer '1' bits
            array<uint32_t,2> c_b_occs = it.second;
            if (color::represent(c_b_color, best_color_B, b_size)) {
                // occurrence tuple has to be swapped if the color was inverted
                std::swap(c_b_occs[0], c_b_occs[1]);
            }

            if (c_b_color == 0b0u) {
                // skip the color if its empty
                continue;
            }
            // check if the representative color was already encountered
            auto it_b = color_table_B.find(c_b_color);
            if (it_b != color_table_B.end()) {
                // occurrences of the current color have to be added to the existing entry
                it_b.value()[0] += c_b_occs[0];
                it_b.value()[1] += c_b_occs[1];
            } else {
                // add the new color to the table
                color_table_B.emplace(c_b_color, c_b_occs);
            }
        }
    }

    if (a_size == 1) {
        // chosen best 'A | B' split is a trivial leaf (A) vs rest (B) split
        // -> only refine the B set further and skip all attempts to refine A...

        // todo refactor, such that everything from a fits in this...
    }


    /*for (const auto& it : color_table) {
        //if (a_size != 1) { // only process A if it is not a leaf
            // erase all '1' bits that don't correspond to the A part
            color_t c_a_color = it.first & best_color_A;
            // check if new A color can be represented with fewer '1' bits
            array<uint32_t,2> c_a_occs = it.second;
            if (color::represent(c_a_color, best_color_A, a_size)) {
                // occurrence tuple has to be swapped if the color was inverted
                std::swap(c_a_occs[0], c_a_occs[1]);
            }
            // ANOTHER TEST
            if (c_a_color != 0b0u) { // skip color if its empty
                // check if the representative color was already encountered
                auto it_a = color_table_A.find(c_a_color);
                if (it_a != color_table_A.end()) {
                    // occurrences of the current color have to be added to the existing entry
                    it_a.value()[0] += c_a_occs[0];
                    it_a.value()[1] += c_a_occs[1];
                } else {
                    // add the new color to the table
                    color_table_A.emplace(c_a_color, c_a_occs);
                }
            }
        //}
        // vice versa for the B part
        // erase all '1' bits that don't correspond to the B part
        color_t c_b_color = it.first & best_color_B;
        // check if new B color can be represented with fewer '1' bits
        array<uint32_t,2> c_b_occs = it.second;
        if (color::represent(c_b_color, best_color_B, b_size)) {
            // occurrence tuple has to be swapped if the color was inverted
            std::swap(c_b_occs[0], c_b_occs[1]);
        }
        // ANOTHER TEST
        if (c_b_color == 0b0u) {
            // skip the color
            continue;
        }
        // check if the representative color was already encountered
        auto it_b = color_table_B.find(c_b_color);
        if (it_b != color_table_B.end()) {
            // occurrences of the current color have to be added to the existing entry
            it_b.value()[0] += c_b_occs[0];
            it_b.value()[1] += c_b_occs[1];
        } else {
            // add the new color to the table
            color_table_B.emplace(c_b_color, c_b_occs);
        }
    }*/ // todo remove, as it was refactored above to be more efficient if one of the sets is a leaf
    // TODO: add verbose output for the split compilation if verbose param is true
    // calculate the new split lists from the corresponding color tables
    // auto a_split_list = multimap_<double, color_t>();
    auto a_split_map = hash_map<color_t, double>();
    pair<double, color_t> best_color_in_a_split_map;
    // auto b_split_list = multimap_<double, color_t>();
    auto b_split_map = hash_map<color_t, double>();
    pair<double, color_t> best_color_in_b_split_map;

    // compile the new split lists from the previously collected color tables
    if (a_size != 1) { // check if the a split list has to be calculated
        // a was not a leaf -> we need the new isolated A split list!
        graph::compile_split_map(util::geometric_mean2, color_table_A, best_color_A, a_split_map, best_color_in_a_split_map);
    }
    // NOTE: B can never be a leaf on its own, because then A would not have been the representative
    // color of the 'A | B' split
    // therefore the isolated B split set has to always be calculated, as the case of A and B being
    // leaves was already handled in the leaf base case!
    graph::compile_split_map(util::geometric_mean2, color_table_B, best_color_B, b_split_map, best_color_in_b_split_map);

    // auto a_split_map_size = a_split_map.size();
    // auto b_split_map_size = b_split_map.size();

    // auto same_size_a = 0; TODO REMOVE
    // auto same_size_b = 0;
    // for (const auto& it : split_map) {
    //     if (b_split_map.find(it.first) != b_split_map.end()) {
    //         ++same_size_b;
    //     }
    //     if (a_split_map.find(it.first) != a_split_map.end()) {
    //         ++same_size_a;
    //     }
    // }


    // ==== RECURSION ====
    // filter both split lists of the best split "A | B" recursively
    auto filtered_a_splits = helper_filter_gdac(best_color_in_a_split_map, a_split_map, color_table_A, best_color_A, verbose);
    auto filtered_b_splits = helper_filter_gdac(best_color_in_b_split_map, b_split_map, color_table_B, best_color_B, verbose);

    // ==== CONQUER / COMBINE STEP ====

    // ---- CUT INFERENCE ----
    // infer the most suited split in A to 'cut' in half / remove
    pair<double, color_t> a_x_vs_y_split;
    if (a_size != 1 && !filtered_a_splits.empty()) {
        a_x_vs_y_split =  graph::query_gdac_cut_split(mask, best_color_A, filtered_a_splits, split_map, verbose);
    } else {
        a_x_vs_y_split = {};
    }
    // vice versa for the best split in B to cut
    pair<double, color_t> b_x_vs_y_split;
    if (b_size != 1 && !filtered_b_splits.empty()) {
        b_x_vs_y_split = graph::query_gdac_cut_split(mask, best_color_B, filtered_b_splits, split_map, verbose);
    } else {
        b_x_vs_y_split = {};
    }

    // TODO remove if refactored method works
    // double max_weight_a {numeric_limits<double>::min()};
    // // iterate over all 'X | Y' splits inside the filtered splits of set A
    // for (const auto& it : filtered_a_splits) {
    //     // query weight of 'X | Y cup B' split
    //
    //     const color_t x_color {it.second};
    //     // check if X is the color representation with less '1' bits (i.e. it is the representative
    //     // key for the split)
    //     auto it_x {split_map.find(x_color)};
    //     double weight_x_vs_yb {};
    //     if (it_x != split_map.end()) {
    //         // X was the representative color, i.e. 'Y cup B' does not need to be calculated
    //         weight_x_vs_yb = it_x->second;
    //     } else {
    //         // 'Y cup B' must be the representative color with less '1' bits
    //         color_t y_cup_b_color {x_color};
    //         // NOTE: complement of X with respect to the mask of 'A cup B' must be 'Y cup B'
    //         // as 'X | Y cup B' is a valid split of the set 'A cup B'
    //         color::complement(y_cup_b_color, mask);
    //         auto it_yb {split_map.find(y_cup_b_color)};
    //         if (it_yb != split_map.end()) {
    //             // should always be true for valid split maps!
    //             weight_x_vs_yb = it_yb->second;
    //         } else {
    //             // only for the sake of defensive programming... should never be reached
    //             cerr << "ERROR: could not find ";
    //             for (size_t i = color::n; i-- > 0;) {
    //                 cerr << ((y_cup_b_color >> i) & 1);
    //             }
    //             cerr << " color (or its inverse) in the given splits!";
    //             return {};
    //         }
    //     }
    //     // query weight of 'Y | X cup B'
    //
    //     color_t y_color {x_color};
    //     color::complement(y_color, best_color_A);
    //     // check if Y is the color representation with less '1' bits (i.e. it is the representative
    //     // key for the split)
    //     auto it_y {split_map.find(y_color)};
    //     double weight_y_vs_xb {};
    //     if (it_y != split_map.end()) {
    //         // Y was the representative color, i.e. 'X cup B' does not need to be calculated
    //         weight_y_vs_xb = it_y->second;
    //     } else {
    //         // 'X cup B' must be the representative color with less '1' bits
    //         color_t x_cup_b_color {x_color};
    //         // NOTE: complement of Y with respect to the mask of 'A cup B' must be 'X cup B'
    //         // as 'Y | X cup B' is a valid split of the set 'A cup B'
    //         color::complement(x_cup_b_color, mask);
    //         auto it_xb {split_map.find(x_cup_b_color)};
    //         if (it_xb != split_map.end()) {
    //             // should always be true for valid split maps!
    //             weight_y_vs_xb = it_xb->second;
    //         } else {
    //             // only for the sake of defensive programming... should never be reached
    //             cerr << "ERROR: could not find ";
    //             for (size_t i = color::n; i-- > 0;) {
    //                 cerr << ((x_cup_b_color >> i) & 1);
    //             }
    //             cerr << " color (or its inverse) in the given splits!";
    //             return {};
    //         }
    //     }
    //     // update the new weight if it is higher than the old max weight
    //     // TODO: try different combination methods, e.g. geom mean, etc
    //     weight_x_vs_yb += weight_y_vs_xb;
    //     if (weight_x_vs_yb > max_weight_a) {
    //         max_weight_a = weight_x_vs_yb;
    //     }
    //
    // }

    // ---- COMBINE RESULTS TO A SPLIT LIST FOR SET 'A CUP B' ----

    // todo remove comment below
    // EDGE CASE: one of the recursions yields an empty split list, i.e. 0b0u as a color


    // edge case: both set A and B did not yield a suitable cut position, i.e. cannot be refined any further
    // -> return the input list from the previous recursion step
    if ((a_x_vs_y_split.second == 0b0u /*|| best_color_A == 0b0u*/) && (b_x_vs_y_split.second == 0b0u /*|| best_color_B == 0b0u*/)) {
        if (verbose) {
            cerr << "cannot refine the split data any further, as the cuts for the filtered A and B parts were empty. returning unfiltered list" << endl;
        }
        // return the input list without any additional filtering in this recursion step
        return split_map;
    }
    /*
    if (a_x_vs_y_split.second == 0b0u || best_color_A == 0b0u) {
        // no suitable cut was found for the A set, but a cut for set B was found
        if (a_x_vs_y_split.second != best_color_A) {
            cerr << "okay that can happen i guess :( A" << endl;
        }
        cerr << "first case hehe" << endl;
        // set A cannot be refined any further, but B can!
        filtered_splits.insert(b_x_vs_y_split);
        return filtered_splits;
    }
    if (b_x_vs_y_split.second == 0b0u || best_color_B == 0b0u) {
        if (a_x_vs_y_split.second != best_color_A) {
            cerr << "okay that can happen i guess :( B" << endl;
        }
        cerr << "second case hehe" << endl;
        filtered_splits.insert(a_x_vs_y_split);
        return filtered_splits;
    }
    */

    // add the best split 'A | B' to the splits list
    // Note: A is always the representative color by construction!
    // get the original value of best_color_A to make a new shared pointer (saving memory)
    //auto best_color_A_it = split_map.find(best_color_A);
    //assert (best_color_A_it != split_map.end());
    filtered_splits.emplace(best_color_A, best_weight_AB);

    // use refined A set, if a cut split was found in A (i.e. 'A^x | A^y' is not empty)
    if (a_x_vs_y_split.second != 0b0u && best_color_A != 0b0u) {

        // add split set of A without the chosen cut split
        auto a_cut_it {a_split_map.find(a_x_vs_y_split.second)}; // retrieve the original weight of the split
        if (a_cut_it == a_split_map.end()) {
            // if it could not be found, try its complement
            a_cut_it = a_split_map.find(~(a_x_vs_y_split.second) & best_color_A);
        }
        assert(a_cut_it != a_split_map.end()); // the key has to exist as it is generated from this map
        // use the original weight of the color to erase the pair from the split list
        auto found_a_remove {filtered_a_splits.find(a_cut_it->first)};
        if (found_a_remove != filtered_a_splits.end()) {
            assert (found_a_remove->first == a_cut_it->first); // todo remove hopefully
        }
        auto num_removed = filtered_a_splits.erase(a_cut_it->first);
        a_split_map.erase(a_cut_it);
        // TODO bring back assert(num_removed == 1);
        if (num_removed != 1 && verbose) {
            cerr << "did not remove the A cut split, as it was not part of the filtered A splits." << endl;
        }
        // add the updated split list of set A to the output
        // NOTE: Insert and move_iterator would be more efficient here, because the pair values already exist and can be cheaply moved, but it is not defined for the tsl hashmap :(
        for (auto& it : filtered_a_splits) {
            filtered_splits.emplace(it.first, it.second);
        }
        //filtered_splits.insert(std::make_move_iterator(filtered_a_splits.begin()), std::make_move_iterator(filtered_a_splits.end()));

        // Add the cut splits individually as 'X | other' splits of set 'A cup B' to the filtered output
        // NOTE: no representative check once found has to be done as these are directly queried from
        // the input split list

        // @ 'Ax | (A cup B) \ Ax'
        auto ax_in_AB_it {split_map.find(a_x_vs_y_split.second)};
        if (ax_in_AB_it != split_map.end()) {
            // Ax was the representative of 'Ax | (A cup B) \ Ax'
            filtered_splits.emplace(ax_in_AB_it->first, ax_in_AB_it->second);
        } else {
            // (A cup B) \ Ax must be the representative of 'Ax | (A cup B) \ Ax'
            // construct via complement
            color_t other_without_ax_color = a_x_vs_y_split.second;
            color::complement(other_without_ax_color, mask);
            auto other_without_ax_it {split_map.find(other_without_ax_color)};
            // has to be a key, as otherwise the split would not be valid
            assert(other_without_ax_it != split_map.end());
            filtered_splits.emplace(other_without_ax_it->first, other_without_ax_it->second);
            //filtered_splits.emplace(other_without_ax_it->second, other_without_ax_it->first);
        }
        // @ 'Ay | (A cup B) \ Ay'
        color_t a_y_vs_x_color {a_x_vs_y_split.second};
        color::complement(a_y_vs_x_color, best_color_A);
        auto ay_in_AB_it {split_map.find(a_y_vs_x_color)};
        if (ay_in_AB_it != split_map.end()) {
            // Ay was the representative of 'Ay | (A cup B) \ Ay'
            filtered_splits.emplace(ay_in_AB_it->first, ay_in_AB_it->second);
            //filtered_splits.emplace(ay_in_AB_it->second, ay_in_AB_it->first);
        } else {
            // (A cup B) \ Ay must be the representative of 'Ay | (A cup B) \ Ay'
            // construct via complement
            color_t other_without_ay_color = a_y_vs_x_color;
            color::complement(other_without_ay_color, mask);
            auto other_without_ay_it {split_map.find(other_without_ay_color)};
            // has to be a key, as otherwise the split would not be valid
            assert(other_without_ay_it != split_map.end());
            filtered_splits.emplace(other_without_ay_it->first, other_without_ay_it->second);
            //filtered_splits.emplace(other_without_ay_it->second, other_without_ay_it->first);
        }
    }

    if (b_x_vs_y_split.second != 0b0u && best_color_B != 0b0u) {
        // vice versa for B split list
        auto b_cut_it {b_split_map.find(b_x_vs_y_split.second)}; // retrieve the original weight of the split
        if (b_cut_it == b_split_map.end()) {
            // if it could not be found, try its complement
            b_cut_it = b_split_map.find(~b_x_vs_y_split.second & best_color_B);
        }
        assert(b_cut_it != b_split_map.end()); // the key has to exist as it is generated from this map
        // use the original weight of the color to erase the pair from the split list
        //auto num_removed = filtered_b_splits.erase(make_pair(b_cut_it->second, b_cut_it->first));
        auto num_removed = filtered_b_splits.erase(b_cut_it->first);
        if (num_removed != 1 && verbose) {
            cerr << "did not remove the B cut split, as it was not part of the filtered B splits." << endl;
        }
        // assert(num_removed == 1);
        // add the updated split list of set A to the output
        // NOTE: Insert and move_iterator is more efficient here, because the pair values already exist and can be cheaply moved

         for (auto& it : filtered_b_splits) {
             filtered_splits.emplace(it.first, it.second);
        }

        // NOT WORKING WITH TSL:SPARSE_MAPS
        //filtered_splits.emplace(std::make_move_iterator(filtered_b_splits.begin()), std::make_move_iterator(filtered_b_splits.end()));

        // @ 'Bx | (A cup B) \ Bx'
        auto bx_in_AB_it {split_map.find(b_x_vs_y_split.second)};
        if (bx_in_AB_it != split_map.end()) {
            // Bx was the representative of 'Bx | (A cup B) \ Bx'
            filtered_splits.emplace(bx_in_AB_it->first, bx_in_AB_it->second);
            //filtered_splits.emplace(bx_in_AB_it->second, bx_in_AB_it->first);
        } else {
            // (A cup B) \ Bx must be the representative of 'Bx | (A cup B) \ Bx'
            // construct via complement
            color_t other_without_bx_color = b_x_vs_y_split.second;
            color::complement(other_without_bx_color, mask);
            auto other_without_bx_it {split_map.find(other_without_bx_color)};
            // has to be a key, as otherwise the split would not be valid
            assert(other_without_bx_it != split_map.end());
            filtered_splits.emplace(other_without_bx_it->first, other_without_bx_it->second);
            //filtered_splits.emplace(other_without_bx_it->second, other_without_bx_it->first);
        }
        // @ 'By | (A cup B) \ By'
        color_t b_y_vs_x_color {b_x_vs_y_split.second};
        color::complement(b_y_vs_x_color, best_color_B);
        auto by_in_AB_it {split_map.find(b_y_vs_x_color)};
        if (by_in_AB_it != split_map.end()) {
            // Ay was the representative of 'Ay | (A cup B) \ Ay'
            filtered_splits.emplace(by_in_AB_it->first, by_in_AB_it->second);
            //filtered_splits.emplace(by_in_AB_it->second, by_in_AB_it->first);
        } else {
            // (A cup B) \ By must be the representative of 'By | (A cup B) \ By'
            // construct via complement
            color_t other_without_by_color = b_y_vs_x_color;
            color::complement(other_without_by_color, mask);
            auto other_without_by_it {split_map.find(other_without_by_color)};
            // has to be a key, as otherwise the split would not be valid
            assert(other_without_by_it != split_map.end());
            filtered_splits.emplace(other_without_by_it->first, other_without_by_it->second);
            //filtered_splits.emplace(other_without_by_it->second, other_without_by_it->first);
        }
    }
    // return the finished (i.e. refined) split set for set 'A cup B'
    return filtered_splits;
}

/**
 * This function calculated the highest weighted cut point split 'X | Y' of set A where the
 * weight is calculated as the following formula:
 *
 * weight(X | Y cup B) + weight(Y | X cup B)
 *
 * given weights of splits for the set 'A cup B'
 * @param whole_mask the bitmask for the splits of A cup B
 * @param a_mask the bitmask for the splits of A
 * @param split_map_A the *gdac filtered* split_map of set A
 * @param whole_split_map the split map (color -> weight of color) of set A cup B
 * @param verbose will print more error information, if true.
 * @return the best suited cut point split in the supplied split_list for usage with
 * graph::helper_filter_gdac
 */
pair<double, color_t> graph::query_gdac_cut_split( color_t& whole_mask,  color_t& a_mask,  hash_map<color_t, double>& split_map_A,  hash_map<color_t, double>& whole_split_map,  bool& verbose){
    // get the number of '1' bits of the a_mask for usage with color::represent
    auto num_a_colors {a_mask.popcnt()};
    // infer the most suited split in A to 'cut' in half / remove
    double max_weight_a {numeric_limits<double>::min()};
    color_t best_X_vs_Y_color {};
    // iterate over all 'X | Y' splits inside the filtered splits of set A
    for (const auto& it : split_map_A) {
        // query weight of 'X | Y cup B' split
        const color_t x_color {it.first};
        // check if X is the color representation with less '1' bits (i.e. it is the representative
        // key for the split)
        auto it_x {whole_split_map.find(x_color)};
        double weight_x_vs_yb {};
        if (it_x != whole_split_map.end()) {
            // X was the representative color, i.e. 'Y cup B' does not need to be calculated
            weight_x_vs_yb = it_x->second;
        } else {
            // 'Y cup B' must be the representative color with less '1' bits
            color_t y_cup_b_color {x_color};
            // NOTE: complement of X with respect to the mask of 'A cup B' must be 'Y cup B'
            // as 'X | Y cup B' is a valid split of the set 'A cup B'
            color::complement(y_cup_b_color, whole_mask);
            auto it_yb {whole_split_map.find(y_cup_b_color)};
            if (it_yb != whole_split_map.end()) {
                // should always be true for valid split maps!
                weight_x_vs_yb = it_yb->second;
            } else {
                if (verbose) {
                    // only for the sake of defensive programming... should never be reached
                    /*cerr << "ERROR: could not find a weight for X | Y cup B";
                    for (size_t i = color::n; i-- > 0;) {
                        cerr << ((y_cup_b_color >> i) & 1);
                    }
                    cerr << " color in the previous splits! Skipping this color" << endl ;*/
                }
                continue;
                //return {};
            }
        }
        // query weight of 'Y | X cup B'

        color_t y_color {x_color};
        color::complement(y_color, a_mask);
        // check if Y is the color representation with less '1' bits (i.e. it is the representative
        // key for the split)
        auto it_y {whole_split_map.find(y_color)};
        double weight_y_vs_xb {};
        if (it_y != whole_split_map.end()) {
            // Y was the representative color, i.e. 'X cup B' does not need to be calculated
            weight_y_vs_xb = it_y->second;
        } else {
            // 'X cup B' must be the representative color with less '1' bits
            color_t x_cup_b_color {y_color};
            // NOTE: complement of Y with respect to the mask of 'A cup B' must be 'X cup B'
            // as 'Y | X cup B' is a valid split of the set 'A cup B'
            color::complement(x_cup_b_color, whole_mask);
            auto it_xb {whole_split_map.find(x_cup_b_color)};
            if (it_xb != whole_split_map.end()) {
                // should always be true for valid split maps!
                weight_y_vs_xb = it_xb->second;
            } else {
                // only for the sake of defensive programming... should never be reached
                if (verbose) {
                    /*cerr << "ERROR: could not find a weight for Y | X cup B";
                    for (size_t i = color::n; i-- > 0;) {
                        cerr << ((x_cup_b_color >> i) & 1);
                    }
                    cerr << " color in the given splits! Skipping this color" << endl ;*/
                }
                continue;
            }
        }
        // update the new weight if it is higher than the old max weight
        //weight_x_vs_yb += weight_y_vs_xb;
        // alternative weight calculation using the geometric mean with pseudocounts.
        // This is not used, as it did not yield different results but is ever so slightly less efficient..
        weight_x_vs_yb = weight_x_vs_yb + weight_y_vs_xb;
        if (weight_x_vs_yb > max_weight_a) {
            max_weight_a = weight_x_vs_yb;
            best_X_vs_Y_color = x_color;
            // ensure that the color with the least '1' bits is stored
            color::represent(best_X_vs_Y_color, a_mask, num_a_colors);
        }
    }
    // retrieve the shared pointer to make a new one pointing to the same color :)
    //auto best_color_it = split_map_A.find(best_X_vs_Y_color);
    //assert (best_color_it != split_map_A.end());

    return make_pair(max_weight_a, best_X_vs_Y_color);
}

// multimap_<double, color_t> graph::helper_filter_gdac1(multimap_<double, color_t>& split_list, bool& verbose) {
//     // filtered result splits for the currently handled genomes
//     auto filtered_splits = multimap_<double, color_t>();
//
//     // split list contains splits of genome set T
//     // get the highest weighted split of the current split set
//     auto whole_split_set_iter = *split_list.begin();
//
//     // BASE CASES
//     if (split_list.empty()) {
//         return filtered_splits; // return empty split_list
//     }
//     if (split_list.size() == 1) {
//         // there is only one split, i.e. is automatically the best available pick
//         filtered_splits.emplace(whole_split_set_iter.first, whole_split_set_iter.second);
//         return filtered_splits;
//     }
//
//     // best split 'A | B'
//     color_t best_split_A_B = whole_split_set_iter.second;
//     // add the best split to the splits list
//     add_split(whole_split_set_iter.first, whole_split_set_iter.second, filtered_splits);
//
//     // for easier reading purposes only. this is always equivalent to best_split_A_B
//     // TODO: maybe rename best_split_color to not waste unnecessary memory
//     color_t a_color = best_split_A_B;
//     // compute the complementary color representing set B of the 'A | B' split
//     color_t b_color = best_split_A_B;
//     color::complement(b_color);
//
//     // collect all kmers corresponding to A and B respectively
//     auto a_kmers = vector<kmer_t>();
//     auto b_kmers = vector<kmer_t>();
//     for (const auto& kmer_map : kmer_table) {
//         // iterate over the current bin map
//         for (const auto& kmer_entry : kmer_map) {
//             // sort the current kmer into the a list if the color matches
//             if (kmer_entry.second == a_color) {
//                 a_kmers.emplace_back(kmer_entry.first);
//             } else {
//                 // kmer does not correspond to the A set of the split, i.e. it must correspond to
//                 // the b set...
//                 b_kmers.emplace_back(kmer_entry.first);
//             }
//         }
//     }
//
//     // RECURSION
//
//     // TODO: construct new graphs from the kmer sets
//     // need q_table for graph::init? recounting kmer occurrence??
//     // main.cpp refrences q_table as the occurrence threshold per color
//     // -> kmer occ in color >= threshold -> kmer is in color? :(
//     // and more importantly, this would absolutely nuke the parallelization as the number of
//     // spinlocks would skyrocket...
//     graph a_graph {}; // dummy stub
//     graph b_graph {}; // dummy stub
//
//     // call helper_filter_gdac on these graphs recursively
//     auto a_split_set = helper_filter_gdac(a_graph.split_list, verbose);
//     hash_map<color_t, double> b_split_colors {};
//     auto b_split_set = helper_filter_gdac(b_graph.split_list, verbose);
//
//     // infer the most suited split in A to cut
//     double max_weight_a {DBL_MIN};
//     double original_weight_a {DBL_MIN};
//     color_t best_split_A {};
//     bool initialized = false;
//     for (auto& a_split : a_split_set) {
//         double a_x_weight {0};
//         double a_y_weight {0};
//
//         color_t a_x_color {a_split.second};
//         color_t a_y_color {a_x_color};
//         color::complement(a_y_color);
//         color_t a_y_B_color {a_y_color & b_color};
//         color_t a_x_B_color {a_x_color & b_color};
//         // find 'A^x | A^y cup B' and 'A^y | A^x cup B' splits in the split_list
//         // where A^x and A^y represent the current a_split
//         // NOTE: Because each color bitmask is forced to be of the same size, A^x | A^y is
//         // equivalent to A^x | A^y cup B ?!
//         auto a_x_iter = split_list_colors.find(a_x_color);
//         if (a_x_iter != split_list_colors.end()) {
//             // a_x_color (A^x) has less 1's than 'A^y cup B' (inverse of A^x)
//             a_x_weight = a_x_iter->second;
//         } else {
//             // A^y cup B has less 1's in its mask, i.e. is used as a key instead
//             auto a_y_B_iter = split_list_colors.find(a_y_B_color);
//             if (a_y_B_iter != split_list_colors.end()) {
//                 a_x_weight = a_y_B_iter->second;
//             } else {
//                 // color was not found
//                 cout << "ERROR: could not find ";
//                 for (size_t i = color::n; i-- > 0;) {
//                     cout << ((a_x_color >> i) & 1);
//                 }
//                 cout << " color (or its inverse) in split list!";
//                 return {};
//             }
//         }
//
//         // vice versa for A^y | A^x cup B
//         auto a_y_iter = split_list_colors.find(a_y_color);
//         if (a_y_iter != split_list_colors.end()) {
//             // a_x_color (A^x) has less 1's than 'A^y cup B' (inverse of A^x)
//             a_y_weight = a_y_iter->second;
//         } else {
//             // A^x cup B has less 1's in its mask, i.e. is used as a key instead
//             auto a_x_B_iter = split_list_colors.find(a_x_B_color);
//             if (a_x_B_iter != split_list_colors.end()) {
//                 a_y_weight = a_x_B_iter->second;
//             } else {
//                 // color was not found
//                 cout << "ERROR: could not find ";
//                 for (size_t i = color::n; i-- > 0;) {
//                     cout << ((a_x_B_color >> i) & 1);
//                 }
//                 cout << " color (or its inverse) in split list!";
//                 return {};
//             }
//         }
//
//         if (!initialized) {
//             best_split_A =  a_split.second;
//             max_weight_a = a_y_weight + a_x_weight;
//             original_weight_a = a_split.first;
//             initialized = true;
//         // TODO: calculation of split weight is subject to change
//         // e.g. try geometric mean, etc.
//         } else if (a_y_weight + a_x_weight > max_weight_a) {
//             // found new candidate for best a_split to cut!
//             max_weight_a = a_y_weight + a_x_weight;
//             color::represent(a_x_color);
//             best_split_A = a_x_color;
//             original_weight_a = a_split.first;
//         }
//     }
//     // vice verse for optimal B cut...
//     color_t best_split_B {};
//     double max_weight_b {DBL_MIN};
//     double original_weight_b {DBL_MIN};
//     // reset counter
//     initialized = false;
//     for (auto& b_split : b_split_set) {
//         double b_x_weight {0};
//         double b_y_weight {0};
//
//         color_t b_x_color {b_split.second};
//         color_t b_y_color {b_x_color};
//         color::complement(b_y_color);
//         color_t b_y_A_color {b_y_color & a_color};
//         color_t b_x_A_color {b_x_color & a_color};
//
//         // find 'B^x | B^y cup A' and 'B^y | B^x cup A' splits in the split_list
//         // where B^x and B^y represent the current b_split
//         // NOTE: Because each color bitmask is forced to be of the same size, B^x | B^y is
//         // equivalent to B^x | B^y cup A ?!
//         auto b_x_iter = split_list_colors.find(b_x_color);
//         if (b_x_iter != split_list_colors.end()) {
//             // a_x_color (A^x) has less 1's than 'A^y cup B' (inverse of A^x)
//             b_x_weight = b_x_iter->second;
//         } else {
//             // A^y cup B has less 1's in its mask, i.e. is used as a key instead
//             auto b_y_A_iter = split_list_colors.find(b_y_A_color);
//             if (b_y_A_iter != split_list_colors.end()) {
//                 b_x_weight = b_y_A_iter->second;
//             } else {
//                 // color was not found
//                 cout << "ERROR: could not find ";
//                 for (size_t i = color::n; i-- > 0;) {
//                     cout << ((b_x_color >> i) & 1);
//                 }
//                 cout << " color (or its inverse) in split list!";
//                 return {};
//             }
//         }
//
//         // vice versa for A^y | A^x cup B
//         auto b_y_iter = split_list_colors.find(b_y_color);
//         if (b_y_iter != split_list_colors.end()) {
//             // a_x_color (A^x) has less 1's than 'A^y cup B' (inverse of A^x)
//             b_y_weight = b_y_iter->second;
//         } else {
//             // B^y cup B has less 1's in its mask, i.e. is used as a key instead
//             auto b_x_A_iter = split_list_colors.find(b_x_A_color);
//             if (b_x_A_iter != split_list_colors.end()) {
//                 b_y_weight = b_x_A_iter->second;
//             } else {
//                 // color was not found
//                 cout << "ERROR: could not find ";
//                 for (size_t i = color::n; i-- > 0;) {
//                     cout << ((b_x_A_color >> i) & 1);
//                 }
//                 cout << " color (or its inverse) in split list!";
//                 return {};
//             }
//         }
//
//         if (!initialized) {
//             best_split_B =  b_split.second;
//             max_weight_b = b_y_weight + b_x_weight;
//             original_weight_b = b_split.first;
//             initialized = true;
//         // TODO: calculation of split weight is subject to change
//         // e.g. try geometric mean, etc.
//         } else if (b_y_weight + b_x_weight > max_weight_b) {
//             // found new candidate for best a_split to cut!
//             max_weight_b = b_y_weight + b_x_weight;
//             color::represent(b_x_color);
//             best_split_B = b_x_color;
//             original_weight_b = b_split.first;
//         }
//     }
//
//     // TODO: Check if insert is ok for ordered sets...
//
//     // add split_listA without split A^x, A^y to the split list
//     a_split_set.erase(std::make_pair(original_weight_a, best_split_A));
//     filtered_splits.insert(a_split_set.begin(), a_split_set.end());
//     // add split_listB without split B^x, B^y to the split list
//     b_split_set.erase(std::make_pair(original_weight_b, best_split_B));
//     filtered_splits.insert(b_split_set.begin(), b_split_set.end());
//
//     // add all the cut splits, i.e. (A^x | rest) (A^y | rest) (B^x | rest) (B^y | rest)...
//     auto it = split_list_colors.find(best_split_A);
//     if (it != split_list_colors.end()) {
//         // best_split_A was representative of split -> add it to the filtered list output
//         filtered_splits.emplace(it->second, it->first);
//     } else {
//         // try to find rest by removing its bits from a & b
//         color_t best_split_A_inv = ~best_split_A;
//         color_t all_without_A = (a_color & b_color) & best_split_A_inv;
//         auto it_all = split_list_colors.find(all_without_A);
//         if (it_all != split_list_colors.end()) {
//             // found corresponding color
//             filtered_splits.emplace(it_all->second, it_all->first);
//         } else {
//             // color was not found
//             // error, should have found this split !!
//             cout << "ERROR: could not find ";
//             for (size_t i = color::n; i-- > 0;) {
//                 cout << ((all_without_A >> i) & 1);
//             }
//             cout << " color (or its inverse) in split list!";
//             return {};
//         }
//
//     }
//     color_t a_inverse = best_split_A;
//     color::complement(a_inverse);
//     it = split_list_colors.find(a_inverse);
//     if (it != split_list_colors.end()) {
//         // best_split_A was representative of split -> add it to the filtered list output
//         filtered_splits.emplace(it->second, it->first);
//     } else {
//         // try to find rest by removing its bits from a & b
//         color_t best_split_A_inv = ~a_inverse;
//         color_t all_without_A = (a_color & b_color) & best_split_A_inv;
//         auto it_all = split_list_colors.find(all_without_A);
//         if (it_all != split_list_colors.end()) {
//             // found corresponding color
//             filtered_splits.emplace(it_all->second, it_all->first);
//         } else {
//             // color was not found
//             // error, should have found this split !!
//             cout << "ERROR: could not find ";
//             for (size_t i = color::n; i-- > 0;) {
//                 cout << ((all_without_A >> i) & 1);
//             }
//             cout << " color (or its inverse) in split list!";
//             return {};
//         }
//     }
//
//     // vice versa for B^x and B^y
//
//     it = split_list_colors.find(best_split_B);
//     if (it != split_list_colors.end()) {
//         // best_split_B was representative of split -> add it to the filtered list output
//         filtered_splits.emplace(it->second, it->first);
//     } else {
//         // try to find rest by removing its bits from a & b
//         color_t best_split_B_inv = ~best_split_B;
//         color_t all_without_B = (a_color & b_color) & best_split_B_inv;
//         auto it_all = split_list_colors.find(all_without_B);
//         if (it_all != split_list_colors.end()) {
//             // found corresponding color
//             filtered_splits.emplace(it_all->second, it_all->first);
//         } else {
//             // color was not found
//             // error, should have found this split !!
//             cout << "ERROR: could not find ";
//             for (size_t i = color::n; i-- > 0;) {
//                 cout << ((all_without_B >> i) & 1);
//             }
//             cout << " color (or its inverse) in split list!";
//             return {};
//         }
//
//     }
//     color_t b_inverse = best_split_B;
//     color::complement(b_inverse);
//     it = split_list_colors.find(b_inverse);
//     if (it != split_list_colors.end()) {
//         // best_split_A was representative of split -> add it to the filtered list output
//         filtered_splits.emplace(it->second, it->first);
//     } else {
//         // try to find rest by removing its bits from a & b
//         color_t best_split_B_inv = ~b_inverse;
//         color_t all_without_B = (a_color & b_color) & best_split_B_inv;
//         auto it_all = split_list_colors.find(all_without_B);
//         if (it_all != split_list_colors.end()) {
//             // found corresponding color
//             filtered_splits.emplace(it_all->second, it_all->first);
//         } else {
//             // color was not found
//             // error, should have found this split !!
//             cout << "ERROR: could not find ";
//             for (size_t i = color::n; i-- > 0;) {
//                 cout << ((all_without_B >> i) & 1);
//             }
//             cout << " color (or its inverse) in split list!";
//             return {};
//         }
//     }
//     return filtered_splits;
// }

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


