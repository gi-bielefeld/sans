#include <iostream>
#include <limits>
#include <algorithm>
#include <functional>
#include <utility>
#include <vector>

#include <thread>
#include <atomic>

#include <iomanip>
#include <string>
#include <random>



#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include "tsl/sparse_map.h"
#include "tsl/sparse_set.h"

using namespace std;

template <typename K, typename V>
    // using hash_map = unordered_map<K,V>;
    using hash_map = tsl::sparse_pg_map<K,V>;
template <typename T>
    // using hash_set = unordered_set<T>;
    using hash_set = tsl::sparse_pg_set<T>;

//stable sorting of split weights
template <typename K, typename V>
  struct compare: public function<bool(pair<K,V>, pair<K,V>)> {
      constexpr bool operator()(const pair<K,V>& x, const pair<K,V>& y) const noexcept {
          return x.first > y.first || x.first == y.first && x.second < y.second;
      }
  };
template <typename K, typename V>
  using multimap_ = set<pair<K,V>, compare<K,V>>;



#include "kmer.h"
#include "kmerAmino.h"



#include "color.h"
  
#pragma once

/**
 * A tree structure that is needed for generating a NEWICK string.
 */
struct node {
    color_t taxa;
    double weight;
    vector<node*> subsets;
};

/**
* A spinlock implementation
* source: https://rigtorp.se/spinlock/
*/
struct spinlock {
  atomic<bool> lock_ = {0};

  void lock() noexcept {
    for (;;) {
      // Optimistically assume the lock is free on the first try
      if (!lock_.exchange(true, std::memory_order_acquire)) {
        return;
      }
      // Wait for lock to be released without generating cache misses
      while (lock_.load(std::memory_order_relaxed)) {
        // Issue X86 PAUSE or ARM YIELD instruction to reduce contention between
        // hyper-threads
        #if defined(__i386__) || defined(__x86_64__)
          __builtin_ia32_pause();
        #else
          std::this_thread::yield();
        #endif
      }
    }
  }

  bool try_lock() noexcept {
    // First do a relaxed load to check if lock is free in order to prevent
    // unnecessary cache misses if someone does while(!try_lock())
    return !lock_.load(std::memory_order_relaxed) &&
           !lock_.exchange(true, std::memory_order_acquire);
  }

  void unlock() noexcept {
    lock_.store(false, std::memory_order_release);
  }
};



/**
 * This class manages the k-mer/color hash tables and split list.
 */
class graph {

private:

    /**
     * This is the size of the top list.
     */
    static uint64_t t;

    /**
     * This is the min. coverage threshold for k-mers.
     */
    static vector<int> q_table;
	static int quality;

    static bool isAmino;
    
	/**
	* Look-up set for k-mers that are ignored, i.e., not stored, counted etc.
	*/
	static hash_set<kmer_t> blacklist;
	static hash_set<kmerAmino_t> blacklist_amino;	
	
    /**
     * This int indicates the number of tables to use for hashing
     */
    static uint64_t table_count;
    
    /**
     * This vector holds the carries of 2**i % table_count for fast distribution of bitset represented kmers
     */
    static vector<uint_fast32_t> period;
    /**
     * This is a vector of hash tables mapping k-mers to colors [O(1)].
     */
    static vector<hash_map<kmer_t, color_t>> kmer_table;

    /**
     * This is a vector of spinlocks protecting the hash tables.
     */
    static vector<spinlock> lock;

    /**
     * This is a hash table mapping k-mers to colors [O(1)].
     */
    static vector<hash_map<kmerAmino_t, color_t>> kmer_tableAmino;


    /**
     * This is a hash set used to filter k-mers for coverage (q > 1).
     */
    static vector<hash_set<kmer_t>> quality_set;
    static vector<hash_set<kmerAmino_t>> quality_setAmino;

	static vector<hash_map<kmer_t, uint16_t>> singleton_kmer_table;
	static vector<hash_map<kmerAmino_t, uint16_t>> singleton_kmer_tableAmino;
	static uint64_t singleton_counters[];
	static spinlock singleton_counters_locks[];

	/**
     * Count unique kmers per genome
     */
    static uint64_t kmer_counters[];

	
    /**
     * This is a hash map used to filter k-mers for coverage (q > 2).
     */
    static vector<hash_map<kmer_t, uint16_t>> quality_map;
    static vector<hash_map<kmerAmino_t, uint16_t>> quality_mapAmino;

public:

	/**
	* This function generates a bootstrap replicate. We mimic drawing n k-mers at random with replacement from all n observed k-mers. Say a k-mer would be drawn x times. Instead, we calculate x for each k-mer (in each split in color_table) from a binomial distribution (n repetitions, 1/n success rate) and calculate a new split weight according to the new number of k-mers.
	* @param mean weight function
	* @return the new list of splits of length at least t ordered by weight as usual
	*/
	static multimap_<double, color_t> bootstrap(double mean(uint32_t&, uint32_t&));

    /**
     * This is an ordered tree collecting the splits [O(log n)].
     */
    static multimap_<double, color_t> split_list;

    /**
     * This is a hash table mapping colors to weights [O(1)].
     */
    static hash_map<color_t, array<uint32_t,2>> color_table;

	/**
    * These are the allowed chars.
    */
    static vector<char> allowedChars;

    /**
     * This function initializes the top list size, coverage threshold, and allowed characters.
     *
     * @param top list size
     * @param isAmino use amino processing
     * @param count_kmers count k_mers per genome
     * @param q_table coverage threshold
	 * @param quality global q or maximum among all q values
	 * @param blacklist k-mers to be ignored
	 * @param blacklist_amino amino k-mers to be ignored
     * @param bins hash_tables to use for parallel processing
     * @param thread_count the number of threads used for processing
     */
    template<bool count_kmers>
    static void init(uint64_t& top_size, bool amino, vector<int>& q_table, int& quality, hash_set<kmer_t>& blacklist, hash_set<kmerAmino_t>& blacklist_amino, uint64_t& thread_count) {
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
                hash_kmer<count_kmers>(bin, kmer, color);
            };
            emplace_kmer_amino_tmp = [&] (const uint64_t& T, uint_fast32_t& bin, const kmerAmino_t& kmer, const uint16_t& color) {
                hash_kmer_amino<count_kmers>(bin, kmer, color);
            };
            break;

        case 2:
            isAmino ? quality_setAmino.resize(thread_count) : quality_set.resize(thread_count);
            if (q_table.size()>0){
                emplace_kmer_tmp = [&] (const uint64_t& T, uint_fast32_t& bin, const kmer_t& kmer, const uint16_t& color) {
                    if (q_table[color]==1){
                        hash_kmer<count_kmers>(bin, kmer, color);
                    } else if (quality_set[T].find(kmer) == quality_set[T].end()) {
                        quality_set[T].emplace(kmer);
                    } else {
                        quality_set[T].erase(kmer);
                        hash_kmer<count_kmers>(bin, kmer, color);
                    }
                };
                emplace_kmer_amino_tmp = [&] (const uint64_t& T, uint_fast32_t& bin, const kmerAmino_t& kmer, const uint16_t& color) {
                    if (q_table[color]==1){
                        hash_kmer_amino<count_kmers>(bin, kmer, color);
                    } else if (quality_setAmino[T].find(kmer) == quality_setAmino[T].end()) {
                        quality_setAmino[T].emplace(kmer);
                    } else {
                        quality_setAmino[T].erase(kmer);
                        hash_kmer_amino<count_kmers>(bin, kmer, color);
                    }
                };
            } else { // global quality value (one if-clause fewer)
                emplace_kmer_tmp = [&] (const uint64_t& T, uint_fast32_t& bin, const kmer_t& kmer, const uint16_t& color) {
                    if (quality_set[T].find(kmer) == quality_set[T].end()) {
                        quality_set[T].emplace(kmer);
                    } else {
                        quality_set[T].erase(kmer);
                        hash_kmer<count_kmers>(bin, kmer, color);
                    }
                };
                emplace_kmer_amino_tmp = [&] (const uint64_t& T, uint_fast32_t& bin, const kmerAmino_t& kmer, const uint16_t& color) {
                    if (quality_setAmino[T].find(kmer) == quality_setAmino[T].end()) {
                        quality_setAmino[T].emplace(kmer);
                    } else {
                        quality_setAmino[T].erase(kmer);
                        hash_kmer_amino<count_kmers>(bin, kmer, color);
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
                        hash_kmer<count_kmers>(bin, kmer, color);
                    }
                };
                emplace_kmer_amino_tmp = [&] (const uint64_t& T, uint_fast32_t& bin, const kmerAmino_t& kmer, const uint16_t& color) {
                    if (quality_mapAmino[T][kmer] < q_table[color]-1) {
                        quality_mapAmino[T][kmer]++;
                    } else {
                        quality_mapAmino[T].erase(kmer);
                        hash_kmer_amino<count_kmers>(bin, kmer, color);
                    }
                };
            }else { // global quality value
                emplace_kmer_tmp = [&] (const uint64_t& T, uint_fast32_t& bin, const kmer_t& kmer, const uint16_t& color) {
                    if (quality_map[T][kmer] < quality-1) {
                        quality_map[T][kmer]++;
                    } else {
                        quality_map[T].erase(kmer);
                        hash_kmer<count_kmers>(bin, kmer, color);
                    }
                };
                emplace_kmer_amino_tmp = [&] (const uint64_t& T, uint_fast32_t& bin, const kmerAmino_t& kmer, const uint16_t& color) {
                    if (quality_mapAmino[T][kmer] < quality-1) {
                        quality_mapAmino[T][kmer]++;
                    } else {
                        quality_mapAmino[T].erase(kmer);
                        hash_kmer_amino<count_kmers>(bin, kmer, color);
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

    static void init_count(uint64_t& top_size, bool amino, vector<int>& q_table, int& quality, hash_set<kmer_t>& blacklist, hash_set<kmerAmino_t>& blacklist_amino, uint64_t& thread_count) {
      init<true>(top_size, amino, q_table, quality, blacklist, blacklist_amino, thread_count);
    }

    static void init_noCount(uint64_t& top_size, bool amino, vector<int>& q_table, int& quality, hash_set<kmer_t>& blacklist, hash_set<kmerAmino_t>& blacklist_amino, uint64_t& thread_count) {
      init<false>(top_size, amino, q_table, quality, blacklist, blacklist_amino, thread_count);
    }


    /**
    Hash map access
    */
    
    /**
     * This function shift updates the bin of a kmer
    */
    static uint_fast32_t shift_update_bin(uint_fast32_t& bin, uint_fast8_t& left, uint_fast8_t& right);
    
    /**
     * This method shift updates the reverse complement bin for a kmer
    */
    static uint_fast32_t shift_update_rc_bin(uint_fast32_t& rc_bin, uint_fast8_t& left, uint_fast8_t& right);

    /**
    * This function shift updates a bin for an amino kmer
    */
    static uint_fast32_t shift_update_amino_bin(uint_fast32_t& bin, kmerAmino_t& kmer, uint_fast8_t& right);

    /**
     *  This method computes the bin of a given kmer(slower than shift update)
     * @param kmer The target kmer
     * @return uint64_t The bin
     */
     
    static uint_fast32_t compute_bin(const kmer_t& kmer);

    /**
     *  This function computes the bin of a given amino kmer(slower than shift update)
     * @param kmer The target kmer
     * @return uint64_t The bin
     */
    static uint_fast32_t compute_amino_bin(const kmerAmino_t& kmer);

    /**
    * This function hashes a k-mer and stores it in the correstponding hash table.
    * The corresponding table is chosen by the carry of the encoded k-mer given the number of tables as module.
    * @param bin The bin, the kmer is stored in
    * @param kmer The kmer to store
    * @param color The color to store 
    */
    template<bool count_kmers>
    static void hash_kmer(uint_fast32_t& bin, const kmer_t& kmer, const uint16_t& color)
    {
        lock[bin].lock();
        hash_map<kmer_t,color_t>::iterator entry=kmer_table[bin].find(kmer); 
        // already in the kmer table?
        if(entry != kmer_table[bin].end()){
            if(count_kmers){
              // check if not seen in this genome before, i.e., count a new (unique) kmer
              if(!entry.value().test(color)){
                // count
                kmer_counters[color]++;
              }
            }
            // add
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
                    if(count_kmers){
                      // count
                      kmer_counters[color]++;
                    }
                }
            }
            // not seen before -> add to singleton_table
            else{
                singleton_kmer_table[bin][kmer]=color;
                singleton_counters_locks[color].lock();
                singleton_counters[color]++;
                singleton_counters_locks[color].unlock();
                if(count_kmers){
                  // count
                  kmer_counters[color]++;
                }
            }
        }
        lock[bin].unlock();
    }



    /**
    * This function hashes an amino k-mer and stores it in the corresponding hash table.
    * The correspontind table is chosen by the carry of the encoded k-mer bitset by the bit-module function.
    * @param bin The bin, the kmer is stored in
    * @param kmer The kmer to store
    * @param color The color to store
    */
    template<bool count_kmers>
    static void hash_kmer_amino(uint_fast32_t& bin, const kmerAmino_t& kmer, const uint16_t& color)
    {
        lock[bin].lock();
        hash_map<kmerAmino_t,color_t>::iterator entry=kmer_tableAmino[bin].find(kmer); 
        // already in the kmer table? -> add
        if(entry != kmer_tableAmino[bin].end()){
            // check if not seen in this genome before, i.e., count a new (unique) kmer
            if(!entry.value().test(color)){
                // add
                entry.value().set(color);
                // count
                kmer_counters[color]++;
            }
        }
        // not yet in the kmer table?
        else{
            // count
            kmer_counters[color]++;
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
     * This function searches the bit-wise corresponding hash table for the given kmer
     * @param kmer The kmer to search
     * @return True if the kmer is stored.
     */
    static bool search_kmer(const kmer_t& kmer);

    /**
     * This function searches the bit-wise corresponding hash table for the given amnio kmer
     * @param kmer The kmer to search
     * @return Ture if the kmer is stored
     */
    static bool search_kmer_amino(const kmerAmino_t& kmer);


    /**
     * This function returns the stored colores of the given kmer
     * @param kmer The target kmer
     * @return color_t The stored colores
     */
    static color_t get_color(const kmer_t& kmer, bool reversed);

    /**
     * This function returns the stored color of the given kmer
     * @param kmer The target kmer
     * @return color_t The stored color
     */
    static color_t get_color_amino(const kmerAmino_t& kmer);


    /**
     * This function removes the kmer entry from the hash map
     */
    static void remove_kmer(const kmer_t& kmer, bool reversed);

    /**
     * This function removes the amino kmer entry from the corresponding hash table
     */
    static void remove_kmer_amino(const kmerAmino_t& kmer);

	/**
	* This function extracts k-mers from a sequence and adds them to the black list.
	*/
	static void fill_blacklist(string& str, bool& reverse);
	
	/**
	* This function tells how many k-mers are in the black list.
	*/
	static uint64_t size_blacklist();
	
	/**
	 * This function activates using the blacklist while inserting kmers.
	 */
	static void activate_blacklist();
	
    /**
     * This function extracts k-mers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     * @param reverse merge complements
     */
    static void add_kmers(uint64_t& T, string& str, uint16_t& color, bool& reverse);

    /**
     * This function extracts k-mer minimizers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     * @param reverse merge complements
     * @param m number of k-mers to minimize
     */
    static void add_minimizers(uint64_t& T, string& str, uint16_t& color, bool& reverse, uint64_t& m);

    /**
     * This function extracts k-mers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     * @param reverse merge complements
     * @param max_iupac allowed number of ambiguous k-mers per position
     */
    static void add_kmers(uint64_t& T, string& str, uint16_t& color, bool& reverse, uint64_t& max_iupac);

    /**
     * This function extracts k-mer minimizers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     * @param reverse merge complements
     * @param m number of k-mers to minimize
     * @param max_iupac allowed number of ambiguous k-mers per position
     */
    static void add_minimizers(uint64_t& T, string& str, uint16_t& color, bool& reverse, uint64_t& m, uint64_t& max_iupac);

	/**
	* This function calculates the weight for all splits and puts them into the split_ölist
	* @param mean weight function
	* @param min_value the minimal weight represented in the top list
	*/
	static void compile_split_list(double mean(uint32_t&, uint32_t&), double min_value);

	/**
	* This function determines the core k-mers, i.e., all k-mers present in all genomes.
	* Core k-mers are output to given file in fasta format, one k-mer per entry
	* @param file output file stream
	* @param verbose print progess
	*/
	static void output_core(ostream& file, bool& verbose);
	
	
	/**
	* Get the number of k-mers in all tables.
	* @return number of k-mers in all tables.
	*/
	static uint64_t number_kmers();
	
	/**
	* Get the number of k-mers in all singleton tables.
	* @return number of k-mers in all singleton tables.
	*/
	static uint64_t number_singleton_kmers();
    
    
    /**
    * Get the number of singleton k-mers in a given genome.
    * @param genome the id (rf. denom_names) of the genome
    * @return number of singleton k-mers read in that genome.
    */
    static uint64_t number_singleton_kmers(uint16_t& genome);
    


    /**
    * Get the number of unique k-mers in a given genome.
    * @param genome the id (rf. denom_names) of the genome
    * @return number of unique k-mers read in that genome.
    */
    static uint64_t number_kmers(uint16_t& genome);

    /**
     * @return mean value of number of k-mers per genome
     */
    static double mean_number_kmers();

    /**
     * @return standard deviation of number of k-mers per genome
     */
    static double stdev_number_kmers();

	/**
     * This function iterates over the hash table and calculates the split weights.
     *
     * @param mean weight function
     * @param verbose print progress
     * @param min_value the minimal weight currently represented in the top list
     */
    static void add_weights(double mean(uint32_t&, uint32_t&), double min_value, bool& verbose);
	
	
	/**
	* This function iterates over the singleton tables and adds the split weights.
	* 
	* @param mean weight function
	* @param min_value the minimal weight represented in the top list
	* @param verbose print progess
	*/
	static void add_singleton_weights(double mean(uint32_t&, uint32_t&), double min_value, bool& verbose);


    /**
     * This function adds a single split (weight and colors) to the output list.
     *
     * @param weight split weight
     * @param color split colors
	 * (@param split_list list of splits to add the split to)
     */
    static void add_split(double& weight, color_t& color);
	static void add_split(double& weight, color_t& color, multimap_<double, color_t>& split_list);


    /**
     * This funtion adds a sigle split from a cdbg to the output list.
     * 
     * @param mean mean function
     * @param kmer_seq kmer
     * @param kmer_color split colors
     * @param min_value the minimal weight currently represented in the top list
     */
     static void add_cdbg_colored_kmer(string kmer_seq, const uint16_t& kmer_color);       

    /**
     * This function clears color-related temporary files.

     * @param T thread to clear
     */
    static void clear_thread(uint64_t& T);

    /**
     * This function filters a greedy maximum weight tree compatible subset.
     *
     * @param split_list list of splits to be filtered
     * @param verbose print progress
     * @return the new minimal weight represented in the top list
     */
    static void filter_strict(multimap_<double, color_t>& split_list, bool& verbose);

    /**
     * This function filters a greedy maximum weight tree compatible subset and returns a newick string.
     *
     * @param map function that maps an integer to the original id, or null
     * @param split_list list of splits to be filtered
	 * @param support_values a hash map storing the absolut support values for each color set
	 * @param bootstrap_no the number of bootstrap replicates for computing the per centage support
     * @param verbose print progress
     */
	static string filter_strict(std::function<string(const uint16_t&)> map, multimap_<double, color_t>& split_list, hash_map<color_t, uint32_t>* support_values, const uint32_t& bootstrap_no, bool& verbose);

    /**
     * This function filters a greedy maximum weight weakly compatible subset.
     *
     * @param split_list list of splits to be filtered
     * @param verbose print progress
     */
    static void filter_weakly(multimap_<double, color_t>& split_list, bool& verbose);

    /**
     * This function filters a planar graph compatible subset.
     * @param split_list list of splits to be filtered
     * @param verbose print progress
     */
    static void filter_planar(multimap_<double, color_t>& split_list, bool& verbose, uint64_t& num);


    /**
     * This function filters a greedy maximum weight n-tree compatible subset.
     *
     * @param n number of trees
     * @param split_list list of splits to be filtered
     * @param verbose print progress
     */
    static void filter_n_tree(uint64_t n, multimap_<double, color_t>& split_list, bool& verbose);

    /**
     * This function filters a greedy maximum weight n-tree compatible subset and returns a string with all trees in newick format.
     *
     * @param n number of trees
     * @param map function that maps an integer to the original id, or null
     * @param split_list list of splits to be filtered
	 * @param support_values a hash map storing the absolut support values for each color set
	 * @param bootstrap_no the number of bootstrap replicates for computing the per centage support
     * @param verbose print progress
     */
    static string filter_n_tree(uint64_t n, std::function<string(const uint16_t&)> map, multimap_<double, color_t>& split_list, hash_map<color_t, uint32_t>* support_values, const uint32_t& bootstrap_no, bool& verbose);
	
	


protected:

    /**
     * This function qualifies a k-mer and places it into the hash table.
     *
     * @param T number of thread processing the kmer
     * @param bin The bin, the kmer is stored in
     * @param kmer bit sequence
     * @param color color flag
     */
    static function<void(const uint64_t& T, uint_fast32_t& bin, const kmer_t&, const uint16_t&)> emplace_kmer;
    static function<void(const uint64_t& T, uint_fast32_t& bin, const kmer_t&, const uint16_t&)> emplace_kmer_tmp;

    /**
     * This function qualifies a k-mer and places it into the hash table.
     *
     * @param T number of thread processing the kmer
     * @param bin The bin, the kmer is stored in
     * @param kmer bit sequence
     * @param color color flag
     */
    static function<void(const uint64_t& T, uint_fast32_t& bin, const kmerAmino_t&, const uint16_t&)> emplace_kmer_amino;
    static function<void(const uint64_t& T, uint_fast32_t& bin, const kmerAmino_t&, const uint16_t&)> emplace_kmer_amino_tmp;

    /**
     * This function tests if a split is compatible with an existing set of splits.
     *
     * @param color new split
     * @param color_set set of splits
     * @return true, if compatible
     */
    static bool test_strict(color_t& color, vector<color_t>& color_set);

    /**
     * This function tests if a split is weakly compatible with an existing set of splits.
     *
     * @param color new split
     * @param color_set set of splits
     * @return true, if weakly compatible
     */
    static bool test_weakly(color_t& color, vector<color_t>& color_set);

    /**
     * This function calculates the multiplicity of iupac k-mers.
     *
     * @param product overall multiplicity
     * @param factors per base multiplicity
     * @param input iupac character
     */
    static void iupac_calc(long double& product, vector<uint8_t>& factors, char& input);

    /**
     * This function shifts a base into a set of ambiguous iupac k-mers.
     *
     * @param prev set of k-mers
     * @param next set of k-mers
     * @param input iupac character
     */
    static void iupac_shift(hash_set<kmer_t>& prev, hash_set<kmer_t>& next, char& input);

    /**
   * This function shifts a base into a set of ambiguous iupac k-mers.
   *
   * @param prev set of k-mers
   * @param next set of k-mers
   * @param input iupac character
   */
    static void iupac_shift_amino(hash_set<kmerAmino_t>& prev, hash_set<kmerAmino_t>& next, char& input);

    /**
     * This function returns a tree structure (struct node) generated from the given list of color sets.
     *
     * @param color_set list of color sets
     * @return tree structure (struct node)
     */
    static node* build_tree(vector<color_t>& color_set);

    /**
     * This function recursively refines a given set/tree structure by a given split.
     *
     * @param current_set node of currently considered (sub-)set/tree structure
     * @param split color set to refine by
     * @return whether or not the given split is compatible with the set/tree structure
     */
    static bool refine_tree(node* current_set, color_t& split, color_t& allTaxa);

    /**
     * This function returns a newick string generated from the given tree structure (set).
     *
     * @param root root of the tree/set structure
     * @param map function that maps an integer to the original id, or null
	 * (@param support_values a hash map storing the absolut support values for each color set)
	 * (@param bootstrap_no the number of bootstrap replicates for computing the per centage support)
     * @return newick string
     */
    static string print_tree(node* root, std::function<string(const uint16_t&)> map, hash_map<color_t, uint32_t>* support_values, const uint32_t& bootstrap_no);
    static string print_tree(node* root, std::function<string(const uint16_t&)> map);

    /**
     * This function checks if the character at the given position is allowed.
     * @param pos position in str
     * @param str the current part of the sequence
     * @return true if allowed, false otherwise
     */
    static bool isAllowedChar(uint64_t pos, string &str);
};
