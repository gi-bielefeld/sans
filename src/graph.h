#include <iostream>
#include <limits>
#include <algorithm>
#include <functional>
#include <utility>
#include <vector>
#include <chrono>

#include <thread>
#include <atomic>
#include <mutex>
#include <shared_mutex>
#include <condition_variable>

#include <iomanip>
#include <string>
#include <random>



#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include "tsl/sparse_map.h"
#include "tsl/sparse_set.h"
#include <bitset>

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

// For Fingerprinting: A bitset 128 bits long should be more than sufficient
#ifndef FL
#define FL 128
#endif

// due to undefined behaviour of negative biteshifting the table_number must be shorter than the fingerprint.
constexpr static uint_fast16_t Fprint_length = FL < 16 ? 16 : FL;      // this could be made a function of maxN

using fingerprint = bitset<Fprint_length>;

// This is a struct holding the pair<counter, color_set_vector>
struct Fprint_table_entry {
    // the counter is ATOMIC for multithreading
    std::atomic<uint_fast32_t> counter{0};
    const color_t color_set;        // the color set is never actually changed in the entry, only between them. So it can be const.
    
    Fprint_table_entry(uint_fast32_t count, color_t c) : counter(count), color_set(c) {}
};
// a std::unique_pointer in the hash_map will point to the Fprint_table_entry stored on the heap!
// an entry will never be copied or moved - that is impossible with an atomic counter. 

class Measure_time {
    // With this you can measure the total execution time of any line of code.
private:
    chrono::high_resolution_clock::duration execution_time{};
    chrono::high_resolution_clock::time_point begin;
    chrono::high_resolution_clock::time_point end;

public:
    Measure_time() : execution_time(chrono::high_resolution_clock::duration::zero()) {}

    void start_clock() {
        begin = chrono::high_resolution_clock::now();
    }

    void stop_clock() {
        end = chrono::high_resolution_clock::now();
        execution_time += end - begin;
    }

    void reset() {
        execution_time = chrono::high_resolution_clock::duration::zero();
    }

    chrono::high_resolution_clock::duration elapsed() const {
        return execution_time;
    }

    void show_time() const {
        cout << util::format_time(execution_time) << endl;
    }
};

// Debug variables:
inline Measure_time opti, caut;
inline uint_fast32_t same_kmer_counter = 0;
inline uint_fast32_t same_kmer_in_singletons_counter = 0;

// Poor multithreading performance troubleshooting 
inline uint_fast32_t lock_counter = 0;      
inline uint_fast64_t optimistic_counter = 0;
inline uint_fast64_t cautious_counter = 0;

/*
 *  This class is for the purpose of synchronising threads during cleanup.
 *  Courtesy of Gemini
*/
class ThreadBarrier {
private:
    std::mutex mtx;
    std::condition_variable cv;
    size_t total_threads;
    size_t waiting_threads;
    size_t generation; // Tracks cycles to allow safe reuse

public:
    // ThreadBarrier needs to be default-constructible
    ThreadBarrier() 
        : total_threads(0), waiting_threads(0), generation(0) {}

    explicit ThreadBarrier(size_t threads) 
        : total_threads(threads), waiting_threads(threads), generation(0) {}

    void init(size_t threads) {
        std::unique_lock<std::mutex> lock(mtx);
        total_threads = threads;
        waiting_threads = threads;
    }

    void arrive_and_wait() {
        std::unique_lock<std::mutex> lock(mtx);
        size_t current_generation = generation;

        if (--waiting_threads == 0) {
            std::cout << "--- BARRIER TRIPPED BY THREAD " << " ---" << std::endl;
            // The last thread has arrived! 
            generation++;               // Advance to the next cycle
            waiting_threads = total_threads; // Reset counter for next time
            cv.notify_all();            // Wake up everyone
        } else {
            // Wait until the generation changes (woken up by the last thread)
            cv.wait(lock, [this, current_generation] {
                return current_generation != generation;
            });
        }
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


    static uint_fast32_t F_print_table_count; 

    /**
     * This is a threshold for efficient cleaning of the Fprint_tables from 0-entries.
     */
    inline static uint_fast32_t cleaning_threshold;

    /**
     * The cleaning is best to do after ALL threads have finished their "round".
     * In the next round, zero-value entries will not be needed for sure.
     */
    inline static uint_fast32_t cleaning_round = 1;

    /**
     * Limits the maximum number of cleaning rounds that can take place.
     */
    inline static uint_fast32_t max_cleaning_rounds;

    /**
     * This is the related counter.
     */
    inline static uint_fast32_t inserts_counter;

    /**
     * This variable is used to check if all threads are active for cleaning.
     */
    inline static atomic<uint_fast32_t> threads_active;

    /**
     *  This barrier makes all threads stop at at sync point to perform the cleaning in a synchronised manner.
     */
    inline static ThreadBarrier sync_point;
    
    /**
     * This vector holds the carries of 2**i % table_count for fast distribution of bitset represented kmers
     */
    static vector<uint_fast32_t> period;

    /**
     *   FINGERPRINTING
     *  The collision probability is bounded above by N^2 / 2^l
     *  where N is the number of kmers and l length of the fingerprint.
     *    - added Fprint_table to store [key : value] = [fingerprint : count, color_set_vector]
     *    - count is the number of kmers that are present in the same set of genomes.
    */    

    /**
     * This vector holds the fingerprints assigned to each genome.
     */
    static vector<fingerprint> genome_fingerprints;

    /**
    * changed kmer_table to store binary fingerprints of length l
    */
    static vector<hash_map<kmer_t, fingerprint>> kmer_table;
    
    /**
    * This table holds unique pointers to pairs of (count, color set) for each fingerprint.  
    */
    static vector<hash_map<fingerprint, std::shared_ptr<Fprint_table_entry>>> Fprint_table;

    /**
     *  Used for the calculation of the fingerprint bin.
     */
    static uint_fast32_t shift_bits_by;

    /**
     * This is a vector of spinlocks protecting the hash tables.
     */
    static vector<spinlock> lock;

    /**
     * This is a vector of shared locks protecting the fingerprint table.
     */
    static vector<std::shared_mutex> F_lock;

    inline static uint64_t threads;    // number of threads

    /**
     * This is a hash table mapping k-mers to fingerprints [O(1)].
     */
    static vector<hash_map<kmerAmino_t, fingerprint>> kmer_tableAmino;

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

    /**     // For now the fingerprinting is not completed in this function for the AMINO case.
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
        threads = thread_count;
        t = top_size;
        isAmino = amino;
        if(!isAmino){

            // Automatic table count
            //    table_count = 45 * thread_count - 33; // Estimated scaling
            //    table_count = table_count % 2 ? table_count : table_count + 1; // Ensure the table count is odd
            table_count = (0b1u << 14) + 1;
            // table_count as a number will be a 1 followed by fourteen 0
            // for the computation of fingerprint table bin:
            // so that after shifting, the fingerprint is 14 valid digits long:

            // default setting:
            short n_of_digits = 14; //14;
            F_print_table_count = (0b1u << n_of_digits);
            shift_bits_by = Fprint_length - n_of_digits;  
            // cout << shift_bits_by;   
            // cout << table_count << "\n";
             cout << "F_print_table_count: " << F_print_table_count << "\n";

            // Create random binary fingerprints
            genome_fingerprints = vector<fingerprint> (color::n);
            random_device rd; // Seed source
            mt19937 gen(rd());
            uniform_int_distribution<uint32_t> dist;

            int8_t takes = Fprint_length / 32;
            int8_t rest  = Fprint_length % 32;
            string zeros_ones; 
            bool random_fingerprints = true;

            if (random_fingerprints){
                for (uint1N_t i = 0; i < color::n; i++){
                    for (int8_t j = 0; j < takes; j++){
                        zeros_ones += bitset<32>(dist(gen)).to_string();
                    }
                    if (rest != 0){
                        zeros_ones += bitset<32>(dist(gen)).to_string().substr(0, rest);
                    }
                    genome_fingerprints[i] = bitset<Fprint_length>(zeros_ones);
                    zeros_ones.clear();

                    // cout << "Fingerprint preview: " << genome_fingerprints[i] << endl;
                }
                //cout << "Fingerprint preview: " << genome_fingerprints[0] << endl;
            }
            else { // each fingerprint is just one color. 
                   // This was for the purpose of testing / debugging of the effect of collisions on results.
                for (uint1N_t i = 0; i < color::n; i++){
                    fingerprint f;
                    f.set(i % Fprint_length);
                    genome_fingerprints[i] = f;
                    //cout << "Fingerprint preview: " << genome_fingerprints[i] << endl;
                }
            }

            // Init base tables
            kmer_table = vector<hash_map<kmer_t, fingerprint>> (table_count);
            singleton_kmer_table = vector<hash_map<kmer_t, uint16_t>> (table_count);
            Fprint_table = vector<hash_map<fingerprint, std::shared_ptr<Fprint_table_entry>>> (F_print_table_count);

            // For multithreading:
            // Init the lock vector
            lock = vector<spinlock> (table_count);
            F_lock = vector<std::shared_mutex> (F_print_table_count);

            // For cleaning zero-value entries:
            cleaning_threshold = F_print_table_count;
            max_cleaning_rounds = color::n / threads;
            // initialise the ThreadBarrier sync_point:
            sync_point.init(threads);

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
            kmer_tableAmino = vector<hash_map<kmerAmino_t, fingerprint>> (table_count);
            singleton_kmer_tableAmino = vector<hash_map<kmerAmino_t, uint16_t>> (table_count);

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
     * This method computes the bin of a given fingerprint
     * @param f The target fingerprint
     * @return uint64_t The bin
     */
    static uint_fast32_t compute_Fprint_bin(const fingerprint& f);

    /**
     *  This function computes the bin of a given amino kmer(slower than shift update)
     * @param kmer The target kmer
     * @return uint64_t The bin
     */
    static uint_fast32_t compute_amino_bin(const kmerAmino_t& kmer);

    // /**     Last version:
    // * This function hashes a k-mer and stores it in the correstponding hash table.
    // * The corresponding table is chosen by the carry of the encoded k-mer given the number of tables as module.
    // * @param bin The bin, the kmer is stored in
    // * @param kmer The kmer to store
    // * @param color The color to store 
    // */
    // template<bool count_kmers>
    // static void hash_kmer(uint_fast32_t& bin, const kmer_t& kmer, const uint16_t& color)
    // {
    //     lock[bin].lock();
    //     hash_map<kmer_t,color_t>::iterator entry=kmer_table[bin].find(kmer); 
    //     // already in the kmer table?
    //     if(entry != kmer_table[bin].end()){
    //         if(count_kmers){
    //           // check if not seen in this genome before, i.e., count a new (unique) kmer
    //           if(!entry.value().test(color)){
    //             // count
    //             kmer_counters[color]++;
    //           }
    //         }
    //         // add
    //         entry.value().set(color);
    //     }
    //     // not yet in the kmer table?
    //     else{
    //         hash_map<kmer_t,uint16_t>::iterator s_entry = singleton_kmer_table[bin].find(kmer);
    //         //seen once before? -> add to kmer table / remove from singleton table
    //         if(s_entry != singleton_kmer_table[bin].end()){
    //             if(s_entry.value() != color){
    //                 kmer_table[bin][kmer].set(s_entry.value());
    //                 kmer_table[bin][kmer].set(color);
    //                 singleton_counters_locks[s_entry.value()].lock();
    //                 singleton_counters[s_entry.value()]--;
    //                 singleton_counters_locks[s_entry.value()].unlock();
    //                 singleton_kmer_table[bin].erase(s_entry);
    //                 if(count_kmers){
    //                   // count
    //                   kmer_counters[color]++;
    //                 }
    //             }
    //         }
    //         // not seen before -> add to singleton_table
    //         else{
    //             singleton_kmer_table[bin][kmer]=color;
    //             singleton_counters_locks[color].lock();
    //             singleton_counters[color]++;
    //             singleton_counters_locks[color].unlock();
    //             if(count_kmers){
    //               // count
    //               kmer_counters[color]++;
    //             }
    //         }
    //     }
    //     lock[bin].unlock();
    // }

    /**     
    * This function hashes a k-mer and stores it in the corresponding hash table.
    * The corresponding table is chosen by the carry of the encoded k-mer given the number of tables as module.
    * Kmers that share a color set (and thus their fingerprint) are written into the Fprint_table.
    * @param bin The bin, the kmer is stored in
    * @param kmer The kmer to store
    * @param color The color to store 
    */
    template<bool count_kmers>
    static void hash_kmer(uint_fast32_t& bin, const kmer_t& kmer, const uint16_t& color)
    {
        lock[bin].lock();
        hash_map<kmer_t,fingerprint>::iterator entry = kmer_table[bin].find(kmer);
        // already in the kmer table?
        if(entry != kmer_table[bin].end()){
            fingerprint F_old = entry.value();  
            uint_fast32_t bin_F_old = compute_Fprint_bin(F_old);

            // // To see if many same fingerprints are seeked in waves at the same time
            // cout << F_old << "\n";   // - not the case!

            // Note to Multi-threading  
            // Although we have 32 000 bins, for some data there might be much fewer color sets, or the majority of
            // kmers concentrated in few color sets, which causes bin contention.
            // To avoid bin contention, find()-ing entries is allowed to take place
            // concurrently. However, an insertion, even if the table doesn't grow and rehash, may invalidate any existing iterators.  
            // Therefore, three levels of protection are implemented:
            // 1. shared_lock       - for concurrent finding of entries in the same table
            // 2. atomic counter    - for incrementing the counter thread-safe
            // 3. unique_lock       - for inserts - the most dangerous operation.
            // 
            // ...This approach is good because adding a new color set (that is, inserting a new entry)
            // is a much less frequent operation than incrementing the counter in an existing entry
            // of the Fprint_table. 

            // Find, if the kmer has been seen in the current genome or not.

            // Create an instance of shared_mutex ( is destructed automatically when it goes out of scope)
            std::shared_lock<std::shared_mutex> shared_F_lock_old(F_lock[bin_F_old]);
            // The bin is now shared-locked

            hash_map<fingerprint, std::shared_ptr<Fprint_table_entry>>::iterator it1 = Fprint_table[bin_F_old].find(F_old);

            // Sanity check, remove after testing is successful: 
            if (it1 == Fprint_table[bin_F_old].end()){
                // this must not happen
                int raise(404);         // page not found error.
            }

            // Get entry stored on the heap through the unique pointer
            Fprint_table_entry& Fpt_entry_old = *(it1->second); 

            // we don't want to xor fingerprints of the same genome twice.
            // check if not seen in this genome before, i.e., count a new (unique) kmer
            if(0 == Fpt_entry_old.color_set.test(color)){
                // remove the kmer from the color-set it belonged to previously:
                // thanks to atomic we need no additional locks
                // Why does Adrain assume this here:
                assert(Fpt_entry_old.counter > 0);
                Fpt_entry_old.counter.fetch_sub(1, std::memory_order_relaxed);    // the fastest way
                // erasing of 0-value entries is done globally once in a time
            
                // no need to hold the lock any longer:
                shared_F_lock_old.unlock();

                if (count_kmers){
                    // count
                    kmer_counters[color]++;
                }

                // do the xor - combine fingerprints
                fingerprint F_new = F_old ^ genome_fingerprints[color];
                uint_fast32_t bin_F_new = compute_Fprint_bin(F_new);

                // update the fingerprint table - increment an existing color set or create a new one (i.e. insert new entry)

                std::shared_lock<std::shared_mutex> shared_F_lock_new(F_lock[bin_F_new]);
                hash_map<fingerprint, std::shared_ptr<Fprint_table_entry>>::iterator it2 = Fprint_table[bin_F_new].find(F_new);
                // We cannot derefence if the iterator is end()!
                
                // cout << "bin_F_new: " << bin_F_new << endl;  // works
                // cout << "bin_F_old: " << bin_F_old << endl;
                    
                if (it2 != Fprint_table[bin_F_new].end()){
                    // Get entry stored on the heap through the unique pointer
                    Fprint_table_entry& Fpt_entry_new = *(it2->second);
                    // Here was the last segfault. Let's print the iterator.
                    // cout << "Printing Fpt_entry_new (color = " << color << ") " << static_cast<const void*>(&(*Fpt_entry_new)) << endl;
                    Fpt_entry_new.counter.fetch_add(1, std::memory_order_relaxed);
                    shared_F_lock_new.unlock();
//                    cout  << bin_F_new << "\n" << flush;
                } else {
                    // add a new fingerprint
                    // WARNING: Fpt_entry_old on the heap is not protected by lock!
                    // the cleaning method could erase it, if all threads were not synchronised.
                    color_t color_set_vector = Fpt_entry_old.color_set;
                    color_set_vector.set(color);
                    // CRITICAL SECTION - insert entry
                    // Before acquiring unique_lock, relinquish the shared lock.
                    shared_F_lock_new.unlock();
                    // in this moment another thread might slip into the queue, acquire unique access to the bin, 
                    // and perform the insertion of the very same element! *
                    // Before a unique lock is acquired, all shared_locks are released. 
                    std::unique_lock<std::shared_mutex> unique_F_lock(F_lock[bin_F_new]);
                    // * therefore, we have to re-check if this didn't happen.
                    hash_map<fingerprint, std::shared_ptr<Fprint_table_entry>>::iterator it2 = Fprint_table[bin_F_new].find(F_new);
                    // if the element is still not in the table:
                    if (it2 == Fprint_table[bin_F_new].end()){
                        // Create the shared_ptr on a new Fprint_table_entry object on the heap
                        Fprint_table[bin_F_new][F_new] = std::make_unique<Fprint_table_entry>(1, color_set_vector);
                        inserts_counter++;
//                        cout  << bin_F_new << "\n" << flush;
                    } else {
                        // we should just increment
                        Fprint_table_entry& Fpt_entry_new = *(it2->second); 
                        Fpt_entry_new.counter.fetch_add(1, std::memory_order_relaxed);      
                    }
                }   // the unique lock is here released automatically.
                // update the kmer_table
                entry.value() = F_new;
                // The kmer table does not need multi-threading optimisation,
                // since with 32 000 bins the likelihood of two threads waiting for each other is very small.
            }
        }
        // not yet in the kmer table?
        else{
            hash_map<kmer_t, uint16_t>::iterator s_entry = singleton_kmer_table[bin].find(kmer);
            //seen once before? -> add to kmer table / remove from singleton table
            if(s_entry != singleton_kmer_table[bin].end()){
                if(s_entry.value() != color){

                    // printout_tables("kmer_table.txt", "Fprint_table.txt");  

                    uint16_t &color2 = s_entry.value();
                    fingerprint F_new = genome_fingerprints[color2] ^ genome_fingerprints[color];
                    kmer_table[bin][kmer] = F_new;  
                                 
                    singleton_counters_locks[s_entry.value()].lock();
                    singleton_counters[s_entry.value()]--;
                    singleton_counters_locks[s_entry.value()].unlock();
                    if(count_kmers){
                      // count
                      kmer_counters[color]++;
                    }

                    // update the Fprint_table
                    uint_fast32_t bin_F_new = compute_Fprint_bin(F_new);
                    std::shared_lock<std::shared_mutex> shared_F_lock_new(F_lock[bin_F_new]);
                    hash_map<fingerprint, std::shared_ptr<Fprint_table_entry>>::iterator it = Fprint_table[bin_F_new].find(F_new);

                    // increment or create new fingerprint
                    if (it != Fprint_table[bin_F_new].end()){
                        // fingerprint already exists, add 1 kmer to the count
                        Fprint_table_entry& Fpt_entry = *(it->second);
                        Fpt_entry.counter.fetch_add(1, std::memory_order_relaxed);
                        shared_F_lock_new.unlock();
//                        cout  << bin_F_new << "\n" << flush;
                    } else {
                        // Initialise an empty color set
                        color_t color_set_vector(0);    
                        color_set_vector.set(color); color_set_vector.set(color2);
                        // cout << "color_set_vector: " << color_set_vector << endl; 
                        // CRITICAL SECTION - insert
                        shared_F_lock_new.unlock();
                        std::unique_lock<std::shared_mutex> unique_F_lock(F_lock[bin_F_new]);
                        // As before, we need to be very careful.
                        hash_map<fingerprint, std::shared_ptr<Fprint_table_entry>>::iterator it = Fprint_table[bin_F_new].find(F_new);
                        // if the element is still not in the table:
                        if (it == Fprint_table[bin_F_new].end()){
                            // Create the shared_ptr on a new Fprint_table_entry object on the heap
                            Fprint_table[bin_F_new][F_new] = std::make_unique<Fprint_table_entry>(1, color_set_vector);
                            inserts_counter++;
//                            cout << bin_F_new << "\n" << flush;
                        } else {
                            // we should just increment
                            Fprint_table_entry& Fpt_entry = *(it->second); 
                            Fpt_entry.counter.fetch_add(1, std::memory_order_relaxed);      
                        }
                    }   // and unique_lock goes out of scope.
                    
                    singleton_kmer_table[bin].erase(s_entry);   

                    // debug
                    // printout_tables("kmer_table.txt", "Fprint_table.txt"); 
                }//else {
                //same_kmer_in_singletons_counter++;
                //}
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

        //##### DO NOT ERASE ANYTHING!
 /*
        //
        // ========== CLEANING zero entries for saving memory ===========
        // 

        // thread_local is a special variable initialised only once for each thread:
        bool thread_local completed_my_round = false;
        // this creates a deterministic order, so that each thread gets one number from 1 to |threads|
        uint_fast16_t thread_local my_segment; 

        // How this works:
        // no thread even thinks about cleaning unless the threshold is reached:
        if (inserts_counter > cleaning_threshold){
            // now, some threads might be behind while others already processing the next genomes.
            // we should delay the cleaning until all threads complete the processing of their last genome. 
            if (completed_my_round == false){
                if (color >= cleaning_round * threads &&
                    cleaning_round < max_cleaning_rounds) {
                    completed_my_round = true;
                    my_segment = threads_active.fetch_add(1, std::memory_order_relaxed);
                    }
                // still not completed? - continue.
                else return;
            }
            // if a thread arrives here it means it is ready for cleaning

            // if (color == color::n - 1){
            //     // trigger a global cleaning for the last time:
            //     threads_active = threads;
            //     cleaning_round = max_cleaning_rounds;       // no more cleaning after this.
            // }
            
            if (threads_active >= threads){
                // means that all threads are ready for cleaning.

                // Distribute the bins across threads evenly using math:
                // If you have 10 threads and 1000 bins, Thread 0 cleans 0-99, Thread 1 cleans 100-199, etc.

                size_t bins_per_thread = F_print_table_count / threads; 
                size_t start_bin = (my_segment % threads) * bins_per_thread;
                size_t end_bin = start_bin + bins_per_thread;
                uint_fast32_t n_of_nonzero_entries = 0;
                
                // Because of that one unprotected reference to the color set vector
                // We have to synchronise all threads here so that no reads are concurrent. (It's a reference to the color set vector to avoid copying it, because it can be huge. )
                sync_point.arrive_and_wait(); 

                // Loop through and clean the designated segment of bins
                for (size_t i = start_bin; i < end_bin; ++i) {
                    cleanup_zero_entries(i);
                    n_of_nonzero_entries += Fprint_table[i].size();
                }

                // to take care of the division rest: this is the last thread's task
                if (color % threads == threads-1 && end_bin < F_print_table_count){
                    for (size_t i = end_bin; i < F_print_table_count; ++i) {
                        cleanup_zero_entries(i);
                        n_of_nonzero_entries += Fprint_table[i].size();
                    }
                }

                // IMPORTANT: reset active threads and inserts_counter
                if (my_segment == 0){
                    // Update the cleaning limit depending on the proportion of non-zero entries in the tables.
                    cleaning_threshold = max(cleaning_threshold, 2 * n_of_nonzero_entries);
                    cleaning_round++;
                    threads_active = 0;     
                    inserts_counter = 0; 
                }
                // all threads:
                completed_my_round = false;

                sync_point.arrive_and_wait(); 

                // So for example, if we have 32 000 bins and 32 threads,
                // each thread cleans 1000 bins. If one finds, that after the last cleaning 
                // the average number of nonzero entries is close to zero, the threshold stays the same.
                // However, if it finds that after cleaning there are 10 nonzero entries, 
                // the cleaning limit will be set to 20, so that the tables are not cleaned too often
                // with little effect. This keeps the cost of cleaning per zero entry bound by a constant.
                // (The same principle as extending a vector to double size.)
            }

        }*/
    }


    static void cleanup_all_zero_entries(){
        // Loop through all bins
        for (size_t i = 0; i < F_print_table_count; ++i) {
            cleanup_zero_entries(i);
        }
    }

    static void cleanup_zero_entries(uint_fast32_t bin) {
    // remove Fprint_table_entries that have count == 0, because they are no longer needed.
    // 1. Lock EXCLUSIVELY. No one else can read/write to this bin during cleanup.
    std::unique_lock<std::shared_mutex> lock(F_lock[bin]);
    
    auto& table = Fprint_table[bin];
    auto it = table.begin();
    
    // 2. Safely iterate and erase
    while (it != table.end()) {
        // it->second is the std::shared_ptr<Fprint_table_entry>
        if (it->second->counter.load(std::memory_order_relaxed) == 0) {
            
            // this calls also the destructor of the fingerprint_entry on the heap
            it = table.erase(it); // table.erase returns the next valid iterator automatically
            
        } else {
            // Move to the next element manually if we didn't delete anything
            ++it;
        }
    }
    
    // force rehash to a minimal size.
    table.rehash(0);    // this invalidates all existing iterators!
    // see sparse_hash.h 
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
        hash_map<kmerAmino_t, fingerprint>::iterator entry=kmer_tableAmino[bin].find(kmer); 
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
     * @param n number of genomes
     * @return mean value of number of k-mers per genome
     */
    static double mean_number_kmers(int n);

    /**
     * @param mu mean value
     * @param n number of genomes
     * @return standard deviation of number of k-mers per genome
     */
    static double stdev_number_kmers(double mu, int n);

    /**
     *  Clears the kmer_table or the kmerAmino_table.
     */
    static void clear_kmer_table();

    /* 
    *  Prints the entire contents of the kmer and fingerprint table to a log file for debugging purposes.
    */
    static void printout_tables(string kmer_filename, string fprint_filename);

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
