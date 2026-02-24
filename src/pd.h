using namespace std;

#include <iostream>

#include "graph.h"

class pd {
   

 protected:
   
   // Used to construct maps having pairs of ints as keys.
   struct PairHash {
      std::size_t operator()(const std::pair<int,int>& p) const noexcept {
         return std::hash<int>{}(p.first) ^ (std::hash<int>{}(p.second) << 1);
      }
   };


   
 private:
 
   //number of taxa
   int n;
   //given split set is filtered for planarity
   multimap_<double, color_t> planar_splits;
   //cyclic order of taxa along a planar split graph
   vector<int> cycle;
   //table to store PD values for pairs of taxa (precomputation for function pd_set)
   vector<vector<double>> pd_pair_vals;
   //table to store PD values for intervals of the cycle
   unordered_map<std::pair<int,int>,double,PairHash> pd_cycle_vals;
   //minimum observed weight of leaf edge; used to truncate leaf edges when calculating PD values
   double min;

   
   
 public:

    /**
     * Simple getter method.
     */
    multimap_<double, color_t> get_planar_splits();

    /**
    * Contructor
    * A given set of splits is copied and filtered to a planar split system.
    * The resulting cycle (cyclic taxa order along the planar graph) is stored as well.
    * The minimum split weight is determined for later use in pd::pd_value. 
    * 
    * @param split_list split list with weights and colors
    * @param n number of taxa
    * 
    */
   explicit pd(const multimap_<double, color_t> split_list, const int n);


   /**
     * Compute PD score, i.e., total weight of all splits that separate a given set of taxa.
     * Trivial splits (leaf edges) are fixed to a minimum weight.
     * 
     * @param taxa list of taxa to compute score of, IDs w.r.t. denom_names (not pos in cycle)
     * @return pd value of above list of taxa
     * 
     */
    double pd_value(const vector<int>& taxa);
    double pd_value_printing(const vector<int>& taxa);
    
    /**
     * Determine set of k taxa of maximum pd value.
     * 
     * DP algo in DOI:10.1093/sysbio/syp058
     * 
     * @param k size of set
     * @return list of IDs (w.r.t. denom_names) of k taxa
     * 
     */
    vector<int> pd_set(const int k);
    
    /**
     * Partition the set of all taxa into subsets containing one representative/seed taxon each 
     * such as to minimize the sum of PD values of all partitions.
     * 
     * @param representatives list of represenative/seed taxa (w.r.t. denom_names)
     * @param score varibale to store result: optimal total PD value (sum over all partitions)
     * @param min_score varibale to store result: minimum PD value among partitions
     * @param max_score varibale to store result: maximum PD value among partitions
     * @return vector assigning a partition ID to each taxon: result[tax_id]=part_id
     * 
     */
    vector<int> partition(vector<int> representatives, double& score, double& min_score, double& max_score);


    
 protected:

    /**
     * Get the PD value for an interval [l,r[ of taxa (along the cycle).
     * If the value has been calculated before, it is looked up in pd_cycle_vals.
     * If not, it is calculated and stored.
     * 
     * @param l left interval boundary (inclusive) in the cycle
     * @param r right interval boundary (exclusive)  in the cycle
     * @return PD value of all taxa in the interval [l,r[ of the cycle.
     * 
     */
    double pd_value_lookup(int l, int r);
    
    
    /**
     * Helper function for partition.
     * Find a combination of further partition boundaries such as:
     * - first partition starts on given left boundary,
     * - each partition includes exactly one represenative/seed taxon,
     * - sum of PD values over all partitions is minimal.
     * Dynamic programming.
     * 
     * @param start a coordinate on the cycle that is fixed as the first interval boundary
     * @param pos list of positions of represenatives/seeds along the cycle
     * @param local_min_pd optimal PD value for given start coordinate is stored in this variable
     * @param local_min_seps interval/partition boundaries (left coordinate on cycle each) leading to optimal PD value
     */
    void partition_dp(int start, vector<int> pos, double& local_min_pd, vector<int>& local_min_seps);
    
};
