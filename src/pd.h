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

   /**
    * Convert a given list of partition boundaries (left boundaries w.r.t. cycle) into a mapping assigning a partition ID to each taxon: result[tax_id]=part_id
    * @param seps
    * @return
    */
   vector<int> seps2map(vector<int>& seps);
   
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
   //last computes partitioning
   vector<int> partition_boundaries;
   
   
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
     * @return vector assigning a partition ID to each taxon: result[tax_id]=part_id
     * 
     */
    vector<int> partition(vector<int> representatives);

    /**
     * Split the partition with largest PD score into two subsets with minimum sum of PD scores.
     * Input (so to say) is the previously computed partitioning.
     * If no partitioning has been computed yet, an initial bipartition is computed.
     *
     * @return vector assigning a partition ID to each taxon: result[tax_id]=part_id
     *
     */
    vector<int> greedily_split();

    /**
     * Determine cluster statistics of previously computed partitioning.
     * 
     * @param pd varibale to store result: optimal total PD value (sum over all partitions)
     * @param pd_list varibale to store result: PD value per cluster
     * @param min_pd varibale to store result: minimum PD value among partitions
     * @param max_pd varibale to store result: maximum PD value among partitions
     * @param min_pd_normalized varibale to store result of min_pd_normalized
     * @param max_pd_normalized varibale to store result of max_pd_normalized
     * @param intra_cluster variable to store result of intra_cluster
     * @param inter_cluster variable to store result of inter_cluster
     * 
     */
    void partition_statistics(double* pd, vector<double>* pd_list, double* min_pd, double* max_pd, double* min_pd_normalized, double* max_pd_normalized, double* intra_cluster, double* inter_cluster);


    /**
     * Minimum PD value among all subsets normalized by subset size.
     * 
     * @param partition_boundaries the left boundaries of a partitioning
     */
    double min_pd_normalized(vector<int> partition_boundaries);

    /**
     * Maximum PD value among all subsets normalized by subset size.
     * (Might be interpreted as an intra-cluster distance.)
     * 
     * @param partition_boundaries the left boundaries of a partitioning
     */
    double max_pd_normalized(vector<int> partition_boundaries);


    /**
     * Maximum pairwise PD value within any subset.
     * 
     * @param partition_boundaries the left boundaries of a partitioning
     */
    double intra_cluster(vector<int> partition_boundaries);



    /**
     * Minimum pairwise PD value between any subsets.
     * 
     * @param partition_boundaries the left boundaries of a partitioning
     */
    double inter_cluster(vector<int> partition_boundaries);

    
    
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
