using namespace std;

#include <iostream>

#include "graph.h"

class pd {
   

 private:
    
  	struct PairHash {
		std::size_t operator()(const std::pair<int,int>& p) const noexcept {
			return std::hash<int>{}(p.first) ^ (std::hash<int>{}(p.second) << 1);
		}
	};
 
    multimap_<double, color_t> planar_splits;
   int n;
   vector<int> cycle;
   vector<vector<double>> pd_pair_vals;
   unordered_map<std::pair<int,int>,double,PairHash> pd_cycle_vals;
   
   double pd_value_lookup(int l, int r);
   void partition_dp(int start, vector<int> pos, double& local_min_pd, vector<int>& local_min_seps);

   
 
   
 public:

   /**
    * Contructor
    * @param split_list split list with weights and colors
    * @param n number of taxa
    * 
    */
   explicit pd(const multimap_<double, color_t> split_list, const int n);


   /**
     * Compute PD score, i.e., total weight of all splits that separate a given set of taxa
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
     * @param k size of set
     * @return list of IDs (w.r.t. denom_names) of k taxa
     * 
     */
    vector<int> pd_set(const int k);
    
    /**
     * Partition the set of all taxa into subsets containing one representative each 
     * such as to minimize the sum of PD values of all partitions.
     * 
     * DP algo in DOI:10.1093/sysbio/syp058
     * 
     * @param result of pd_set
     * @return vector assigning a partition ID to each taxon: result[tax_id]=part_id
     * 
     */
    vector<int> partition(vector<int> representatives, double& score);
    
 protected:

};
