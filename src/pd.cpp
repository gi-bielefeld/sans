#include "pd.h"


#include <utility>
#include <iomanip> // für std::setw und std::right
#include <algorithm>


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
pd::pd(const multimap_<double, color_t> split_list, const int n): n(n){

	// determine cycle
	for (auto& it : split_list){
		planar_splits.insert({it.first,it.second});
	}
	pc_tree::PCTree tree = graph::filter_planar(planar_splits,false,color::n);
	pc_tree::IntrusiveList<pc_tree::PCNode> leaves=tree.getLeaves();
	for (auto* leaf : tree.currentLeafOrder()) {
		cycle.push_back(leaf->index()-1);
    }

	//minimum split weight for truncating leaf edges
	auto it = planar_splits.begin();
	min=-1;
	while (it != planar_splits.end()) {
		if(min=-1 or it->first<min) min=it->first;
		it++;
	}
	
// 	cerr<<"pd.cpp, l39: truncate leaf edges to " << 0.5*min << endl;
	
	//PD value for a pair of taxa on the cycle (not for an interval!)
	//needed for the DP algo in pd::pd_set
	vector<vector<double>> x(n,vector<double>(n));
	pd_pair_vals=x;
	for (int i=0;i<n;i++){
		for (int j=i+1;j<n;j++){
			vector<int> taxa = {cycle[i], cycle[j]};
			pd_pair_vals[i][j]=pd_value(taxa);
		}
	}
}


/*
 * Simple getter method.
 */
multimap_<double, color_t> pd::get_planar_splits(){
	return planar_splits;
}


/**
* Compute PD score, i.e., total weight of all splits that separate a given set of taxa.
* Trivial splits (leaf edges) are fixed to a minimum weight.
* 
* @param taxa list of taxa to compute score of, IDs w.r.t. denom_names (not pos in cycle)
* @return pd value of above list of taxa
* 
*/
double pd::pd_value(const vector<int>& taxa) {
	double val=0.0;
	if (taxa.size()<=1) return val;
	
	// represent taxa as colors
	color_t tax_col = 0;
	for (int t : taxa) {
		tax_col.set(t);
		// Trivial splits (leaf edges) are fixed to a minimum weight and ignored below.
		val+=0.5*min;
	}
	
	// Iterating over all splits and check if disjoint with taxa
	// Trivial splits (leaf edges) are excluded here / considered above
	auto it = planar_splits.begin();
	while (it != planar_splits.end()) {

		double weight = it->first;
		color_t colors = it->second;
		it++;
				
		//separating split? (and non-trivial)
		if( !color::is_singleton(colors) && ((tax_col & colors) != tax_col) && ((tax_col & colors) != 0)){
			val+=weight;
		}
	}

	return val;
}



/**
* Determine set of k taxa of maximum pd value.
* 
* DP algo in DOI:10.1093/sysbio/syp058
* 
* @param k size of set
* @return list of IDs (w.r.t. denom_names) of k taxa
* 
*/
vector<int> pd::pd_set(const int k) {

	
	// DP according to DOI:10.1093/sysbio/syp058
	
	// L[k][{u,v}]
	std::vector<std::unordered_map<std::pair<int,int>,double,PairHash>> L(k+1);
 	std::vector<std::unordered_map<std::pair<int,int>,int,PairHash>> alpha(k+1);
	
	//k=2
	for (int u=0; u<n; u++){
		for (int v=u+1; v<n; v++){
			L[2][{u,v}]=pd_pair_vals[u][v];//pre-computed
		}
	}
	
	//larger k
	for (int i=3;i<=k;i++) {
		for (int u=0; u<=n-i+1; u++){
			for (int v=u+i-1; v<n; v++){
				double max=0;
				int arg_max=-1;
				for (int s=u+1;s<v;s++) {
					double l=L[i-1][{u,s}];
					double d=pd_pair_vals[s][v]; //pre-computed
					if(l+d>max){
						max=l+d;
						arg_max=s;
					}
				}
				L[i][{u,v}]=max;
				alpha[i][{u,v}]=arg_max;
			}
		}
	}

	//final optimum
	//max of L[k](u,v)+d(u,v)
	double max=0;
	int max_u=-1;
	int max_v=-1;
	for (int u=0; u<=n-k+1; u++){
		for (int v=u+k-1; v<n; v++){
			double l_k=L[k][{u,v}];
			double d=L[2][{u,v}];
			if(l_k+d>max){
				max=l_k+d;
				max_u=u;
				max_v=v;
			}
		}
	}
	//PD score: 
	//max=max/2;
	
	//trace back
	int s;
	vector<int> max_set;
	max_set.push_back(cycle[max_v]);
	for (int i=k;i>2;i--) {
		s=alpha[i][{max_u,max_v}];
		max_set.push_back(cycle[s]);
		max_v=s;
	}
	max_set.push_back(cycle[max_u]);
	
	return max_set;
}

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
double pd::pd_value_lookup(int l, int r){
	//stored in the table already?
	auto it=pd_cycle_vals.find({l,r});
	// Yes: look up
	if(it!=pd_cycle_vals.end()){
		return it->second;
	} else { //No: compute from scratch
		//collect taxa
		vector<int> taxa;
		for(int i=l%n;i!=r;i=(i+1)%n){
			taxa.push_back(cycle[i]);
		}
		//compute, store, return
		double v = pd_value(taxa);
		pd_cycle_vals[{l,r}]=v;
		return v;
	}
}

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
void pd::partition_dp(int start, vector<int> pos, double& local_min_pd, vector<int>& local_min_seps){
	
	int k=pos.size();
	vector<unordered_map<int, double>> min_tab(k);
	vector<unordered_map<int, int>> argmin_tab(k);
	
	//first column
	unordered_map<int, double> c_min_tab;
	unordered_map<int, int> c_argmin_tab;
	for (int sep=pos[0]+1;sep<=pos[1];sep++){
		c_min_tab[sep]=pd_value_lookup(start,sep);
		c_argmin_tab[sep]=start;
	}
	min_tab.push_back(c_min_tab);
	argmin_tab.push_back(c_argmin_tab);
		
	//other columns
	unordered_map<int, double> prev_col = c_min_tab;
	for (int col=1;col<k-1;col++){
		c_min_tab.clear();
		c_argmin_tab.clear();

		//rows = different seperation points
		for (int sep=pos[col]+1;sep<=pos[col+1];sep++){
			double pd_min=-1;
			int pd_argmin=-1;
			// minimize over all previous seperation points
			for (auto it = prev_col.begin(); it != prev_col.end(); ++it) {
				int prev_sep = it->first;
				double prev_pd = it->second;
				double additional_pd =  pd_value_lookup(prev_sep,sep);
				double sum = prev_pd + additional_pd;
				if (pd_min==-1 or sum < pd_min){ // new min
					pd_min=sum;
					pd_argmin=prev_sep;
				}
			}
			// extend current column
			c_min_tab[sep]=pd_min;
			c_argmin_tab[sep]=pd_argmin;
		}
		// done with this column
		min_tab[col]=c_min_tab;
		argmin_tab[col]=c_argmin_tab;
		prev_col = c_min_tab;

	}
	
	// closing the cycle
	double pd_min=-1;
	int pd_argmin=-1;
	// minimize over all previous seperation points
	for (auto it = prev_col.begin(); it != prev_col.end(); ++it) {
		int prev_sep = it->first;
		double prev_pd = it->second;
		double additional_pd =  pd_value_lookup(prev_sep,start);
		double sum = prev_pd + additional_pd;
		if (pd_min==-1 or sum < pd_min){ // new min
			pd_min=sum;
			pd_argmin=prev_sep;
		}
	}
	
	// compose optimum (backtracing from sep to start)
 	local_min_pd=pd_min;
	int sep=pd_argmin;
	local_min_seps[k-1]=sep;
	for(int i=k-2;i>0;i--){
		sep=argmin_tab[i][sep];
		local_min_seps[i]=sep;
	}
	local_min_seps[0]=start;
}



/**
* Partition the set of all taxa into subsets containing one representative/seed taxon each 
* such as to minimize the sum of PD values of all partitions.
* 
* @param representatives list of represenative/seed taxa (w.r.t. denom_names)
* @return vector assigning a partition ID to each taxon: result[tax_id]=part_id
* 
*/
vector<int> pd::partition(vector<int> representatives){
	return partition(representatives, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
}


/**
* Same as above storing some cluster statistics.
* 
* @param representatives list of represenative/seed taxa (w.r.t. denom_names)
* @param pd varibale to store result: optimal total PD value (sum over all partitions)
* @param min_pd varibale to store result: minimum PD value among partitions
* @param max_pd varibale to store result: maximum PD value among partitions
* @param min_pd_normalized varibale to store result of min_pd_normalized
* @param max_pd_normalized varibale to store result of max_pd_normalized
* @param intra_cluster varibale to store result of intra_cluster
* @param inter_cluster varibale to store result of inter_cluster
* @return vector assigning a partition ID to each taxon: result[tax_id]=part_id
* 
*/
vector<int> pd::partition(vector<int> representatives, double* pd, double* min_pd, double* max_pd, double* min_pd_normalized, double* max_pd_normalized, double* intra_cluster, double* inter_cluster){
	
	int k=representatives.size();
	
	// where in the cycle is taxon i?
	vector<int> inv_cycle(n);
	for (int i=0;i<n;i++){
		inv_cycle[cycle[i]]=i;
	}

	// where in the cycle are the representatives?
	vector<int> pos(k);
	for (int i=0;i<k;i++){
		pos[i]=inv_cycle[representatives[i]];
	}

	// sorted positions
	std::sort(pos.begin(), pos.end());
	
	// call DP algo for any possible starting point
	double global_min_pd=-1;
	vector<int> global_min_seps(k);
	// consider interval that contains the end/start of cycle: (rep_k,rep_0]
	for (int start=(pos[k-1]+1)%n; start!=(pos[0]+1)%n; start=(start+1)%n) {
		double local_min_pd=-1;
		vector<int> local_min_seps(k);
		partition_dp(start,pos,local_min_pd,local_min_seps);
		if (global_min_pd==-1 or local_min_pd < global_min_pd) {
                global_min_pd=local_min_pd;
                global_min_seps=local_min_seps;
		}
	}
	
	//compose mapping tax_id->cluster_id
	vector<int> mapping(n);
	for(int cluster=0;cluster<k;cluster++){
		int sep=global_min_seps[cluster];
		for(int t=sep%n;t!=global_min_seps[(cluster+1)%k];t=(t+1)%n){
			mapping[cycle[t]]=cluster;
		}
	}
	
	//statistics if asked for
	if(pd){*pd=global_min_pd;}
	if(min_pd or max_pd){
		if(min_pd){*min_pd=-1;}
		if(max_pd){*max_pd=0;}
		for(int cluster=0;cluster<k;cluster++){
			int sep=global_min_seps[cluster];
			double s=pd_value_lookup(sep,global_min_seps[(cluster+1)%k]);
			if(min_pd and (*min_pd==-1 or s<*min_pd)){*min_pd=s;}
			if(max_pd and s>*max_pd){*max_pd=s;}
		}
	}
	if(min_pd_normalized){*min_pd_normalized=pd::min_pd_normalized(global_min_seps);}
	if(max_pd_normalized){*max_pd_normalized=pd::max_pd_normalized(global_min_seps);}
	if(intra_cluster){*intra_cluster=pd::intra_cluster(global_min_seps);}
	if(inter_cluster){*inter_cluster=pd::inter_cluster(global_min_seps);}
	
	return mapping;
	
}

/**
* Minimum PD value among all subsets normalized by subset size.
* 
* @param partition_boundaries the left boundaries of a partitioning
*/
double pd::min_pd_normalized(vector<int> partition_boundaries){
	int k=partition_boundaries.size();
	double min=-1;
	//maximize over all clusters
	for(int cluster=0;cluster<k;cluster++){
		int curr=partition_boundaries[cluster];
		int next=partition_boundaries[(cluster+1)%k];
		//PD value for cluster
		double s=pd_value_lookup(curr,next);
		//cluster size
		int size=(next>curr)?(next-curr):(next+n-curr);
		//normalize
		s/=size;
		if(min==-1 or s<min){min=s;}
	}
	return min;
}

/**
* Maximum PD value among all subsets normalized by subset size.
* (Might be interpreted as an intra-cluster distance.)
* 
* @param partition_boundaries the left boundaries of a partitioning
*/
double pd::max_pd_normalized(vector<int> partition_boundaries){
	int k=partition_boundaries.size();
	double max=0;
	//maximize over all clusters
	for(int cluster=0;cluster<k;cluster++){
		int curr=partition_boundaries[cluster];
		int next=partition_boundaries[(cluster+1)%k];
		//PD value for cluster
		double s=pd_value_lookup(curr,next);
		//cluster size
		int size=(next>curr)?(next-curr):(next+n-curr);
		//normalize
		s/=size;
		if(s>max){max=s;}
	}
	return max;
}


/**
* Maximum pairwise PD value within any subset.
* 
* @param partition_boundaries the left boundaries of a partitioning
*/
double pd::intra_cluster(vector<int> partition_boundaries){
	int k=partition_boundaries.size();
	double max=0;
	int argmax_i, argmax_j;
	//maximize over all clusters
	for(int cluster=0;cluster<k;cluster++){
		int curr=partition_boundaries[cluster];
		int next=partition_boundaries[(cluster+1)%k];
		//maximize over all pairs
		for(int i=curr; (i+1)%n!=next; i=(i+1)%n){
			for(int j=(i+1)%n; j!=next; j=(j+1)%n){
				//PD value (considered as "distance" here)
				int l = std::min(i,j);
				int r = std::max(i,j) ;
				double s=pd_pair_vals[l][r];
				//max?
				if(s>max){
					max=s;
					argmax_i=l;
					argmax_j=r;
				}
			}
		}
	}
	return max;
}



/**
* Minimum pairwise PD value between any subsets.
* 
* @param partition_boundaries the left boundaries of a partitioning
*/
double pd::inter_cluster(vector<int> partition_boundaries){
	int k=partition_boundaries.size();
	double min=-1;
	//minimize over all neighboring clusters
	for(int cluster=0;cluster<k;cluster++){
		//PD of gap between clusters
		int sep=partition_boundaries[cluster];
		int l = std::min(sep,(sep-1));
		int r = std::max(sep,(sep-1));
		if(l<0){l=0;r=n-1;}
		double s=pd_pair_vals[l][r];
		if(min==-1 or s<min){
			min=s;
		}
	}	
	return min;
}
