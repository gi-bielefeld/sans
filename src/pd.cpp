#include "pd.h"


#include <utility>
#include <iomanip> // für std::setw und std::right


pd::pd(const multimap_<double, color_t> split_list, const int n): n(n){

	// determine cycle
	for (auto& it : split_list){
		planar_splits.insert({it.first,it.second});
	}
	
	pc_tree::PCTree tree = graph::filter_planar(planar_splits,false,color::n);
	pc_tree::IntrusiveList<pc_tree::PCNode> leaves=tree.getLeaves();
	
	//minimum split weight for truncating leaf edges
	auto it = planar_splits.begin();
	min=-1;
	while (it != planar_splits.end()) {
		if(min=-1 or it->first<min) min=it->first;
		it++;
	}

	for (auto* leaf : tree.currentLeafOrder()) {
		cycle.push_back(leaf->index()-1);
    }

    vector<vector<double>> x(n,vector<double>(n));
	pd_pair_vals=x;
	for (int i=0;i<n;i++){
		for (int j=i+1;j<n;j++){
			vector<int> taxa = {cycle[i], cycle[j]};
			pd_pair_vals[i][j]=pd_value(taxa);
		}
	}

//  std::cerr << "tree.currentLeafOrder(): ";
// 	for(int i : cycle){std::cerr<<i<<" ";}
// 	std::cerr << endl;

}

multimap_<double, color_t> pd::get_planar_splits(){
	return planar_splits;
}

double pd::pd_value(const vector<int>& taxa) {
	double val=0.0;
	if (taxa.size()<=1) return val;
	
	// represent taxa as colors
	color_t tax_col = 0;
	for (int t : taxa) {
		tax_col.set(t);
		val+=min;
	}
	
	// Iterating over all splits and check if disjoint with taxa
	auto it = planar_splits.begin();
	while (it != planar_splits.end()) {

		double weight = it->first;
		color_t colors = it->second;
		it++;
		
		//if (color::is_singleton(colors)) {/*val+=min; */continue;}
		
		//separating split?
		if( !color::is_singleton(colors) && ((tax_col & colors) != tax_col) && ((tax_col & colors) != 0)){
			val+=weight;
		}
	}

	return val;
}



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


// including l, excluding r: [l,r)
// l and r are positions in cycle
double pd::pd_value_lookup(int l, int r){
	auto it=pd_cycle_vals.find({l,r});
	if(it!=pd_cycle_vals.end()){
		return it->second;
	} else {
		vector<int> taxa;
		for(int i=l%n;i!=r;i=(i+1)%n){
			taxa.push_back(cycle[i]);
		}
		double v = pd_value(taxa);
		pd_cycle_vals[{l,r}]=v;
		return v;
	}
}


void pd::partition_dp(int start, vector<int> pos, double& local_min_pd, vector<int>& local_min_seps){
	
	int k=pos.size();
	vector<unordered_map<int, double>> min_tab(k);
	vector<unordered_map<int, int>> argmin_tab(k);
	
// 	for(int i=0;i<k;i++){cerr<<" "<<pos[i];}cerr<<endl;
	
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
// 		cerr << "prevcol " << endl;
// 		for (auto it = prev_col.begin(); it != prev_col.end(); ++it) {
// 			cerr << it->first << "->" << it->second << " ";
// 		}
// 		cerr << endl;

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
// 		cerr << "col " << col << endl;
// 		for (auto it = argmin_tab[col].begin(); it != argmin_tab[col].end(); ++it) {
// 			cerr << it->first << "->" << it->second << " ";
// 		}
// 		cerr << endl;

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
	
// 	for (int i=0;i<k;i++){
// 		cerr << i << ": ";
// 	}
	
	
	// compose optimum (backtracing from sep to start)
 	local_min_pd=pd_min;
	int sep=pd_argmin;
	local_min_seps[k-1]=sep;
	for(int i=k-2;i>0;i--){
// 		cerr << "SEP: " << sep << endl;
// 		for (auto it = argmin_tab[i].begin(); it != argmin_tab[i].end(); ++it) {
// 			cerr << it->first << "->" << it->second << " ";
// 		}
// 		cerr << endl;
		sep=argmin_tab[i][sep];
// 		cerr << "NEW SEP: " << sep << endl;
		local_min_seps[i]=sep;
	}
	local_min_seps[0]=start;
}

vector<int> pd::partition(vector<int> representatives, double& score){
	
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
	for (int start=(pos[k-1]+1)%n; start!=pos[0]; start=(start+1)%n) {
		double local_min_pd=-1;
		vector<int> local_min_seps(k);
		partition_dp(start,pos,local_min_pd,local_min_seps);
		if (global_min_pd==-1 or local_min_pd < global_min_pd) {
                global_min_pd=local_min_pd;
                global_min_seps=local_min_seps;
		}
// 		cerr << start << ": " << local_min_pd << endl;
	}
	
	score=global_min_pd;

// 	for (int i=0;i<k;i++){
// 		cerr << "["<<global_min_seps[i]<<","<<global_min_seps[(i+1)%k]<<"): " << pd_value_lookup(global_min_seps[i],global_min_seps[(i+1)%k]) << endl;
// 	}
// 	cerr<<endl;

	//compose mapping tax_id->cluster_id
	vector<int> mapping(n);
	for(int cluster=0;cluster<k;cluster++){
		int sep=global_min_seps[cluster];
		for(int t=sep%n;t!=global_min_seps[(cluster+1)%k];t=(t+1)%n){
			mapping[cycle[t]]=cluster;
		}
	}
	return mapping;
	
}

