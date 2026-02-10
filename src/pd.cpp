#include "pd.h"







double pd::pd_value(const multimap_<double, color_t> split_list, const vector<int>& taxa) {
	double val=0.0;
	// represent taxa as colors
	color_t tax_col = 0;
	for (int t : taxa) {
		tax_col.set(t);
	}
	// Iterating over all splits and check if disjoint with taxa
	auto it = split_list.begin();
	while (it != split_list.end()) {

		double weight = it->first;
		color_t colors = it->second;
		
		//separating split?
		if( ((tax_col & colors) != tax_col) && ((tax_col & colors) != 0)){
			val+=weight;
		}
		it++;
	}

	return val;
}


