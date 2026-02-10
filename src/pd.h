using namespace std;

#include <iostream>

#include "graph.h"

class pd {

 private:


 public:

    /**
     * This is a template.
     */
    //static size1N_t n;

    /**
     * This is a template
     *
     * @param number color number
     */
    //static void init(const size1N_t& number);


    /**
     * Compute PD score, i.e., total weight of all splits that separate a given set of taxa
     * 
     * @param split_list split list with weights and colors
     * @param taxa list of taxa
     * 
     */
    static double pd_value(const multimap_<double, color_t> split_list, const vector<int>& taxa);
    
 protected:

};
