#include <algorithm>
#include <functional>
#include <cmath>
#include <vector>
#include <set>

#include "color.h"
using namespace std;

template <typename K, typename V>
  struct compare: public function<bool(pair<K,V>, pair<K,V>)> {
      constexpr bool operator()(const pair<K,V>& x, const pair<K,V>& y) const noexcept {
          return x.first > y.first || x.first == y.first && x.second < y.second;
      }
  };
template <typename K, typename V>
  using multimap_ = set<pair<K,V>, compare<K,V>>;

/**
 * This class manages the split/tree filters and Newick output.
 */
class tree {

 private:

    /**
     * This is a tree structure used for generating a Newick string.
     */
    struct node;

 public:

    /**
     * This is the max. size of the splits list.
     */
    static uint64_t t;

    /**
     * This is a list collecting all the splits ordered by weight.
     */
    static multimap_<double, color_t> splits;

    /**
     * This function initializes the max. size of the splits list.
     *
     * @param top_size list size
     */
    static void init(const uint64_t& top_size);

    /**
     * This function adds a single split (weight and colors) to the output list.
     *
     * @param weight split weight
     * @param color split colors
     */
    static void insert_split(const double& weight, const color_t& color);

    /**
     * This function builds/refines/prints trees and generates a Newick string.
     *
     * @return Newick string
     */
    static string build_string();

    /**
     * This function filters a greedy maximum weight tree compatible subset.
     *
     * @param verbose print progress
     */
    static void filter_strict(const bool& verbose);

    /**
     * This function filters a greedy maximum weight weakly compatible subset.
     *
     * @param verbose print progress
     */
    static void filter_weakly(const bool& verbose);

    /**
     * This function filters a greedy maximum weight n-tree compatible subset.
     *
     * @param n number of trees
     * @param verbose print progress
     */
    static void filter_n_tree(const uint64_t& n, const bool& verbose);

    /**
     * This function maps a color position to a file name.
     *
     * @param index color position
     */
    static function<string(const size1N_t&)> color_name;

    /**
     * This function maps a file name to a color position.
     *
     * @param name file name
     */
    static function<size1N_t(const string&)> color_index;

 protected:

    /**
     * This is a collection of (weakly) compatible splits trees.
     */
    static vector<vector<color_t>> forest;

    /**
     * This function returns a tree structure generated from a list of color sets.
     *
     * @param color_sets list of color sets
     * @return tree structure (struct node)
     */
    static node* build_tree(vector<color_t>& color_sets);

    /**
     * This function recursively refines a set/tree structure by a given split.
     *
     * @param root root node of the current (sub-)set/tree structure
     * @param split color set to refine by
     * @param all_taxa color set of the top node of the set/tree structure
     * @return whether the split is compatible with the set/tree structure
     */
    static bool refine_tree(node* root, color_t& split, const color_t& all_taxa);

    /**
     * This function returns a Newick string generated from a given tree structure.
     *
     * @param root root node of the current (sub-)set/tree structure
     * @return Newick string
     */
    static string print_tree(const node* root);

    /**
     * This function tests if a split is compatible with an existing set of splits.
     *
     * @param color new split
     * @param color_set set of splits
     * @return true, if compatible
     */
    static bool test_strict(const color_t& color, const vector<color_t>& color_set);

    /**
     * This function tests if a split is weakly compatible with an existing set of splits.
     *
     * @param color new split
     * @param color_set set of splits
     * @return true, if weakly compatible
     */
    static bool test_weakly(const color_t& color, const vector<color_t>& color_set);

};
