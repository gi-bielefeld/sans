#include "tree.h"
#include "ansi.h"

/*
 * This class manages the split/tree filters and Newick output.
 */
uint64_t                             tree::t;        // max. size of the splits list
multimap<double, color_t, greater<>> tree::splits;   // list collecting all the splits ordered by weight
vector<vector<color_t>>              tree::forest;   // collection of (weakly) compatible splits trees

/**
 * This is a tree structure used for generating a Newick string.
 */
struct tree::node {
    color_t taxa;    // color_t coding all taxa beneath this node
    double weight;    // weight of the edge going to this node
    vector<node*> subsets;    // list of subsets
};

/**
 * This function initializes the max. size of the splits list.
 *
 * @param top_size list size
 */
void tree::init(const uint64_t& top_size) {
    t = top_size;
}

/**
 * This function adds a single split (weight and colors) to the output list.
 *
 * @param weight split weight
 * @param color split colors
 */
void tree::insert_split(const double& weight, const color_t& color) {
    splits.emplace(weight, color);    // insert it at the correct position ordered by weight
    if (splits.size() > t) {      // if the splits list exceeds its limit, erase the last entry
        splits.erase(--splits.end());
    }
}

/**
 * This function builds/refines/prints trees and generates a Newick string.
 *
 * @param color_name function to map a color bit to a file name
 * @return Newick string
 */
string tree::build_string(const function<string(const size1N_t&)>& color_name) {
    string str;
    for (auto& tree : forest) {
        node* root = build_tree(tree);
        str += print_tree(root, color_name);
        str += ";\n";
    }
    return str;
}

/**
 * This function returns a tree structure generated from a list of color sets.
 *
 * @param color_sets list of color sets
 * @return tree structure (struct node)
 */
tree::node* tree::build_tree(vector<color_t>& color_sets) {
    vector<node*> subsets = {};
    color_t root_taxa = 0b0u;

    for (size1N_t i = 0; i < color::n; ++i) { // initialize set of trivial splits
        color_t leaf_taxa = 0b0u;
        color::set(root_taxa, i);
        color::set(leaf_taxa, i);

        double weight = 0; // get weight of split
        for (auto& it : splits) {
            if (leaf_taxa == it.second) {
                weight = it.first;
                break;
            }
        }
        node* leaf = new node { leaf_taxa, weight, {} };
        subsets.emplace_back(leaf);
    }
    node* root = new node { root_taxa, 0, subsets };

    for (color_t& split : color_sets) { // split if possible
        if (!refine_tree(root, split, root_taxa)) {
            $err << "Error: splits are incompatible" << _end$;
            exit(1);
        }
    }
    return root;
}

/**
 * This function recursively refines a set/tree structure by a given split.
 *
 * @param root root node of the current (sub-)set/tree structure
 * @param split color set to refine by
 * @param all_taxa color set of the top node of the set/tree structure
 * @return whether the split is compatible with the set/tree structure
 */
bool tree::refine_tree(node* root, color_t& split, const color_t& all_taxa) {
    // possible cases:
    //  - split size < 2: nothing has to be done
    //  - split equals one subset -> warning: split twice
    //  - split is fully contained in one subset -> recurse
    //  - inverse split ... (i.e. split covers one subset partially) -> recurse with inverse
    //  - split covers several subsets completely -> introduce new split

    if (color::count(split) < 2 || color::count(all_taxa) - color::count(split) < 2)
        return true;

    vector<node*>* current_subsets = &root->subsets;
    vector<node*> fully_covered_subsets = {};
    node* partially_covered_subset = nullptr;

    for (node* subset : *current_subsets) {
        if (split == subset->taxa)
            return true;
        else if ((split & subset->taxa) == split) // split.issubset(subtaxa)
            return refine_tree(subset, split, all_taxa);
        else if ((split & subset->taxa) == subset->taxa) // subtaxa.issubset(split)
            fully_covered_subsets.emplace_back(subset);
        else if ((split & subset->taxa) != 0b0u) { // subtaxa.intersects(split)
            if (partially_covered_subset != nullptr)
                return false; // there cannot be more than one
            else
                partially_covered_subset = subset;
        }
    }

    if (partially_covered_subset != nullptr) {
        if (fully_covered_subsets.size() == current_subsets->size()-1) {
            color_t inverse_split = split; // recurse into this subset with inverse split
            color::complement(inverse_split);
            if ((inverse_split & partially_covered_subset->taxa) == inverse_split)
                return refine_tree(partially_covered_subset, inverse_split, all_taxa);
            else
                return false;
        } else {
            return false;
        }
    }
    else if (fully_covered_subsets.size() > 1) {
        color_t combined_taxa = 0b0u; // introduce new split
        for (node* subset : fully_covered_subsets)
            combined_taxa |= subset->taxa;

        double weight = 0; // get weight of split
        for (auto& it : splits) {
            if (split == it.second) {
                weight = it.first;
                break;
            }
        }
        node* combined_set = new node { combined_taxa, weight, fully_covered_subsets };
        for (node* subset : fully_covered_subsets) // remove old sets
            current_subsets->erase(remove(current_subsets->begin(), current_subsets->end(), subset), current_subsets->end());
        current_subsets->emplace_back(combined_set); // add new set
        return true;
    }
    else {
        $err << "Error: one fully covered subset and nothing else!?" << _end$;
        exit(1);
    }
}

/**
 * This function returns a Newick string generated from a given tree structure.
 *
 * @param root root node of the current (sub-)set/tree structure
 * @param color_name function to map a color bit to a file name
 * @return Newick string
 */
string tree::print_tree(const node* root, const function<string(const size1N_t&)>& color_name) {
    if (root->subsets.empty()) {    // leaf set
        if (color::count(root->taxa) == 0) {
            $err << "Error: child with no taxon!?" << _end$;
            exit(1);
        } else if (color::count(root->taxa) == 1) {
            return color_name(color::index(root->taxa)) + ":" + to_string(root->weight);
        } else {
            $err << "Error: child with more than one taxon!?" << _end$;
            exit(1);
        }
    } else {
        string str = "(";
        for (node* subset : root->subsets) {
            str += print_tree(subset, color_name);
            if (subset != root->subsets.back())
                str += ",";
        }
        str += "):" + to_string(root->weight);
        return str;
    }
}

/**
 * This function filters a greedy maximum weight tree compatible subset.
 *
 * @param verbose print progress
 */
void tree::filter_strict(const bool& verbose) {
    forest = vector<vector<color_t>>(1);    // create a set for compatible splits
    auto& tree = forest[0];
    auto it = splits.begin();

    uint64_t cur = 0, prog = 0, next;
    uint64_t max = splits.size();
loop:
    while (it != splits.end()) {
        if (verbose) {
            next = 100 * cur / max;
             if (prog < next)  $log $_ << "Filtering splits... " << next << "%" << $;
            prog = next; cur++;
        }
        if (test_strict(it->second, tree)) {
            tree.emplace_back(it->second);
            ++it; goto loop;    // if compatible, add the new split to the set
        }
        it = splits.erase(it);    // otherwise, remove split
    }
}

/**
 * This function filters a greedy maximum weight weakly compatible subset.
 *
 * @param verbose print progress
 */
void tree::filter_weakly(const bool& verbose) {
    forest = vector<vector<color_t>>(1);    // create a set for compatible splits
    auto& network = forest[0];
    auto it = splits.begin();

    uint64_t cur = 0, prog = 0, next;
    uint64_t max = splits.size();
loop:
    while (it != splits.end()) {
        if (verbose) {
            next = 100 * (cur * sqrt(cur)) / (max * sqrt(max));
             if (prog < next)  $log $_ << "Filtering splits... " << next << "%" << $;
            prog = next; cur++;
        }
        if (test_weakly(it->second, network)) {
            network.emplace_back(it->second);
            ++it; goto loop;    // if compatible, add the new split to the set
        }
        it = splits.erase(it);    // otherwise, remove split
    }
}

/**
 * This function filters a greedy maximum weight n-tree compatible subset.
 *
 * @param n number of trees
 * @param verbose print progress
 */
void tree::filter_n_tree(const uint64_t& n, const bool& verbose) {
    forest = vector<vector<color_t>>(n);    // create a set for compatible splits
    auto it = splits.begin();

    uint64_t cur = 0, prog = 0, next;
    uint64_t max = splits.size();
loop:
    while (it != splits.end()) {
        if (verbose) {
            next = 100 * cur / max;
             if (prog < next)  $log $_ << "Filtering splits... " << next << "%" << $;
            prog = next; cur++;
        }
        for (auto& tree : forest)
            if (test_strict(it->second, tree)) {
                tree.emplace_back(it->second);
                ++it; goto loop;    // if compatible, add the new split to the set
            }
        it = splits.erase(it);    // otherwise, remove split
    }
}

/**
 * This function tests if a split is compatible with an existing set of splits.
 *
 * @param color new split
 * @param color_set set of splits
 * @return true, if compatible
 */
bool tree::test_strict(const color_t& color, const vector<color_t>& color_set) {
    for (auto& elem : color_set) {
        if (!color::is_compatible(elem, color)) {
            return false;    // compare to each split in the set
        }
    }
    return true;
}

/**
 * This function tests if a split is weakly compatible with an existing set of splits.
 *
 * @param color new split
 * @param color_set set of splits
 * @return true, if weakly compatible
 */
bool tree::test_weakly(const color_t& color, const vector<color_t>& color_set) {
    for (auto& elem1 : color_set) {
        if (!color::is_compatible(elem1, color)) {
            for (auto& elem2 : color_set) {
                if (!color::is_weakly_compatible(elem1, elem2, color)) {
                    return false;    // compare to each pair of splits in the set
                }
            }
        }
    }
    return true;
}
