#ifndef BIFROST_TOPSPLITS_H
#define BIFROST_TOPSPLITS_H

#include "TopTrie.h"

using namespace std;

void searchGraph(ColoredCDBG<> &graph, const SANS_opt& opt);

void addColors(TopTrie *trie, ColoredCDBG<> &graph, const SANS_opt& opt,
               UnitigColors::const_iterator it, UnitigColors::const_iterator end, unsigned int value);

void printSplits(TopTrie *trie, ostream &out);

#endif //BIFROST_TOPSPLITS_H
