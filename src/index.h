#include <iostream>
#include <fstream>
#include <map>
using namespace std;

enum class type {
    NONE,
    INPUT_FILE,  // list of sequence files, one per line
    GRAPH_FILE,  // Bifrost graph, .gfa.gz + .color.bfg
    INDEX_OR_SPLITS_FILE,
        INDEX_FILE,   // k-mer index, e.g. counts table
        SPLITS_FILE,  // list of splits, sorted by weight
    NEWICK_FILE,  // tree topology in Newick format
};

enum class mode {
    NONE,
    READ,
    WRITE,
};

/**
 * This class provides an interface for handling different file types.
 */
class file {
 private:
    string _path;
    type   _type;
    mode   _mode;

 protected:
    map<string, string> properties;

 public:
    file() { _path.clear(); _type = type::NONE; _mode = mode::NONE; properties.clear(); }
   ~file() { _path.clear(); _type = type::NONE; _mode = mode::NONE; properties.clear(); }

    file(const string& path, const type& type, const mode& mode);

    operator  string()                  { return _path; }
    inline    bool empty()              { return _path.empty(); }
    constexpr bool is(const type& type) { return _type == type; }
    constexpr bool is(const mode& mode) { return _mode == mode; }

    inline string&  operator[](const string& key)             { return properties[key]; }
    friend ostream& operator<<(ostream& os, const file& file) { return os << file._path; }
};
