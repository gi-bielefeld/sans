#include <iostream>
using namespace std;

enum class type {
    NONE,
    INPUT_FILE,  // list of sequence files, one per line
    GRAPH_FILE,  // Bifrost graph
    INDEX_OR_SPLITS_FILE,
        INDEX_FILE,   // k-mer index, e.g. counts table
        SPLITS_FILE,  // list of splits, sorted by weight
    NEWICK_FILE,  // tree topology
};

enum class mode {
    NONE,
    READ,
    WRITE,
};

class file {
 private:
    string _path;
    type   _type;
    mode   _mode;

 public:
    file() { _path.clear(); _type = type::NONE; _mode = mode::NONE; }
   ~file() { _path.clear(); _type = type::NONE; _mode = mode::NONE; }

    file(const string& path, const type& type, const mode& mode);

    operator  string()                  { return _path; }
    inline    bool empty()              { return _path.empty(); }
    constexpr bool is(const type& type) { return _type == type; }
    constexpr bool is(const mode& mode) { return _mode == mode; }

    friend ostream& operator<<(ostream& os, const file& file) { return os << file._path; }
};
