#include "index.h"
#include "ansi.h"

/**
 * This class provides an interface for handling different file types.
 */
file::file(const string& path, const type& type, const mode& mode) {
    _path = path; _type = type; _mode = mode; properties.clear();

    if (_path.empty())       $err << "Error: file path not specified"            << _end$$;
    if (_type == type::NONE) $err << "Error: file type not specified: " << _path << _end$$;
    if (_mode == mode::NONE) $err << "Error: file mode not specified: " << _path << _end$$;

    if (_type == type::GRAPH_FILE) {
        #ifndef useBF
            $err << "Error: requires compiler flag -DuseBF, please see makefile" << _end$$;
        #endif
        if (_mode == mode::WRITE) $err << "Error: file mode not supported: " << _path << _end$$;

        size_t pos; string ext;
         if ((pos = _path.rfind(".gfa.gz"))     != string::npos && (ext = _path.substr(pos)) == ".gfa.gz"
         ||  (pos = _path.rfind(".gfa"))        != string::npos && (ext = _path.substr(pos)) == ".gfa"
         ||  (pos = _path.rfind(".color.bfg"))  != string::npos && (ext = _path.substr(pos)) == ".color.bfg"
         ||  (pos = _path.rfind(".bfg_colors")) != string::npos && (ext = _path.substr(pos)) == ".bfg_colors");
        string prefix = _path.substr(0, pos);

         if (ext != ".gfa" && ifstream(properties["graph_fn"] = prefix + ".gfa.gz").good()
         ||  ext != ".gfa.gz" && ifstream(properties["graph_fn"] = prefix + ".gfa").good());
        else $err << "Error: could not find Bifrost graph file (.gfa.gz)" << _end$$;

         if (ext != ".bfg_colors" && ifstream(properties["color_fn"] = prefix + ".color.bfg").good()
         ||  ext != ".color.bfg" && ifstream(properties["color_fn"] = prefix + ".bfg_colors").good());
        else $err << "Error: could not find Bifrost color file (.color.bfg)" << _end$$;
    }
}
