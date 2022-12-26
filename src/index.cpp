#include "index.h"
#include "ansi.h"

file::file(const string& path, const type& type, const mode& mode) {
    _path = path; _type = type; _mode = mode;

    if (_path.empty())       $err << "Error: file path not specified"            << _end$$;
    if (_type == type::NONE) $err << "Error: file type not specified: " << _path << _end$$;
    if (_mode == mode::NONE) $err << "Error: file mode not specified: " << _path << _end$$;

    if (_type == type::GRAPH_FILE) {
        #ifndef useBF
            $err << "Error: requires compiler flag -DuseBF, please see makefile" << _end$$;
        #endif
    }
}
