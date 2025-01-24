#include "nexus_color.h"

struct rgb_color {
    int r, g, b;

    bool is_white(){ return ((r == 255) && (g == 255) && (b == 255));}
    bool is_default(){ return ((r == -1) && (g == -1) && (b == -1));} // none color default
    bool is_equal(rgb_color color){
        return ((r == color.r) && (g == color.g) && (b == color.b));
    }

    void set_white(){ r = g = b = 255;}
    void set_black(){ r = g = b = 0;}
    void set_default(){ r = g = b = -1;}

    void print(){ cout << r << " " << g << " " << b;}
};


bool nexus_color::program_in_path(const string& programName) {
    // suppress output of this command depending on system
    #ifdef _WIN32
    #define DEV_NULL "nul"
    #else
    #define DEV_NULL "/dev/null"
    #endif

    string command = "which " + programName + " > " + string(DEV_NULL);
    // Run command and check result
    int result = system(command.c_str());
    // If result is 0, the program is in PATH
    return result == 0;
}


string nexus_color::modify_filename(string& old_name, string front_extension){

    // TODO handle filepath for windows and unix \\ / and such

    size_t lastSlash = old_name.find_last_of("/\\");
    string new_name;
    if (lastSlash != string::npos) {
        string path = old_name.substr(0, lastSlash + 1);
        new_name = old_name.substr(lastSlash + 1);
        new_name = path + front_extension + new_name;
    } else {
        new_name = front_extension + old_name;
    }
    return new_name;
}

string nexus_color::remove_extensions(string& orig_filename) {

    vector<string> extensions_fst = {".fa",".fas",".fastq",".mfasta",".fasta",".fsa",".fna"};
    vector<string> extensions_snd {".gz",".gzip",".zip"};

    string filename = orig_filename;

    // remove folder prefix
    size_t lastSlash = filename.find_last_of("/\\");
    if (lastSlash != string::npos) {
        filename = filename.substr(lastSlash + 1);
    }

    size_t lastDotPos = filename.find_last_of('.'); // last occurrence of the dot
    // removing gz, gzip, zip
    if (lastDotPos != std::string::npos) {
        string ext = filename.substr(lastDotPos); // extract extension
        // check if one of the above
        auto it = find(extensions_snd.begin(), extensions_snd.end(), ext);
        if (it != extensions_snd.end()) {
            filename = filename.substr(0, lastDotPos); // remove extension
        }
    }
    // removing rest
    lastDotPos = filename.find_last_of('.');
    if (lastDotPos != std::string::npos) {
        string ext = filename.substr(lastDotPos); // extract extension
        // check if one of the above
        auto it = find(extensions_fst.begin(), extensions_fst.end(), ext);
        if (it != extensions_fst.end()) {
            filename = filename.substr(0, lastDotPos); // remove extension
        }
    }
    return filename;
}

/**
 * Creates a map of a given tab separated file with two fields (key,value) per line.
 * Also creates another map with the given values as keys for further usage.
 * Here: Key = Taxon, Value = Group
 * @param map The dictionary consistent of the provided (key,value) pairs.
 * @param filename The input file to read from.
 * @param value_map The dictionary with the provided values as keys.
 */
void create_maps(unordered_map<string, string>& map, const string& filename, unordered_map<string, rgb_color>& value_map){

    string line;
    ifstream infile(filename);
    if(!infile.is_open()){
        cerr << "Error: Could not read input file " << filename << endl;
        return;
    }
    while (getline(infile, line)){
        // String stream to split and read line
        istringstream iss(line);
        string key, value; // two fields, tab separated (taxon and group)
        if (getline(iss, key, '\t') && getline(iss, value)) {
            key = nexus_color::remove_extensions(key);

            if(map.find(key) != map.end()){ // taxon has already been assigned to group
                if(map[key] != value){ // the groups are not the same
                    cerr << "Warning: Several groups have been assigned to " << key << endl;
                    cerr << " Using latest assigned group: " << value << endl;
                }
            }
            // Store fields in dictionary
            map[key] = value;
            rgb_color clr; clr.set_default();
            value_map[value] = clr; // default color = -1 (no color)
        }
    }
    infile.close();
}

// Using the golden ratio to create distinct hsv colors and then convert them to rgb
vector<rgb_color> create_colors(int n) { // !not completely my function!
    vector<rgb_color> colors;

    // Golden ratio to ensure visual distinctiveness
    const double goldenRatioConjugate = 0.618033988749895; // used to create different colors
    double hue = 0.01; // probably adjust for better colors

    double saturation = 0.6; // changed later for distinct colors
    double value = 0.9; // = brightness

    for (int i = 0; i < n; ++i) {
        // Convert HSV to RGB

        int hi = static_cast<int>(std::floor(hue * 6.0));
        double f = hue * 6.0 - hi;
        double p = value * (1.0 - saturation);
        double q = value * (1.0 - f * saturation);
        double t = value * (1.0 - (1.0 - f) * saturation);

        int r, g, b;
        switch (hi % 6) {
            case 0: r = value * 255; g = t * 255; b = p * 255; break;
            case 1: r = q * 255; g = value * 255; b = p * 255; break;
            case 2: r = p * 255; g = value * 255; b = t * 255; break;
            case 3: r = p * 255; g = q * 255; b = value * 255; break;
            case 4: r = t * 255; g = p * 255; b = value * 255; break;
            case 5: r = value * 255; g = p * 255; b = q * 255; break;
        }

        rgb_color color;
        color.r = r;
        color.g = g;
        color.b = b;
        colors.push_back(color);

        // Increment hue using golden ratio
        hue += goldenRatioConjugate;
        hue = fmod(hue, 1.0);
        // Gradually decrease saturation and value
        saturation -= 0.2 / double(n);
        saturation = std::max(0.2, saturation);
        //value -= 0.2 / double(n);
        //value = std::max(0.6, value);
    }

    return colors;
}


void reading_grp_clr_file(const string& grp_clr_file, unordered_map<string, rgb_color>& grp_clr_map, int& no_grps){

    no_grps = grp_clr_map.size(); // number of groups

    if(!grp_clr_file.empty()){
        // add colors to grp_clr_map
        string grp_clr;
        ifstream infile(grp_clr_file);
        if(!infile.is_open()){
            cerr << "Error: Could not read input file " << grp_clr_file << endl;
            return;
        }
        while (getline(infile, grp_clr)){
            // String stream to split and read line
            istringstream iss(grp_clr);
            string grp, clr; // two fields, tab separated
            if (getline(iss, grp, '\t') && getline(iss, clr)) {

                istringstream rgb(clr);
                int r, g, b;
                // check if clr is proper rgb value: only ints and in range 0-255
                if(rgb >> r >> g >> b && r >= 0 && r <= 255 && g >= 0 && g <= 255 && b >= 0 && b <= 255){
                    rgb_color color;
                    color.r = r;
                    color.g = g;
                    color.b = b;

                    // check if group already has a color
                    if(!grp_clr_map[grp].is_default() && !grp_clr_map[grp].is_equal(color)){
                        cerr << "Warning: Several colors given for group " << grp << ". Using latest: ";
                        color.print(); cerr << endl;
                    }
                    grp_clr_map[grp] = color;
                } else {
                    cerr << "Warning: Invalid color (rgb value) given: " << grp << "  " << clr << endl;
                }
            }
        }
        infile.close();

    } else { // create colors if they are not given
        vector<rgb_color> colors = create_colors(no_grps);
        int i = 0;
        for(auto& vertex : grp_clr_map){
            vertex.second = colors[i];
            ++i;
        }
    }
}

/* Unused for now TODO saving pdf doesnt work (the created file itself is good) **/
void create_legend_pdf(unordered_map<string, rgb_color>& grp_clr_map){

    const char* temp_file = "./temp_legend.nex";
    ofstream temp_nexus(temp_file);
    if(!temp_nexus.is_open()){
        cerr << "Error creating temporary file  "<< temp_file <<" for legend creation.\n";
        return;
    }

    ostringstream labels;
    string pdf = "./example_data/legend_test.pdf"; // TODO saving doesn't work
    int nvertices = grp_clr_map.size()+1;
    double i = 0;
    rgb_color clr;

    // Nexus format stuff
    temp_nexus << "#nexus\n\nBEGIN Taxa;\nDIMENSIONS ntax=0;\nTAXLABELS\n;\nEND; [Taxa]\n";
    temp_nexus << "\nBEGIN Network;\nDIMENSIONS ntax=0 nvertices=" << nvertices << " nedges=0;\nDRAW to_scale;\nTRANSLATE\n;\n";
    temp_nexus << "VERTICES\n";
    for (pair<string, rgb_color> grp_clr : grp_clr_map) { // adding nodes for legend
        ++i;
        clr = grp_clr.second;
        if(clr.is_white()){
            temp_nexus << i << " 0.0 " << i*0.1 << " w=10 h=10 bg="<<clr.r<<" "<<clr.g<<" "<< clr.b<<",\n";
        } else {
            temp_nexus << i << " 0.0 " << i*0.1 << " w=10 h=10 fg="<<clr.r<<" "<<clr.g<<" "<< clr.b<<" bg="<< clr.r<<" "<<clr.g<<" "<< clr.b<<",\n";
        }
        labels << i << " " << grp_clr.first << " l=9 x=15 y=5 f='Dialog-PLAIN-12',\n";
    }
    // invisible node to keep good zoom
    temp_nexus << ++i << " 0.0 " << (2.0*double(nvertices)*0.1)+2.0 << " w=0 h=0,\n";

    temp_nexus << ";\nVLABELS\n" << labels.str() << ";\nEDGES\n;\nEND; [Network]\n";
    temp_nexus << "BEGIN st_Assumptions;\nuptodate;\n exclude  no missing;\nautolayoutnodelabels;\nEND; [st_Assumptions]";

    nexus_color::open_in_splitstree(temp_file, "", "", true, false, temp_file);
    // TODO weird error messages (: must start with #nexus , tut sie aber doch??!?)

    //cout << "Legend saved in " << pdf << endl;
    //remove(temp_file);
}


void nexus_color::scale_nexus(const string& unopened_nexus_file, bool verbose, bool scale_notification){

    if(verbose && scale_notification) cout << "Scaling nexus" << endl;

    bool matrix = false;
    int max_weight = -1;
    string line;
    string temp_file = "./temp_scaled";

    ifstream plain_nexus(unopened_nexus_file); // nexus input file
    ofstream scaled_nexus(temp_file.c_str());
    if(!plain_nexus.is_open()){
        cerr << "Error opening nexus file "<< unopened_nexus_file <<" to modify with color.\n";
        return;
    }
    if(!scaled_nexus.is_open()){
        cerr << "Error creating scaled temporary nexus file " << temp_file << endl;
        return;
    }

    // reading nexus file
    while (getline(plain_nexus, line)) {
        // end of block
        if(line.find(";") != string::npos){
            matrix = false;
        }
        // read matrix block and scale length/weights
        if(matrix){
            istringstream iss(line);
            string token;
            string split;
            string split_members;
            double weight;

            if (std::getline(iss, token, '\t')) split = token; // extract first field
            if (std::getline(iss, token, '\t')) { // modify weight
                weight = stod(token);
                if(max_weight == -1) max_weight = weight; // first entry is the max weight
                weight /= max_weight; // normalize using first/max weight
            }
            if (std::getline(iss, token)) split_members = token; // extract third field

            // save to write to output
            ostringstream oss;
            oss << split << "\t" << weight << "\t " << split_members;
            line = oss.str();
        }

        // identify matrix block
        if(line.find("MATRIX") != string::npos || line.find("Matrix") != string::npos || line.find("matrix") != string::npos){
            matrix = true;
        }

        scaled_nexus << line << endl; // write to output
    }

    plain_nexus.close();
    scaled_nexus.close();

    //either rename or return name of the scaled nexus file
    std::rename(temp_file.c_str(), unopened_nexus_file.c_str());

    if(scale_notification) {
        cerr << "Attention: For a reliable visualization using SplitsTree, split weights in the Nexus file have been scaled "
            "to the range 0 to 1 by division by the maximum split weight " << max_weight << endl;
    }
}


void nexus_color::open_in_splitstree(const string& nexus_file, const string& pdf, const string& svg, bool verbose, bool update, const string& save_as, const string splitstree_path){
    // = "../splitstree4/SplitsTree"
    if (!program_in_path(splitstree_path)) {
        cerr << splitstree_path << " is not in the PATH." << endl;
        return;
    }
    // temporary file for splitstree commands
    const char* temp = "./temp_splitstree_commands";
    ofstream temp_file(temp);
    if (temp_file.is_open()) {

        // Writing commands for splitstree
        temp_file << "begin SplitsTree;\nEXECUTE FILE=" << nexus_file << endl;
        if(update) temp_file << "UPDATE\n";// add network in SplitsTree graphic
        if(!save_as.empty()) temp_file << "SAVE FILE=" << save_as << " REPLACE=YES\n"; // save the network to the file (only needed if color should be added)
        if(!pdf.empty()) temp_file << "EXPORTGRAPHICS format=PDF file=" << pdf << " REPLACE=yes\n";
        if(!svg.empty()) temp_file << "EXPORTGRAPHICS format=SVG file=" << svg << " REPLACE=yes\n";
        temp_file << "QUIT\nend;";
        temp_file.close();

        // Running splitstree
        if(verbose) cout  << "Running SplitsTree\n";
        /// .../SplitsTree -g -S -c run_splitstree
        string command = splitstree_path + " -g -c " + temp + " > splitstree.log";
        const char* splitstree_command = command.c_str();
        int result = system(splitstree_command); // executing command
        if(result != 0){
            cerr << "Error while running Splitstree " << splitstree_path << endl;
        }
        remove(temp);
    } else {
        cerr << "Unable to create temp command file for splitstree: " << temp << "\nTherefor unable to add network to nexus file" << endl;
    }
}


void nexus_color::color_nexus(const string& nexus_file, const string& tax_grp_file, const string& grp_clr_file){

    bool translate = false;
    bool vertices = false;
    bool vlabels = false;
    int no_vertices = 0; // number of all vertices in network
    int no_grps; // number of groups given
    int size = 10; // size of colored vertex
    int textsize = 12; // size of label text
    double max_x = 0; // max x value of graph nodes
    double max_y = 0; // max y value of graph nodes
    double min_y = std::numeric_limits<double>::max(); // min y value of graph nodes
    double pos_x; // value to shift legend by
    double pos_y; // value to shift legend by
    ostringstream legend; // to save info for legend later
    string line = "";

    unordered_map<string, string> tax_grp_map; // map containing taxname -> group
    unordered_map<string, rgb_color> grp_clr_map; // map conatining group -> color
    unordered_map<int, rgb_color> no_clr_map; // map containing number -> color to color vertices

    // Reading and saving taxa -> grp mapping
    create_maps(tax_grp_map, tax_grp_file, grp_clr_map);

    // (Reading &) Saving grp -> col mapping
    reading_grp_clr_file(grp_clr_file, grp_clr_map, no_grps);

    // temporary naming output file
    //string filename = modify_filename(nexus_file, "clrd_");
    string filename = "./temp_nexus";

    ifstream plain_nexus(nexus_file); // nexus input file
    ofstream colored_nexus(filename.c_str()); // colored nexus output file
    if(!plain_nexus.is_open()){
        cerr << "Error opening nexus file "<< nexus_file <<" to modify with color.\n";
        return;
    }
    if(!colored_nexus.is_open()){
        cerr << "Error creating colored nexus file " << filename << endl;
        return;
    }


    // read nexus data
    while (getline(plain_nexus, line)) {

        // end of any block
        if(line.find(";") != string::npos){

            // if end of vertices (#vertices has been assigned) add nodes at the end for legend
            if(vertices && no_vertices>0){
                int i = 0;
                //create_legend_pdf(grp_clr_map); // Extra pdf with only legend, doesn't work yet
                for (auto grp_clr : grp_clr_map) { // add node for each group
                    ostringstream oss;
                    ++no_vertices;
                    ++i;
                    rgb_color clr = grp_clr.second; // get color of group
                    // no_clr_map[no_vertices] = clr; // save new node with color with all other clrd nodes

                    // legend positioning TODO improve
                    pos_x = 0.3; //max((max_x+2)/2.0, 50.0); // +100
                    double h = abs(max_y-min_y);
                    pos_y = double(i) * 0.12; //min(h/no_grps, 1.0/no_grps); //double(i) * abs(max_y-min_y)/(double(no_grps)*3.0); // min(abs((min_y+1)/double(grp_clr_map.size()+1)), 150.0)
                    
                    if(clr.is_white()){ // if color is white, turn outline black (don't give values for fg, they are treated as 0 then)
                        oss << no_vertices << " " << (max_x+pos_x) <<" "<< (min_y+pos_y) <<" w="<<size<<" h="<<size<<" bg="<< clr.r<<" "<<clr.g<<" "<< clr.b<<",";
                    } else if (clr.is_default()){ // if no color given change the shape to rectangle (s=r)
                        cerr << "Invalid/No color given for group " << grp_clr.first << ". Using white rectangle\n";
                        clr.set_white();
                        oss << no_vertices << " " << (max_x+pos_x) <<" "<< min_y+(pos_y) <<" w="<<size<<" h="<<size<<" s=r bg="<< clr.r<<" "<<clr.g<<" "<< clr.b<<",";
                    } else { // normal, shape is circle, color is the given one
                        oss << no_vertices << " " << (max_x+pos_x) <<" "<< min_y+(pos_y) <<" w="<<size<<" h="<<size<<" fg="<<clr.r<<" "<<clr.g<<" "<< clr.b<<" bg="<< clr.r<<" "<<clr.g<<" "<< clr.b<<",";
                    }
                    colored_nexus << oss.str() << endl;

                    // save info for vlabels
                    legend << no_vertices << " '" << grp_clr.first << "' l=9 x=15 y=5 f='Dialog-PLAIN-" << textsize <<"',\n";

                }
            }
            // if end of vlabels add labels for groups
            if(vlabels && !vertices){
                colored_nexus << legend.str();
            }

            translate = false;
            vertices = false;
            vlabels = false;
        }

        // work with found blocks
        if(translate){ // map the vertex number to its taxname and then color
            istringstream iss(line);
            string no, taxname; // number, taxname of node
            if (getline(iss, no, ' ') && getline(iss, taxname)) {

                // removing quotes and comma
                taxname.erase(remove_if(taxname.begin(), taxname.end(),
                                        [](char c) { return (c == '\'' || c == ','); }), taxname.end());

                // get color via taxname -> group -> color
                if(tax_grp_map.find(taxname) != tax_grp_map.end()){ // check if taxa has been assigned to a group
                    rgb_color clr = grp_clr_map[tax_grp_map[taxname]];
                    no_clr_map[stoi(no)] = clr; // save node number -> color
                } else {
                    // test if several taxa joined at node
                    size_t space_pos = taxname.find(' ');
                    size_t pos = 0;
                    string name;
                    if(space_pos != string::npos){
                        name = taxname.substr(pos, space_pos - pos);
                        rgb_color current_clr= grp_clr_map[tax_grp_map[name]];
                        string current_group = tax_grp_map[name];

                        // search through node name and extract all actual taxnames
                        while(pos < taxname.size() && space_pos != string::npos){

                            space_pos = taxname.find(' ', pos);
                            name = taxname.substr(pos, space_pos - pos);
                            pos = space_pos+1;
                            // check if all the same group (= same color)
                            if( !(current_group == tax_grp_map[name])){ // !current_clr.is_equal(grp_clr_map[tax_grp_map[name]])
                                // if non equal colors: several taxa groups merged at one node
                                current_clr.set_black();
                                legend << no << " '" << taxname << "' l=9 x=15 y=5 f='Dialog-PLAIN-" << textsize <<"',\n";
                                //grp_clr_map[taxname] = current_clr; //problems bc of nvertices TODO instead of taxname all group names

                                cerr << "Warning: Several taxa of different groups have been joined at one node\n";
                                cerr << " Node will be colored black" << endl;
                                break;
                            }
                        }
                        no_clr_map[stoi(no)] = current_clr;

                    } else { // message if no group is assigned
                        cerr << "Warning: No group assigned to " << taxname << endl;
                    }
                }

            } else {
                cerr << "Error reading TRANSLATE block in nexus file.\n";
                cerr << "Problematic line: " << line << endl;
            }

        }

        if(vertices) { // color certain vertices
            // variables to save info contained for the vertex
            istringstream iss(line);
            int vertex_no;
            double x, y;
            char c_w, c_h, c_s, s, eq;
            int w, h;
            string rest;

            // splitting line into its components: [no.] [x] [y] w=[..] h=[..] s=n,
            if (iss >> vertex_no >> x >> y >> c_w >> eq >> w >> c_h >> eq >> h >> c_s >> eq >> s >> rest) {

                if (no_clr_map.find(vertex_no) != no_clr_map.end()){ // only color vertex if leaf

                    if(rest == ","){ // if last char is ',' add color
                        rgb_color clr = no_clr_map[vertex_no]; // get defined color for vertex
                        // add color information to line
                        ostringstream oss;
                        if(clr.is_white()){ // if color is white, turn outline black
                            oss <<vertex_no<<" "<<x<<" "<<y<<" w="<<size<<" h="<<size<<" bg="<< clr.r<<" "<<clr.g<<" "<< clr.b<<",";
                        } else if (clr.is_default()){ // if no color given change the shape to rectangle (s=r)
                            clr.set_white();
                            oss <<vertex_no<<" "<<x<<" "<<y<<" w="<<size<<" h="<<size<<" s=r bg="<< clr.r<<" "<<clr.g<<" "<< clr.b<<",";
                        } else {
                            oss <<vertex_no<<" "<<x<<" "<<y<<" w="<<size<<" h="<<size<<" fg="<<clr.r<<" "<<clr.g<<" "<< clr.b<<" bg="<< clr.r<<" "<<clr.g<<" "<< clr.b<<",";
                        }
                        line = oss.str();

                    } else { /* if already colored do not change */ }

                    // find 'corner' of graph to put legend later
                    if(x > max_x) max_x = x;
                    if(y > max_y) max_y = y;
                    if(y < min_y) min_y = y;
                }
            } else {
                cerr << "Problems reading vertex: " << line << endl;
            }
            no_vertices = vertex_no;
        }

        if(vlabels){ // label nodes without a group

            istringstream iss(line);
            string no, taxname; // number, taxname of node, additional info
            if ((getline(iss, no, ' ') && getline(iss, taxname))) {

                // if additional info is given for label, extract only the name
                size_t pos = taxname.find_first_of(' ');
                string name, more_info = "";
                if (pos != string::npos) {
                    more_info = taxname.substr(pos + 1); // extract everything after that
                    more_info.erase(remove_if(more_info.begin(), more_info.end(),
                                              [](char c) { return (c == ','); }), more_info.end());
                    taxname = taxname.substr(0, pos); // extract name
                }
                // removing quotes and comma
                taxname.erase(remove_if(taxname.begin(), taxname.end(),
                                        [](char c) { return (c == '\'' || c == ','); }), taxname.end());

                if(tax_grp_map.find(taxname) != tax_grp_map.end()){ // skip these
                    continue;
                } else { // if no group is assigned to taxa, keep its name at the node
                    ostringstream oss;
                    oss << no << " \'" << taxname << "\' " << more_info << ",";
                    line = oss.str();
                }

            }
        }

        // identify block
        if(line.find("TRANSLATE") != string::npos || line.find("Translate") != string::npos || line.find("translate") != string::npos){// starting point of: [number] [taxname]
            translate = true;
        }
        if(line.find("DIMENSIONS") != string::npos){
            // change nvertices=xx in DIMENSIONS
            size_t pos_str = line.find("nvertices=");
            if(pos_str != string::npos){

                pos_str += string("nvertices=").length();
                size_t end = line.find(" ", pos_str);
                // get old number of vertices and add number of labels
                int nvert = stoi(line.substr(pos_str, end - pos_str));
                nvert += grp_clr_map.size();
                // replace with the new number of vertices
                line.replace(pos_str, end - pos_str, to_string(nvert));
            }
        }
        if (line.find("VERTICES") != string::npos || line.find("Vertices") != string::npos || line.find("vertices") != string::npos) {// starting point of vertices
            vertices = true;
        }
        if(line.find("VLABELS") != string::npos){
            vlabels = true;
        }
        colored_nexus << line << endl;
    }

    plain_nexus.close();
    colored_nexus.close();

    std::rename(filename.c_str(), nexus_file.c_str());

    // special treatment if nodes would have multiple colors
    //      (draw two dots? other shape? different inner/outer color?)
    // nice-to-have: color for colorblind, or diff shape
}
