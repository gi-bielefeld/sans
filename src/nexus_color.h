#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cctype>
#include <algorithm>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <limits>

using namespace std;

#ifndef SRC_NEXUS_COLOR_H
#define SRC_NEXUS_COLOR_H

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
    
    string toString(){return "rgb("+to_string(r)+","+to_string(g)+","+to_string(b)+")";}
};


class nexus_color{
	
 private:
	 
	static unordered_map<string, string> tax_grp_map; // map containing taxname -> group
	static unordered_map<string, rgb_color> grp_clr_map; // map conatining group -> color
	static unordered_map<int, rgb_color> no_clr_map; // map containing number -> color to color vertices
	static int no_grps; // number of groups given


 public:

    /**
     * Tests whether a program can be found and executed.
     * @param programName Program to be tested.
     * @return bool if program is executable
     */
    static bool program_in_path(const string& programName);

    /**
     * Modifies a filename with the given extension at the front.
     * @param file THe filename to modify
     * @param front_extension The extension to add
     * @return The new name
     */
    static string modify_filename(string& file, string front_extension);


    static string remove_extensions(string& filename);

    /**
     * This function adds a network to the initially generated nexus file via SplitsTree.
     * @param nexus_file Path to the initial nexus file
     * @param verbose If info should be printed
     * @param splitstree_path Path to SplitsTree
     * @param update If the given nexus file needs to be updated (if a network needs to be added)
     * @param save If the network should be saved to the given nexus file
     *             (chance of SplitsTree saving a not openable network but needed to add color)
     */
    static void open_in_splitstree(const string& nexus_file, const string& pdf, bool verbose = false, bool update = true, const string& save_as = "", const string splitstree_path = "SplitsTree");

    /**
     * This function adds color values to the nodes of a given nexus file. Already colored nodes
     * will not be recolored.
     * @param nexus_file The nexus file containing the network to be colored.
     * @param tax_grp_file A tab seperated file containing the taxa name and their respective group.
     * @param grp_clr_file A tab separated file containing the group and 3 integers for the rgb value of the color.
     */
    static void color_nexus(const string& nexus_file, const string& tax_grp_file, const string& grp_clr_file = "");



    static void scale_nexus(const string& unopened_nexus_file, bool verbose = false, bool scale_notification = false);
	
	static void prepare_marking(const string& tax_grp_file, const string& grp_clr_file);
	
	static string get_mark(string id);
	
	static string get_grp(string id);
};

#endif //SRC_NEXUS_COLOR_H
