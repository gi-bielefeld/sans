#include "index.h"


hash_map<uint64_t, color_t> Index::colors;
hash_map<color_t, uint64_t> Index::color_address;
hash_map<kmer_t, uint64_t> Index::colored_kmer;
hash_map<uint64_t, uint64_t> Index::colored_support;

uint64_t Index::next_address = 0;

void Index::add_colored_kmer(const kmer_t& kmer, const uint64_t& color)
{
    auto show_color = [] (color_t c) 
    {cout << c.test(0b0u) << c.test(0b1u) << c.test(0b10u) << endl;};
    
    /*
    cout << "\nADD_TO-INDEX" << endl;
    cout << "CURRENT STATE" << endl;
    for (auto it : color_address)
    {
        show_color(it.first);
    }
    */

    // If the kmer is new it has not occured in any other color
    if (!colored_kmer.contains(kmer))
    {
        // cout << "NEW KMER DETECTED" << endl;
        // Create the corresponding color
        color_t bit_color;
        bit_color.set(color);

        // If the corresponding color is new
        if (!color_address.contains(bit_color))
        {
            // cout << "\t LINKING NEW COLOR : ";
            // show_color(bit_color);

            // Create a new entry for the new color
            uint64_t bit_color_address = next_address;
            next_address++;
            colors[bit_color_address] = bit_color;
            color_address[bit_color] = bit_color_address;

            // Create the color support entry
            colored_support.emplace(bit_color_address, 0b1u);
            // Create the kmer link
            colored_kmer.emplace(kmer, bit_color_address);
        }

        // If the corresponding color already occured
        else
        {
            // cout << "\t SUPPORTING KNOWN COLOR : ";
            // show_color(bit_color);
            // Increase the colored support of the target color
            colored_support[color_address[bit_color]]+=1;

        }
    }

    // If the kmer is not new, it has occured before
    else
    {
        // cout << "UPDATING REOCCURING KMER" << endl;
        // if (!colored_kmer.contains(kmer)){cout << "KMER ERROR" << endl;}


        // Get the current color of the kmer
        uint64_t stored_kmer_color_address = colored_kmer[kmer];
        color_t stored_color = colors[stored_kmer_color_address];

        // If the kmer reoccured in this color, the color is already set -> nothing to do
        if (stored_color.test(color))
        {
            // cout << "\t TARGET COLOR " << color << " already set ... done" << endl;
            // goto end;
            return;
        }

        // Remove the kmer support from the current color
        colored_support[stored_kmer_color_address]--;
        // If there is no support, remove color entry
        if (colored_support[stored_kmer_color_address] == 0)
        {
            color_t old_color = colors[stored_kmer_color_address];
            colors.erase(stored_kmer_color_address);
            color_address.erase(old_color);
            colored_support.erase(stored_kmer_color_address);
        }

        // Create the new correct color of the kmer
        color_t update_color;
        update_color = stored_color;
        update_color.set(color);
        
        // cout << "\t UPDATE COLOR : ";
        // show_color(update_color);

        // If the update color has not uccured before 
        if (!color_address.contains(update_color))
        {
            //cout << "\t SUPPORTING NEW UPDATED COLOR : ";
            // show_color(update_color);

            // Create the update color entries
            uint64_t update_color_address = next_address;
            next_address++;
            colors[update_color_address] = update_color;
            color_address[update_color] = update_color_address;

            // Create the color support entry
            colored_support.emplace(update_color_address, 0b1u);
            // Create the kmer link
            colored_kmer[kmer] = update_color_address;
        }

        // If the update color has occured before, update the stored support kmer_color and color support
        else
        {
            // Get the stored color address
            // cout << "\t SUPPORTING KNOWN UPDATED COLOR : ";
            // show_color(update_color);
            uint64_t stored_update_color_address = color_address[update_color];
            colored_kmer[kmer] = stored_update_color_address;
            colored_support[stored_update_color_address]+=1;
        }
    }
    /*
    end : 
    cout << "UPDATED INDEX" << endl;
    for (auto it : color_address)
    {
        cout << it.first.test(0b0u) << it.first.test(0b1u) << (it).first.test(0b10u) << endl;
    }
    cout << "\n" << endl;
    */
}
