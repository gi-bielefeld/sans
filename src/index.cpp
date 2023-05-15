#include "index.h"


hash_map<uint64_t, color_t> Index::color_by_id;
hash_map<color_t, uint64_t> Index::id_by_color;
hash_map<kmer_t, uint64_t> Index::kmer_color;

uint64_t Index::next_address = 0b1u;

void Index::add_colored_kmer(const kmer_t& kmer, const uint64_t& color)
{

    auto show_color = [](color_t c)
    {
            for (int i=0; i<maxN; i++)
            {cout << c.test(i);}
            cout << endl;
    };


    color_t bit_color  = (0b1u << color);

    // If the kmer is new it has not occured in any other color
    if (!kmer_color[kmer])
    {
        // If the corresponding color is new
        if (!id_by_color[bit_color])
        {
            // cout << "\nNew color " << color << endl;
            // show_color(bit_color);
            // Create a new entry for the new color
            color_by_id[next_address].set(color);
            id_by_color[bit_color] = next_address;
            // Create the kmer link
            kmer_color[kmer] = next_address;
            next_address++;
        }

        // If the corresponding color is not new
        else
        {
            kmer_color[kmer] = id_by_color[bit_color];
        }
    }

    // If the kmer is not new, it has occured before
    else
    {
        // Create the new correct color of the kmer
        bit_color |= color_by_id[kmer_color[kmer]];
        // If the update color has not uccured before 
        if (!id_by_color[bit_color])
        {
            // cout << "New Color : ";
            // show_color(bit_color);
            // Create the update color entries
            color_by_id[next_address] = bit_color;
            id_by_color[bit_color] = next_address;
            // Create the kmer link
            kmer_color[kmer] = next_address;
            next_address++;
        }

        // If the update color has occured before, update the stored kmer_color 
        else
        {
            kmer_color[kmer] = id_by_color[bit_color];
        }
    }

    /*
    cout << "INDEX" << endl;
    int i = 0;
    for (auto it : kmer_color)
    {
        i++;
        cout << "KMER : " << i << " COLOR_ADDRESS : " << it.second << endl;
    }
    */
}
