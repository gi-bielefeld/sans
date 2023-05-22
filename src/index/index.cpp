#include "index.h"


uint64_t Index::bins;

vector<hash_map<kmer_t, uint64_t>> Index::kmerMatrix;

// The mask for tree-id and color-id separation
uint32_t Index::color_mask;

subtree<colorQ0_t> Index::q0Tree;
subtree<colorQ2_t> Index::q2Tree;




void Index::init()
{
    // The number of hash_tables in the
    bins = (0b1u << MOD_POWER) + 0b1u;

    kmerMatrix = vector<hash_map<kmer_t, uint64_t>> (bins);
    uint32_t color_mask_bits = ~0b0u;
    color_mask = 0b0u | color_mask_bits;
}

/*
* Creating a compressed kmer entry
*/
kmer_t Index::compress_kmer(const kmer_t& kmer, uint64_t& bin)
{
    auto show_bits = [&] (kmer_t kmer)
    {
        for (int i = 0; i < 2*kmer::k; i++)
        {cout << kmer.test(i);}
        cout << endl;
    };

    // This lambda function performs binary subtraction
    auto sub = [&] (kmer_t minuend, kmer_t subtrahend)
    {
        kmer_t result = minuend;
        for (int i = 0; i < MOD_POWER + 1; i++)
        {
            int j = i;
            if (subtrahend.test(i))
            {

                if (result.test(i)){result.reset(i);}
                else
                {
                    result.set(i);
                    j++;
                    while (!result.test(j))
                    {
                        result.set(j);
                        j++;
                    }
                    result.reset(j);
                }
            }
        }
        return result;
    };

    // This lambda function performs binary division
    auto div = [&] (kmer_t dividend, kmer_t divisor)
    {
        kmer_t temp = 0b0u;
        kmer_t result = 0b0u;
        
        for (int i = 2*kmer::k - 1; i >= 0; i--)
        {
            temp <<= 1;
            result <<=1;
            if (dividend.test(i)){temp |= 0b1u;}
            if (temp >= divisor)
            {
                temp = sub(temp, divisor);
                result |= 0b1u;
            }
        }
        return result;
    };

    // Compute the kmer multiplier c, such that kmer = c * bins + bin
    // Compute kmer - bin
    kmer_t bit_bin = bin;
    kmer_t mod_multiple = sub(kmer, bit_bin);
    // Compute the multiplier (kmer - bin) / #bins
    kmer_t bit_module = bins;
    kmer_t mod_multipler = div(mod_multiple, bit_module);

    return mod_multipler;
}

// --- COLOR HANDLING ---
void Index::add_colored_kmer(const kmer_t& kmer, uint64_t& bin, const uint64_t& color)
{
    // Create the compressed kmer key
    // kmer_t kmer_entry = compress_kmer(kmer, bin);
    kmer_t kmer_entry = kmer;
    // If the kmer has not occured before
    if (!kmerMatrix[bin][kmer_entry])
    {
        // Create the quantile id and initialize the identifier 
        uint32_t color_id;
        color_id = q0Tree.incrementColorSupport(color);
        kmerMatrix[bin][kmer_entry] = color_id;
    }

    // If the kmer occurred before
    else{
        // Get the current concat identifier
        uint32_t color_id = kmerMatrix[bin][kmer_entry];
        uint32_t new_color_id = q0Tree.incrementColorSupport(color, color_id);
        kmerMatrix[bin][kmer_entry] = new_color_id;
    }   
}
