#include "index.h"


uint64_t Index::bins;

// The idQueues link colors and kmers
vector<IDQueue> Index::idQueue;

vector<mutex> Index::kmer_lock;
vector<hash_map<kmer_t, uint32_t>> Index::kmerMatrix;          // The binned hash table

vector<mutex> Index::color_lock;
vector<uint32_t> Index::color_period;
vector<hash_map<color_t, uint32_t>> Index::id_by_color; 
vector<hash_map<uint32_t, color_t>> Index::color_by_id;                              // The set of current colors
vector<hash_map<uint32_t, uint32_t>> Index::support;


void Index::init()
{
    // The number of hash_tables in the
    bins = (0b1u << MOD_POWER) + 0b1u;

    // Queue init
    idQueue = vector<IDQueue> (bins);
    for (int i = 0; i < bins; i++){idQueue[i].bins = bins; idQueue[i].offset = i;}

    // Kmer data init
    kmer_lock = vector<mutex> (bins);
    kmerMatrix = vector<hash_map<kmer_t, uint32_t>> (bins);
    
    // Color data init
    color_lock = vector<mutex> (bins);
    if (bins == 1){color_period.push_back(0);}

    else{
        uint32_t carry = 1;
        for (int i = 0; i <= maxN; i++)
        {
            color_period.push_back(carry);
            carry = (carry * 2) % bins;
        }
    }
    
    id_by_color = vector<hash_map<color_t, uint32_t>> (bins);
    color_by_id = vector<hash_map<uint32_t, color_t>> (bins);
    support = vector<hash_map<uint32_t, uint32_t>> (bins);
    cout << "INIT DONE" << endl;
}

/*
* Creating a compressed kmer entry
BUGGED
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

    // If the kmer has not occured before
    std::lock_guard<mutex> kl(kmer_lock[bin]); 
    if (!kmerMatrix[bin][kmer])
    {
        // Create the target color
        uint32_t color_id;
        color_t bit_color = (0b1u << color);
        uint32_t color_bin = color_period[color];
        std::lock_guard<mutex> bcl(color_lock[color_bin]); // Lock the bit_bin color column
        // If the color has not occured before
        if (!id_by_color[color_bin][bit_color]){
            color_id = idQueue[color_bin].pop();
            // Create the new color entry
            color_by_id[color_bin][color_id] = bit_color;
            id_by_color[color_bin][bit_color] = color_id;
            support[color_bin][color_id] = 0b1u;
        }
        // If the color has occured before
        else {
            color_id = id_by_color[color_bin][bit_color];
            support[color_bin][color_id]++;
        }
        kmerMatrix[bin][kmer] = color_id;
    }

    // If the kmer occurred before
    else{
        
        // Get the current concat identifier
        uint32_t color_id = kmerMatrix[bin][kmer];
        uint32_t color_bin = color_id % bins;

        color_t current_color;
        color_t update_color;

        // Start of current color lock section
        {
        std::lock_guard<mutex> ccl(color_lock[color_bin]); // Lock the current_color column
        current_color = color_by_id[color_bin][color_id];
        // Decrement the support value
        // If the bit is already set -> The kmer reoccurred in the same sequence -> nothing to do
        if (current_color.test(color)){return;}

        // Otherwise the color changes
        support[color_bin][color_id]--;
        if (support[color_bin][color_id] == 0)
        {
            id_by_color[color_bin].erase(current_color);
            color_by_id[color_bin].erase(color_id);
            support[color_bin].erase(color_id);
            idQueue[color_bin].push(color_id);
        } 
        } // End of current color lock section

        // Start of update color lock section
        {
        // Create the updated color
        update_color = current_color;
        update_color.set(color);
        uint32_t update_color_bin = (color_bin + color_period[color]) % bins;
        std::lock_guard<mutex> ucl(color_lock[update_color_bin]); // Lock the update color bin column
        // If the updated color already exists
        if (id_by_color[update_color_bin][update_color])
        {
            support[update_color_bin][id_by_color[update_color_bin][update_color]]++;
            kmerMatrix[bin][kmer] = id_by_color[update_color_bin][update_color];
        }

        // Create the update color
        else
        {
            uint32_t update_color_id = idQueue[update_color_bin].pop();
            color_by_id[update_color_bin][update_color_id] = update_color;
            id_by_color[update_color_bin][update_color] = update_color_id;
            support[update_color_bin][update_color_id] = 0b1u;
            kmerMatrix[bin][kmer] = update_color_id;
        }
        }
        // End of the update color lock section
    }
    return;
}
