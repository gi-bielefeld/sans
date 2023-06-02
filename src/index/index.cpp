#include "index.h"


uint64_t Index::bins;

// The idQueues link colors and kmers
vector<IDQueue> Index::idQueue;

vector<mutex> Index::kmer_lock;
vector<hash_map<kmer_t, uint64_t>> Index::kmerMatrix;          // The binned hash table

vector<mutex> Index::color_lock;

uint64_t Index::relevant_bits;
vector<uint64_t> Index::color_period;
uint64_t Index::color_period_sum;

vector<hash_map<color_t, uint64_t>> Index::id_by_color; 
vector<hash_map<uint64_t, color_t>> Index::color_by_id;                              // The set of current colors
vector<hash_map<uint64_t, array<uint32_t, 2>>> Index::support;


void Index::init(uint64_t num)
{
    // The number of hash_tables in the
    bins = (0b1u << MOD_POWER) + 0b1u;

    // Queue init
    idQueue = vector<IDQueue> (bins);
    for (int i = 0; i < bins; i++){idQueue[i].bins = bins; idQueue[i].offset = i;}

    // Kmer data init
    kmer_lock = vector<mutex> (bins);
    kmerMatrix = vector<hash_map<kmer_t, uint64_t>> (bins);
    
    // Color data init
    color_lock = vector<mutex> (bins); // vector storing the binning locks
    
    color_period_sum = 0;

    relevant_bits = num;
    if (bins == 1){color_period.push_back(0);}
    else{
        uint64_t carry = 1;
        for (int i = 0; i < num; i++)
        {
            color_period.push_back(carry);
            color_period_sum += carry;
            carry = (carry * 2) % bins;
        }
    }
    
    id_by_color = vector<hash_map<color_t, uint64_t>> (bins);
    color_by_id = vector<hash_map<uint64_t, color_t>> (bins);
    support = vector<hash_map<uint64_t, array<uint32_t, 2>>> (bins);
    cout << "INIT DONE" << endl;
}

kmer_t Index::compress_kmer(const kmer_t& kmer, uint64_t& kmer_bin)
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

    // Compute the kmer multiplier c, such that kmer = c * bins + kmer_bin
    // Compute kmer - kmer_bin
    kmer_t bit_bin = kmer_bin;
    kmer_t mod_multiple = sub(kmer, bit_bin);
    // Compute the multiplier (kmer - kmer_bin) / #bins
    kmer_t bit_module = bins;
    kmer_t mod_multipler = div(mod_multiple, bit_module);
    return mod_multipler;
}

// --- COLOR HANDLING ---
void Index::add_colored_kmer(const kmer_t& kmer, uint64_t& kmer_bin, const uint64_t& color)
{

    // Create the bit color and compute its bin
    color_t bit_color;
    bit_color.set(color);
    bool bit_color_parity = color::represent(bit_color);
    uint64_t bit_color_bin = 0;
    for (uint64_t it = 0b0u; it < maxN; it++)
    {
        if(bit_color.test(it)){bit_color_bin += color_period[it];}
    }
    bit_color_bin %= bins;

    kmer_t  kmer_multiplier = compress_kmer(kmer, kmer_bin);  // The representative c, where kmer = c * M + kmer_bin (Uses 2k-floor(log2(M)) bits)
    kmer_t  kmer_parity_probe = kmer_multiplier << 1;
    kmer_t  kmer_key;
    bool    kmer_parity;

    // Lock the kmer_bin column
    std::lock_guard<mutex> kl(kmer_lock[kmer_bin]);

    // Probe if the kmer has occured before (with either parity)
    bool found = kmerMatrix[kmer_bin][kmer_parity_probe];
    if(!found){kmer_parity_probe |= 0b1u; found = kmerMatrix[kmer_bin][kmer_parity_probe];}

    // If the kmer has not occured before
    if (!found)
    {
        kmer_parity = bit_color_parity; // Set the kmer_parity to the color bit_color_parity
        kmer_key = (kmer_multiplier << 1) | kmer_parity; // Store the kmer bit_color_parity in the kmer key 

        // Lock the bit color bin
        std::lock_guard<mutex> bcl(color_lock[bit_color_bin]); 

        uint64_t bit_color_id; // The identifier of the bit color

        // If the color has not occured before
        if (!id_by_color[bit_color_bin][bit_color])
        {
            bit_color_id = idQueue[bit_color_bin].pop(); // create a new ID for the new color
            
            // Create the new color entry
            color_by_id[bit_color_bin][bit_color_id] = bit_color;
            id_by_color[bit_color_bin][bit_color] = bit_color_id;
            
            // Create the support entry for the bit color
            support[bit_color_bin][bit_color_id][kmer_parity]++;
        }

        // If the bit_color has occured before
        else
        {
            bit_color_id = id_by_color[bit_color_bin][bit_color]; // Get the id of the bit-color
            support[bit_color_bin][bit_color_id][kmer_parity]++;
        }
        kmerMatrix[kmer_bin][kmer_key] = bit_color_id; // Set the color id entry in the kmer_matrix
    }

    // If the kmer occurred before
    else
    {
        // Get the kmer parity from the probe
        kmer_key = kmer_parity_probe;
        kmer_parity = kmer_parity_probe & 0b1u;

        // Fetch current color and identifier
        uint64_t current_color_id = kmerMatrix[kmer_bin][kmer_key];  // Get the current color identifier
        uint64_t current_color_bin = current_color_id % bins;        // Compute the current color bin

        color_t  update_color;          // The updated color of the kmer 
        uint64_t update_color_bin;      // The bin of the updated color of the kmer

        // Start of current color lock section
        {
            std::lock_guard<mutex> ccl(color_lock[current_color_bin]); // Lock the current_colors hash table

            color_t current_color = color_by_id[current_color_bin][current_color_id];

            // If the color is already represented -> The kmer reoccurred in the same sequence -> nothing to do
            if (kmer_parity ? !current_color.test(color) : current_color.test(color)){return;}

            // Otherwise the color changes
            support[current_color_bin][current_color_id][kmer_parity]--;
        

            // Delete the color if there is no support left
            if (support[current_color_bin][current_color_id][0] + support[current_color_bin][current_color_id][1] == 0)
            {
                id_by_color[current_color_bin].erase(current_color);
                color_by_id[current_color_bin].erase(current_color_id);
                support[current_color_bin].erase(current_color_id);
                idQueue[current_color_bin].push(current_color_id);
            }

            // Remove the old kmer_key
            kmerMatrix[kmer_bin].erase(kmer_key);

            update_color = current_color;
            if (kmer_parity)
            {
                color::complement(update_color);
            }

            // Create the new color
            update_color.set(color);
            kmer_parity = color::represent(update_color);
            kmer_key = (kmer_multiplier << 1) | kmer_parity;
            

            update_color_bin = 0;
            for (uint64_t it = 0b0u; it < maxN; it++)
            {
                    if(update_color.test(it)){update_color_bin += color_period[it];}
            }
            update_color_bin %= bins;
        } // End of current color lock section

        // Start of update color lock section
        {
            // Create the updated color
            std::lock_guard<mutex> ucl(color_lock[update_color_bin]); // Lock the update color bin

            // If the updated color already exists
            if (id_by_color[update_color_bin][update_color])
            {
                uint64_t update_color_id = id_by_color[update_color_bin][update_color];
                support[update_color_bin][update_color_id][kmer_parity]++;
                kmerMatrix[kmer_bin][kmer_key] = id_by_color[update_color_bin][update_color];
            }
            // Create the update color
            else
            {
                uint64_t update_color_id = idQueue[update_color_bin].pop();
                color_by_id[update_color_bin][update_color_id] = update_color;
                id_by_color[update_color_bin][update_color] = update_color_id;
                support[update_color_bin][update_color_id][kmer_parity] = 0b1u;
                kmerMatrix[kmer_bin][kmer_key] = update_color_id;
            }
        }
        // End of the update color lock section
    }
    return;
}
