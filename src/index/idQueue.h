#include <vector>
#include <thread>
#include <mutex>

using namespace std;

class IDQueue
{

    protected:
        uint32_t max_multiplier = 0b0u;
        vector<uint32_t> queue;

    public:
        uint32_t bins;
        uint32_t offset;
        uint32_t pop();
        void push(uint32_t id);
};