#include <vector>
#include <thread>
#include <mutex>

using namespace std;

class IDQueue
{

    protected:
        uint64_t max_multiplier = 0b0u;
        vector<uint64_t> queue;

    public:
        uint64_t bins;
        uint64_t offset;
        uint64_t pop();
        void push(uint64_t id);
};