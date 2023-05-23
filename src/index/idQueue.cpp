#include "idQueue.h"


uint64_t IDQueue::pop()
{
    if (queue.size()==0)
    {
        max_multiplier++;
        queue.push_back(max_multiplier * bins + offset);
    }

    uint64_t id = *queue.begin();
    queue.erase(queue.begin());
    return id;
}

void IDQueue::push(uint64_t id)
{
    queue.push_back(id);
}




