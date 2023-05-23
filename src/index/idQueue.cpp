#include "idQueue.h"


uint32_t IDQueue::pop()
{
    if (queue.size()==0)
    {
        max_multiplier++;
        queue.push_back(max_multiplier * bins + offset);
    }

    uint32_t id = *queue.begin();
    queue.erase(queue.begin());
    return id;
}

void IDQueue::push(uint32_t id)
{
    queue.push_back(id);
}




