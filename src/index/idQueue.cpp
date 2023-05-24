#include "idQueue.h"
#include <iostream>

using namespace std;

uint64_t IDQueue::pop()
{
    if (queue.size()==0)
    {
        max_multiplier++;
        uint64_t new_id = (max_multiplier * bins + offset);
        if (new_id > numeric_limits<uint64_t>::max()){cout << "IDENTIFIER EXCEEDED MAXIMUM VALUE" << endl;}
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




