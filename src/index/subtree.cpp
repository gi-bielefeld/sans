#include "subtree.h"



template<typename T>
uint32_t subtree<T>::incrementColorSupport(const uint64_t& new_color)
{
    uint32_t color_id;
    T bit_color = (0b1u << new_color);
    // If the color has not occured before
    if (!id_by_color[bit_color]){
        color_id = idQueue.pop();
        // Create the new color entry
        color_by_id[color_id] = bit_color;
        id_by_color[bit_color] = color_id;
    }
    // If the color has occured before
    else {
        color_id = id_by_color[bit_color];
    }
    return color_id;
}

template<typename T>
uint32_t subtree<T>::incrementColorSupport(const uint64_t& new_color, uint32_t& current_color_id)
{
    T current_color = color_by_id[current_color_id];
    
    // If the bit is already set -> The kmer reoccurred in the same sequence -> nothing to do
    if (current_color.test(new_color))
    {
        return current_color_id;
    }

    // If the bit is not set
    else
    {
        // Create the updated color
        T update_color;
        update_color = current_color;
        update_color.set(new_color);
        // If the updated color already exists
        if (id_by_color[update_color])
        {
                return id_by_color[update_color];
        }
        // Create the update color
        else
        {
            uint32_t update_color_id = idQueue.pop();
            color_by_id[update_color_id] = update_color;
            id_by_color[update_color] = update_color_id;
            return update_color_id;
        }
    }
}