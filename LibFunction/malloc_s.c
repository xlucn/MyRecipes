#include <stdlib.h>
#include <stdio.h>
#include "NR.h"


#undef malloc
/**
 * @brief safer malloc function with bland pointer check
 */
void* malloc_s(size_t size)
{
    if(size < 0)
    {
        fprintf(stderr, "malloc: allocated space size is negative\n");
        exit(1);
    }
    void* new_mem = malloc(size);
    if (new_mem == NULL)
    {
        fprintf(stderr, "malloc: Allocation of memory failed!\n");
        exit(1);
    }
    return new_mem;
}
