#include <stdio.h>
#include "NR.h"


#undef realloc
/**
 * @brief safer realloc function with null pointer check
 */
void* realloc_s(void* ptr, size_t size)
{
    if(size < 0)
    {
        fprintf(stderr, "realloc: allocated space size is negative\n");
        exit(1);
    }
    void* new_mem = realloc(ptr, size);
    if (new_mem == NULL)
    {
        fprintf(stderr, "realloc: Allocation of memory failed!\n");
        exit(1);
    }
    return new_mem;
}
