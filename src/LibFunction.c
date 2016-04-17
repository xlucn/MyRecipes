#include <stdlib.h>
#include <stdio.h>
#include "LibFunction.h"
#undef malloc
void* malloc_s(size_t size)
{
    if(size < 0)
    {
        fprintf(stderr, "allocated space size is negative\n");
        exit(1);
    }
    void* new_mem = malloc(size);
    if (new_mem == NULL)
    {
        fprintf(stderr, "Allocation of memory failed!\n");
        exit(1);
    }
    return new_mem;
}

#undef realloc
void* realloc_s(void* ptr, size_t size)
{
    if(size < 0)
    {
        fprintf(stderr, "allocated space size is negative\n");
        exit(1);
    }
    void* new_mem = realloc(ptr, size);
    if (new_mem == NULL)
    {
        fprintf(stderr, "Allocation of memory failed!\n");
        exit(1);
    }
    return new_mem;
}
