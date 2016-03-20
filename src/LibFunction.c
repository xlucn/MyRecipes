#include <stdlib.h>
#include <stdio.h>
#include "LibFunction.h"
#undef malloc
void* malloc_s(size_t size)
{
    if(size < 0)
    {
        printf("allocated space size is negative\n");
        exit(0);
    }
    void* new_mem = malloc(size);
    if (new_mem == NULL)
    {
        printf("Out of memory!\n");
        exit(1);
    }
    return new_mem;
}
