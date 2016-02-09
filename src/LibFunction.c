#include <stdlib.h>
#include <stdio.h>
#include <LibFunction.h>
#undef malloc
void* malloc_s(size_t size)
{
    void* new_mem = malloc(size);
    if (new_mem == NULL)
    {
        printf("Out of memory!\n");
        exit(1);
    }
    return new_mem;
}
