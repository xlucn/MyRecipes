#include <stdlib.h>
#ifndef _LF_H_
#define _LF_H_
#define malloc
#define realloc
void* malloc_s(size_t size);
void* realloc_s(void* ptr, size_t size);
#endif
