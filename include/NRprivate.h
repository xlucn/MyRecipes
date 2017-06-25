/** @file NRprivate.h */
#ifndef _NR_PRIVATE_H_
#define _NR_PRIVATE_H_

#include <stdlib.h>

#define malloc
#define realloc
void* malloc_s(size_t size);
void* realloc_s(void* ptr, size_t size);

#endif
