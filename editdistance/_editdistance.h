#ifndef ___EDITDISTANCE__H__
#define ___EDITDISTANCE__H__

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

unsigned int edit_distance_c(const char *a, const size_t asize, const char *b, const size_t bsize);

void edit_distance_matrix(char **words, const size_t *sizes, size_t count);

#ifdef __cplusplus
}
#endif

#endif
