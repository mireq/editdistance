# distutils: language = c++
# distutils: sources = editdistance/_editdistance.cpp

from libc.stdlib cimport malloc, free
from libc.stdint cimport int64_t


cdef extern from "./_editdistance.h":
    unsigned int edit_distance_c(char *a, const size_t asize, char *b, size_t bsize)

cpdef unsigned int eval(object a, object b):
    return edit_distance_c(a, len(a), b, len(b))
