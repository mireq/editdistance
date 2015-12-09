# distutils: language = c++
# distutils: sources = editdistance/_editdistance.cpp

from libc.stdlib cimport malloc, free
from libc.stdint cimport int64_t
from cython.view cimport array as cvarray
from cpython.string cimport PyString_AsString

cdef char ** to_cstring_array(list_str):
    cdef char **ret = <char **>malloc(len(list_str) * sizeof(char *))
    for i in xrange(len(list_str)):
        ret[i] = PyString_AsString(list_str[i])
    return ret

cdef extern from "./_editdistance.h":
    unsigned int edit_distance_c(char *a, const size_t asize, char *b, size_t bsize)

cpdef unsigned int eval(object a, object b):
    return edit_distance_c(a, len(a), b, len(b))

cpdef unsigned int[:,:] distance_matrix(object words):
    cdef unsigned int[:,:] result = cvarray(shape=(len(words), len(words)), itemsize=sizeof(unsigned int), format="I")
    cdef unsigned int i, j, dist;
    cdef unsigned int[:] lengths = cvarray(shape=(len(words),), itemsize=sizeof(unsigned int), format="I")
    cdef list words_b = []

    for i in xrange(len(words)):
        lengths[i] = len(words[i])
        words_b.append(words[i].encode('utf-8'))

    cdef char **c_words = to_cstring_array(words_b)

    for i in xrange(len(words)):
        result[i, i] = 0
        for j in xrange(i + 1, len(words)):
            dist = edit_distance_c(c_words[i], lengths[i], c_words[j], lengths[j])
            result[i, j] = dist
            result[j, i] = dist

    free(c_words)
    return result
