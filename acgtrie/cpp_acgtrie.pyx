# distutils: language = c++
# distutils: sources = src/cpp_acgtrie.cpp


### ## # # # cython: c_string_type=unicode, c_string_encoding=ascii


from libc.stdint cimport uint32_t, uint64_t
from libcpp.string cimport string

from .base_acgtrie import ACGTrieBase, Row


cdef extern from "../src/cpp_acgtrie.hpp":
    cdef cppclass VectorStack[T]:
        void clear()
        T& operator[](int)
        int size()
        void push_back(T&)


    cdef cppclass UpTo29:
        uint64_t value

        UpTo29(char *seq, int start, int end)
        string to_string()

    cdef cppclass LLRow "Row":
        uint32_t count
        uint32_t warps[4]
        UpTo29 seq

    cdef cppclass LLCPPACGTrie "CPPACGTrie":
        VectorStack[LLRow] rows

        void add_sequence(const char *seq, int start, int end, int count)
        void add_subsequence(const char *seq, int start, int end, int count)

        uint64_t allocated()
        uint32_t lookup(const char *seq, int start, int end)

        uint32_t get_row_count()
        LLRow get_row(uint32_t idx)



cdef class CPPACGTrieMixin:
    cdef LLCPPACGTrie *_trie

    def __init__(self):
        self._trie = new LLCPPACGTrie()

    def __del__(self):
        del self._trie

    def __len__(self):
        return self._trie.get_row_count()

    @property
    def rows(self):
        cdef uint32_t i
        cdef LLRow row_struct
        for i in xrange(self._trie.get_row_count()):
            row = self._trie.get_row(i)
            yield Row(
                row.count,
                row.warps[0],
                row.warps[1],
                row.warps[2],
                row.warps[3],
                row.seq.to_string().c_str().decode('ASCII'),
            )

    def add_sequence(self, seq, count):
        seq = seq.encode('ASCII')
        self._trie.add_sequence(seq, 0, len(seq), count)

    def add_subsequence(self, seq, start, end, count):
        seq = seq.encode('ASCII')
        self._trie.add_subsequence(seq, 0, len(seq), count)

    def lookup(self, seq):
        seq = seq.encode('ASCII')
        cdef int row_idx = self._trie.lookup(seq, 0, len(seq))
        if row_idx < 0:
            return None
        cdef LLRow row = self._trie.get_row(row_idx)
        return Row(
            row.count,
            row.warps[0],
            row.warps[1],
            row.warps[2],
            row.warps[3],
            str(row.seq.to_string().c_str()),
        )

    def allocated(self):
        return self._trie.allocated()

    def save(self, stream):
        cdef uint32_t i
        cdef LLRow *row
        cdef char *row_char_p
        cdef int size = self.row_struct.size
        for i in xrange(self._trie.get_row_count()):
            row = &self._trie.rows[i]
            stream.write(<bytes>(<char *>row)[:size])

    def clear(self):
        self._trie.rows.clear()

    def add_row(self, row):
        cdef LLRow ll_row
        ll_row.count = row.count
        ll_row.warps[0] = row.a
        ll_row.warps[1] = row.c
        ll_row.warps[2] = row.g
        ll_row.warps[3] = row.t
        ll_row.seq = UpTo29(row.seq.encode('ASCII'), 0, len(row.seq))
        self._trie.rows.push_back(ll_row)


class CPPACGTrie(CPPACGTrieMixin, ACGTrieBase):
    pass
