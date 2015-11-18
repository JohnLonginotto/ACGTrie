import ctypes
from os.path import dirname, join

from .base_acgtrie import ACGTrieBase, Row


# ctypes bindings for the library with a c interface
# The library is expected to be located in the same directory as
# this module.
dll = ctypes.CDLL(join(
    dirname(__file__),
    'libcpp_acgtrie.so',
))

cpp_acg_trie_add_subsequence = dll.cpp_acg_trie_add_subsequence
cpp_acg_trie_add_subsequence.argtypes = (
    ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_int,
)

cpp_acg_trie_add_sequence = dll.cpp_acg_trie_add_sequence
cpp_acg_trie_add_sequence.argtypes = (
    ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_int,
)


# There can be a maximum of 31 letters in a 64bit up2bit field,
# and the string should also be NULL terminated.
SEQ_BUF_LEN = 32


class RowStruct(ctypes.Structure):
    _fields_ = [
        ("count", ctypes.c_uint32),
        ("a", ctypes.c_uint32),
        ("c", ctypes.c_uint32),
        ("g", ctypes.c_uint32),
        ("t", ctypes.c_uint32),
        ("seq", ctypes.c_uint64),
    ]


class CPPACGTrie(ACGTrieBase):
    def __init__(self):
        self._trie = dll.cpp_acg_trie_new()

    def __del__(self):
        dll.cpp_acg_trie_delete(self._trie)

    @property
    def rows(self):
        i = 0
        l = dll.cpp_acg_trie_get_row_count(self._trie)
        row_struct = RowStruct()
        seq_buf = ctypes.create_string_buffer(SEQ_BUF_LEN)
        while i < l:
            dll.cpp_acg_trie_get_row(
                self._trie, i,
                ctypes.byref(row_struct),
                seq_buf,
            )
            yield Row(
                row_struct.count,
                row_struct.a, row_struct.c, row_struct.g, row_struct.t,
                seq_buf.value,
            )
            i += 1

    def add_sequence(self, seq, count):
        cpp_acg_trie_add_sequence(
            self._trie, seq, 0, len(seq), count,
        )

    def add_subsequence(self, seq, start, end, count):
        cpp_acg_trie_add_subsequence(
            self._trie, seq, start, end, count,
        )

    def lookup(self, seq):
        row_idx = dll.cpp_acg_trie_lookup(
            self._trie, seq, 0, len(seq),
        )
        if row_idx == -1:
            return None

        row_struct = RowStruct()
        seq_buf = ctypes.create_string_buffer(SEQ_BUF_LEN)

        dll.cpp_acg_trie_get_row(
            self._trie, row_idx,
            ctypes.byref(row_struct),
            seq_buf,
        )
        return Row(
            row_struct.count,
            row_struct.a, row_struct.c, row_struct.g, row_struct.t,
            seq_buf.value,
        )

    def allocated(self):
        return dll.cpp_acg_trie_allocated(self._trie)

    def __len__(self):
        return dll.cpp_acg_trie_get_row_count(self._trie)
