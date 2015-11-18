'''
This implementation is written in Cython.

It uses structs and low level types to operate in a known amount
of memory.

TODO

It does not currently support Up2Bit overflow. Sequences bigger than 63
letters should be broken up into multiple rows but this has not yet been
implemented.

It currently uses a fixed size table so may overflow.
'''


from libc.stdint cimport uint64_t

from .base_acgtrie import ACGTrieBase


#
# Up2Bit utils
#

cdef uint64_t up2bit_from_string(char *seq, int start, int end):
    cdef uint64_t value = 1
    cdef int bits
    cdef char c
    while start < end:
        c = seq[start]
        if c == 'A':
            bits = 0
        elif c == 'C':
            bits = 1
        elif c == 'G':
            bits = 2
        else:
            bits = 3
        value = (value << 2) | bits
        start += 1
    return value


cdef int up2bit_start_shift(uint64_t value):
    cdef int shift = 62
    while value >> shift == 0:
        shift -= 2
        if shift < 0:
            return -2
    return shift - 2


cdef object up2bit_to_string(uint64_t value):
    cdef int shift = up2bit_start_shift(value)
    result = ''
    while shift >= 0:
        bits = (value >> shift) & 0x3
        result += 'ACGT'[bits]
        shift -= 2
    return result


cdef int up2bit_len(uint64_t value):
    return ((up2bit_start_shift(value)) >> 1) + 1


cdef class Up2Bit:
    cdef uint64_t value

    def __init__(self, seq):
        self.value = up2bit_from_string(seq, 0, len(seq))

    def to_str(self):
        return up2bit_to_string(self.value)

    def to_int(self):
        return self.value

    def __len__(self):
        return up2bit_len(self.value)

    def __repr__(self):
        return 'Up2Bit(0x{:x} = {})'.format(self.value, self.to_str())


#
# Low level row structure - actually stored and manipulated.
#


cdef struct LLRow:
    int count
    int warps[4]
    uint64_t seq


cdef void ll_row_clear(LLRow *row):
    row.count = 0
    row.warps[0] = 0
    row.warps[1] = 0
    row.warps[2] = 0
    row.warps[3] = 0
    row.seq = 1


cdef object ll_row_make_wrapper(LLRow *row):
    cdef Row wrapper
    wrapper = Row()
    wrapper.ll_row = row[0]
    return wrapper


cdef int ascii_to_warp_idx(char c):
    if c == 'A':
        return 0
    elif c == 'C':
        return 1
    elif c == 'G':
        return 2
    else:
        return 3


#
# High level row wrapper, a temporary interface around a low level row.
#

cdef class Row(object):
    '''
    This row implementation is interface compatible with base_acgtrie's.
    '''
    cdef LLRow ll_row

    @property
    def count(self):
        return self.ll_row.count

    @property
    def a(self):
        return self.ll_row.warps[0]

    @property
    def c(self):
        return self.ll_row.warps[1]

    @property
    def g(self):
        return self.ll_row.warps[2]

    @property
    def t(self):
        return self.ll_row.warps[3]

    @property
    def seq(self):
        return up2bit_to_string(self.ll_row.seq)


#
# ScanResults capture enough information about a sequence lookup that it
# is possible to subsequently add a new sequence.
#


cdef struct ScanResults:
    bint matched
    int row_idx
    int start
    int seq_match


cdef ScanResults new_scan_results(
    bint matched, int row_idx, int start, int seq_match,
):
    cdef ScanResults result
    result.matched = matched
    result.row_idx = row_idx
    result.start = start
    result.seq_match = seq_match
    return result


#
# Low level trie implementation.
#
# Since Cython does not allow extension types to subclass Python types
# the low level code is implemented separately and combined as a mixin.
#


cdef class LLCythonACGTrie:
    cdef LLRow _rows[1000000]
    cdef int _row_count

    def __init__(self):
        ll_row_clear(&self._rows[0])
        self._row_count = 1

    def __len__(self):
        return self._row_count

    @property
    def rows(self):
        cdef int i
        for i in range(self._row_count):
            yield ll_row_make_wrapper(&self._rows[i])

    def lookup(self, seq):
        cdef ScanResults results = self.scan(seq, 0, len(seq), 0)
        if results.matched:
            return ll_row_make_wrapper(&self._rows[results.row_idx])
        return None

    cpdef add_subsequence(self, seq, int start, int end, int count):
        cdef ScanResults results = self.scan(seq, start, end, count)
        cdef LLRow *row
        cdef LLRow *split_row
        cdef int seq_idx, row_idx, split_row_idx, warp_idx
        cdef char sub
        row = &self._rows[results.row_idx]

        if results.seq_match >= 0:
            # The sequence did not match exactly, split this row.
            seq_idx = results.seq_match
            row_seq = up2bit_to_string(row.seq)
            sub = ord(row_seq[seq_idx])

            split_row_idx = self._row_count
            self._row_count += 1
            split_row = &self._rows[split_row_idx]
            split_row.count = row.count
            split_row.warps[0] = row.warps[0]
            split_row.warps[1] = row.warps[1]
            split_row.warps[2] = row.warps[2]
            split_row.warps[3] = row.warps[3]
            split_row.seq = up2bit_from_string(
                row_seq, seq_idx+1, len(row_seq),
            )

            row.seq = up2bit_from_string(row_seq, 0, seq_idx)
            row.warps[0] = 0
            row.warps[1] = 0
            row.warps[2] = 0
            row.warps[3] = 0
            row.warps[ascii_to_warp_idx(sub)] = split_row_idx

        row.count += count

        start = results.start
        if start < end:
            # Add a new row for the remaining sequence.
            warp_idx = ascii_to_warp_idx(ord(seq[start]))

            # Wire in the warp to a new row.
            row_idx = self._row_count
            self._row_count += 1
            row.warps[warp_idx] = row_idx

            # Create a new row with the remaining seq.
            row = &self._rows[row_idx]
            ll_row_clear(row)
            row.count = count
            row.seq = up2bit_from_string(seq, start+1, end)

    cdef ScanResults scan(self, seq, int start, int end, int add_count):
        cdef LLRow *row
        cdef int sub_idx, row_idx, warp_idx, new_row_idx
        row_idx = 0
        while True:
            row = &self._rows[row_idx]

            for sub_idx, sub in enumerate(up2bit_to_string(row.seq)):
                if start >= end:
                    return new_scan_results(True, row_idx, start, sub_idx)
                elif sub == seq[start]:
                    start += 1
                else:
                    return new_scan_results(False, row_idx, start, sub_idx)

            if start >= end:
                return new_scan_results(True, row_idx, start, -1)

            warp_idx = ascii_to_warp_idx(ord(seq[start]))
            new_row_idx = row.warps[warp_idx]
            if new_row_idx == 0:
                return new_scan_results(False, row_idx, start, -1)
            row.count += add_count
            start += 1
            row_idx = new_row_idx


class CythonACGTrie(LLCythonACGTrie, ACGTrieBase):
    pass
