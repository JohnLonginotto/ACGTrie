'''
This implementation tries to be as simple as possible.

It is idiomatic Python and uses an unreasonable amount of
memory, but can serve as a interesting starting point for
trie construction.
'''

from .base_acgtrie import ACGTrieBase


class Row(object):
    '''
    This row implementation is interface compatible with base_acgtrie's.
    '''
    def __init__(self):
        self.count = 0
        self.warps = [0, 0, 0, 0]
        self.seq = ''

    @property
    def a(self):
        return self.warps[0]

    @property
    def c(self):
        return self.warps[1]

    @property
    def g(self):
        return self.warps[2]

    @property
    def t(self):
        return self.warps[3]


class ScanResults(object):
    def __init__(self, matched, row_idx, start, seq_match):
        self.matched = matched
        self.row_idx = row_idx
        self.start = start
        self.seq_match = seq_match

    def __repr__(self):
        return 'ScanResults({}, {}, {}, {})'.format(
            self.matched, self.row_idx, self.start, self.seq_match,
        )


class SimpleACGTrie(ACGTrieBase):
    def __init__(self):
        self._rows = [Row()]

    def __len__(self):
        return len(self._rows)

    @property
    def rows(self):
        return self._rows

    def lookup(self, seq):
        results = self.scan(seq, 0, len(seq), 0)
        if results.matched:
            return self._rows[results.row_idx]
        return None

    def add_subsequence(self, seq, start, end, count):
        results = self.scan(seq, start, end, count)
        row = self._rows[results.row_idx]

        if results.seq_match >= 0:
            # The sequence did not match exactly, split this row.
            seq_idx = results.seq_match
            sub = row.seq[seq_idx]

            split_row_idx = len(self._rows)
            split_row = Row()
            split_row.count = row.count
            split_row.warps = row.warps
            split_row.seq = row.seq[seq_idx+1:]
            self._rows.append(split_row)

            row.seq = row.seq[:seq_idx]
            row.warps = [0, 0, 0, 0]
            row.warps[self.ascii_to_warp_idx(sub)] = split_row_idx

        row.count += count

        start = results.start
        if start < end:
            # Add a new row for the remaining sequence.
            warp_idx = self.ascii_to_warp_idx(seq[start])

            # Wire in the warp to a new row.
            row_idx = len(self._rows)
            row.warps[warp_idx] = row_idx

            # Create a new row with the remaining seq.
            row = Row()
            self._rows.append(row)
            row.count += count
            row.seq = seq[start+1:end]

    def scan(self, seq, start, end, add_count):
        row_idx = 0
        while True:
            row = self._rows[row_idx]

            for sub_idx, sub in enumerate(row.seq):
                if start >= end:
                    return ScanResults(True, row_idx, start, sub_idx)
                elif sub == seq[start]:
                    start += 1
                else:
                    return ScanResults(False, row_idx, start, sub_idx)

            if start >= end:
                return ScanResults(True, row_idx, start, -1)

            warp_idx = self.ascii_to_warp_idx(seq[start])
            new_row_idx = row.warps[warp_idx]
            if new_row_idx == 0:
                return ScanResults(False, row_idx, start, -1)
            row.count += add_count
            start += 1
            row_idx = new_row_idx

    def ascii_to_warp_idx(self, c):
        if c == 'A':
            return 0
        elif c == 'C':
            return 1
        elif c == 'G':
            return 2
        else:
            return 3
