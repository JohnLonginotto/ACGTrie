'''
This base implementation defines the interface of all ACGTrie
implementations as well as some utility methods.
'''


from struct import Struct


class Row(object):
    '''
    A row interface exposing the high level details
    of each row.

    The seq is a string rather than any low level memory representation, as
    such Rows are typically a view over low level data.
    '''
    def __init__(self, count=0, a=0, c=0, g=0, t=0, seq=''):
        self.count = count
        self.a = a
        self.c = c
        self.g = g
        self.t = t
        self.seq = seq


class ACGTrieBase(object):
    def __len__(self):
        '''
        Returns the number of rows in the trie.
        '''
        raise NotImplementedError()

    @property
    def rows(self):
        '''
        A generator yielding a Row object for every row in the
        trie.
        '''
        raise NotImplementedError()

    def add_sequence(self, seq, count):
        '''
        Adds a sequence and all its subsequences to the trie.
        '''
        end = len(seq)
        for start in range(end):
            self.add_subsequence(seq, start, end, count)

    def lookup(self, seq):
        '''
        Returns the row matching the given seq, or None
        if it does not exist.
        '''
        raise NotImplementedError()

    def get_count(self, seq):
        '''
        Returns the count of occurences of seq found.
        '''
        row = self.lookup(seq)
        if row is None:
            return 0
        return row.count

    def add_subsequence(self, seq, start, end, count):
        raise NotImplementedError()

    def to_table(self, max_rows=None):
        fmt = '{:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:<8}'
        output = [fmt.format('row', 'A', 'C', 'T', 'G', '#', 'Seq')]
        for idx, row in enumerate(self.rows):
            if max_rows is not None and idx >= max_rows:
                output.append('...')
                break
            output.append(fmt.format(
                idx,
                row.a, row.c, row.t, row.g,
                row.count, row.seq,
            ))
        return '\n'.join(output)

    def allocated(self):
        raise NotImplementedError()

    def __repr__(self):
        return self.to_table(10)

    row_struct = Struct('=i4iq')

    def save(self, stream):
        for row in self.rows:
            stream.write(self.row_struct.pack(
                row.count,
                row.a,
                row.c,
                row.g,
                row.t,
                seq_to_up_to_29(row.seq),
            ))

    def load(self, stream):
        self.clear()

        while True:
            data = stream.read(self.row_struct.size)
            if len(data) != self.row_struct.size:
                break

            count, a, c, g, t, seq = self.row_struct.unpack(data)
            self.add_row(Row(count, a, c, g, t, up_to_29_to_seq(seq)))

    def clear(self):
        raise NotImplementedError()

    def add_row(self, row):
        raise NotImplementedError()


def seq_to_up_to_29(seq):
    shift = 64 - 6
    result = len(seq) << shift
    for letter in seq:
        if letter == 'A':
            bits = 0
        elif letter == 'C':
            bits = 1
        elif letter == 'G':
            bits = 2
        else:
            bits = 3
        shift -= 2
        result |= bits << shift
    return result


def up_to_29_to_seq(up_to):
    seq = ''
    l = up_to >> (64 - 6)
    shift = 64 - 8
    for i in range(l):
        seq += 'ACGT'[(up_to >> shift) & 3]
        shift -= 2
    return seq
