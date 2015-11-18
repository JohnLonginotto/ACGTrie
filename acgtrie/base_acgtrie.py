'''
This base implementation defines the interface of all ACGTrie
implementations as well as some utility methods.
'''


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
        for start in xrange(end):
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
