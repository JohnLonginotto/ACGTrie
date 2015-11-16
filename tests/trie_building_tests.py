from unittest import TestCase

from acgtrie import ACGTrie


def get_rows(trie):
    '''
    The Trie class has a rows property with a generator compatible
    interface. Gather each row into a list.
    '''
    return [row for row in trie.rows]


class Tests(TestCase):
    def test_blank_trie(self):
        trie = ACGTrie()
        rows = get_rows(trie)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].count, 0)
        self.assertEqual(rows[0].a, 0)
        self.assertEqual(rows[0].c, 0)
        self.assertEqual(rows[0].g, 0)
        self.assertEqual(rows[0].t, 0)
        self.assertEqual(rows[0].seq, '')

    def test_single_seq_adds_single_row(self):
        trie = ACGTrie()
        trie.add_subsequence('ACG', 0, 3, 1)
        rows = get_rows(trie)
        self.assertEqual(len(rows), 2)
        self.assertRowEqual(rows[0], 1, 1, 0, 0, 0, '')
        self.assertRowEqual(rows[1], 1, 0, 0, 0, 0, 'CG')

    def test_full_match_does_not_add_row(self):
        trie = ACGTrie()
        trie.add_subsequence('ACG', 0, 3, 1)
        rows = get_rows(trie)
        self.assertEqual(len(rows), 2)

        trie.add_subsequence('ACG', 0, 3, 1)
        rows = get_rows(trie)
        self.assertEqual(len(rows), 2)

        self.assertRowEqual(rows[0], 2, 1, 0, 0, 0, '')
        self.assertRowEqual(rows[1], 2, 0, 0, 0, 0, 'CG')

    def test_superseq_match_does_adds_row(self):
        trie = ACGTrie()
        trie.add_subsequence('ACG', 0, 3, 1)
        rows = get_rows(trie)
        self.assertEqual(len(rows), 2)

        trie.add_subsequence('ACGT', 0, 4, 1)
        rows = get_rows(trie)
        self.assertEqual(len(rows), 3)

        self.assertRowEqual(rows[0], 2, 1, 0, 0, 0, '')
        self.assertRowEqual(rows[1], 2, 0, 0, 0, 2, 'CG')
        self.assertRowEqual(rows[2], 1, 0, 0, 0, 0, '')

    def test_partial_seq_match_causes_split(self):
        trie = ACGTrie()
        trie.add_subsequence('ACG', 0, 3, 1)
        rows = get_rows(trie)
        self.assertEqual(len(rows), 2)

        trie.add_subsequence('AC', 0, 2, 1)
        rows = get_rows(trie)
        self.assertEqual(len(rows), 3)

        self.assertRowEqual(rows[0], 2, 1, 0, 0, 0, '')
        self.assertRowEqual(rows[1], 2, 0, 0, 2, 0, 'C')
        self.assertRowEqual(rows[2], 1, 0, 0, 0, 0, '')

    def test_empty_partial_seq_match_causes_split(self):
        trie = ACGTrie()
        trie.add_subsequence('ACG', 0, 3, 1)
        rows = get_rows(trie)
        self.assertEqual(len(rows), 2)

        trie.add_subsequence('A', 0, 1, 1)
        rows = get_rows(trie)
        self.assertEqual(len(rows), 3)

        self.assertRowEqual(rows[0], 2, 1, 0, 0, 0, '')
        self.assertRowEqual(rows[1], 2, 0, 2, 0, 0, '')
        self.assertRowEqual(rows[2], 1, 0, 0, 0, 0, 'G')

    def test_sequence_adds_all_subsequences(self):
        trie = ACGTrie()
        trie.add_sequence('ACG', 1)
        self.assertEqual(trie.lookup('A').count, 1)
        self.assertEqual(trie.lookup('C').count, 1)
        self.assertEqual(trie.lookup('G').count, 1)
        self.assertEqual(trie.lookup('AC').count, 1)
        self.assertEqual(trie.lookup('CG').count, 1)

        self.assertIsNone(trie.lookup('AG'))

    def assertRowEqual(self, row, count, a, c, g, t, seq):
        self.assertEqual(row.count, count)
        self.assertEqual(row.a, a)
        self.assertEqual(row.c, c)
        self.assertEqual(row.g, g)
        self.assertEqual(row.t, t)
        self.assertEqual(row.seq, seq)
