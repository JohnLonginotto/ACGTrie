from unittest import TestCase

from acgtrie.cpp_acgtrie import CPPACGTrie

from trie_building_mixin import TrieBuildingMixin, get_rows


class TrieBuildingTests(TrieBuildingMixin, TestCase):
    ACGTrie = CPPACGTrie
    seq_limit = 31

    def test_seq_overflows_to_next_row(self):
        seq = 'ACGT' * (self.seq_limit / 3)
        trie = self.ACGTrie()
        trie.add_subsequence(seq, 0, len(seq), 1)
        rows = get_rows(trie)
        print rows
        self.assertEqual(len(rows), 3)
        self.assertRowEqual(rows[0], 1, 1, 0, 0, 0, '')

        sub_seq = seq[1:]
        a = 2 if sub_seq[self.seq_limit] == 'A' else 0
        c = 2 if sub_seq[self.seq_limit] == 'C' else 0
        g = 2 if sub_seq[self.seq_limit] == 'G' else 0
        t = 2 if sub_seq[self.seq_limit] == 'T' else 0
        self.assertRowEqual(rows[1], 1, a, c, g, t, sub_seq[:self.seq_limit])
        self.assertRowEqual(rows[2], 1, 0, 0, 0, 0, sub_seq[self.seq_limit+1:])
