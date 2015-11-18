try:
    from CStringIO import StringIO
except ImportError:
    from io import BytesIO as StringIO


def get_rows(trie):
    '''
    The Trie class has a rows property with a generator compatible
    interface. Gather each row into a list.
    '''
    return [row for row in trie.rows]


class TrieBuildingMixin(object):
    ACGTrie = None

    def test_blank_trie(self):
        trie = self.ACGTrie()
        rows = get_rows(trie)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].count, 0)
        self.assertEqual(rows[0].a, 0)
        self.assertEqual(rows[0].c, 0)
        self.assertEqual(rows[0].g, 0)
        self.assertEqual(rows[0].t, 0)
        self.assertEqual(rows[0].seq, '')

    def test_single_seq_adds_single_row(self):
        trie = self.ACGTrie()
        trie.add_subsequence('ACG', 0, 3, 1)
        rows = get_rows(trie)
        self.assertEqual(len(rows), 2)
        self.assertRowEqual(rows[0], 1, 1, 0, 0, 0, '')
        self.assertRowEqual(rows[1], 1, 0, 0, 0, 0, 'CG')

    def test_full_match_does_not_add_row(self):
        trie = self.ACGTrie()
        trie.add_subsequence('ACG', 0, 3, 1)
        rows = get_rows(trie)
        self.assertEqual(len(rows), 2)

        trie.add_subsequence('ACG', 0, 3, 1)
        rows = get_rows(trie)
        self.assertEqual(len(rows), 2)

        self.assertRowEqual(rows[0], 2, 1, 0, 0, 0, '')
        self.assertRowEqual(rows[1], 2, 0, 0, 0, 0, 'CG')

    def test_superseq_match_does_adds_row(self):
        trie = self.ACGTrie()
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
        trie = self.ACGTrie()
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
        trie = self.ACGTrie()
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
        trie = self.ACGTrie()
        trie.add_sequence('ACG', 1)
        self.assertEqual(trie.get_count('A'), 1)
        self.assertEqual(trie.get_count('C'), 1)
        self.assertEqual(trie.get_count('G'), 1)
        self.assertEqual(trie.get_count('AC'), 1)
        self.assertEqual(trie.get_count('CG'), 1)

        self.assertEqual(trie.get_count('AG'), 0)

    def test_save_load(self):
        stream = StringIO()
        src = self.ACGTrie()
        src.add_sequence('ACG', 1)
        src.save(stream)

        dst = self.ACGTrie()
        stream.seek(0, 0)
        dst.load(stream)

        self.assertEqual(len(src), len(dst))
        for src_row, dst_row in zip(src.rows, dst.rows):
            self.assertRowEqual(
                src_row,
                dst_row.count,
                dst_row.a,
                dst_row.c,
                dst_row.g,
                dst_row.t,
                dst_row.seq,
            )

    def assertRowEqual(self, row, count, a, c, g, t, seq):
        self.assertEqual(row.count, count)
        self.assertEqual(row.a, a)
        self.assertEqual(row.c, c)
        self.assertEqual(row.g, g)
        self.assertEqual(row.t, t)
        self.assertEqual(row.seq, seq)
