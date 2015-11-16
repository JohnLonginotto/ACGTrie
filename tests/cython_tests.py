import pyximport
pyximport.install()

from unittest import TestCase

from acgtrie.cython_acgtrie import CythonACGTrie, Up2Bit

from trie_building_mixin import TrieBuildingMixin


class TrieBuildingTests(TrieBuildingMixin, TestCase):
    ACGTrie = CythonACGTrie


class Up2BitTests(TestCase):
    def test_storage(self):
        bseq = Up2Bit('')
        self.assertEqual(len(bseq), 0)
        self.assertEqual(bseq.to_str(), '')

        bseq = Up2Bit('A')
        self.assertEqual(len(bseq), 1)
        self.assertEqual(bseq.to_str(), 'A')

        bseq = Up2Bit('AC')
        self.assertEqual(len(bseq), 2)
        self.assertEqual(bseq.to_str(), 'AC')

    def test_encoding(self):
        bseq = Up2Bit('')
        self.assertEqual(bseq.to_int(), 1)

        bseq = Up2Bit('A')
        self.assertEqual(bseq.to_int(), 4)

        bseq = Up2Bit('AC')
        self.assertEqual(bseq.to_int(), 0x11)
