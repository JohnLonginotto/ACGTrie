from unittest import TestCase

from acgtrie.simple_acgtrie import SimpleACGTrie

from trie_building_mixin import TrieBuildingMixin


class TrieBuildingTests(TrieBuildingMixin, TestCase):
    ACGTrie = SimpleACGTrie
