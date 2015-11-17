from unittest import TestCase

from acgtrie.cpp_acgtrie import CPPACGTrie

from trie_building_mixin import TrieBuildingMixin


class TrieBuildingTests(TrieBuildingMixin, TestCase):
    ACGTrie = CPPACGTrie
