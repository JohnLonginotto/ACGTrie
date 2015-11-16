'''
This module exports one of the ACGTrie implementations.

There are a few different implementations to chose from, by default
the one considered fastest is chosen, but this can be overridden using
the ACGTRIE environment variable.
'''

from simple_acgtrie import ACGTrie


__all__ = ['ACGTrie']
