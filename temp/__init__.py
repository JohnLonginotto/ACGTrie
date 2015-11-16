'''
This module exports one of the ACGTrie implementations.

There are a few different implementations to chose from, by default
the one considered fastest is chosen, but this can be overridden using
the ACGTRIE environment variable.
'''

import logging
from os import environ


options = dict()


# Load CythonACGTrie.
try:
    import pyximport
    pyximport.install()
except:
    pass
else:
    from .cython_acgtrie import CythonACGTrie
    options['CythonACGTrie'] = CythonACGTrie


# Load SimpleACGTrie.
from simple_acgtrie import SimpleACGTrie
options['SimpleACGTrie'] = SimpleACGTrie


priorities = [
    'CythonACGTrie',
    'SimpleACGTrie',
]


if 'ACGTrie' in environ:
    ACGTrie = options[environ['ACGTrie']]
else:
    for option in priorities:
        if option in options:
            ACGTrie = options[option]
            break
    else:
        raise ValueError('No ACGTrie implementation selected')


logging.info('Using {}.'.format(ACGTrie))
