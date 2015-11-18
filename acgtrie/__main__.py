import csv
import logging
import sys


logging.basicConfig(level=logging.INFO)


from . import ACGTrie


def human_mem(bytes):
    if bytes > 1000000000:
        return '{:.2f}GB'.format(bytes / 1000000000.0)
    if bytes > 1000000:
        return '{:.2f}MB'.format(bytes / 1000000.0)
    if bytes > 1000:
        return '{:.2f}KB'.format(bytes / 1000.0)
    return '{:s}B'.format(bytes)


def main():
    logging.info('Using {}.'.format(ACGTrie))

    reader = csv.reader(sys.stdin, delimiter=',')
    trie = ACGTrie()

    for sequence, text_count in reader:
        trie.add_sequence(sequence.rstrip(), int(text_count))

    print(trie)
    try:
        allocated = trie.allocated()
    except NotImplementedError:
        pass
    else:
        print('Allocated {}.'.format(human_mem(allocated)))
        print('{} bytes per row.'.format(allocated / len(trie)))


if __name__ == '__main__':
    main()
