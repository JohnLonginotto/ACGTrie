import argparse
import csv
import struct
import sys
import time


from . import ACGTrie


def human_mem(bytes):
    if bytes > 1000000000:
        return '{:.2f}GB'.format(bytes / 1000000000.0)
    if bytes > 1000000:
        return '{:.2f}MB'.format(bytes / 1000000.0)
    if bytes > 1000:
        return '{:.2f}KB'.format(bytes / 1000.0)
    return '{:s}B'.format(bytes)


def up_to_29(seq):
    shift = 64 - 6
    result = len(seq) << shift
    for letter in seq:
        if letter == 'A':
            bits = 0
        elif letter == 'C':
            bits = 1
        elif letter == 'G':
            bits = 2
        else:
            bits = 3
        shift -= 2
        result |= bits << shift
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--out', dest='out', action='store',
        help='Save the trie to the specified file.',
    )
    args = parser.parse_args()

    print('Using {}.'.format(ACGTrie))

    start = time.time()
    reader = csv.reader(sys.stdin, delimiter=',')
    trie = ACGTrie()

    for sequence, text_count in reader:
        trie.add_sequence(sequence.rstrip(), int(text_count))

    print(trie)
    print('{} rows.'.format(len(trie)))

    try:
        allocated = trie.allocated()
    except NotImplementedError:
        pass
    else:
        print('Allocated {}.'.format(human_mem(allocated)))
        print('{} bytes per row.'.format(allocated / len(trie)))

    print('Trie built in {}s.'.format(time.time() - start))

    if args.out:
        start = time.time()
        fmt = struct.Struct('=i4iq')
        with open(args.out, 'wb') as stream:
            for row in trie.rows:
                stream.write(fmt.pack(
                    row.count,
                    row.a,
                    row.c,
                    row.g,
                    row.t,
                    up_to_29(row.seq),
                ))
        print('Output {} written in {}s.'.format(
            repr(args.out), time.time() - start,
        ))


if __name__ == '__main__':
    main()
