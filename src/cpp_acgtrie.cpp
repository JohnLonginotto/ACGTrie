#include <algorithm>
#include <vector>
#include <string.h>

#include "cpp_acgtrie.hpp"


using namespace std;


inline int letter_to_code(char c) {
    switch (c) {
    case 'A':
        return 0;
    case 'C':
        return 1;
    case 'G':
        return 2;
    default:
        return 3;
    }
}


//
// UpTo29 implementation.
//


UpTo29::UpTo29() {
    value = 0;
}


UpTo29::UpTo29(const char *seq, int start, int end) {
    value = end - start;
    value <<= (64 - 6);

    int shift = 64 - 8;
    uint64_t bits;
    char c;

    while (start < end) {
        c = seq[start];
        bits = letter_to_code(c);
        value |= bits << shift;
        ++start;
        shift -= 2;
    }
}


string UpTo29::to_string() {
    string result;
    int shift = 64 - 8;
    int l = len();
    result.reserve(l);
    for (int i=0; i<l; ++i) {
        int bits = (value >> shift) & 0x3;
        result += "ACGT"[bits];
        shift -= 2;
    }
    return result;
}


UpTo29 UpTo29::subseq(int start, int end) {
    UpTo29 result;
    uint64_t l = end-start;
    uint64_t data = value;
    // Shift data to the right to lose the tail.
    data >>= (29 - end) * 2;
    // Shift data to the left to lose the head.
    data <<= 64 - (l * 2);

    result.value = l << (64 - 6);
    result.value |= data >> 6;

    return result;
}


int UpTo29::bits(int idx) {
    int shift = 64 - 8 - (idx * 2);
    return (value >> shift) & 0x3;
}


char UpTo29::operator[](int idx) {
    int shift = 64 - 8 - (idx<<1);
    int bits = (value >> shift) & 0x3;
    return "ACGT"[bits];
}


int UpTo29::len() {
    return value >> (64 - 6);
}


//
// Row implementation.
//


Row::Row() {
    count = 0;
    warps[0] = 0;
    warps[1] = 0;
    warps[2] = 0;
    warps[3] = 0;
    seq = UpTo29();
}


Row::Row(const Row &other) {
    count = other.count;
    warps[0] = other.warps[0];
    warps[1] = other.warps[1];
    warps[2] = other.warps[2];
    warps[3] = other.warps[3];
    seq = other.seq;
}


CPPACGTrie::CPPACGTrie() {
    rows.push_back(Row());
}


void CPPACGTrie::add_sequence(const char *seq, int start, int end, int count) {
    for (; start<end; ++start) {
        add_subsequence(seq, start, end, count);
    }
}


void CPPACGTrie::add_subsequence(
    const char *seq, int start, int end, int count
) {
    ScanResults results = scan(seq, start, end, count);
    Row *row, split_row;
    int seq_idx, row_idx, split_row_idx, warp_idx;

    row = &rows[results.row_idx];

    if (results.seq_match >= 0) {
        // The sequence did not match exactly, split this row.
        seq_idx = results.seq_match;
        warp_idx = row->seq.bits(seq_idx);

        split_row_idx = rows.size();
        split_row = Row(*row);
        split_row.seq = row->seq.subseq(seq_idx+1, row->seq.len());
        rows.push_back(split_row);

        row->seq = row->seq.subseq(0, seq_idx);
        row->warps[0] = 0;
        row->warps[1] = 0;
        row->warps[2] = 0;
        row->warps[3] = 0;
        row->warps[warp_idx] = split_row_idx;
    }

    row->count += count;

    start = results.start;
    while (start < end) {
        // Add a new row for the remaining sequence.
        // If the sequence is more than 29 letters long
        // it is split into more rows, as that is the capacity
        // of a single UpTo29 field.
        warp_idx = letter_to_code(seq[start]);

        // Wire in the warp to a new row.
        row_idx = rows.size();
        row->warps[warp_idx] = row_idx;

        // Create a new row with the remaining seq.
        Row new_row;
        new_row.count = count;
        new_row.seq = UpTo29(seq, start+1, min(start+1+29, end));
        rows.push_back(new_row);

        row = &rows[row_idx];
        start += 30;
    }
}


uint64_t CPPACGTrie::allocated() {
    return rows.allocated();
}


uint32_t CPPACGTrie::lookup(const char *seq, int start, int end) {
    ScanResults results = scan(seq, start, end, 0);
    if (results.matched) {
        return results.row_idx;
    }
    return -1;
}


uint32_t CPPACGTrie::get_row_count() {
    return rows.size();
}


Row CPPACGTrie::get_row(uint32_t idx) {
    return rows[idx];
}


CPPACGTrie::ScanResults CPPACGTrie::scan(
    const char *seq, int start, int end, int add_count
) {
    int row_idx = 0;
    for (;;) {
        Row &row = rows[row_idx];
        int l = row.seq.len();

        // Consume characters from this row's seq.
        for (int i=0; i<l; ++i) {
            if (start >= end) {
                // The scan matched, but was not a whole row match.
                // The row's seq goes on to define a longer seq.
                // For insertion this will require a split.
                return {true, row_idx, start, i};
            } else if (row.seq[i] == seq[start]) {
                // Still matching.
                ++start;
            } else {
                // Match failed, there is no such seq in the trie.
                // The returned details are enough to split and
                // continue.
                return {false, row_idx, start, i};
            }
        }

        if (start >= end) {
            // The seq has been matched exactly.
            return {true, row_idx, start, -1};
        }

        int warp_idx = letter_to_code(seq[start]);
        int new_row_idx = row.warps[warp_idx];
        if (new_row_idx == 0) {
            // The seq is not in the trie, a new warp would be
            // required to add it.
            return {false, row_idx, start, -1};
        }

        // Since insertion needs to add the count to every visited
        // row in the trie the count is added to this row before moving
        // on.
        row.count += add_count;

        // Continue with the next warp.
        ++start;
        row_idx = new_row_idx;
    }
}


extern "C" {
    CPPACGTrie *cpp_acg_trie_new() {
        return new CPPACGTrie();
    }

    void cpp_acg_trie_delete(CPPACGTrie *trie) {
        delete trie;
    }

    void cpp_acg_trie_add_sequence(
        CPPACGTrie *trie, const char *seq, int start, int end, int count
    ) {
        trie->add_sequence(seq, start, end, count);
    }

    void cpp_acg_trie_add_subsequence(
        CPPACGTrie *trie, const char *seq, int start, int end, int count
    ) {
        trie->add_subsequence(seq, start, end, count);
    }

    int cpp_acg_trie_lookup(
        CPPACGTrie *trie, const char *seq, int start, int end
    ) {
        return trie->lookup(seq, start, end);
    }

    int cpp_acg_trie_get_row_count(CPPACGTrie *trie) {
        return trie->rows.size();
    }

    void cpp_acg_trie_get_row(CPPACGTrie *trie, int idx, Row *row, char *seq) {
        *row = trie->rows[idx];
        strcpy(seq, row->seq.to_string().c_str());
    }

    size_t cpp_acg_trie_allocated(CPPACGTrie *trie) {
        return trie->allocated();
    }
}
