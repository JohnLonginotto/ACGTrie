#include <deque>
#include <stdint.h>
#include <string>
#include <string.h>


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


static const char LogTable256[256] =
{
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
    -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
    LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
    LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
#undef LT
};


int highest_bit(uint64_t v) {
    uint64_t t;
    int add = 0;

    if ((t = v >> 32)) {
        add = 32;
        v = t;
    }
    if ((t = v >> 16)) {
        add += 16;
        v = t;
    }
    if ((t = v >> 8)) {
        add += 8;
        v = t;
    }

    return add + LogTable256[v];
}


struct Up2Bit {
    uint64_t value;

    Up2Bit() {
        value = 1;
    }

    Up2Bit(const char *seq, int start, int end) {
        value = 1;

        int bits;
        char c;

        while (start < end) {
            c = seq[start];
            bits = letter_to_code(c);
            value = (value << 2) | bits;
            ++start;
        }
    }

    inline int start_shift() {
        return highest_bit(value) - 2;
    }

    string to_string() {
        string result;
        int shift = start_shift();
        while (shift >= 0) {
            int bits = (value >> shift) & 0x3;
            result += "ACGT"[bits];
            shift -= 2;
        }
        return result;
    }

    int len() {
        return (start_shift() >> 1) + 1;
    }
};


struct Row {
    uint32_t count;
    uint32_t warps[4];
    Up2Bit seq;

    Row() {
        clear();
    }

    void clear() {
        count = 0;
        warps[0] = 0;
        warps[1] = 0;
        warps[2] = 0;
        warps[3] = 0;
        seq = Up2Bit();
    }
};


struct ScanResults {
    bool matched;
    int row_idx;
    int start;
    int seq_match;
};


struct CPPACGTrie {
    deque<Row> rows;

    CPPACGTrie() {
        rows.push_back(Row());
    }

    void add_sequence(const char *seq, int start, int end, int count) {
        for (; start<end; ++start) {
            add_subsequence(seq, start, end, count);
        }
    }

    void add_subsequence(const char *seq, int start, int end, int count) {
        ScanResults results = scan(seq, start, end, count);
        Row *row, split_row;
        int seq_idx, row_idx, split_row_idx, warp_idx;
        string row_seq;
        char sub;

        row = &rows[results.row_idx];

        if (results.seq_match >= 0) {
            // The sequence did not match exactly, split this row.
            seq_idx = results.seq_match;
            row_seq = row->seq.to_string();
            sub = row_seq[seq_idx];

            split_row_idx = rows.size();
            split_row = Row();
            split_row.count = row->count;
            split_row.warps[0] = row->warps[0];
            split_row.warps[1] = row->warps[1];
            split_row.warps[2] = row->warps[2];
            split_row.warps[3] = row->warps[3];
            split_row.seq = Up2Bit(row_seq.c_str(), seq_idx+1, row_seq.size());
            rows.push_back(split_row);

            row->seq = Up2Bit(row_seq.c_str(), 0, seq_idx);
            row->warps[0] = 0;
            row->warps[1] = 0;
            row->warps[2] = 0;
            row->warps[3] = 0;
            row->warps[letter_to_code(sub)] = split_row_idx;
        }

        row->count += count;

        start = results.start;
        if (start < end) {
            // Add a new row for the remaining sequence.
            warp_idx = letter_to_code(seq[start]);

            // Wire in the warp to a new row.
            row_idx = rows.size();
            row->warps[warp_idx] = row_idx;

            // Create a new row with the remaining seq.
            Row new_row;
            new_row.count = count;
            new_row.seq = Up2Bit(seq, start+1, end);
            rows.push_back(new_row);
        }
    }

    ScanResults scan(const char *seq, int start, int end, int add_count) {
        int row_idx = 0;
        for (;;) {
            Row &row = rows[row_idx];
            if (row.seq.value != 1) {
                string row_seq = row.seq.to_string();

                for (int i=0, j=row_seq.size(); i<j; ++i) {
                    if (start >= end) {
                        return {true, row_idx, start, i};
                    } else if (row_seq[i] == seq[start]) {
                        ++start;
                    } else {
                        return {false, row_idx, start, i};
                    }
                }
            }

            if (start >= end) {
                return {true, row_idx, start, -1};
            }

            int warp_idx = letter_to_code(seq[start]);
            int new_row_idx = row.warps[warp_idx];
            if (new_row_idx == 0) {
                return {false, row_idx, start, -1};
            }
            row.count += add_count;
            ++start;
            row_idx = new_row_idx;
        }

        return {};
    }
};


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
        ScanResults results = trie->scan(seq, start, end, 0);
        if (results.matched) {
            return results.row_idx;
        }
        return -1;
    }

    int cpp_acg_trie_get_row_count(CPPACGTrie *trie) {
        return trie->rows.size();
    }

    void cpp_acg_trie_get_row(CPPACGTrie *trie, int idx, Row *row, char *seq) {
        *row = trie->rows[idx];
        strcpy(seq, row->seq.to_string().c_str());
    }
}
