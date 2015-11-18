#include <stdint.h>
#include <string>
#include <vector>


template<typename T, int bucket_size_power=10>
struct VectorStack {
    // This class is very similar in operation to std::deque, but
    // seems to perform better than some implementations.
    // Also its memory usage can be calculated as an extra benefit.
    std::vector<T*> buckets;
    size_t _size;

    static const int bucket_size = 1 << bucket_size_power;
    static const int bucket_mask = bucket_size - 1;

    VectorStack() {
        _size = 0;
        buckets.reserve(1<<10);
    }

    ~VectorStack() {
        clear();
    }

    void clear() {
        for (T* bucket : buckets) {
            delete[] bucket;
        }
        _size = 0;
        buckets.clear();
        buckets.reserve(1<<10);
    }

    void push_back(const T &value) {
        size_t index = _size++;
        size_t bucket_index = index >> bucket_size_power;
        int item_index = index & bucket_mask;

        if (bucket_index >= buckets.size()) {
            buckets.push_back(new T[bucket_size]);
        }
        buckets[bucket_index][item_index] = value;
    }

    T &operator[](size_t index) {
        size_t bucket_index = index >> bucket_size_power;
        int item_index = index & bucket_mask;
        return buckets[bucket_index][item_index];
    }

    size_t size() {
        return _size;
    }

    size_t allocated() {
        return buckets.size() * bucket_size * sizeof(T);
    }
};



struct UpTo29 {
    // This struct can represent up to 29 letters in a 64bit field.
    // The first (top) 6 bits give the count of letters.
    // The next 58 bits are the letters themselves, left aligned, two
    // bits per letter for A C G and T.
    // Bits to the right of the letters are always set to 0.
    //
    // e.g. An empty sequence is length 0b000000 followed by 0 letters and
    // 58 bits of 0. This is 0.
    // The sequence C is length 0b000001 followed by 1 letter 0b01 and 56
    // bits of 0. This is
    // 0b0000010100000000000000000000000000000000000000000000000000000000
    //
    // Because only 5 bits are required for the length (0-29) the top bit
    // is always 0 and could be used for something else with some careful
    // masking.
    uint64_t value;

    UpTo29();
    UpTo29(const char *seq, int start, int end);

    std::string to_string();
    UpTo29 subseq(int start, int end);

    int bits(int idx);
    char operator[](int idx);

    int len();
};


#pragma pack(push, 1)
struct Row {
    uint32_t count;
    uint32_t warps[4];
    UpTo29 seq;

    Row();
    Row(const Row &other);
};
#pragma pack(pop)


struct CPPACGTrie {
    struct ScanResults {
        bool matched;
        int row_idx;
        int start;
        int seq_match;
    };

    VectorStack<Row> rows;

    CPPACGTrie();

    void add_sequence(const char *seq, int start, int end, int count);
    void add_subsequence(const char *seq, int start, int end, int count);

    uint64_t allocated();
    uint32_t lookup(const char *seq, int start, int end);

    uint32_t get_row_count();
    Row get_row(uint32_t idx);

    ScanResults scan(const char *seq, int start, int end, int add_count);
};
