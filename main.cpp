#include "debug_utilities.h"
#include "io_utilities.h"

#include <iostream>
#include <vector>
#include <numeric>
#include <queue>
#include <unordered_map>
#include <map>

const size_t BYTE_SIZE = 256;

struct BTree {
    unsigned char value;

    BTree *left;
    BTree *right;

    [[nodiscard]] bool has_value() const {
        return !this->left && !this->right;
    }

    explicit BTree(unsigned char byte_value = 0, BTree *left = nullptr, BTree *right = nullptr) :
            left(left), right(right), value(byte_value) {}

};

class bwt_cmp_reverse {
    const std::vector<unsigned char> &data;
public:
    explicit bwt_cmp_reverse(const std::vector<unsigned char> &bwt_data) : data(bwt_data) {}

    bool operator()(size_t left, size_t right) {
        return data[left] < data[right];
    }
};

size_t cyclic_index(
        const size_t &start_index,
        const size_t &offset,
        const size_t &n
) {
    return start_index + offset < n ? start_index + offset : start_index + offset - n;
}

class bwt_cmp_straight {
    const std::vector<unsigned char> &data;
public:
    explicit bwt_cmp_straight(const std::vector<unsigned char> &bwt_data) : data(bwt_data) {}

    bool operator()(size_t left, size_t right) {
        size_t i = 0;
        while (data[cyclic_index(left, i, data.size())] == data[cyclic_index(right, i, data.size())] &&
               i < data.size()) {
            ++i;
        }
        return data[cyclic_index(left, i, data.size())] < data[cyclic_index(right, i, data.size())];
    }
};

std::vector<unsigned char> bwt_reverse(
        const std::vector<unsigned char> &bwt_data,
        size_t row_index
) {
    std::vector<size_t> l_shift(bwt_data.size());
    std::iota(l_shift.begin(), l_shift.end(), 0);
    std::stable_sort(l_shift.begin(), l_shift.end(), bwt_cmp_reverse(bwt_data));

    std::vector<unsigned char> initial_data(bwt_data.size(), '0');
    for (size_t i = 0; i < bwt_data.size(); ++i) {
        initial_data[i] = bwt_data[l_shift[row_index]];
        row_index = l_shift[row_index];
    }
    return initial_data;
}

std::pair<size_t, std::vector<unsigned char>> bwt(
        std::vector<unsigned char> data
) {
    std::vector<size_t> shift_order(data.size());
    std::iota(shift_order.begin(), shift_order.end(), 0);
    std::stable_sort(shift_order.begin(), shift_order.end(), bwt_cmp_straight(data));

    std::vector<unsigned char> encoded(data.size());
    size_t shift_position;
    for (size_t i = 0; i < data.size(); ++i) {
        encoded[i] = data[cyclic_index(shift_order[i], data.size() - 1, data.size())];
        if (shift_order[i] == 0) shift_position = i;
    }
    return std::make_pair(shift_position, encoded);
}


std::vector<unsigned char> rle_encode(
        const std::vector<unsigned char> &data
) {
    std::vector<unsigned char> encoded_data;

    unsigned char sequence_length = 0;
    unsigned char previous_byte = data[0];
    encoded_data.push_back(previous_byte);

    for (size_t i = 1; i < data.size(); ++i) {
        if (previous_byte == data[i] && sequence_length < 255) {
            ++sequence_length;
            if (sequence_length == 1) {
                encoded_data.push_back(previous_byte);
            }
            if (i == data.size() - 1) {
                encoded_data.push_back(sequence_length - 1);
            }
        } else if (sequence_length > 0) {
            encoded_data.push_back(sequence_length - 1);
            sequence_length = 0;
            previous_byte = data[i];
            encoded_data.push_back(previous_byte);
        } else {
            encoded_data.push_back(data[i]);
            previous_byte = data[i];
        }
    }
    return encoded_data;
}

std::vector<unsigned char> rle_decode(
        const std::vector<unsigned char> &data
) {
    std::vector<unsigned char> decoded_data;
    unsigned char previous_byte = data[0];
    bool is_previous_byte_read = true;
    decoded_data.push_back(previous_byte);
    for (size_t i = 1; i < data.size(); ++i) {
        if (data[i] == previous_byte && is_previous_byte_read) {
            decoded_data.push_back(previous_byte);
            for (unsigned char j = 0; j < data[i + 1]; ++j) {
                decoded_data.push_back(previous_byte);
            }
            ++i;
            is_previous_byte_read = false;
        } else {
            previous_byte = data[i];
            is_previous_byte_read = true;
            decoded_data.push_back(data[i]);
        }
    }

    return decoded_data;
}

std::vector<unsigned char> move_to_front(
        std::vector<unsigned char> data
) {
    std::vector<unsigned char> alphabet(BYTE_SIZE);
    std::iota(alphabet.begin(), alphabet.end(), 0);
    std::vector<unsigned char> encoded_data(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        auto current_symbol = data[i];
        auto found_index_it = std::find_if(alphabet.begin(), alphabet.end(),
                                           [&current_symbol](const unsigned char &c) -> bool {
                                               return c == current_symbol;
                                           });
        auto found_index = found_index_it - alphabet.begin();
        encoded_data[i] = found_index;
        if (found_index != 0 && found_index_it != alphabet.end()) {
            std::rotate(alphabet.begin(), found_index_it, found_index_it + 1);
        }
    }
    return encoded_data;
}

std::vector<unsigned char> move_to_front_reverse(
        std::vector<unsigned char> data
) {
    std::vector<unsigned char> alphabet(BYTE_SIZE);
    std::iota(alphabet.begin(), alphabet.end(), 0);
    std::vector<unsigned char> decoded_data(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        auto current_index = data[i];
        decoded_data[i] = alphabet[current_index];

        auto found_index_it = alphabet.begin() + current_index;
        if (found_index_it != alphabet.begin() && found_index_it != alphabet.end()) {
            std::rotate(alphabet.begin(), found_index_it, found_index_it + 1);
        }
    }
    return decoded_data;
}

void traverse(
        BTree *root,
        const std::vector<bool> &code_word,
        std::unordered_map<unsigned char, std::vector<bool>> &code_words
) {
    if (root->has_value()) {
        code_words[root->value] = code_word;
        return;
    }
    std::vector<bool> left_codeword = code_word;
    std::vector<bool> right_codeword = code_word;
    left_codeword.push_back(false);
    right_codeword.push_back(true);
    traverse(root->left, left_codeword, code_words);
    traverse(root->right, right_codeword, code_words);
}

std::unordered_map<unsigned char, std::vector<bool>> build_hashmap(
        BTree *root
) {
    std::unordered_map<unsigned char, std::vector<bool>> code_words;
    std::vector<bool> initial_codeword;
    traverse(root, initial_codeword, code_words);
    return code_words;
}

std::vector<unsigned char> encode_with_huffman(
        const std::vector<unsigned char> &data,
        BTree *huffman_tree_root
) {
    auto encoded_data = std::vector<unsigned char>(1);
    size_t first_free_pos = 7;
    auto code_words = build_hashmap(huffman_tree_root);
    for (const auto &byte: data) {
        const auto &code_word = code_words[byte];
        for (bool bit: code_word) {
            append_bit(encoded_data, first_free_pos, bit);
        }
    }
    return encoded_data;
}

void dfs(
        BTree *node,
        std::vector<unsigned char> &result,
        size_t &bit_index
) {
    if (node->has_value()) {
        append_bit(result, bit_index, 0);
        append_byte(result, bit_index, node->value);
        return;
    }
    append_bit(result, bit_index, 1);
    dfs(node->left, result, bit_index);
    dfs(node->right, result, bit_index);
}

std::vector<unsigned char> tree_to_bytes(
        BTree *huffman_tree_root
) {
    std::vector<unsigned char> result(1);
    size_t bit_index = 7;
    dfs(huffman_tree_root, result, bit_index);
    return result;
}

BTree *bytes_to_tree_dfs(
        const std::vector<unsigned char> &encoded_tree,
        size_t &byte_index,
        size_t &bit_index
) {
    bool bit = read_bit(encoded_tree, byte_index, bit_index);
    if (!bit) {
        unsigned char value = read_byte(encoded_tree, byte_index, bit_index);

        auto leaf = new BTree(value);
        return leaf;
    }

    auto node = new BTree();
    auto left = bytes_to_tree_dfs(encoded_tree, byte_index, bit_index);
    node->left = new BTree(left->value, left->left, left->right);

    auto right = bytes_to_tree_dfs(encoded_tree, byte_index, bit_index);
    node->right = new BTree(right->value, right->left, right->right);

    return node;
}

BTree *bytes_to_tree(
        const std::vector<unsigned char> &encoded_tree
) {
    size_t byte_index = 0;
    size_t bit_index = 0;
    return bytes_to_tree_dfs(encoded_tree, byte_index, bit_index);
}

std::pair<std::vector<unsigned char>, BTree *>
huffman(const std::vector<unsigned char> &data) {
    std::vector<long> frequencies(BYTE_SIZE, 0);
    std::priority_queue<std::pair<long, BTree *>> vertex_pool;

    std::vector<bool> already_in_queue(BYTE_SIZE);
    for (const unsigned char &byte: data) {
        ++frequencies[byte];
    }
    for (const unsigned char &byte: data) {
        if (!already_in_queue[byte] && frequencies[byte] != 0) {
            auto new_vertex = new BTree(byte);
            vertex_pool.push(std::make_pair(-frequencies[byte], new_vertex));
            already_in_queue[byte] = true;
        }
    }
    while (vertex_pool.size() != 1) {
        auto left = vertex_pool.top();
        vertex_pool.pop();
        auto right = vertex_pool.top();
        vertex_pool.pop();

        long new_freq = left.first + right.first;
        auto new_node = new BTree(0, left.second, right.second);
        vertex_pool.push(std::make_pair(new_freq, new_node));
    }
    auto encoded_data = encode_with_huffman(data, vertex_pool.top().second);
    return std::make_pair(encoded_data, vertex_pool.top().second);
}

std::vector<unsigned char> huffman_reverse(
        const std::vector<unsigned char> &data,
        std::unordered_map<std::vector<bool>, unsigned char> &code_words,
        const size_t &initial_data_size
) {
    std::vector<unsigned char> decoded_data;
    size_t byte_index = 0;
    size_t bit_index = 0;
    size_t bytes_decoded = 0;
    while (bytes_decoded < initial_data_size) {

        std::vector<bool> cur_codeword;
        while (!code_words.contains(cur_codeword)) {
            bool bit = read_bit(data, byte_index, bit_index);
            cur_codeword.push_back(bit);
        }
        unsigned char decoded_word = code_words[cur_codeword];

        decoded_data.push_back(decoded_word);
        ++bytes_decoded;
    }
    return decoded_data;
}

template<class K, class V>
std::unordered_map<V, K> reverse_map(
        const std::unordered_map<K, V> &map
) {
    std::unordered_map<V, K> new_map;
    for (const auto &pair: map) {
        new_map[pair.second] = pair.first;
    }
    return new_map;
}

void print_metrics(
        const std::string &output_file,
        const size_t &initial_data,
        const size_t &encoded_data_size
);

const uint16_t END_OF_FILE = BYTE_SIZE + 1;
const uint16_t MAX_CODE_VALUE = (1 << 16) - 1;
const uint16_t QUARTER = MAX_CODE_VALUE / 4 + 1;

std::vector<uint16_t> byte_to_index;
std::vector<unsigned char> index_to_byte;
std::vector<uint16_t> freq;
std::vector<uint16_t> cumu;

uint16_t left, right, code_value = 0;
size_t byte_index_io, bit_index, first_free_pos = 0;
uint8_t free_bits = 0;

void update(uint16_t byte_index) {
    if (cumu[0] == QUARTER - 1) {
        uint16_t current_cumu = 0;
        for (int i = BYTE_SIZE + 1; i >= 0; --i) {
            cumu[i] = current_cumu;
            ++freq[i];
            freq[i] /= 2;
            current_cumu += freq[i];
        }
    }
    while (freq[byte_index] == freq[byte_index - 1]) {
        --byte_index;
    }
    unsigned char new_byte = index_to_byte[byte_index];
    unsigned char old_byte = index_to_byte[byte_index];
    index_to_byte[byte_index] = old_byte;
    index_to_byte[byte_index] = new_byte;
    byte_to_index[new_byte] = byte_index;
    byte_to_index[old_byte] = byte_index;
    ++freq[byte_index];
    while (byte_index) {
        --byte_index;
        ++cumu[byte_index];
    }
}

void reset_state() {
    left = 0;
    right = MAX_CODE_VALUE;
    code_value = 0;
    byte_index_io = 0;
    bit_index = 0;
    free_bits = 0;
    first_free_pos = 7;

    byte_to_index = std::vector<uint16_t>(BYTE_SIZE);
    index_to_byte = std::vector<unsigned char>(BYTE_SIZE + 2);
    freq = std::vector<uint16_t>(BYTE_SIZE + 2);
    cumu = std::vector<uint16_t>(BYTE_SIZE + 2);

    for (size_t i = 0; i < BYTE_SIZE; ++i) {
        byte_to_index[i] = i + 1;
        index_to_byte[i + 1] = i;

        freq[i] = i > 0 ? 1 : 0;
        cumu[i] = BYTE_SIZE + 1 - i;
    }
    freq[BYTE_SIZE] = 1;
    freq[BYTE_SIZE + 1] = 1;
    cumu[BYTE_SIZE] = 1;
}

void write_bit_with_free_bits(
        std::vector<unsigned char> &encoded_data,
        const bool &bit
) {
    append_bit(encoded_data, first_free_pos, bit);
    while (free_bits) {
        --free_bits;
        append_bit(encoded_data, first_free_pos, !bit);
    }
}

void encode_next_byte(
        const uint16_t &byte_index,
        std::vector<unsigned char> &encoded_data
) {
    uint32_t range = 1 + right - left;
    right = left + range * cumu[byte_index - 1] / cumu[0] - 1;
    left = left + range * cumu[byte_index] / cumu[0];
    for (;;) {
        if (right < 2 * QUARTER) {
            write_bit_with_free_bits(encoded_data, false);
        } else if (left >= 2 * QUARTER) {
            write_bit_with_free_bits(encoded_data, true);
            left -= 2 * QUARTER;
            right -= 2 * QUARTER;
        } else if (left >= QUARTER && right < 3 * QUARTER) {
            ++free_bits;
            left -= QUARTER;
            right -= QUARTER;
        } else {
            break;
        }
        left *= 2; right *= 2; ++right;
    }
}

std::vector<unsigned char> aak_encode(
        const std::vector<unsigned char> &initial_data
) {
    reset_state();
    std::vector<unsigned char> encoded_data = std::vector<unsigned char>(1);
    for (unsigned char byte: initial_data) {
        uint16_t byte_index = byte_to_index[byte];
        encode_next_byte(byte_index, encoded_data);
        update(byte_index);
    }
    encode_next_byte(END_OF_FILE, encoded_data);
    ++free_bits;
    bool last_bit = left < QUARTER ? last_bit = false : last_bit = true;
    write_bit_with_free_bits(encoded_data, last_bit);
    return encoded_data;
}

uint16_t next_byte_index(
        const std::vector<unsigned char> &encoded_data,
        size_t &bit_index_io
) {
    uint32_t range = 1 + right - left;
    uint16_t byte_index = 1;
    while (cumu[byte_index] > ((code_value - left + 1) * cumu[0] - 1) / range) {
        ++byte_index;
    }
    right = left + range * cumu[byte_index - 1] / cumu[0] - 1;
    left += range * cumu[byte_index] / cumu[0];
    for (;;) {
        if (right < 2 * QUARTER) {
        } else if (left >= 2 * QUARTER) {
            code_value -= 2 * QUARTER;
            left -= 2 * QUARTER;
            right -= 2 * QUARTER;
        } else if (left >= QUARTER && right < 3 * QUARTER) {
            code_value -= QUARTER;
            left -= QUARTER;
            right -= QUARTER;
        } else {
            break;
        }
        left *= 2;
        right *= 2;
        ++right;
        bool next_bit = read_bit(encoded_data, byte_index_io, bit_index_io);
        if (byte_index_io >= encoded_data.size()) {
            break;
        }
        code_value *= 2;
        code_value += next_bit;
    }
    return byte_index;
}

std::vector<unsigned char> aak_decode(
        const std::vector<unsigned char> &encoded_data
) {
    reset_state();
    size_t bit_index_io = 0;
    std::vector<unsigned char> initial_data = std::vector<unsigned char>(1);
    for (size_t i = 0; i < 16; ++i) {
        bool next_bit = read_bit(encoded_data, byte_index_io, bit_index_io);
        code_value *= 2;
        code_value += next_bit;
    }
    for (;;) {
        uint16_t byte_index = next_byte_index(encoded_data, bit_index_io);
        if (byte_index == END_OF_FILE) {
            break;
        }
        unsigned char byte = index_to_byte[byte_index];
        append_byte(initial_data, first_free_pos, byte);
        update(byte_index);
    }
    return initial_data;
}

void compress(
        const std::string &initial_file_name,
        const std::string &encoded_file_name
) {
    const auto &[bytes_input, dummy1, dummy2, dummy3] = read_bytes(initial_file_name);
    auto bwt_result = bwt(bytes_input);
    auto bwt_data = bwt_result.second;
    auto bwt_shift_position = bwt_result.first;

    auto mtf_data = move_to_front(bwt_data);
//    size_t mtf_data_size = mtf_data.size();
    auto encoded_rle = rle_encode(mtf_data);
    size_t mtf_data_size = encoded_rle.size();

    auto huffman_result = huffman(encoded_rle);
    auto huffman_data = huffman_result.first;
    auto huffman_tree_root = huffman_result.second;

    auto huffman_tree_encoded = tree_to_bytes(huffman_tree_root);
    auto size_of_tree = huffman_tree_encoded.size() * sizeof(unsigned char);

    size_t encoded_data_size = sizeof(unsigned char) * huffman_tree_encoded.size();
    encoded_data_size += 2 * sizeof(size_t) + sizeof(unsigned long);
    std::cout << "header size: " << double(encoded_data_size) << " $$ ";
    encoded_data_size += sizeof(unsigned char) * huffman_data.size();
    if (encoded_data_size + 1 < bytes_input.size()) {
        write_bytes(encoded_file_name, huffman_data, bwt_shift_position, mtf_data_size, size_of_tree,
                    huffman_tree_encoded);
        print_metrics(encoded_file_name, bytes_input.size(), encoded_data_size);
    } else {
        auto encoded2 = bytes_input;
        encoded2.push_back(255);
        print_metrics(encoded_file_name, bytes_input.size(), encoded2.size());
        std::ofstream fout(encoded_file_name, std::ios::binary);
        fout.write(reinterpret_cast<const char *>(encoded2.data()), static_cast<long>(encoded2.size()));
    }
}

void decompress(
        const std::string &encoded_file_name,
        const std::string &decoded_file_name
) {
    std::ifstream fin(encoded_file_name, std::ios::binary);
    std::vector<unsigned char> bytes((std::istreambuf_iterator<char>(fin)), {});
    if (bytes[bytes.size() - 1] == 255) {
        auto result = bytes;
        result.pop_back();
        write_bytes(decoded_file_name, result);
    } else {
        const auto &[encoded_huffman, bwt_shift_position, initial_data_size, huffman_tree_encoded] = read_bytes(
                encoded_file_name, true);
        auto huffman_tree_decoded = bytes_to_tree(huffman_tree_encoded);
        auto code_words = build_hashmap(huffman_tree_decoded);
        auto reversed_code_words = reverse_map(code_words);

        auto decoded_huffman = huffman_reverse(encoded_huffman, reversed_code_words, initial_data_size);
        auto decoded_rle = rle_decode(decoded_huffman);
        auto decoded_mtf = move_to_front_reverse(decoded_rle);

        auto decoded_data = bwt_reverse(decoded_mtf, bwt_shift_position);

        write_bytes(decoded_file_name, decoded_data);
    }
}

void compress_ak(
        const std::string &initial_file_name,
        const std::string &encoded_file_name
) {
    const auto &[bytes_input, dummy1] = read_bytes_ak(initial_file_name);
    auto bwt_result = bwt(bytes_input);
    auto bwt_data = bwt_result.second;
    auto bwt_shift_position = bwt_result.first;
    auto mtf_data = move_to_front(bwt_data);
    auto encoded_rle = rle_encode(mtf_data);

    auto ak_data = aak_encode(encoded_rle);

    size_t encoded_data_size = sizeof(size_t);
    std::cout << "header size: " << double(encoded_data_size) << " $$ ";
    encoded_data_size += sizeof(unsigned char) * ak_data.size();

    if (encoded_data_size + 1 < bytes_input.size()) {
        write_bytes_ak(encoded_file_name, ak_data, bwt_shift_position);
        print_metrics(encoded_file_name, bytes_input.size(), encoded_data_size);
    } else {
        auto encoded2 = bytes_input;
        encoded2.push_back(255);
        print_metrics(encoded_file_name, bytes_input.size(), encoded2.size());
        std::ofstream fout(encoded_file_name, std::ios::binary);
        fout.write(reinterpret_cast<const char *>(encoded2.data()), static_cast<long>(encoded2.size()));
    }
}

void decompress_ak(
        const std::string &encoded_file_name,
        const std::string &decoded_file_name
) {
    std::ifstream fin(encoded_file_name, std::ios::binary);
    std::vector<unsigned char> bytes((std::istreambuf_iterator<char>(fin)), {});
    if (bytes[bytes.size() - 1] == 255) {
        auto result = bytes;
        result.pop_back();
        write_bytes(decoded_file_name, result);
    } else {
        const auto &[encoded_ak, bwt_shift_position] = read_bytes_ak(
                encoded_file_name, true);
        auto decoded_ak = aak_decode(encoded_ak);
        auto decoded_rle = rle_decode(decoded_ak);
        auto decoded_mtf = move_to_front_reverse(decoded_rle);
        auto decoded_data = bwt_reverse(decoded_mtf, bwt_shift_position);
        write_bytes(decoded_file_name, decoded_data);
    }
}

const int COUNT = 10;

void print2DUtil(BTree *root, int space) {
    if (root == nullptr) {
        return;
    }
    space += COUNT;
    print2DUtil(root->right, space);
    std::cout << std::endl;
    for (int i = COUNT; i < space; ++i) {
        std::cout << " ";
    }
    std::cout << root->value << std::endl;
    print2DUtil(root->left, space);
}

void print2D(BTree *root) {
    print2DUtil(root, 0);
}

double compression_rate = 0;
uint32_t initial_data_size_sum = 0;
uint32_t encoded_data_size_sum = 0;

void print_metrics(
        const std::string &output_file,
        const size_t &initial_data_size,
        const size_t &encoded_data_size
) {
    initial_data_size_sum += initial_data_size;
    encoded_data_size_sum += encoded_data_size;

    compression_rate += (double(encoded_data_size)) / (double(initial_data_size));
    std::cout << "file_name: " << output_file
              << " $$ initial_data_size: " << initial_data_size
              << " $$ encoded_file_size: " << encoded_data_size
              << " $$ bits_avg: " << (8 * double(encoded_data_size)) / (double(initial_data_size))
              << " $$ compress_rate = " << (double(encoded_data_size)) / (double(initial_data_size))
            << std::endl;
}

int main(int argc, char *argv[]) {
#ifdef FULL_PIPELINE
//    std::string dir = "calgarycorpus/";
//    std::vector<std::string> file_list = {"bib", "book1", "book2", "geo", "news", "obj1", "obj2", "paper1", "paper2",
//                                          "pic", "progc", "progl", "progp", "trans"};
//    std::string dir = "jpeg30/";
//    std::vector<std::string> file_list = {"airplane30.jpg", "arctichare30.jpg", "baboon30.jpg",
//                                          "cat30.jpg", "fruits30.jpg", "frymire30.jpg",
//                                          "girl30.jpg", "lena30.jpg", "monarch30.jpg",
//                                          "peppers30.jpg", "pool30.jpg", "sails30.jpg", "serrano30.jpg",
//                                          "tulips30.jpg", "watch30.jpg"};

    std::string dir = "jpeg80/";
    std::vector<std::string> file_list = {"airplane80.jpg", "arctichare80.jpg", "baboon80.jpg",
                                          "cat80.jpg", "fruits80.jpg", "frymire80.jpg",
                                          "girl80.jpg", "lena80.jpg", "monarch80.jpg",
                                          "peppers80.jpg", "pool80.jpg", "sails80.jpg", "serrano80.jpg",
                                          "tulips80.jpg", "watch80.jpg"};

    std::string encoding_suffix = ".bzap";
    std::string decoding_suffix = ".decoded";
    size_t counter = 1;
//    testing
//    std::string test_string = "abacabaaaaaaaaaaaaaaaaaacabaaaaaaaaaaaaaaaaaacabaaaaaaaaaaaaaaaaaacabaaaaaaaaaaaaaaaaaaaaaaaaaaa";
//    std::vector<unsigned char> bytes_input;
//    for (unsigned char c: test_string) {
//        bytes_input.push_back(c);
//    }
//    const auto &[bytes_input, dummy1, dummy2, dummy3] = read_bytes("jpeg30/pool30.jpg");
//    print_vector(bytes_input);
//    std::cout << std::endl;
//    auto encoded = aak_encode(bytes_input);
//    print_vector(encoded);
//    std::cout << std::endl;
//    auto decoded = aak_decode(encoded);
//    print_vector(decoded);
//    std::cout << std::endl;
//    std::cout << bytes_input.size() << " " << encoded.size() << " " << decoded.size() << std::endl;
//    return 0;
    bool all_tests_passed = true;
    for (auto &current_file: file_list) {
        std::cout << counter << "/" << file_list.size() << ' ';
        ++counter;

        std::string initial_file_name, encoded_file_name, decoded_file_name;
        initial_file_name.append(dir + current_file);
        encoded_file_name.append(dir).append(current_file + encoding_suffix);
        decoded_file_name.append(dir).append(current_file + decoding_suffix);

        compress_ak(initial_file_name, encoded_file_name);
        decompress_ak(encoded_file_name, decoded_file_name);

        bool operation_succeeded = compare_files(initial_file_name, decoded_file_name);
        all_tests_passed &= operation_succeeded;
        std::cout << (operation_succeeded ? "success" : "fail") << std::endl;
    }
    std::cout << (all_tests_passed ? "done, all tests passed " : "!!️ tests failed !!") << " $$ compression rate = " << compression_rate / 15 << '\n';
    std::cout << "encoded_data_size_sum = " << encoded_data_size_sum << '\n';
    std::cout << "initial_data_size_sum = " << initial_data_size_sum << '\n';
    std::cout << "rate = " << (1 - (double) encoded_data_size_sum / (double) initial_data_size_sum) * 100;
#endif
#ifdef COMPRESS
    if (argc != 3) {
        std::cout << "Wrong arguments. Pass only input and output file as parameters";
        return 1;
    }
    std::string input_file = argv[1];
    std::string output_file = argv[2];
    compress_ak(input_file, output_file);
#endif
#ifdef DECOMPRESS
    if (argc != 3) {
        std::cout << "Wrong arguments. Pass only input and output file as parameters";
        return 1;
    }
    std::string input_file = argv[1];
    std::string output_file = argv[2];
    decompress_ak(input_file, output_file);
#endif
}
