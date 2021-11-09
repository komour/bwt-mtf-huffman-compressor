#include "debug_utilities.h"
#include "io_utilities.h"

#include <iostream>
#include <vector>
#include <numeric>
#include <queue>
#include <unordered_map>

const size_t BYTE_SIZE = 256;

struct BTree {
    unsigned char value;

    BTree *left;
    BTree *right;

    [[nodiscard]] bool has_value() const {
        return !this->left && !this->right;
    }

    BTree() : left(nullptr), right(nullptr), value(0) {}

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
    for (std::size_t i = 0; i < bwt_data.size(); ++i) {
        initial_data[i] = bwt_data[l_shift[row_index]];
        row_index = l_shift[row_index];
    }
    return initial_data;
}

std::pair<size_t, std::vector<unsigned char>> bwt(
        std::vector<unsigned char> data
) {
    std::vector<std::size_t> shift_order(data.size());
    std::iota(shift_order.begin(), shift_order.end(), 0);
    std::stable_sort(shift_order.begin(), shift_order.end(), bwt_cmp_straight(data));

    std::vector<unsigned char> encoded(data.size());
    std::size_t shift_position;
    for (std::size_t i = 0; i < data.size(); ++i) {
        encoded[i] = data[cyclic_index(shift_order[i], data.size() - 1, data.size())];
        if (shift_order[i] == 0) shift_position = i;
    }
    return std::make_pair(shift_position, encoded);
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
    left_codeword.push_back(0);
    right_codeword.push_back(1);
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
        std::unordered_map<unsigned char, std::vector<bool>> &code_words
) {
    auto encoded_data = std::vector<unsigned char>(1);
    size_t first_free_pos = 7;

    for (const auto &byte: data) {
        const auto &code_word = code_words[byte];
        for (bool bit: code_word) {
            append_bit(encoded_data, first_free_pos, bit);
        }
    }
    return encoded_data;
}

std::pair<std::vector<unsigned char>, std::unordered_map<unsigned char, std::vector<bool>>>
huffman(const std::vector<unsigned char> &data) {
    std::vector<long> frequencies(BYTE_SIZE, 0);
    std::priority_queue<std::pair<long, BTree *>> vertex_pool;

    std::vector<bool> already_in_queue(BYTE_SIZE);
    for (unsigned char c: data) {
        ++frequencies[c];
    }

    for (unsigned char byte: data) {
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
    auto code_words = build_hashmap(vertex_pool.top().second);
    auto encoded_data = encode_with_huffman(data, code_words);
    return std::make_pair(encoded_data, code_words);
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
        std::unordered_map<K, V> &map
) {
    std::unordered_map<V, K> new_map;
    for (const auto &pair: map) {
        new_map[pair.second] = pair.first;
    }
    return new_map;
}

void compress(
        const std::string &initial_file_name,
        const std::string &encoded_file_name
) {
    auto bytes_input = read_bytes(initial_file_name).second;
    auto bwt_result = bwt(bytes_input);
    auto bwt_data = bwt_result.second;
    auto bwt_shift_position = bwt_result.first;

    auto mtf_data = move_to_front(bwt_data);

    write_bytes(encoded_file_name, mtf_data, bwt_shift_position);
}

void decompress(
        const std::string &encoded_file_name,
        const std::string &decoded_file_name
) {
    auto result_input = read_bytes(encoded_file_name, true);
    auto encoded_bytes_input = result_input.second;
    size_t bwt_shift_position = result_input.first;

    auto decoded_mtf = move_to_front_reverse(encoded_bytes_input);
    auto decoded_data = bwt_reverse(decoded_mtf, bwt_shift_position);

    write_bytes(decoded_file_name, decoded_data);
}

void full_pipeline(
        const std::string &initial_file_name,
        const std::string &encoded_file_name,
        const std::string &decoded_file_name
) {
    auto bytes_input = read_bytes(initial_file_name).second;
    auto bwt_result = bwt(bytes_input);
    auto bwt_data = bwt_result.second;
    auto bwt_shift_position = bwt_result.first;

    auto mtf_data = move_to_front(bwt_data);
    const size_t initial_data_size = mtf_data.size();

    auto huffman_result = huffman(mtf_data);
    auto code_words = huffman_result.second;
    const auto encoded_huffman = huffman_result.first;
    write_bytes(encoded_file_name, encoded_huffman, bwt_shift_position);

    auto reversed_code_words = reverse_map(code_words);
    auto decoded_huffman = huffman_reverse(encoded_huffman, reversed_code_words, initial_data_size);


    auto decoded_mtf = move_to_front_reverse(decoded_huffman);
    auto decoded_data = bwt_reverse(decoded_mtf, bwt_shift_position);
    write_bytes(decoded_file_name, decoded_data);
}

int main() {
    std::string dir = "calgarycorpus/";
    std::vector<std::string> file_list = {"bib", "book1", "book2", "geo", "news", "obj1", "obj2", "paper1", "paper2",
                                          "pic", "progc", "progl", "progp", "trans"};
    std::string encoding_suffix = ".bzap";
    std::string decoding_suffix = ".decoded";
    size_t counter = 1;

    for (auto &current_file: file_list) {
        std::cout << counter << "/" << file_list.size() << ' ';
        ++counter;

        std::string initial_file_name, encoded_file_name, decoded_file_name;
        initial_file_name.append(dir + current_file);
        encoded_file_name.append(dir).append(current_file + encoding_suffix);
        decoded_file_name.append(dir).append(current_file + decoding_suffix);

//        compress(initial_file_name, encoded_file_name);
//        decompress(encoded_file_name, decoded_file_name);
        full_pipeline(initial_file_name, encoded_file_name, decoded_file_name);

        std::cout << compare_files(initial_file_name, decoded_file_name) << std::endl;
    }
    return 0;
}
