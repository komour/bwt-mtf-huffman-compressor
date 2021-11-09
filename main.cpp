#include <iostream>
#include <vector>
#include <numeric>
#include <fstream>
#include <queue>
#include <unordered_map>
#include <climits>

std::string char_to_bin(unsigned char c) {
    static char bin[CHAR_BIT + 1] = {0};
    int i;
    for (i = CHAR_BIT - 1; i >= 0; i--) {
        bin[i] = (c % 2) + '0';
        c /= 2;
    }
    return bin;
}

const size_t BYTE_SIZE = 256;

struct BTree {
    long freq;
    unsigned char value;

    BTree *left;
    BTree *right;

    bool operator<(const BTree &other) const {
        return this->freq < other.freq;
    }

    [[nodiscard]] bool has_value() const {
        return !this->left && !this->right;
    }

    BTree() : left(nullptr), right(nullptr), freq(-1), value(0) {}

    explicit BTree(long freq, unsigned char byte_value = 0, BTree *left = nullptr, BTree *right = nullptr) :
            freq(freq), left(left), right(right), value(byte_value) {}

};

class bwt_cmp_reverse {
    const std::vector<unsigned char> &data;
public:
    explicit bwt_cmp_reverse(const std::vector<unsigned char> &bwt_data) : data(bwt_data) {}

    bool operator()(size_t left, size_t right) {
        return data[left] < data[right];
    }
};

class bwt_cmp_straight {
    const std::vector<unsigned char> &data;
public:
    explicit bwt_cmp_straight(const std::vector<unsigned char> &bwt_data) : data(bwt_data) {}

    bool operator()(size_t left, size_t right) {
        size_t i = 0;
        while (data[left + i] == data[right + i] && i < data.size()) {  // TODO optimise memory
            ++i;
        }
        return data[left + i] < data[right + i];
    }
};

std::vector<unsigned char> bwt_reverse(const std::vector<unsigned char> &bwt_data, size_t row_index) {
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

void write_bytes(const std::string &file_name, const std::vector<unsigned char> &data,
                 const size_t bwt_shift_position = SIZE_MAX) {
    std::ofstream fout(file_name, std::ios::binary);
    if (bwt_shift_position != SIZE_MAX) {
        fout.write(reinterpret_cast<const char *>(&bwt_shift_position), sizeof(size_t));
    }
    fout.write(reinterpret_cast<const char *>(data.data()), static_cast<long>(data.size()));
    fout.close();
}

std::pair<size_t, std::vector<unsigned char>> read_bytes(const std::string &file_name, const bool read_index = false) {
    std::ifstream fin(file_name, std::ios::binary);
    std::pair<size_t, std::vector<unsigned char>> input_result;
    std::vector<unsigned char> bytes((std::istreambuf_iterator<char>(fin)), {});
    if (read_index) {
        input_result.first = *((std::size_t *) bytes.data());
        input_result.second = std::vector(bytes.begin() + sizeof(size_t), bytes.end());
    } else {
        input_result.first = SIZE_MAX;
        input_result.second = bytes;
    }
    fin.close();
    return input_result;
}

std::pair<size_t, std::vector<unsigned char>> bwt(std::vector<unsigned char> data) {
    std::vector<std::size_t> shift_order(data.size());
    std::iota(shift_order.begin(), shift_order.end(), 0);
    data.insert(data.end(), data.begin(), data.end());
    std::stable_sort(shift_order.begin(), shift_order.end(), bwt_cmp_straight(data));

    std::vector<unsigned char> encoded(data.size() / 2);
    std::size_t shift_position;
    for (std::size_t i = 0; i < data.size() / 2; ++i) {
        encoded[i] = data[shift_order[i] + data.size() / 2 - 1];
        if (shift_order[i] == 0) shift_position = i;
    }
    return std::make_pair(shift_position, encoded);
}

// https://stackoverflow.com/a/37575457/12658435
bool compare_files(const std::string &p1, const std::string &p2) {
    std::ifstream f1(p1, std::ifstream::binary | std::ifstream::ate);
    std::ifstream f2(p2, std::ifstream::binary | std::ifstream::ate);

    if (f1.fail() || f2.fail()) {
        return false; //file problem
    }

    if (f1.tellg() != f2.tellg()) {
        return false; //size mismatch
    }

    //seek back to beginning and use std::equal to compare contents
    f1.seekg(0, std::ifstream::beg);
    f2.seekg(0, std::ifstream::beg);
    return std::equal(std::istreambuf_iterator<char>(f1.rdbuf()),
                      std::istreambuf_iterator<char>(),
                      std::istreambuf_iterator<char>(f2.rdbuf()));
}

std::vector<unsigned char> move_to_front(std::vector<unsigned char> data) {
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

std::vector<unsigned char> move_to_front_reverse(std::vector<unsigned char> data) {
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

std::unordered_map<unsigned char, std::vector<bool>> build_hashmap(BTree *root) {
    std::unordered_map<unsigned char, std::vector<bool>> code_words;
    std::vector<bool> initial_codeword;
    traverse(root, initial_codeword, code_words);
    return code_words;
}

template<class T>
void print_vector(const std::vector<T> &vec) {
    for (const auto &el: vec) {
        std::cout << el;
    }
}

void print_binary_vector(const std::vector<unsigned char> &vec) {
    for (const auto &el: vec) {
        std::cout << char_to_bin(el);
    }
}

void print_map(const std::unordered_map<unsigned char, std::vector<bool>> &map) {
    std::cout << "##################\n";
    for (const auto &pair: map) {
        std::cout << pair.first << ' ';
        print_vector(pair.second);
        std::cout << std::endl;
    }
    std::cout << "##################\n";
}

void print_map(const std::unordered_map<std::vector<bool>, unsigned char> &map) {
    std::cout << "##################\n";
    for (const auto &pair: map) {
        print_vector(pair.first);
        std::cout << ' ' << pair.second;
        std::cout << std::endl;
    }
    std::cout << "##################\n";
}

void append_bit(std::vector<unsigned char> &output_data, size_t &first_free_pos, const unsigned char &bit) {
    if (first_free_pos == -1) {
        output_data.push_back(0u);
        first_free_pos = 7;
    }
    output_data[output_data.size() - 1] |= bit << first_free_pos;
    --first_free_pos;
}

std::vector<unsigned char> encode_with_huffman(
        const std::vector<unsigned char> &data,
        std::unordered_map<unsigned char, std::vector<bool>> &code_words
) {
    auto encoded_data = std::vector<unsigned char>(1);
    size_t first_free_pos = 7;

    for (const auto &byte: data) {
        const auto &code_word = code_words[byte];
        for (const auto &bit: code_word) {
            append_bit(encoded_data, first_free_pos, bit);
        }
    }
    return encoded_data;
}

std::pair<std::vector<unsigned char>, std::unordered_map<unsigned char, std::vector<bool>>>
huffman(const std::vector<unsigned char> &data) {
    std::vector<long> frequencies(BYTE_SIZE, 0);
    std::priority_queue<BTree *> vertex_pool;
    std::vector<bool> already_in_set(BYTE_SIZE);
    for (unsigned char c: data) {
        ++frequencies[c];
    }
    for (unsigned char byte: data) {
        if (!already_in_set[byte] && frequencies[byte] != 0) {
            auto new_vertex = new BTree(frequencies[byte], byte);
            vertex_pool.push(new_vertex);
            already_in_set[byte] = true;
        }
    }
    while (vertex_pool.size() != 1) {
        auto left = vertex_pool.top();
        vertex_pool.pop();
        auto right = vertex_pool.top();
        vertex_pool.pop();

        long new_freq = left->freq + right->freq;
        auto new_node = new BTree(new_freq, 0, left, right);
        vertex_pool.push(new_node);
    }
    auto code_words = build_hashmap(vertex_pool.top());
    auto encoded_data = encode_with_huffman(data, code_words);
    return std::make_pair(encoded_data, code_words);
}

bool get_bit(unsigned char byte, size_t bit_number) {
    return byte & (1 << (bit_number - 1));
}

bool read_bit(const std::vector<unsigned char> &data, size_t &byte_index, size_t &bit_index) {
    if (bit_index == 8) {
        ++byte_index;
        bit_index = 0;
    }
    auto byte = data[byte_index];
    return get_bit(byte, 8 - bit_index++);
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
std::unordered_map<V, K> reverse_map(std::unordered_map<K, V> &map) {
    std::unordered_map<V, K> new_map;
    for (const auto &pair: map) {
        new_map[pair.second] = pair.first;
    }
    return new_map;
}

void compress(const std::string &initial_file_name, const std::string &encoded_file_name) {
    auto bytes_input = read_bytes(initial_file_name).second;
    auto bwt_result = bwt(bytes_input);
    auto bwt_data = bwt_result.second;
    auto bwt_shift_position = bwt_result.first;

    auto mtf_data = move_to_front(bwt_data);

    write_bytes(encoded_file_name, mtf_data, bwt_shift_position);
}

void decompress(const std::string &encoded_file_name, const std::string &decoded_file_name) {
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
//    auto huffman_result = huffman(bytes_input);
//    write_bytes(encoded_file_name, mtf_data, bwt_shift_position);
    auto code_words = huffman_result.second;
    const auto encoded_huffman = huffman_result.first;
    write_bytes(encoded_file_name, encoded_huffman, bwt_shift_position);

    auto reversed_code_words = reverse_map(code_words);
    auto decoded_huffman = huffman_reverse(encoded_huffman, reversed_code_words, initial_data_size);


//    auto decoded_mtf = move_to_front_reverse(mtf_data);
    auto decoded_mtf = move_to_front_reverse(decoded_huffman);
    auto decoded_data = bwt_reverse(decoded_mtf, bwt_shift_position);
    write_bytes(decoded_file_name, decoded_data);
//    write_bytes(decoded_file_name, decoded_huffman);
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
