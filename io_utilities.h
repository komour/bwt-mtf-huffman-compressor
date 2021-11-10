#pragma once

#include <iostream>
#include <vector>
#include <fstream>

void write_bytes(
        const std::string &file_name,
        const std::vector<unsigned char> &data,
        const size_t bwt_shift_position = SIZE_MAX,
        const size_t initial_data_size = SIZE_MAX,
        const unsigned long size_of_tree = SIZE_MAX,
        const std::vector<unsigned char> &huffman_tree_encoded = std::vector<unsigned char>()
) {
    std::ofstream fout(file_name, std::ios::binary);
    if (bwt_shift_position != SIZE_MAX) {
        fout.write(reinterpret_cast<const char *>(&bwt_shift_position), sizeof(size_t));
        fout.write(reinterpret_cast<const char *>(&initial_data_size), sizeof(size_t));
        fout.write(reinterpret_cast<const char *>(&size_of_tree), sizeof(unsigned long));
        fout.write(
                reinterpret_cast<const char *>(huffman_tree_encoded.data()),
                static_cast<long>(huffman_tree_encoded.size())
        );
    }
    fout.write(reinterpret_cast<const char *>(data.data()), static_cast<long>(data.size()));
    fout.close();
}

std::tuple<std::vector<unsigned char>, size_t, size_t, std::vector<unsigned char>>
read_bytes(
        const std::string &file_name,
        const bool read_meta = false
) {
    size_t initial_data_size = SIZE_MAX;
    size_t bwt_shift_position = SIZE_MAX;
    long size_of_tree = SIZE_MAX;
    std::vector<unsigned char> data;
    std::vector<unsigned char> huffman_tree_encoded;

    std::ifstream fin(file_name, std::ios::binary);
    std::vector<unsigned char> bytes((std::istreambuf_iterator<char>(fin)), {});
    fin.close();

    if (read_meta) {
        bwt_shift_position = *((size_t *) bytes.data());
        initial_data_size = *((size_t *) bytes.data() + 1);
        size_of_tree = *((long *) bytes.data() + 2);
        auto begin_tree_it = bytes.begin() + 2 * sizeof(size_t) + sizeof(unsigned long);
        huffman_tree_encoded = std::vector(begin_tree_it, begin_tree_it + size_of_tree);
        data = std::vector(begin_tree_it + size_of_tree, bytes.end());
    } else {
        data = bytes;
    }
    return {data, bwt_shift_position, initial_data_size, huffman_tree_encoded};
}

bool get_bit(unsigned char byte, size_t bit_number) {
    return byte & (1 << (bit_number - 1));
}

bool read_bit(
        const std::vector<unsigned char> &data,
        size_t &byte_index,
        size_t &bit_index
) {
    if (bit_index == 8) {
        ++byte_index;
        bit_index = 0;
    }
    auto byte = data[byte_index];
    return get_bit(byte, 8 - bit_index++);
}

unsigned char read_byte(
        const std::vector<unsigned char> &data,
        size_t &byte_index,
        size_t &bit_index
) {
    unsigned char byte = 0u;
    for (size_t i = 0; i < 8; ++i) {
        bool bit = read_bit(data, byte_index, bit_index);
        byte |= bit << (7 - i);
    }
    return byte;
}

void append_bit(std::vector<unsigned char> &output_data, size_t &bit_index, const bool &bit) {
    if (bit_index == -1) {
        output_data.push_back(0);
        bit_index = 7;
    }
    output_data[output_data.size() - 1] |= bit << bit_index;
    --bit_index;
}

void append_byte(std::vector<unsigned char> &output_data, size_t &bit_index, const unsigned char &byte) {
    for (size_t i = 8; i > 0; --i) {
        bool bit = get_bit(byte, i);
        append_bit(output_data, bit_index, bit);
    }
}

void append_byte(std::vector<unsigned char> &output_data, size_t &bit_index, const std::vector<bool> &code_word) {
    for (const auto &bit: code_word) {
        append_bit(output_data, bit_index, bit);
    }
}