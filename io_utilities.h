#pragma once

#include <iostream>
#include <vector>
#include <fstream>

void write_bytes(
        const std::string &file_name,
        const std::vector<unsigned char> &data,
        const size_t bwt_shift_position = SIZE_MAX
) {
    std::ofstream fout(file_name, std::ios::binary);
    if (bwt_shift_position != SIZE_MAX) {
        fout.write(reinterpret_cast<const char *>(&bwt_shift_position), sizeof(size_t));
    }
    fout.write(reinterpret_cast<const char *>(data.data()), static_cast<long>(data.size()));
    fout.close();
}

std::pair<size_t, std::vector<unsigned char>> read_bytes(
        const std::string &file_name,
        const bool read_meta = false
) {
    std::ifstream fin(file_name, std::ios::binary);
    std::pair<size_t, std::vector<unsigned char>> input_result;
    std::vector<unsigned char> bytes((std::istreambuf_iterator<char>(fin)), {});
    if (read_meta) {
        input_result.first = *((std::size_t *) bytes.data());
        input_result.second = std::vector(bytes.begin() + sizeof(size_t), bytes.end());
    } else {
        input_result.first = SIZE_MAX;
        input_result.second = bytes;
    }
    fin.close();
    return input_result;
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

void append_bit(std::vector<unsigned char> &output_data, size_t &first_free_pos, const bool &bit) {
    if (first_free_pos == -1) {
        output_data.push_back(0);
        first_free_pos = 7;
    }
    output_data[output_data.size() - 1] |= bit << first_free_pos;
    --first_free_pos;
}
