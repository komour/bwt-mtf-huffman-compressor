//
// Created by Andrey Komarov on 12/13/21.
//

#pragma once

#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "io_utilities.h"

const uint16_t NUMBER_OF_CHARS = 256;
const uint16_t END_OF_FILE = NUMBER_OF_CHARS + 1;
const uint16_t MAX_CODE_VALUE = (1 << 16) - 1;
const uint16_t QUARTER = MAX_CODE_VALUE / 4 + 1;

std::vector<uint16_t> byte_to_index;
std::vector<unsigned char> index_to_byte;
std::vector<uint16_t> freq;
std::vector<uint16_t> cumu;

void update(uint16_t byte_index) {
    if (cumu[0] == QUARTER - 1) {
        uint16_t current_cum = 0;
        for (int i = NUMBER_OF_CHARS + 1; i >= 0; --i) {
            freq[i] = (freq[i] + 1) / 2;
            cumu[i] = current_cum;
            current_cum += freq[i];
        }
    }
    uint16_t new_byte_index = byte_index;
    while (freq[new_byte_index] == freq[new_byte_index - 1]) {
        --new_byte_index;
    }

    if (new_byte_index < byte_index) {
        unsigned char new_byte = index_to_byte[new_byte_index];
        unsigned char old_byte = index_to_byte[byte_index];

        index_to_byte[new_byte_index] = old_byte;
        index_to_byte[byte_index] = new_byte;

        byte_to_index[new_byte] = byte_index;
        byte_to_index[old_byte] = new_byte_index;
    }

    ++freq[new_byte_index];
    uint16_t index_to_update = new_byte_index;
    while (index_to_update) {
        --index_to_update;
        cumu[index_to_update] += 1;
    }
}

void init_frequency() {
    byte_to_index = std::vector<uint16_t>(NUMBER_OF_CHARS);
    index_to_byte = std::vector<unsigned char>(NUMBER_OF_CHARS + 1 + 1);
    freq = std::vector<uint16_t>(NUMBER_OF_CHARS + 1 + 1);
    cumu = std::vector<uint16_t>(NUMBER_OF_CHARS + 1 + 1);


    for (size_t i = 0; i < NUMBER_OF_CHARS; ++i) {
        index_to_byte[i + 1] = i;
        byte_to_index[i] = i + 1;

        freq[i] = 1;
        cumu[i] = NUMBER_OF_CHARS + 1 - i;
    }
    freq[NUMBER_OF_CHARS] = 1;
    freq[NUMBER_OF_CHARS + 1] = 1;
    cumu[NUMBER_OF_CHARS] = 1;
    cumu[NUMBER_OF_CHARS + 1] = 0;

    freq[0] = 0;
}

void write_bit_with_free_bits(
        std::vector<unsigned char> &encoded_data,
        const bool &bit,
        size_t &first_free_pos,
        uint16_t &free_bits
) {
    append_bit(encoded_data, first_free_pos, bit);
    bool oppositeBit = !bit;
    while (free_bits) {
        append_bit(encoded_data, first_free_pos, oppositeBit);
        --free_bits;
    }
}

void encode_next_byte(
        const uint16_t &byte_index,
        uint16_t &right,
        uint16_t &left,
        std::vector<unsigned char> &encoded_data,
        size_t &first_free_pos,
        uint16_t &free_bits
) {
    uint32_t range = 1 + right - left;
    uint16_t total = cumu[0];
    uint16_t byte_left = cumu[byte_index];
    uint16_t byte_right = cumu[byte_index - 1];
    right = left + range * byte_right / total - 1;
    left = left + range * byte_left / total;

    for (;;) {
        if (right < 2 * QUARTER) {
            write_bit_with_free_bits(encoded_data, false, first_free_pos, free_bits);
        } else if (left >= 2 * QUARTER) {
            write_bit_with_free_bits(encoded_data, true, first_free_pos, free_bits);
            left -= 2 * QUARTER;
            right -= 2 * QUARTER;
        } else if (left >= QUARTER && right < 3 * QUARTER) {
            ++free_bits;
            left -= QUARTER;
            right -= QUARTER;
        } else break;
        left = 2 * left;
        right = 2 * right + 1;
    }
}

std::vector<unsigned char> aak_encode(
        const std::vector<unsigned char> &initial_data
) {
    init_frequency();
    std::vector<unsigned char> encoded_data = std::vector<unsigned char>(1);
    uint16_t left = 0;
    uint16_t right = MAX_CODE_VALUE;
    uint16_t free_bits = 0;
    size_t first_free_pos = 7;

    for (unsigned char byte: initial_data) {
        uint16_t byte_index = byte_to_index[byte];
        encode_next_byte(byte_index, right, left, encoded_data, first_free_pos, free_bits);
        update(byte_index);
    }
    encode_next_byte(END_OF_FILE, right, left, encoded_data, first_free_pos, free_bits);
    ++free_bits;

    bool last_bit;
    left < QUARTER ? last_bit = false : last_bit = true;
    write_bit_with_free_bits(encoded_data, last_bit, first_free_pos, free_bits);
    return encoded_data;
}

uint16_t next_byte_index(
        uint16_t &left,
        uint16_t &right,
        uint16_t &code_value,
        size_t &byte_index_io,
        size_t &bit_index,
        const std::vector<unsigned char> &encoded_data
) {
    uint32_t range = right - left + 1;
    uint16_t total = cumu[0];

    uint16_t byte_index = 1;
    while (cumu[byte_index] > ((code_value - left + 1) * total - 1) / range) {
        ++byte_index;
    }

    uint16_t byte_left = cumu[byte_index];
    uint16_t byte_right = cumu[byte_index - 1];
    right = left + range * byte_right / total - 1;
    left = left + range * byte_left / total;

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
        bool next_bit = read_bit(encoded_data, byte_index_io, bit_index);
        code_value *= 2;
        code_value += next_bit;
    }

    return byte_index;
}

std::vector<unsigned char> aak_decode(
        const std::vector<unsigned char> &encoded_data
) {
    init_frequency();
    std::vector<unsigned char> initial_data = std::vector<unsigned char>(1);
    uint16_t left = 0;
    uint16_t right = MAX_CODE_VALUE;
    uint16_t code_value = 0;

    size_t byte_index_io = 0;
    size_t bit_index = 0;
    size_t first_free_pos = 7;

    for (size_t i = 0; i < 16; ++i) {
        bool next_bit = read_bit(encoded_data, byte_index_io, bit_index);
        code_value *= 2;
        code_value += next_bit;
    }

    for (;;) {
        uint16_t byte_index = next_byte_index(left, right, code_value, byte_index_io, bit_index, encoded_data);
        if (byte_index == END_OF_FILE) {
            break;
        }
        unsigned char byte = index_to_byte[byte_index];
        append_byte(initial_data, first_free_pos, byte);
        update(byte_index);
    }
    return initial_data;
}