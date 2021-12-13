//
// Created by Andrey Komarov on 12/13/21.
//

#pragma once

#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "io_utilities.h"

using namespace std;

const uint16_t NUMBER_OF_CHARS = 256;
const uint16_t CODE_VALUE_BITS = 16;
const uint16_t CODE_VALUE_MAX = (1 << CODE_VALUE_BITS) - 1;
const uint16_t CODE_VALUE_FIRST_QUARTER = CODE_VALUE_MAX / 4 + 1;
const uint16_t CODE_VALUE_HALF = 2 * CODE_VALUE_FIRST_QUARTER;
const uint16_t CODE_VALUE_THIRD_QUARTER = 3 * CODE_VALUE_FIRST_QUARTER;
const uint16_t MAX_FREQUENCY = CODE_VALUE_FIRST_QUARTER - 1;
const uint16_t EOF_SYMBOL = NUMBER_OF_CHARS + 1;
const uint16_t NUMBER_OF_SYMBOLS = NUMBER_OF_CHARS + 1;


class A0FrequencyModel {
public:
    vector<uint16_t> symbol_to_index = vector<uint16_t>(NUMBER_OF_CHARS);
    vector<uint16_t> index_to_symbol = vector<uint16_t>(NUMBER_OF_SYMBOLS + 1);
    vector<uint16_t> freq = vector<uint16_t>(NUMBER_OF_SYMBOLS + 1);
    vector<uint16_t> cum = vector<uint16_t>(NUMBER_OF_SYMBOLS + 1);

    A0FrequencyModel() {
        for (size_t i = 0; i < NUMBER_OF_CHARS; ++i) {
            symbol_to_index[i] = i + 1;
            index_to_symbol[i + 1] = i;
        }
        for (size_t i = 0; i <= NUMBER_OF_SYMBOLS; ++i) {
            freq[i] = 1;
            cum[i] = NUMBER_OF_SYMBOLS - i;
        }
        freq[0] = 0;
    }

    void halve_all_frequencies_if_exceeded() {
        if (cum[0] == MAX_FREQUENCY) {
            uint16_t current_cum = 0;
            for (int i = NUMBER_OF_SYMBOLS; i >= 0; --i) {
                freq[i] = (freq[i] + 1) / 2;
                cum[i] = current_cum;
                current_cum += freq[i];
            }
        }
    }

    void update(uint16_t symbol_index) {
        halve_all_frequencies_if_exceeded();
        uint16_t new_symbol_index = update_symbol_index_if_needed(symbol_index);
        update_frequencies(new_symbol_index);
    }

    uint16_t update_symbol_index_if_needed(uint16_t symbol_index) {
        uint16_t new_symbol_index = symbol_index;
        while (freq[new_symbol_index] == freq[new_symbol_index - 1]) {
            --new_symbol_index;
        }

        if (new_symbol_index < symbol_index) {
            uint16_t new_symbol = index_to_symbol[new_symbol_index];
            uint16_t old_symbol = index_to_symbol[symbol_index];

            index_to_symbol[new_symbol_index] = old_symbol;
            index_to_symbol[symbol_index] = new_symbol;

            symbol_to_index[new_symbol] = symbol_index;
            symbol_to_index[old_symbol] = new_symbol_index;
        }
        return new_symbol_index;
    }

    void update_frequencies(uint16_t symbol_index) {
        ++freq[symbol_index];
        uint16_t index_to_update = symbol_index;
        while (index_to_update) {
            --index_to_update;
            cum[index_to_update] += 1;
        }
    }
};

class A0CoderWriter {
public:

    A0CoderWriter(std::vector<unsigned char> initial_data) {
        this->initial_data = initial_data;
    }

    std::vector<unsigned char> encoded_data = std::vector<unsigned char>(1);
    std::vector<unsigned char> initial_data;
    A0FrequencyModel frequency_model = A0FrequencyModel();
    uint16_t low = 0;
    uint16_t high = CODE_VALUE_MAX;
    uint16_t bits_to_follow = 0;
    size_t first_free_pos = 7;


    void write_bit_with_follow(bool bit) {
        append_bit(encoded_data, first_free_pos, bit);
        bool oppositeBit = !bit;
        while (bits_to_follow > 0) {
            append_bit(encoded_data, first_free_pos, oppositeBit);
            --bits_to_follow;
        }
    }


    void encode_next(uint16_t symbol_index) {
        uint32_t range = 1 + high - low;
        uint16_t total = frequency_model.cum[0];
        uint16_t symbol_low = frequency_model.cum[symbol_index];
        uint16_t symbol_high = frequency_model.cum[symbol_index - 1];
        high = low + range * symbol_high / total - 1;
        low = low + range * symbol_low / total;

        for (;;) {
            if (high < CODE_VALUE_HALF) {
                write_bit_with_follow(false);
            } else if (low >= CODE_VALUE_HALF) {
                write_bit_with_follow(true);
                low -= CODE_VALUE_HALF;
                high -= CODE_VALUE_HALF;
            } else if (low >= CODE_VALUE_FIRST_QUARTER && high < CODE_VALUE_THIRD_QUARTER) {
                ++bits_to_follow;
                low -= CODE_VALUE_FIRST_QUARTER;
                high -= CODE_VALUE_FIRST_QUARTER;
            } else break;
            low = 2 * low;
            high = 2 * high + 1;
        }
    }

    std::vector<unsigned char> write_encoded() {
        for (unsigned char byte : initial_data) {
            uint16_t symbol_index = frequency_model.symbol_to_index[byte];
            encode_next(symbol_index);
            frequency_model.update(symbol_index);
        }
        encode_next(EOF_SYMBOL);
        ++bits_to_follow;
        if (low < CODE_VALUE_FIRST_QUARTER) {
            write_bit_with_follow(false);
        } else {
            write_bit_with_follow(true);
        }
        return encoded_data;
    }
};

class A0DecoderWriter {
public:
    explicit A0DecoderWriter(std::vector<unsigned char> encoded_data) {
        this->encoded_data = encoded_data;
        for (size_t i = 0; i < CODE_VALUE_BITS; ++i) {
            bool next_bit = read_bit(encoded_data, byte_index, bit_index);
            code_value *= 2;
            code_value += next_bit;
        }
    }

    std::vector<unsigned char> initial_data = std::vector<unsigned char>(1);
    std::vector<unsigned char> encoded_data;
    A0FrequencyModel frequency_model = A0FrequencyModel();
    uint16_t low = 0;
    uint16_t high = CODE_VALUE_MAX;
    uint16_t code_value = 0;

    size_t byte_index = 0;
    size_t bit_index = 0;
    size_t first_free_pos = 7;

    uint16_t next_symbol_index() {
        uint32_t range = high - low + 1;
        uint16_t total = frequency_model.cum[0];
        uint16_t cum = ((code_value - low + 1) * total - 1) / range;

        uint16_t symbol_index = 1;
        while (frequency_model.cum[symbol_index] > cum) {
            ++symbol_index;
        }

        uint16_t symbol_low = frequency_model.cum[symbol_index];
        uint16_t symbol_high = frequency_model.cum[symbol_index - 1];
        high = low + range * symbol_high / total - 1;
        low = low + range * symbol_low / total;

        for (;;) {
            if (high < CODE_VALUE_HALF) {
            } else if (low >= CODE_VALUE_HALF) {
                code_value -= CODE_VALUE_HALF;
                low -= CODE_VALUE_HALF;
                high -= CODE_VALUE_HALF;
            } else if (low >= CODE_VALUE_FIRST_QUARTER && high < CODE_VALUE_THIRD_QUARTER) {
                code_value -= CODE_VALUE_FIRST_QUARTER;
                low -= CODE_VALUE_FIRST_QUARTER;
                high -= CODE_VALUE_FIRST_QUARTER;
            } else {
                break;
            }
            low *= 2;
            high *= 2;
            ++high;
            bool next_bit = read_bit(encoded_data, byte_index, bit_index);
            code_value *= 2;
            code_value += next_bit;
        }

        return symbol_index;
    }

    vector<unsigned char> write_decoded() {
        for (;;) {
            uint16_t symbol_index = next_symbol_index();
            if (symbol_index == EOF_SYMBOL) {
                break;
            }
            uint16_t symbol = frequency_model.index_to_symbol[symbol_index];
            append_byte(initial_data, first_free_pos, symbol);
            frequency_model.update(symbol_index);
        }
        return initial_data;
    }
};

std::vector<unsigned char> ak_encode(
        const std::vector<unsigned char> &data
) {
    A0CoderWriter coder = A0CoderWriter(data);
    std::vector<unsigned char> encoded_data = coder.write_encoded();
    return encoded_data;
}

std::vector<unsigned char> ak_decode(
        const std::vector<unsigned char> &encoded_data
) {
    A0DecoderWriter decoder = A0DecoderWriter(encoded_data);
    std::vector<unsigned char> decoded_data = decoder.write_decoded();
    return decoded_data;
}