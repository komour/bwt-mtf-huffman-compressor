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
    vector<uint16_t> char_to_index = vector<uint16_t>(NUMBER_OF_CHARS);
    vector<uint16_t> index_to_symbol = vector<uint16_t>(NUMBER_OF_SYMBOLS + 1);
    vector<uint16_t> frequencies = vector<uint16_t>(NUMBER_OF_SYMBOLS + 1);
    vector<uint16_t> cumulative = vector<uint16_t>(NUMBER_OF_SYMBOLS + 1);

    A0FrequencyModel() {
        for (size_t i = 0; i < NUMBER_OF_CHARS; ++i) {
            char_to_index[i] = i + 1;
            index_to_symbol[i + 1] = i;
        }
        for (size_t i = 0; i <= NUMBER_OF_SYMBOLS; ++i) {
            frequencies[i] = 1;
            cumulative[i] = NUMBER_OF_SYMBOLS - i;
        }
        frequencies[0] = 0;
    }

    int getTotal() {
        return cumulative[0];
    }

    int getLow(int symbolIndex) {
        return cumulative[symbolIndex];
    }

    int getHigh(int symbolIndex) {
        return cumulative[symbolIndex - 1];
    }


    void halveAllFrequenciesIfExceeded() {
        if (getTotal() == MAX_FREQUENCY) {
            int cum = 0;
            for (int i = NUMBER_OF_SYMBOLS; i >= 0; --i) {
                frequencies[i] = (frequencies[i] + 1) / 2;
                cumulative[i] = cum;
                cum += frequencies[i];
            }
        }
    }

    void update(int symbolIndex) {
        halveAllFrequenciesIfExceeded();
        int newSymbolIndex = updateSymbolIndexIfNeeded(symbolIndex);
        updateFrequencies(newSymbolIndex);
    }

    int updateSymbolIndexIfNeeded(int symbolIndex) {
        int newSymbolIndex = symbolIndex;
        while (frequencies[newSymbolIndex] == frequencies[newSymbolIndex - 1]) {
            newSymbolIndex--;
        }

        if (newSymbolIndex < symbolIndex) {
            int newChar = index_to_symbol[newSymbolIndex];
            int oldChar = index_to_symbol[symbolIndex];

            index_to_symbol[newSymbolIndex] = oldChar;
            index_to_symbol[symbolIndex] = newChar;

            char_to_index[newChar] = symbolIndex;
            char_to_index[oldChar] = newSymbolIndex;
        }
        return newSymbolIndex;
    }

    void updateFrequencies(int symbolIndex) {
        frequencies[symbolIndex]++;
        int indexToUpdate = symbolIndex;
        while (indexToUpdate > 0) {
            indexToUpdate--;
            cumulative[indexToUpdate] += 1;
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
    A0FrequencyModel frequencyModel = A0FrequencyModel();
    int low = 0;
    int high = CODE_VALUE_MAX;
    int bitsToFollow = 0;
    size_t first_free_pos = 7;


    void writeBitWithFollow(bool bit) {
        append_bit(encoded_data, first_free_pos, bit);
        int oppositeBit = bit == 0 ? 1 : 0;
        while (bitsToFollow > 0) {
            append_bit(encoded_data, first_free_pos, oppositeBit);
            bitsToFollow--;
        }
    }


    void encodeNext(int charIndex) {
        int range = high - low + 1;
        int total = frequencyModel.getTotal();
        int charLow = frequencyModel.getLow(charIndex);
        int charHigh = frequencyModel.getHigh(charIndex);
        high = low + range * charHigh / total - 1;
        low = low + range * charLow / total;

        while (true) {
            if (high < CODE_VALUE_HALF) {
                writeBitWithFollow(0);
            } else if (low >= CODE_VALUE_HALF) {
                writeBitWithFollow(1);
                low -= CODE_VALUE_HALF;
                high -= CODE_VALUE_HALF;
            } else if (low >= CODE_VALUE_FIRST_QUARTER && high < CODE_VALUE_THIRD_QUARTER) {
                bitsToFollow++;
                low -= CODE_VALUE_FIRST_QUARTER;
                high -= CODE_VALUE_FIRST_QUARTER;
            } else break;
            low = 2 * low;
            high = 2 * high + 1;
        }
    }

    std::vector<unsigned char> writeEncoded() {
        for (size_t i = 0; i < initial_data.size(); ++i) {
            int charIndex = frequencyModel.char_to_index[initial_data[i]];
            encodeNext(charIndex);
            frequencyModel.update(charIndex);
        }
        encodeNext(EOF_SYMBOL);
        bitsToFollow++;
        if (low < CODE_VALUE_FIRST_QUARTER) {
            writeBitWithFollow(0);
        } else {
            writeBitWithFollow(1);
        }
        return encoded_data;
    }
};

class A0DecoderWriter {
public:
    explicit A0DecoderWriter(std::vector<unsigned char> encoded_data) {
        this->encoded_data = encoded_data;
        for (int i = 0; i < CODE_VALUE_BITS; ++i) {
            moveToTheNextBit();
        }
    }

    std::vector<unsigned char> initial_data = std::vector<unsigned char>(1);
    std::vector<unsigned char> encoded_data;
    A0FrequencyModel frequencyModel = A0FrequencyModel();
    uint16_t low = 0;
    uint16_t high = CODE_VALUE_MAX;
    uint16_t codeValue = 0;

    size_t byte_index = 0;
    size_t bit_index = 0;
    size_t first_free_pos = 7;

    void moveToTheNextBit() {
        int nextBit = read_bit(encoded_data, byte_index, bit_index);
        if (nextBit == -1) {
            std::cout << "LEL\n";
            codeValue = 2 * codeValue;
        } else {
            codeValue = 2 * codeValue + nextBit;
        }
    }

    int findSymbolIndex(int cum) {
        int symbolIndex = 1;
        while (frequencyModel.getLow(symbolIndex) > cum) {
            symbolIndex++;
        }
        return symbolIndex;
    }

    int nextSymbolIndex() {
        int range = high - low + 1;
        int total = frequencyModel.getTotal();
        int cum = ((codeValue - low + 1) * total - 1) / range;

        int symbolIndex = findSymbolIndex(cum);
        int symbolLow = frequencyModel.getLow(symbolIndex);
        int symbolHigh = frequencyModel.getHigh(symbolIndex);
        high = low + range * symbolHigh / total - 1;
        low = low + range * symbolLow / total;

        while (true) {
            if (high < CODE_VALUE_HALF) {
                // do nothing
            } else if (low >= CODE_VALUE_HALF) {
                codeValue -= CODE_VALUE_HALF;
                low -= CODE_VALUE_HALF;
                high -= CODE_VALUE_HALF;
            } else if (low >= CODE_VALUE_FIRST_QUARTER && high < CODE_VALUE_THIRD_QUARTER) {
                codeValue -= CODE_VALUE_FIRST_QUARTER;
                low -= CODE_VALUE_FIRST_QUARTER;
                high -= CODE_VALUE_FIRST_QUARTER;
            } else break;
            low = 2 * low;
            high = 2 * high + 1;
            moveToTheNextBit();
        }

        return symbolIndex;
    }

    vector<unsigned char> writeDecoded() {
        while (true) {
            int symbolIndex = nextSymbolIndex();
            if (symbolIndex == EOF_SYMBOL) break;
            int charr = frequencyModel.index_to_symbol[symbolIndex];

            append_byte(initial_data, first_free_pos, charr);
            frequencyModel.update(symbolIndex);
        }
        return initial_data;
    }
};

std::vector<unsigned char> ak_encode(
        const std::vector<unsigned char> &data
) {

    A0CoderWriter coder = A0CoderWriter(data);
    return coder.writeEncoded();
}

std::vector<unsigned char> ak_decode(
        const std::vector<unsigned char> &data
) {
    std::vector<unsigned char> decoded_data;
    A0DecoderWriter decoder = A0DecoderWriter(data);
    return decoder.writeDecoded();
}