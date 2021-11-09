#pragma once

#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

template<class T>
void print_vector(const std::vector<T> &vec) {
    for (const auto &el: vec) {
        std::cout << el;
    }
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

std::string char_to_bin(unsigned char c) {
    static char bin[CHAR_BIT + 1] = {0};
    int i;
    for (i = CHAR_BIT - 1; i >= 0; i--) {
        bin[i] = (c % 2) + '0';
        c /= 2;
    }
    return bin;
}

void print_binary_vector(const std::vector<unsigned char> &vector) {
    for (const auto &element: vector) {
        std::cout << char_to_bin(element);
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