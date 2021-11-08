#include <iostream>
#include <vector>
#include <numeric>
#include <fstream>

class bwt_cmp_reverse {
    const std::vector<unsigned char> &data;
public:
    bwt_cmp_reverse(const std::vector<unsigned char> &bwt_data) : data(bwt_data) {}

    bool operator()(size_t left, size_t right) {
        return data[left] < data[right];
    }
};

class bwt_cmp_straight {
    const std::vector<unsigned char> &data;
public:
    bwt_cmp_straight(const std::vector<unsigned char> &bwt_data) : data(bwt_data) {}

    bool operator()(size_t left, size_t right) {
        size_t i = 0;
        while (data[left + i] == data[right + i] && i < data.size()) {
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

void write_bytes(const std::string &file_name, std::vector<unsigned char> &data) {
    std::ofstream fout(file_name, std::ios::binary);
    fout.write(reinterpret_cast<const char *>(data.data()), data.size());
}

std::vector<unsigned char> read_bytes(const std::string &file_name) {
    std::ifstream fin_encoded(file_name, std::ios::binary);
    std::vector<unsigned char> bytes((std::istreambuf_iterator<char>(fin_encoded)), {});
    return bytes;
}

std::pair<size_t, std::vector<unsigned char>> bwt(std::vector<unsigned char> data) {
    std::vector<std::size_t> shift_order(data.size());
    std::iota(shift_order.begin(), shift_order.end(), 0);
    data.insert(data.end(), data.begin(), data.end());
    std::stable_sort(shift_order.begin(), shift_order.end(), bwt_cmp_straight(data));

    std::vector<unsigned char> res(data.size() / 2);
    std::size_t shift_position;
    for (std::size_t i = 0; i < data.size() / 2; ++i) {
        if (shift_order[i] == 0) {
            shift_position = i;
        }
        res[i] = data[shift_order[i] + data.size() / 2 - 1];
    }
    return std::make_pair(shift_position, res);
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

int main() {


    std::string dir = "calgarycorpus/";
    std::string work_dir = "calgarycorpus/";

    std::vector<std::string> file_list = {"bib", "book1", "book2", "geo", "news", "obj1", "obj2", "paper1", "paper2",
                                          "pic", "progc", "progl", "progp", "trans"};

    std::string encoding_suffix = ".bzap";
    std::string decoding_suffix = ".decoded";
    size_t counter = 1;


    for (auto &current_file : file_list) {
        std::cout << counter << "/" << file_list.size() << ' ';
        ++counter;
        std::string initial_file_name = dir + current_file;
        std::string encoded_file_name = work_dir + current_file + encoding_suffix;
        std::string decoded_file_name = work_dir + current_file + decoding_suffix;

        auto bytes_input = read_bytes(initial_file_name);
        auto bwt_result = bwt(bytes_input);
        auto bwt_data = bwt_result.second;
        auto bwt_position = bwt_result.first;
        write_bytes(encoded_file_name, bwt_data);

        auto encoded_bytes_input = read_bytes(encoded_file_name);
        auto decoded_data = bwt_reverse(encoded_bytes_input, bwt_position);
        write_bytes(decoded_file_name, decoded_data);

        std::cout << compare_files(initial_file_name, decoded_file_name) << std::endl;
    }
    return 0;
}
