#include <iostream>
#include <vector>
#include <numeric>

bool comp(const std::pair<size_t, char> &left, const std::pair<size_t, char> &right) {
    return left.second < right.second;
}

class bwt_cmp {
    const std::string &data;
public:
    bwt_cmp(const std::string &bwt_string) : data(bwt_string) {}

    bool operator()(size_t left, size_t right) {
        return data[left] < data[right];
    }
};

std::string bwt_reverse(const std::string &bwt_string, size_t row_index) {
    std::vector<size_t> l_shift(bwt_string.size());
    std::iota(l_shift.begin(), l_shift.end(), 0);
    std::stable_sort(l_shift.begin(), l_shift.end(), bwt_cmp(bwt_string));

    std::string initial_string(bwt_string.size(), '0');
    for (std::size_t i = 0; i < bwt_string.size(); ++i) {
        std::cout << "row_index = " << row_index << " | l_shift[" << row_index << "] = " << l_shift[row_index]
                  << " | bwt_string[" << l_shift[row_index] << "] = " << bwt_string[l_shift[row_index]] << '\n';
        initial_string[i] = bwt_string[l_shift[row_index]];
        row_index = l_shift[row_index];
    }
    std::cout << "\n\n\n";
    return initial_string;

}

int main() {
    std::string initial_string = "abrakadabra";
    std::string bwt_string = "rdakraaaabb";
    size_t row_index = 2;

    std::cout << bwt_reverse(bwt_string, row_index) << std::endl;
    return 0;
}
