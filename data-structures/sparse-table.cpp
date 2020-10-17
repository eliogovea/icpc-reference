#include <bits/stdc++.h>

using namespace std;

template <typename T, const std::size_t MaxLogSize>
class SparseTable {
   public:
    template <typename Iterator>
    SparseTable(Iterator begin, Iterator end, std::function<T(T, T)> merge)
        : size_(std::distance(begin, end)), merge_(merge) {
        Build(begin, end);
    }

    SparseTable(std::size_t size, std::function<T(T, T)> merge)
        : size_(size), merge_(merge) {
        std::vector<T> dummy(size);
        Build(dummy.begin(), dummy.end());
    }

    template <typename Container>
    SparseTable(const Container& c, std::function<T(T, T)> merge)
        : size(c.size()), merge_(merge) {
        Build(c.begin(), c.end());
    }

   private:
    template <typename Iterator>
    void Build(Iterator begin, Iterator end) {
        log2_[1] = 0;
        for (std::size_t i = 2; i <= size_; i++) {
            log2_[i] = log2_[i >> 1] + 1;
        }
        for (Iterator i = begin; i != end; i++) {
            values_[0][i - begin] = static_cast<T>(*i);
        }
        for (std::size_t lg = 1; lg <= size_; lg <<= 1) {
            for (std::size_t i = 0; i + (1 << lg) <= size_; i++) {
                values_[lg][i] = merge_(values_[lg - 1][i],
                                        values_[lg - 1][i + (1 << (lg - 1))]);
            }
        }
    }

    T RangeQuery(std::size_t from, std::size_t to) const {
        assert(to <= from);
        assert(0 <= from && to <= size);
        std::size_t index = log2_[to - from];
        return merge_(values_[index][from], values_[index][to - (1 << lg)]);
    }

    std::size_t size_;
    std::function<T(T, T)> merge_;

    static constexpr MaxSize = 1 << MaxLogSize;
    std::array<std::size_t, MaxSize> log2_;
    std::array<std::array<T, MaxSize>, MaxSize> values_;
};