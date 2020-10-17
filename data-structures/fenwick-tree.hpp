#pragma once

#include <cassert>
#include <functional>
#include <vector>

template <typename T>
class FenwickTree {
   public:
    explicit FenwickTree(std::size_t size, std::function<T(T, T)> merge)
        : values_(size, T{}), merge_{merge} {}
    explicit FenwickTree(std::size_t size)
        : FenwickTree(size,
                      [](const T& lhs, const T& rhs) { return lhs + rhs; }) {}
    void Update(std::size_t index, const T& delta) {
        assert(0 <= index && index < values_.size());
        index++;
        while (index <= values_.size()) {
            values_[index - 1] = merge_(values_[index - 1], delta);
            index += (index & -index);
        }
    }
    T PrefixQuery(std::size_t size) const {
        T result{};
        while (size > 0) {
            result = merge_(result, values_[size - 1]);
            size -= (size & -size);
        }
        return result;
    }
    T RangeQuery(std::size_t begin, std::size_t end) const {
        return PrefixQuery(end) - PrefixQuery(begin);  // [begin, end)
    }

   private:
    std::vector<T>                       values_;
    std::function<T(const T&, const T&)> merge_;
};

// example:
// FenwickTree<int> f(10, [] (int a, int b) {
//     return a + b;
// });