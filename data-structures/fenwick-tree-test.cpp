#include "fenwick-tree.hpp"

#include <algorithm>
#include <iostream>

template <typename T>
class SlowRangeQueriesSolver {
   public:
    SlowRangeQueriesSolver(std::size_t size, std::function<T(T, T)> merge)
        : values_(size), merge_(merge) {}
    void Update(std::size_t index, const T& delta) {
        assert(0 <= index && index < values_.size());
        values_[index] = merge_(values_[index], delta);
    }
    T RangeQuery(std::size_t from, std::size_t to) const {
        assert(0 <= from && from < to && to <= values_.size());
        T result{};
        for (std::size_t i = from; i < to; i++) {
            result = merge_(result, values_[i]);
        }
        return result;
    }

   private:
    std::vector<T>         values_;
    std::function<T(T, T)> merge_;
};

void Test() {
    const std::size_t size = 100;

    auto add = [](const auto& lhs, const auto& rhs) { return lhs + rhs; };

    SlowRangeQueriesSolver<std::int64_t> slow(size, add);
    FenwickTree<std::int64_t>            fast(size, add);

    const std::int32_t steps    = 1000;
    const std::int64_t maxDelta = 1000 * 1000;
    for (std::int32_t i = 0; i < steps; i++) {
        if (rand() & 1) {
            std::size_t  index = rand() % size;
            std::int64_t delta = rand() % maxDelta;
            slow.Update(index, delta);
            fast.Update(index, delta);
        } else {
            std::size_t from = rand() % size;
            std::size_t to   = rand() % size;
            if (to <= from) { std::swap(from, to); }
            auto slowAnwer  = slow.RangeQuery(from, to + 1);
            auto fastAnswer = fast.RangeQuery(from, to + 1);
            if (slowAnwer != fastAnswer) {
                std::cout << "ERROR: expected" << slowAnwer << ", found"
                          << fastAnswer << std::endl;
                return;
            }
        }
    }
    std::cout << "Test Passed\n";
}

int main() {
    Test();
    return EXIT_SUCCESS;
}