#include <cassert>
#include <cmath>
#include <cstdint>
#include <vector>

struct primes_counter {
    std::int64_t n;
    std::int64_t sqrtn;
    std::vector<std::int64_t> values;  // all [n / x]
    std::vector<std::int64_t> primes_count;
    std::size_t index_of(std::int64_t v) const {
        std::size_t i = (v > sqrtn) ? n / v - 1 : values.size() - v;
        assert(values[i] == v);
        return i;
    }
    std::int64_t get_count(std::int64_t n) const {
        return primes_count[index_of(n)];
    }
    primes_counter(std::int64_t n) : n{n}, values() {
        assert(n > 1);
        sqrtn = static_cast<std::int64_t>(sqrtl(n));
        while ((sqrtn + 1) * (sqrtn + 1) <= n) {
            sqrtn++;
        }
        for (int64_t x = 1; x <= n; x = n / (n / x) + 1) {
            values.push_back(n / x);
        }
        primes_count.resize(values.size());
        for (auto v : values) {
            primes_count[index_of(v)] = v - 1;
        }
        for (int p = 2; p <= sqrtn; p++) {
            if (primes_count[index_of(p)] == primes_count[index_of(p - 1)]) {
                continue;
            }
            for (auto v : values) {
                if (static_cast<std::int64_t>(p) * p > v) {
                    break;
                }
                auto delta = primes_count[index_of(v / p)] -
                             primes_count[index_of(p - 1)];
                // delta is the ammount of numbers that will be removed when you
                // use p in the sieve
                primes_count[index_of(v)] -= delta;
            }
        }
    }
};
