#pragma once

#include <cstdint>
#include <vector>

struct linear_sieve {
    std::vector<bool> is_prime;
    std::vector<std::int32_t> min_prime;
    std::vector<std::int32_t> primes;

    linear_sieve(const std::size_t n)
        : is_prime(n, true), min_prime(n), primes() {
        is_prime[0] = false;
        is_prime[1] = false;
        // TODO primes.reserve(~aprox)
        for (std::size_t i = 2; i < n; i++) {
            if (is_prime[i]) primes.push_back(i);
            for (auto p : primes) {
                if (p * i >= n) break;
                is_prime[p * i] = false;
                min_prime[p * i] = p;
                if (i % p == 0) { break; }
            }
        }
    }
};
