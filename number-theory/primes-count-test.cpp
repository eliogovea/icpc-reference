#include <bits/stdc++.h>

using namespace std;

template <typename T>
std::map<T, T> primes_count(const T& n) {
    static_assert(std::is_integral<T>::value);
    std::vector<T> values;  // all [n / x]
    for (T x {1}; x <= n; x = n / (n / x) + 1) values.push_back(n / x);
    std::map<T, T> pi;
    for (const auto& v : values) pi[v] = v - 1;
    for (T p {2}; p * p <= n; p++) {
        if (pi[p] == pi[p - 1]) continue;
        for (const auto& v : values) {
            if (p * p > v) break;
            const auto delta = pi[v / p] - pi[p - 1];
            // delta is the ammount of numbers that will
            // be removed when p is used in the sieve
            pi[v] -= delta;
        }
    }
    return pi;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);

    int64_t n;
    cin >> n;

    auto pi = primes_count(n);

    cout << pi[n] << "\n";
}