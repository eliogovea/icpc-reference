#pragma once

#include <cassert>
#include <cstdint>
#include <iostream>

template <typename I>
I egcd(I a, I b, I& x, I& y) {
    if (a == I {0}) {
        x = I {0};
        y = I {1};
        return b;
    }
    I g = egcd(b % a, a, y, x);
    x -= (b / a) * y;
    return g;
}

template <typename T, T modulo_>
struct int_mod_t {
    constexpr T modulo = modulo_;
    static_assert(std::is_integral<T>::value);
    // static_assert(sizeof(T) <= 8);
    static_assert(modulo > 0);

    T value;

    int_mod_t(T value_ = T {0}) : value(value_ % modulo) { fix(); }

    void fix() {
        if (value >= modulo) { value -= modulo; }
        if (value < 0) { value += modulo; }
    }

    int_mod_t inverse() const {
        std::int64_t x, y;
        auto v = static_cast<std::int64_t>(value);
        auto m = static_cast<std::int64_t>(modulo);

        auto g = egcd(v, m, x, y);
        assert(g == 1);

        x %= modulo;
        if (x < 0) { x += modulo; }

        assert((v * x) % m == 1);
        return {static_cast<T>(x)};
    }

    static T multiply(T lhs, T rhs) {
        static constexpr auto limit32 = 2 * 1000 * 1000 * 1000;
        static constexpr auto limit64 =
          2LL * 1000LL * 1000LL * 1000LL * 1000LL * 1000LL * 1000LL;
        if constexpr (sizeof(T) <= 4 || modulo <= limit32) {
            return static_cast<std::int64_t>(lhs) * rhs % modulo;
        }
        if constexpr (modulo < limit64) {
            std::int64_t q = static_cast<std::int64_t>(
              static_cast<long double>(lhs) * rhs / modulo);
            std::int64_t r = lhs * rhs - q * modulo;
            r %= modulo;
            if (r < 0) { r += modulo; }
            return r;
        }
        return static_cast<__int128_t>(lhs) * rhs % modulo;
    }

    friend std::istream& operator>>(std::istream& input, int_mod_t& x) {
        input >> x.value;
        x.fix();
        return input;
    }

    friend std::ostream& operator<<(std::ostream& output, const int_mod_t& x) {
        output << x.value;
        return output;
    }

    int_mod_t& operator+=(const int_mod_t& other) {
        value += other.value;
        fix();
        return *this;
    }

    int_mod_t& operator-=(const int_mod_t& other) {
        value -= other.value;
        fix();
        return *this;
    }

    int_mod_t& operator*=(const int_mod_t& other) {
        value = multiply(value, other.value);
    }

    int_mod_t& operator/=(const int_mod_t& other) {
        assert(other.value != 0);
        *this *= other.inverse();
        return *this;
    }

    int_mod_t operator+(const int_mod_t& other) {
        return {value + other.value};
    }

    int_mod_t operator-(const int_mod_t& other) {
        return {value - other.value};
    }

    int_mod_t operator*(const int_mod_t& other) {
        return {multiply(value, other.value)};
    }

    int_mod_t operator/(const int_mod_t& other) {
        assert(other.value != 0);
        return *this * other.inverse();
    }

    static int_mod_t power(int_mod_t x, std::int64_t n) {
        int_mod_t y {1};
        while (n > 0) {
            if (n & 1) { y *= x; }
            x *= x;
            n >>= 1;
        }
        return y;
    }
};

template <std::int32_t modulo>
using int32_mod_t = int_mod_t<std::int32_t, modulo>;

template <std::int64_t modulo>
using int64_mod_t = int_mod_t<std::int64_t, modulo>;
