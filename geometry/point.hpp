#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <vector>

template <typename T>
constexpr std::int32_t sign(T x) {
    if constexpr (std::is_floating_point<T>::value) {
        static constexpr double epsilon =
          1e-9;  // TODO change epsilon depending on the problem
        return x < -epsilon ? -1 : (epsilon < x);
    }
    return x < T {0} ? -1 : (T {0} < x);
}

template <typename T>
struct point {
    T x, y;

    point(const T& x_ = {0}, const T& y_ = {0}) : x(x_), y(y_) {}

    template <typename T_>
    explicit point(const point<T_>& P)
        : x(static_cast<T>(P.x)), y(static_cast<T>(P.y)) {}

    template <typename T_>
    explicit operator point<T_>() const {
        return {static_cast<T_>(x), static_cast<T_>(y)};
    }

    template <typename U, typename V>
    explicit point(const std::pair<U, V>& P)
        : x(static_cast<T>(P.first)), y(static_cast<T>(P.second)) {}

    template <typename T_>
    explicit operator std::pair<T_, T_>() const {
        return std::pair<T_, T_> {static_cast<T_>(x), static_cast<T_>(y)};
    }

    friend std::istream& operator>>(std::istream& input, point& P) {
        return input >> P.x >> P.y;
    }

    friend std::ostream& operator<<(std::ostream& output, point& P) {
        output << "(" << P.x << ", " << P.y << ")";
        return output;
    }

    void normalize() {
        if constexpr (std::is_integral<T>::value) {
            auto g = std::__gcd(std::abs(x), std::abs(y));
            x /= g;
            y /= g;
        }
        if constexpr (std::is_floating_point<T>::value) {
            auto l = length();
            assert(sign(l) != 0);
            x /= l;
            y /= l;
        }
    }

    point& operator+=(const point& P) {
        x += P.x;
        y += P.y;
        return *this;
    }

    point& operator-=(const point& P) {
        x -= P.x;
        y -= P.y;
        return *this;
    }

    point& operator*=(const T& k) {
        x *= k;
        y *= k;
        return *this;
    }

    point& operator/=(const T& k) {
        x /= k;
        y /= k;
        return *this;
    }

    point operator+(const point& P) const { return {x + P.x, y + P.y}; }

    point operator-(const point& P) const { return {x - P.x, y - P.y}; }

    point operator*(const T& k) const { return {x * k, y * k}; }

    point operator/(const T& k) const { return {x / k, y / k}; }

    bool operator==(const point& P) const {
        return (sign(x - P.x) == 0 && sign(x - P.y) == 0);
    }

    T norm() const { return dot(*this, *this); }

    auto length() const -> double { return sqrt(static_cast<double>(norm())); }

    static auto dot(const point& P, const point& Q) -> T {
        return P.x * Q.x + P.y * Q.y;
    }

    static auto cross(const point& P, const point& Q) -> T {
        return P.x * Q.y - P.y * Q.x;
    }

    static auto half_plane(const point& P) -> int const {
        if (P.y != 0) { return P.y < 0 ? -1 : +1; }
        return P.x < 0 ? -1 : +1;
    }

    static auto xy_comparator() {
        return [](const point& P, const point& Q) {
            if (sign(P.x - Q.x) != 0) { return sign(P.x - Q.x) == -1; }
            return sign(P.y - Q.y) == -1;
        };
    }

    static auto yx_comparator() {
        return [](const point& P, const point& Q) {
            if (sign(P.y - Q.y) != 0) { return sign(P.y - Q.y) == -1; }
            return sign(P.x - Q.x) == -1;
        };
    }

    static auto angle_comparator(const point& O) {
        return [O](const point& P, const point& Q) {
            auto hP = half_plane(P);
            auto hQ = half_plane(Q);
            if (hP != hQ) { return hP < hQ; }
            auto c = sign(cross(P - O, Q - O));
            if (c != 0) { return (c > 0); }
            return sign(norm(P - O) - norm(Q - O)) == -1;
        };
    }

    point project(const point& P, const point& Q) const {
        static_assert(std::is_floating_point<T>::value);
        return P + (Q - P) * (dot(Q - P, *this - P) / (Q - P).norm());
    }

    bool inside_segment(const point<T>& A, const point<T>& B) const {
        return sign<T>(point<T>::cross(*this - A, B - A)) == 0 &&
               sign<T>(point<T>::dot(A - *this, B - *this)) <= 0;
    }
};
