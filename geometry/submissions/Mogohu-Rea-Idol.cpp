// Problem: https://codeforces.com/problemset/problem/87/E

#include <bits/stdc++.h>

constexpr double epsilon = 1e-9;  // TODO change epsilon depending on the problem

template <typename T>
constexpr std::int32_t sign(T x) {
    return x < T {0} ? -1 : (T {0} < x);
}

template <typename T>
struct point {
    T x, y;

    point(const T& x_ = {0}, const T& y_ = {0}) : x(x_), y(y_) {
    }

    template <typename T_>
    explicit point(const point<T_>& P) : x(static_cast<T>(P.x)), y(static_cast<T>(P.y)) {
    }

    template <typename T_>
    explicit operator point<T_>() const {
        return {static_cast<T_>(x), static_cast<T_>(y)};
    }

    template <typename U, typename V>
    explicit point(const std::pair<U, V>& P)
        : x(static_cast<T>(P.first)), y(static_cast<T>(P.second)) {
    }

    template <typename T_>
    explicit operator std::pair<T_, T_>() const {
        return std::pair<T_, T_> {static_cast<T_>(x), static_cast<T_>(y)};
    }

    friend std::istream& operator>>(std::istream& input, point& P) {
        return input >> P.x >> P.y;
    }

    friend std::ostream& operator<<(std::ostream& output, const point& P) {
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

    point operator+(const point& P) const {
        return {x + P.x, y + P.y};
    }

    point operator-(const point& P) const {
        return {x - P.x, y - P.y};
    }

    point operator*(const T& k) const {
        return {x * k, y * k};
    }

    point operator/(const T& k) const {
        return {x / k, y / k};
    }

    bool operator==(const point& P) const {
        return (sign(x - P.x) == 0 && sign(y - P.y) == 0);
    }

    static auto norm(const point& P) -> T {
        return dot(P, P);
    }

    static auto length(const point& P) -> double {
        return sqrt(static_cast<double>(norm(P)));
    }

    static auto dot(const point& P, const point& Q) -> T {
        return P.x * Q.x + P.y * Q.y;
    }

    static auto cross(const point& P, const point& Q) -> T {
        return P.x * Q.y - P.y * Q.x;
    }

    static auto half_plane(const point& P) -> int const {
        if (P.y != 0) {
            return P.y < 0 ? -1 : +1;
        }
        return P.x < 0 ? -1 : +1;
    }

    static auto xy_comparator() {
        return [](const point& P, const point& Q) {
            if (sign(P.x - Q.x) != 0) {
                return sign(P.x - Q.x) == -1;
            }
            return sign(P.y - Q.y) == -1;
        };
    }

    static auto yx_comparator() {
        return [](const point& P, const point& Q) {
            if (sign(P.y - Q.y) != 0) {
                return sign(P.y - Q.y) == -1;
            }
            return sign(P.x - Q.x) == -1;
        };
    }

    static auto angle_comparator(const point& O) {
        return [O](const point& P, const point& Q) {
            const auto hP = half_plane(P - O);
            const auto hQ = half_plane(Q - O);
            if (hP != hQ) {
                return hP > hQ;
            }
            const auto c = sign(cross(P - O, Q - O));
            if (c != 0) {
                return (c > 0);
            }
            return sign(norm(P - O) - norm(Q - O)) == -1;
        };
    }

    point project(const point& P, const point& Q) const {
        static_assert(std::is_floating_point<T>::value);
        return P + (Q - P) * (dot(Q - P, *this - P) / (Q - P).norm());
    }

    bool inside_segment(const point<T>& A, const point<T>& B) const {
        return sign<T>(point<T>::cross(*this - A, B - A)) == 0
               && sign<T>(point<T>::dot(A - *this, B - *this)) <= 0;
    }
};

enum class relative_position : std::int32_t { in = -1, on = 0, out = 1 };

template <typename T>
relative_position inside_convex(point<T> P, const std::vector<point<T>>& G) {
    if (sign(point<T>::cross(P - G[0], G[1] - G[0])) > 0) {
        return relative_position::out;
    }

    if (sign(point<T>::cross(G.back() - G[0], P - G[0])) > 0) {
        return relative_position::out;
    }

    if (sign(point<T>::cross(P - G[0], G[1] - G[0])) == 0) {
        return P.inside_segment(G[0], G[1]) ? relative_position::on : relative_position::out;
    }

    if (sign(point<T>::cross(G.back() - G[0], P - G[0])) == 0) {
        return P.inside_segment(G[0], G.back()) ? relative_position::on : relative_position::out;
    }

    std::size_t lo = 2;
    std::size_t hi = G.size() - 1;
    std::size_t pos = hi;

    while (lo <= hi) {
        int mid = (lo + hi) >> 1;
        if (sign(point<T>::cross(P - G[0], G[mid] - G[0])) >= 0) {
            pos = mid;
            hi = mid - 1;
        } else {
            lo = mid + 1;
        }
    }

    const auto c = sign(point<T>::cross(G[pos] - G[pos - 1], P - G[pos - 1]));

    if (c > 0) {
        return relative_position::in;
    }

    if (c < 0) {
        return relative_position::out;
    }

    return relative_position::on;
}

template <typename T>
void normalize_convex(std::vector<point<T>>& G) {
    G.erase(std::unique(G.begin(), G.end()), G.end());

    while (G.size() > 1 && G[0] == G.back()) {
        G.pop_back();
    }

    std::rotate(G.begin(), std::min_element(G.begin(), G.end(), point<T>::yx_comparator()),
                G.end());

    std::size_t size = 1;
    for (std::size_t i = 1; i < G.size(); i++) {
        const auto& P = G[i - 1];
        const auto& Q = G[i + 1 == G.size() ? 0 : i + 1];
        if (!G[i].inside_segment(P, Q)) {
            G[size++] = G[i];
        }
    }

    G.resize(size);
}

template <typename T>
std::vector<point<T>> minkowsky_sum(std::vector<point<T>> lhs, std::vector<point<T>> rhs) {
    if (lhs.size() == 0 || rhs.size() == 0) {
        return {};
    }

    std::vector<point<T>> result(std::size(lhs) + std::size(rhs));

    normalize_convex(lhs);
    normalize_convex(rhs);

    const auto B = lhs[0] + rhs[0];  // bottom left point of result

    static const auto to_sides_vectors = [](std::vector<point<T>>& convex) {
        const auto last = convex.back();  // save it first
        std::adjacent_difference(std::begin(convex), std::end(convex), std::begin(convex));
        *std::begin(convex) -= last;
        std::rotate(std::begin(convex), std::next(std::begin(convex)), std::end(convex));
    };

    to_sides_vectors(lhs);
    to_sides_vectors(rhs);

    std::merge(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs), std::cend(rhs),
               std::begin(result), point<T>::angle_comparator({0, 0}));

    std::partial_sum(std::begin(result), std::end(result), std::begin(result));

    std::transform(std::begin(result), std::end(result), std::begin(result),
                   [&B](const point<T>& P) { return P + B; });

    normalize_convex(result);

    return result;
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    static auto read_polygon = [] {
        std::int32_t n;
        std::cin >> n;

        std::vector<point<std::int64_t>> pts(n);
        for (auto& P : pts) {
            std::cin >> P;
        }

        return pts;
    };

    std::vector<point<std::int64_t>> sum;

    for (std::int32_t i = 0; i < 3; i++) {
        auto pts = read_polygon();

        normalize_convex(pts);

        if (i == 0) {
            sum = std::move(pts);
        } else {
            sum = minkowsky_sum(std::move(sum), std::move(pts));
        }
    }

    std::int32_t t;
    std::cin >> t;

    while (t--) {
        point<std::int64_t> P;
        std::cin >> P.x >> P.y;

        std::cout << (inside_convex(P * 3, sum) != relative_position::out ? "YES" : "NO") << "\n";
    }
}