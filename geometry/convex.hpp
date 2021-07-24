#pragma once

#include <numeric>

#include "point.hpp"

template <typename T>
void normalize_convex(const std::vector<point<T>>& G) {
    G.erase(std::unique(G.begin(), G.end()), G.end());
    auto new_end = std::remove_if(std::next(std::begin(G)), std::end(G),
                                  [O = *std::begin(G)](const auto& P) { return P == O; });
    while (G.size() > 1 && G[0] == G.back()) {
        G.pop_back();
    }
    auto c = std::min_element(G.begin(), G.end());
    std::rotate(G.begin(), c, G.end());
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
std::vector<point<T>> convex_hull(std::vector<point<T>> pts) {
    sort(pts.begin(), pts.end(), point<T>::yx_comparator());

    std::vector<point<T>> hull(2 * pts.size());
    std::size_t size = 0;

    hull[size++] = pts[0];

    for (std::size_t step = 0; step < 2; step++) {
        auto prev_size = size;
        for (std::size_t i = 1; i < pts.size(); i++) {
            while (
              size > prev_size
              && sign(point<T>::cross(pts[i] - hull[size - 2], hull[size - 1] - hull[size - 2]))
                   <= 0) {
                size--;
            }
            hull[size++] = pts[i];
        }
        std::reverse(pts.begin(), pts.end());
    }

    if (size > 1 && hull[0] == hull[size - 1]) {
        size--;
    }

    hull.resize(size);

    return hull;
}

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

    std::size_t l = 2;
    std::size_t r = G.size() - 1;
    std::size_t p = r;

    while (l <= r) {
        std::size_t m = (l + r) >> 1;
        if (sign(point<T>::cross(P - G[0], G[m] - G[0])) >= 0) {
            p = m;
            r = m - 1;
        } else {
            l = m + 1;
        }
    }

    const auto c = sign(point<T>::cross(G[p] - G[p - 1], P - G[p - 1]));

    if (c > 0) {
        return relative_position::in;
    }

    if (c < 0) {
        return relative_position::out;
    }

    return relative_position::on;
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
/**
 * TODO:
 * - common_tangents
 */
