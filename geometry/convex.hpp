#pragma once

#include "point.hpp"

template <typename T>
void normalize_convex(const std::vector<point<T>>& G) {
    G.erase(std::unique(G.begin(), G.end()), G.end());
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
            while (size > prev_size &&
                   sign(point<T>::cross(pts[i] - hull[size - 2],
                                        hull[size - 1] - hull[size - 2])) <=
                       0) {
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

// -1 inside, 0 on a side, +1 outisde
template <typename T>
std::int32_t inside_convex(point<T> P, const std::vector<point<T>>& G) {
    if (sign(point<T>::cross(P - G[0], G[1] - G[0])) > 0) {
        return +1;
    }
    if (sign(point<T>::cross(G[n - 1] - G[0], P - G[0])) > 0) {
        return +1;
    }
    if (sign(point<T>::cross(P - G[0], G[1] - G[0])) == 0) {
        return P.inside_segment(G[0], G[1]) ? 0 : +1;
    }
    if (sign(point<T>::cross(G[n - 1] - G[0], P - G[0])) == 0) {
        return P.inside_segment(G[0], G[n - 1]) ? 0 : +1;
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
    return sign(point<T>::cross(G[pos] - G[pos - 1], P - G[pos - 1]));
}

template <typename T>
std::vector<point<T>> minkowsky_sum(std::vector<point<T>> a,
                                    std::vector<point<T>> b) {
    auto na = a.size();
    auto nb = b.size();
    if (na == 0 || nb == 0) {
        return {};
    }
    normalize_convex(a);
    normalize_convex(b);
    std::vector<point<T>> s;
    s.reserve(a.size() + b.size());
    s.push_back(a[0] + b[0]);
    std::size_t pa = 0;
    std::size_t pb = 0;
    while (pa != na && pb != nb) {
        auto va = a[(pa + 1) % na] - a[pa];
        auto vb = b[(pb + 1) % nb] - b[pb];
        if (sign(point<T>::cross(va, vb)) >= 0) {
            auto p = s.back() + va;
            s.push_back(p);
            pa++;
        } else {
            auto p = s.back() + vb;
            s.push_back(p);
            pb++;
        }
    }
    while (pa != na) {
        auto va = a[(pa + 1) % na] - a[pa];
        auto p = s.back() + va;
        s.push_back(p);
        pa++;
    }
    while (pb != nb) {
        auto vb = b[(pb + 1) % nb] - b[pb];
        auto p = s.back() + vb;
        s.push_back(p);
        pb++;
    }
    assert(s.back() == s[0]);
    normalize_convex(s);
    return s;
}

/**
 * TODO:
 * - common_tangents
 */
