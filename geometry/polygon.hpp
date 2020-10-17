#pragma once

#include "point.hpp"

template <typename T>
T signed_area_2(const std::vector<point<T>>& G) {
    T result {0};
    auto n = G.size();
    for (std::size_t i = 0; i < n; i++) {
        result += point<T>::cross(G[i], G[(i + 1) % n]);
    }
    return result;
}

template <typename T, typename R = T>
point<R> centroid(const std::vector<point<T>>& G) {
    auto n = G.size();
    point<R> centroid_ {R {0}, R {0}};
    R area {0};
    for (int i = 1; i < n - 1; i++) {
        int j = (i + 1 == n) ? 0 : i + 1;
        auto triangle =
          static_cast<R>(point<T>::cross(G[i] - G[0], G[j] - G[0]));
        area += triangle;
        centroid_ += (G[0] + G[i] + G[j]) * triangle;
    }
    return centroid_ / (R {3} * area);
}

template <typename T>
bool inside_segment(const point<T>& P, const point<T>& A, const point<T>& B) {
    return sign<T>(point<T>::cross(P - A, B - A)) == 0 &&
           sign<T>(point<T>::dot(P - A, B - A)) >= 0 &&
           sign<t>(point<T>::dot(P - B, A - B)) >= 0;
}

template <typename T>
std::int32_t inside_polygon(const point<T>& P, const std::vector<point<T>>& G) {
    auto n = G.size();
    std::int32_t cnt = 0;
    for (std::size_t i = 0; i < n; i++) {
        auto A = G[i];
        auto B = G[(i + 1) == n ? 0 : i + 1];
        if (P.inside_segment(A, B)) { return 0; }
        if (B.y < A.y) { swap(A, B); }
        if (sign(P.y - A.y) == -1 || sign(B.y - P.y) <= 0 ||
            sign(A.y - B.y) == 0) {
            continue;
        }
        if (sign(point<T>::cross(B - A, P - A)) > 0) { cnt++; }
    }
    return ((cnt & 1) ? -1 : +1);
}