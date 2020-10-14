#include <bits/stdc++.h>

using namespace std;

typedef long long LL;
const double INF = 1e17, EPS = 1e-10, PI = acos(-1);

// be careful with EPS
// !!! implement sign for LL and double
// most of the functions work for both

inline int sign(const LL x) { return (x < 0) ? -1 : (x > 0); }

inline int sign(const double x) { return x < -EPS ? -1 : x > EPS; }

inline bool is_in(double a, double b, double x) {
    if (a > b) {
        swap(a, b);
    }
    return (a - EPS <= x && x <= b + EPS);
}

struct point {
    double x, y;
    point(double _x = 0, double _y = 0) : x(_x), y(_y) {}
};

bool operator<(const point& P, const point& Q) {
    if (sign(P.y - Q.y) != 0) {
        return P.y < Q.y;
    }
    return (sign(P.x - Q.x) == -1);
}

bool operator==(const point& P, const point& Q) { return !(P < Q) && !(Q < P); }

struct compare_x {
    bool operator()(const point& P, const point& Q) {
        if (sign(P.x - Q.x) != 0) return P.x < Q.x;
        return P.y < Q.y;
    }
};

struct compare_y {
    bool operator()(const point& P, const point& Q) {
        if (sign(P.y - Q.y) != 0) return P.y < Q.y;
        return P.x < Q.x;
    }
};

inline void read(point& P) { cin >> P.x >> P.y; }

point operator+(const point& P, const point& Q) {
    return point(P.x + Q.x, P.y + Q.y);
}

point operator-(const point& P, const point& Q) {
    return point(P.x - Q.x, P.y - Q.y);
}

point operator*(const point& P, const double k) {
    return point(P.x * k, P.y * k);
}

point operator/(const point& P, const double k) {
    assert(sign(k) != 0);
    return point(P.x / k, P.y / k);
}

inline int half_plane(const point& P) {
    if (sign(P.y) != 0) {
        return sign(P.y);
    }
    return sign(P.x);
}

inline double dot(const point& P, const point& Q) {
    return P.x * Q.x + P.y * Q.y;
}

inline double cross(const point& P, const point& Q) {
    return P.x * Q.y - P.y * Q.x;
}

inline double norm2(const point& P) { return dot(P, P); }

inline double norm(const point& P) { return sqrt(dot(P, P)); }

// squared euclidean distance
inline double dist2(const point& P, const point& Q) { return norm2(P - Q); }

// euclidean distance
inline double dist(const point& P, const point& Q) {
    return sqrt(dot(P - Q, P - Q));
}

/// returns true if P belongs in segment AB
inline bool is_in(point A, point B, point P) {
    if (sign(cross(B - A, P - A)) != 0) {
        return false;
    }
    return (is_in(A.x, B.x, P.x) && is_in(A.y, B.y, P.y));
}

point rotate_point(point P, double angle) {
    return point(P.x * cos(angle) - P.y * sin(angle),
                 P.y * cos(angle) + P.x * sin(angle));
}

point rotate_90_ccw(point P) { return point(-P.y, P.x); }

point rotate_90_cw(point P) { return point(P.y, -P.x); }

void normalize(point& P) {  // for vectors
    assert(sign(P.x) != 0 || sign(P.y) != 0);
    // LL g = __gcd(abs(P.x), abs(P.y));
    // P.x /= g; P.y /= g;
    if (P.x < 0 || (P.x == 0 && P.y < 0)) {
        P.x = -P.x, P.y = -P.y;
    }
}

struct compare_angle {
    point O;
    compare_angle(point _O = point(0, 0)) { O = _O; }
    bool operator()(const point& P, const point& Q) {
        if (half_plane(P - O) != half_plane(Q - O)) {
            return half_plane(P - O) < half_plane(Q - O);
        }
        int c = sign(cross(P - O, Q - O));
        if (c != 0) {
            return (c > 0);
        }
        return dist2(P, O) < dist2(Q, O);
    }
};

inline double signed_area_2(const vector<point>& G) {
    double res = 0;
    int n = G.size();
    for (int i = 0; i < n; i++) {
        res += cross(G[i], G[(i + 1) % n]);
    }
    return res;
}

inline double abs_area(const vector<point>& G) {
    return abs(0.5 * signed_area_2(G));
}

inline bool is_convex(const vector<point>& G) {
    int n = G.size();
    assert(n >= 3);
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        int k = (i + 2) % n;
        if (sign(cross(G[j] - G[i], G[k] - G[i])) < 0) {
            return false;
        }
    }
    return true;
}

void normalize_convex(vector<point>& G) {
    G.erase(unique(G.begin(), G.end()), G.end());
    while (G.size() > 1 && G[0] == G.back()) {
        G.pop_back();
    }
    rotate(G.begin(), min_element(G.begin(), G.end()), G.end());
    int ptr = 1;
    for (int i = 1; i < G.size(); i++) {
        if (is_in(G[i], G[i - 1], G[i + 1 == G.size() ? 0 : i + 1])) continue;
        G[ptr++] = G[i];
    }
    G.resize(ptr);
}

const int OUT = 0, ON = 1, IN = 2;
int inside_polygon(point P, const vector<point>& G) {
    int n = G.size(), cnt = 0;
    for (int i = 0; i < n; i++) {
        point A = G[i], B = G[(i + 1) == n ? 0 : i + 1];
        if (is_in(P, A, B)) {
            return ON;
        }
        if (B.y < A.y) {
            swap(A, B);
        }
        if (P.y < A.y || B.y <= P.y || A.y == B.y) {
            continue;
        }
        if (sign(cross(B - A, P - A)) > 0) {
            cnt++;
        }
    }
    return ((cnt & 1) ? IN : OUT);
}

/// O(log(n)) !!! apply normalize_convex before !!!
int inside_convex(point P, const vector<point>& G) {
    int n = G.size();
    assert(n >= 3);
    if (sign(cross(P - G[0], G[1] - G[0])) > 0) {
        return OUT;
    }
    if (sign(cross(G[n - 1] - G[0], P - G[0])) > 0) {
        return OUT;
    }
    if (sign(cross(P - G[0], G[1] - G[0])) == 0) {
        return (is_in(P, G[0], G[1]) ? ON : OUT);
    }
    if (sign(cross(G[n - 1] - G[0], P - G[0])) == 0) {
        return (is_in(P, G[0], G[n - 1]) ? ON : OUT);
    }
    int lo = 2, hi = n - 1, pos = hi;
    while (lo <= hi) {
        int mid = (lo + hi) >> 1;
        if (sign(cross(P - G[0], G[mid] - G[0])) >= 0) {
            pos = mid, hi = mid - 1;
        } else {
            lo = mid + 1;
        }
    }
    int s = sign(cross(G[pos] - G[pos - 1], P - G[pos - 1]));
    if (s == 0) {
        return ON;
    }
    return ((s > 0) ? IN : OUT);
}

vector<point> convex_hull(vector<point> pts) {
    if (pts.size() <= 2) {
        return pts;
    }
    sort(pts.begin(), pts.end());
    int n = pts.size(), t = 0;
    vector<point> ch(2 * n);
    for (int i = 0; i < n; i++) {
        while (t > 1 &&
               sign(cross(ch[t - 1] - ch[t - 2], pts[i] - ch[t - 2])) <= 0) {
            t--;
        }
        ch[t++] = pts[i];
    }
    int pt = t;
    for (int i = n - 2; i >= 0; i--) {
        while (t > pt &&
               sign(cross(ch[t - 1] - ch[t - 2], pts[i] - ch[t - 2])) <= 0) {
            t--;
        }
        ch[t++] = pts[i];
    }
    if (ch[0] == ch[t - 1]) t--;
    ch.resize(t);
    return ch;
}

int count_commont_tangents(point C1, LL r1, point C2, LL r2) {
    LL d2 = dist2(C1, C2);
    if (d2 > (r1 + r2) * (r1 + r2)) return 4;
    if (d2 == (r1 + r2) * (r1 + r2)) return 3;
    if (d2 < (r1 - r2) * (r1 - r2)) return 0;
    if (d2 == (r1 - r2) * (r1 - r2)) return 1;
    return 2;
}

inline int get_upper_point(const vector<point>& G) {
    int n = G.size(), upper = 0;
    while (upper + 1 < n && G[upper] < G[upper + 1]) {
        upper++;
    }
    return upper;
}

/// O(log(n)) !!! apply normalize_convex and
/// get_upper_point before !!!
/// !!! TODO test
bool convex_to_segment_intersection(vector<point>& G, int upper, point A,
                                    point B) {
    if (sign(cross(G[0] - A, B - A)) * sign(cross(G[upper] - A, B - A)) <= 0) {
        return true;
    }
    int n = G.size();
    if (B < A) {
        swap(A, B);
    }
    if (cross(B - A, G[0] - A) > 0) {
        int lo = 1, hi = upper, id = 0;
        while (lo <= hi) {
            int mid = (lo + hi) >> 1;
            if (cross(G[mid] - A, B - A) >= cross(G[mid - 1] - A, B - A)) {
                id = mid, lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }
        return (cross(G[id] - A, B - A) >= 0);
    }
    int lo = upper, hi = ((int)G.size()) - 1, id = 0;
    while (lo <= hi) {
        int mid = (lo + hi) >> 1;
        if (cross(B - A, G[mid]) >= cross(B - A, G[(mid + 1) % n])) {
            id = mid, hi = mid - 1;
        } else {
            lo = mid + 1;
        }
    }
    return (cross(B - A, G[id] - A) >= 0);
}

/// O(n)
vector<point> minkowsky_sum(vector<point> a, vector<point> b) {
    int na = a.size(), nb = b.size();
    if (na == 0 || nb == 0) return {};
    normalize_convex(a);
    normalize_convex(b);
    vector<point> s;
    s.push_back(a[0] + b[0]);
    int pa = 0, pb = 0;
    while (pa != na && pb != nb) {
        point va = a[(pa + 1) % na] - a[pa];
        point vb = b[(pb + 1) % nb] - b[pb];
        if (sign(cross(va, vb)) >= 0) {
            point p = s.back() + va;
            s.push_back(p);
            pa++;
        } else {
            point p = s.back() + vb;
            s.push_back(p);
            pb++;
        }
    }
    while (pa != na) {
        point va = a[(pa + 1) % na] - a[pa];
        point p = s.back() + va;
        s.push_back(p);
        pa++;
    }
    while (pb != nb) {
        point vb = b[(pb + 1) % nb] - b[pb];
        point p = s.back() + vb;
        s.push_back(p);
        pb++;
    }
    assert(s.back() == s[0]);
    normalize_convex(s);
    return s;
}

// Rotating calipers to find the further pair
//  of points O(n). Returns the d ^ 2
// !!! apply normalize(G) before use
// for width:
// same idea but using distance from point to line
LL convex_diameter_2(vector<point>& G) {
    int n = G.size(), p0 = 0, p1 = 0;
    for (int i = 1; i < n; i++) {
        // < compare y first then x
        if (G[i] < G[p0]) {
            p0 = i;
        }
        if (G[p1] < G[i]) {
            p1 = i;
        }
    }
    LL res = dist2(G[p0], G[p1]);
    int c0 = p0, c1 = p1;
    do {
        point v1 = G[p0 + 1 == n ? 0 : p0 + 1] - G[p0];
        point v2 = G[p1] - G[p1 + 1 == n ? 0 : p1 + 1];
        int s = sign(cross(v1, v2));
        if (s == 1) {
            p0 = p0 + 1 == n ? 0 : p0 + 1;
        } else if (s == -1) {
            p1 = p1 + 1 == n ? 0 : p1 + 1;
        } else {
            p0 = p0 + 1 == n ? 0 : p0 + 1;
            p1 = p1 + 1 == n ? 0 : p1 + 1;
        }
        res = max(res, (LL)dist2(G[p0], G[p1]));
    } while (c0 != p0 || c1 != p1);
    return res;
}

point project(point P, point P1, point P2) {
    return P1 + (P2 - P1) * (dot(P2 - P1, P - P1) / norm2(P2 - P1));
}

point reflect(point P, point P1, point P2) {
    return project(P, P1, P2) * 2.0 - P;
}

double point_to_line(point P, point A, point B) {
    return abs(cross(B - A, P - A) / norm(B - A));
}

double point_to_segment(point P, point A, point B) {
    if (sign(dot(P - A, B - A)) <= 0) {
        return dist(P, A);
    }
    if (sign(dot(P - B, A - B)) <= 0) {
        return dist(P, B);
    }
    return point_to_line(P, A, B);
}

// lines intersection
point intersect(point A, point B, point C, point D) {
    return A + (B - A) * (cross(C - A, C - D) / cross(B - A, C - D));
}

bool intersectQuery(point A, point B, point C, point D) {
    if (sign(cross(B - A, C - A)) == 0) {
        if (sign(cross(B - A, C - A)) == 0) {
            if (B < A) swap(A, B);
            if (D < C) swap(C, D);
            if (C < A) swap(A, C), swap(B, D);
            return C < B || C == B;
        }
        return false;
    }
    if (sign(cross(C - A, B - A)) * sign(cross(D - A, B - A)) == 1)
        return false;
    if (sign(cross(A - C, D - C)) * sign(cross(B - C, D - C)) == 1)
        return false;
    return true;
}  // segment AB intersects segment CD

double segment_to_segment(point A, point B, point C, point D) {
    if (intersectQuery(A, B, C, D)) return 0.0;
    return min(min(point_to_segment(A, C, D), point_to_segment(B, C, D)),
               min(point_to_segment(C, A, B), point_to_segment(D, A, B)));
}

// TODO test
point circle_center(point A, point B, point C) {
    assert(abs(cross(B - A, C - A)) > EPS);  // no colinear
    return intersect((A + B) / 2.0, (A + B) / 2.0 + rotate_90_ccw(B - A),
                     (B + C) / 2.0, (B + C) / 2.0 + rotate_90_ccw(C - B));
}

pair<point, point> point_circle_tangent(point P, point C, double r) {
    double d = dist(P, C), a = asin(r / d);
    double l = sqrt(d * d - r * r);
    point A = P + rotate_point((C - P) * (l / d), a);
    point B = P + rotate_point((C - P) * (l / d), -a);
    return {A, B};
}

vector<pair<point, point>> common_tangents(point C1, double r1, point C2,
                                           double r2) {
    double d = dist(C1, C2);
    assert(!(d <= EPS && abs(r1 - r2) <= EPS));
    if (r2 > r1) swap(C1, C2), swap(r1, r2);
    if (r1 > d + r2 + EPS) return {};
    if (abs(r1 - d - r1) <= EPS)
        return {{C1 + (C2 - C1) * (r1 / d), C1 + (C2 - C1) * (r1 / d)}};
    vector<pair<point, point>> answer;
    {
        auto t = point_circle_tangent(C2, C1, r1 - r2);
        auto V_first =
            rotate_point((t.first - C2) * (r2 / dist(t.first, C2)), 0.5 * PI);
        point V_second = rotate_point(
            (t.second - C2) * (r2 / dist(t.second, C2)), -0.5 * PI);
        answer.push_back(make_pair(C2 + V_first, t.first + V_first));
        answer.push_back(make_pair(C2 + V_second, t.second + V_second));
    }
    if (abs(d - r1 - r2) <= EPS) {
        answer.push_back(
            make_pair(C1 + (C2 - C1) * (r1 / d), C1 + (C2 - C1) * (r1 / d)));
    } else if (d > r1 + r2 + EPS) {
        auto t = point_circle_tangent(C2, C1, r1 + r2);
        point V_first =
            rotate_point((t.first - C2) * (r2 / dist(t.first, C2)), -0.5 * PI);
        point V_second =
            rotate_point((t.second - C2) * (r2 / dist(t.second, C2)), 0.5 * PI);
        answer.push_back(make_pair(C2 + V_first, t.first + V_first));
        answer.push_back(make_pair(C2 + V_second, t.second + V_second));
    }  // TODO avoid rotate(P, 0.5 PI), use rotate_90_ccw
    return answer;
}

vector<point> line_circle_intersect(point A, point B, point C, double r) {
    point PC = project(C, A, B);
    double d = dist(C, PC);
    if (d > r + EPS) return {};
    if (sign(d - r) == 0) return {PC};
    double l = sqrt(r * r - d * d);
    point v = (B - A) * (l / dist(A, B));
    return {PC + v, PC - v};
}

vector<point> circle_circle_intersect(point C1, double r1, point C2,
                                      double r2) {
    if (r2 > r1) swap(r2, r1), swap(C2, C1);
    double d = dist(C1, C2);
    assert(!(d <= EPS && abs(r1 - r2) <= EPS));
    if (d > r1 + r2 + EPS || r1 > d + r2 + EPS) return {};
    if (abs(d - (r1 + r2)) <= EPS || abs(r1 - (d + r2)) <= EPS)
        return {C1 + (C2 - C1) * (r1 / d)};
    double a = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
    double b = sqrt(r1 * r1 - a * a);
    point P = C1 + (C2 - C1) * (a / d);
    point V = rotate_point(C2 - C1, 0.5 * PI) * (b / d);
    return {P + V, P - V};
}

double closest_pair_of_points(vector<point> pts) {
    sort(pts.begin(), pts.end(), compare_x());
    multiset<point> S;
    int n = pts.size();
    double res = INF;
    for (int i = 0, last = 0; i < n; i++) {
        while (last < i && pts[i].x - pts[last].x >= res + EPS)
            S.erase(S.find(pts[last++]));
        auto lo = S.lower_bound(point(-INF, pts[i].y - res - EPS));
        auto hi = S.upper_bound(point(INF, pts[i].y + res + EPS));
        while (lo != hi) res = min(res, dist(pts[i], *lo)), lo++;
        S.insert(pts[i]);
    }
    return res;
}

/// half planes intersection
struct line {  // a * x + b * y + c = 0
    double a, b, c, angle;
    line() {}
    line(double _a, double _b, double _c) : a(_a), b(_b), c(_c) {
        set_angle(); /* !!! */
    }
    line(point P, point Q) {
        double dx = Q.x - P.x, dy = Q.y - P.y;
        double len = sqrt(dx * dx + dy * dy);
        dx /= len;
        dy /= len;
        a = -dy;
        b = dx;
        // c = -cross(point(Q - P, P));
        c = -(a * P.x + b * P.y);
        set_angle();  /// !!!
    }
    inline int side(const point& P) { return sign(a * P.x + b * P.y + c); }
    inline void set_angle() { angle = atan2(-a, b); }
};  // half plane too

inline point intersect(const line& a, const line& b) {
    double det = a.a * b.b - a.b * b.a;
    // assert(abs(det) > EPS); // !!!
    double det_x = (-a.c) * b.b - a.b * (-b.c);
    double det_y = a.a * (-b.c) - (-a.c) * b.a;
    return point(det_x / det, det_y / det);
}

// Tested: Campamento UCI2017(Gleb) Contest 3 Problem G
vector<point> half_planes_intersection(vector<line> all) {
    const double INF = 1e9;  /// segun el problema
    all.push_back(line(point(-INF, -INF), point(INF, -INF)));
    all.push_back(line(point(INF, -INF), point(INF, INF)));
    all.push_back(line(point(INF, INF), point(-INF, INF)));
    all.push_back(line(point(-INF, INF), point(-INF, -INF)));
    int n = all.size();
    sort(all.begin(), all.end(), [&](const line& a, const line& b) {
        if (sign(a.angle - b.angle) != 0) return a.angle < b.angle;
        return a.c > b.c;
    });
    int ptr = 1;
    for (int i = 1; i < n; i++) {
        if (sign(all[i].angle - all[ptr - 1].angle) == 0) continue;
        // TODO cambiar eps para comparar angulos
        all[ptr++] = all[i];
    }
    if (ptr > 1 && sign(all[0].angle - all[ptr - 1].angle - 2.0 * PI) == 0) {
        if (all[ptr - 1].c < all[0].c) swap(all[ptr - 1], all[0]);
        ptr--;
    }
    all.resize(ptr);
    n = all.size();
    vector<line> Q(n);
    int head = 0, tail = 0;
    for (int i = 0; i < n; i++) {
        while (head + 1 < tail &&
               all[i].side(intersect(Q[tail - 2], Q[tail - 1])) != 1)
            tail--;
        while (head + 1 < tail &&
               all[i].side(intersect(Q[head], Q[head + 1])) != 1)
            head++;
        Q[tail++] = all[i];
    }
    while (head + 1 < tail &&
           Q[head].side(intersect(Q[tail - 1], Q[tail - 2])) != 1)
        tail--;
    while (head + 1 < tail &&
           Q[tail - 1].side(intersect(Q[head], Q[head + 1])) != 1)
        head++;  /// not sure

    vector<point> hull(tail - head);
    for (int i = head; i < tail; i++) {
        int j = (i + 1 == tail) ? head : i + 1;
        hull[i - head] = intersect(Q[i], Q[j]);
    }
    return hull;
}

// Pick Theorem
// A = I + B / 2 - 1

/// cut G, left side of PQ (Q - P)
vector<point> convex_cut(const vector<point>& G, point P, point Q) {
    int n = G.size();
    vector<point> res;
    for (int i = 0; i < n; i++) {
        int si = sign(cross(Q - P, G[i] - P));
        if (si >= 0) res.push_back(G[i]);
        int j = (i + 1 == n) ? 0 : i + 1;
        int sj = sign(cross(Q - P, G[j] - P));
        if (si * sj == -1) res.push_back(intersect(P, Q, G[i], G[j]));
    }
    return res;
}

point centroid(const vector<point>& g) {
    int n = g.size();
    point C(0, 0);
    double area = 0.0;
    for (int i = 1; i < n - 1; i++) {
        int j = (i + 1 == n) ? 0 : i + 1;
        double a = cross(g[i] - g[0], g[j] - g[0]);
        area += a;
        C = C + (g[0] + g[i] + g[j]) * a;
    }
    C = C / (3.0 * area);
    return C;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);
}
