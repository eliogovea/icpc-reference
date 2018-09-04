/// Geometry - 2D - INT

#include <bits/stdc++.h>

using namespace std;

typedef long long LL;

inline int sign(const LL x) {
	if (x < 0) {
		return -1;
	}
	if (x > 0) {
		return 1;
	}
	return 0;
}

struct point {
	LL x, y;
	point() {}
	point(LL _x, LL _y) : x(_x), y(_y) {}
};

bool operator < (const point &P, const point &Q) {
	if (P.y != Q.y) {
		return P.y < Q.y;
	}
	return P.x < Q.x;
}

bool operator == (const point &P, const point &Q) {
	return !(P < Q) && !(Q < P);
}

struct compare_x {
	bool operator () (const point &P, const point &Q) {
		if (P.x != Q.x) {
			return P.x < Q.x;
		}
		return P.y < Q.y;
	}
};

void normalize(point &P) {
	assert(P.x != 0 || P.y != 0);
	LL g = __gcd(abs(P.x), abs(P.y));
	P.x /= g;
	P.y /= g;
	if (P.x < 0 || (P.x == 0 && P.y < 0)) {
		P.x = -P.x;
		P.y = -P.y;
	}
}

inline int half_plane(const point &P) {
	if (P.y != 0) {
		return sign(P.y);
	}
	return sign(P.x);
}

point operator + (const point &P, const point &Q) {
	return point(P.x + Q.x, P.y + Q.y);
}

point operator - (const point &P, const point &Q) {
	return point (P.x - Q.x, P.y - Q.y);
}

point operator * (const point &P, const LL k) {
	return point(P.x * k, P.y * k);
}

point operator / (const point &P, const LL k) {
	assert(k != 0 && P.x % k == 0 && P.y % k == 0);
	return point(P.x / k, P.y / k);
}

inline LL dot(const point &P, const point &Q) {
	return P.x * Q.x + P.y * Q.y;
}

inline LL cross(const point &P, const point &Q) {
	return P.x * Q.y - P.y * Q.x;
}

inline LL dist2(const point &P, const point &Q) {
	LL dx = P.x - Q.x;
	LL dy = P.y - Q.y;
	return dx * dx + dy * dy;
}

struct compare_angle {
	point O;
	compare_angle() {}
	compare_angle(point _O) {
		O = _O;
	}
	bool operator () (const point &P, const point &Q) {
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

inline point rotate_ccw_90(const point &P) {
	return point(-P.y, P.x);
}

inline point rotate_cw_90(const point &P) {
	return point(P.y, -P.x);
}

inline bool is_in(LL x, LL a, LL b) {
	if (a > b) {
		swap(a, b);
	}
	return (a <= x && x <= b);
}

inline bool is_in(const point &P, const point &A, const point &B) {
	if (cross(B - A, P - A) != 0) {
		return false;
	}
	return (is_in(P.x, A.x, B.x) && is_in(P.y, A.y, B.y));
}

/// test ok
inline bool segment_segment_intersect(const point &A, const point &B, const point &C, const point &D) {
	if (cross(B - A, D - C) == 0) { // lines are parallel
		return (is_in(A, C, D) || is_in(B, C, D) || is_in(C, A, B) || is_in(D, A, B));
	}
	if (sign(cross(C - A, B - A)) * sign(cross(D - A, B - A)) > 0) {
		return false;
	}
	if (sign(cross(A - C, D - C)) * sign(cross(B - C, D - C)) > 0) {
		return false;
	}
	return true;
}

/// signed polygon area O(n)
inline LL signed_area_2(const vector <point> &polygon) {
	LL res = 0;
	int n = polygon.size();
	for (int i = 0; i < n; i++) {
		int j = (i + 1 == n) ? 0 : i + 1;
		res += cross(polygon[i], polygon[j]);
	}
	return res;
}

/// O(n)
/// test ok
inline bool is_convex(const vector <point> &polygon) {
	int n = polygon.size();
	assert(n >= 3);
	for (int i = 0; i < n; i++) {
		int j = (i + 1) % n;
		int k = (i + 2) % n;
		if (sign(cross(polygon[j] - polygon[i], polygon[k] - polygon[i])) < 0) {
			return false;
		}
	}
	return true;
}

/// O(n)
const int OUT = 0;
const int ON = 1;
const int IN = 2;
/// 0 outside, 1 boundary, 2 inside
/// test ok
inline int point_inside_polygon(const point &P, const vector <point> &polygon) {
	int n = polygon.size();
	int cnt = 0;
	for (int i = 0; i < n; i++) {
		point A = polygon[i];
		point B = polygon[(i + 1) % n];
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
	if (cnt & 1) {
		return IN;
	}
	return OUT;
}

/// O(n * log(n)) -> n size of pts
/// !!! no se como mantener los puntos colineales del convex hull
/// test ???
vector <point> convex_hull_graham_scan(vector <point> pts) {
	sort(pts.begin(), pts.end());
	int n = pts.size();
	if (n <= 2) {
        return pts;
	}
	sort(++pts.begin(), pts.end(), compare_angle(pts[0]));
	vector <point> ch(n);
	int top = 0;
	for (int i = 0; i < n; i++) {
		while (top >= 2 && cross(ch[top - 1] - ch[top - 2], pts[i] - ch[top - 2]) <= 0) {
			top--;
		}
		ch[top++] = pts[i];
	}
	ch.resize(top);
	return ch;
}

/// O(n * log(n)) -> n size of pts
/// test ok
vector <point> convex_hull_monotone_chain(vector <point> pts) {
	if (pts.size() <= 2) {
        return pts;
	}
	sort(pts.begin(), pts.end());
	int n = pts.size();
	vector <point> ch(2 * n);
	int top = 0;
	for (int i = 0; i < n; i++) {
		while (top >= 2 && cross(ch[top - 1] - ch[top - 2], pts[i] - ch[top - 2]) <= 0) {
			top--;
		}
		ch[top++] = pts[i];
	}
	int size = top;
	for (int i = n - 2; i >= 0; i--) {
		while (top - size >= 1 && cross(ch[top - 1] - ch[top - 2], pts[i] - ch[top - 2]) <= 0) {
			top--;
		}
		ch[top++] = pts[i];
	}
	if (ch[0] == ch[top - 1]) {
		top--;
	}
	ch.resize(top);
	return ch;
}

/// test ok
int count_commont_tangents(point C1, LL r1, point C2, LL r2) {
	LL d2 = dist2(C1, C2);
	// if (r1 > r2) {
		// swap(C1, C2);
		// swap(r1, r2);
	// }
	if (d2 > (r1 + r2) * (r1 + r2)) {
		return 4;
	}
	if (d2 == (r1 + r2) * (r1 + r2)) {
		return 3;
	}
	if (d2 < (r1 - r2) * (r1 - r2)) {
		return 0;
	}
	if (d2 == (r1 - r2) * (r1 - r2)) {
		return 1;
	}
	return 2;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

/// O(n)
void normalize_convex(vector <point> &polygon) {
	polygon.erase(unique(polygon.begin(), polygon.end()), polygon.end());
	while (polygon.size() > 1 && polygon[0] == polygon.back()) {
		polygon.pop_back();
	}
	rotate(polygon.begin(), min_element(polygon.begin(), polygon.end()), polygon.end());
	int ptr = 1;
	for (int i = 1; i < polygon.size(); i++) {
		if (is_in(polygon[i], polygon[i - 1], polygon[i + 1 == polygon.size() ? 0 : i + 1]));
		polygon[ptr++] = polygon[i];
	}
	polygon.resize(ptr);
}

/// O(log(n)) -> n size of the polygon
/// apply normalize_convex before
/// test: OK
inline int inside_convex(const point &P, const vector <point> &polygon) {
	int n = polygon.size();
	assert(n >= 3);
	if (cross(P - polygon[0], polygon[1] - polygon[0]) > 0) {
		return OUT;
	}
	if (cross(polygon[n - 1] - polygon[0], P - polygon[0]) > 0) {
		return OUT;
	}
	if (cross(P - polygon[0], polygon[1] - polygon[0]) == 0) {
		return (is_in(P, polygon[0], polygon[1]) ? ON : OUT);
	}
	if (cross(polygon[n - 1] - polygon[0], P - polygon[0]) == 0) {
		return (is_in(P, polygon[0], polygon[n - 1]) ? ON : OUT);
	}
	int lo = 2;
	int hi = n - 1;
	int pos = hi;
	while (lo <= hi) {
		int mid = (lo + hi) >> 1;
		if (cross(P - polygon[0], polygon[mid] - polygon[0]) >= 0) {
			pos = mid;
			hi = mid - 1;
		} else {
			lo = mid + 1;
		}
	}
	int s = sign(cross(polygon[pos] - polygon[pos - 1], P - polygon[pos - 1]));
	if (s == 0) {
		return ON;
	}
	return ((s > 0) ? IN : OUT);
}

inline int get_upper_point(const vector <point> &polygon) {
	int n = polygon.size();
	int upper = 0;
	while (upper + 1 < n && polygon[upper] < polygon[upper + 1]) {
		upper++;
	}
	return upper;
}

/// O(log(n)) -> n size of the polygon
/// apply normalize_convex and get_upper_point to the polygon before
/// TODO: test
bool convex_to_segment_intersection(const vector <point> &polygon, int upper_point, point A, point B) {
	if (sign(cross(polygon[0] - A, B - A)) * sign(cross(polygon[upper_point] - A, B - A)) <= 0) {
		return true;
	}
	if (B < A) {
		swap(A, B);
	}
	if (cross(B - A, polygon[0] - A) > 0) {
		int lo = 1;
		int hi = upper_point;
		int id = 0;
		while (lo <= hi) {
			int mid = (lo + hi) >> 1;
			if (cross(polygon[mid] - A, B - A) >= cross(polygon[mid - 1] - A, B - A)) {
				id = mid;
				lo = mid + 1;
			} else {
				hi = mid - 1;
			}
		}
		return (cross(polygon[id] - A, B - A) >= 0);
	}
	int lo = upper_point;
	int hi = polygon.size() - 1;
	int id = 0;
	while (lo <= hi) {
		int mid = (lo + hi) >> 1;
		if (cross(B - A, polygon[mid]) >= cross(B - A, polygon[mid + 1 == polygon.size() ? 0 : mid + 1])) {
			id = mid;
			hi = mid - 1;
		} else {
			lo = mid + 1;
		}
	}
	return (cross(B - A, polygon[id] - A) >= 0);
}

/// TODO -> terminar y revisar!!!
/// creo que esta mal
/// si pts es un rectangulo ???
// double minimum_area_bounding_box(vector <point> pts) {
	// pts = convex_hull_monotone_chain(pts);
	// if (pts.size() == 3) {
		// return cross(pts[1] - pts[0], pts[2] - pts[0]);
	// }
	// int n = pts.size();
	// int ind[4] = {0, 0, 0, 0};
	// const int low = 0, right = 1, up = 2, left = 3;
	// for (int i = 1; i < n; i++) {
		// if (pts[i] < pts[low]) {
			// ind[low] = i;
		// }
		// if (pts[up] < pts[i]) {
			// ind[up] = i;
		// }
		// if (compare_x(pts[i], pts[left])) {
			// ind[left] = i;
		// }
		// if (compare_x(pts[right], pts[i])) {
			// ind[right] = i;
		// }
	// }
	// double min_area = (double)(pts[up].y - pts[low].y) * (double)(pts[right].x - pts[left].x);
	// for (int i = 0; i < n; i++) { /// i < n ???
		// point v[4];
		// v[low] = pts[ind[low] + 1 == n ? 0 : ind[low] + 1] - pts[ind[low]];
		// v[up] = pts[ind[up]] - pts[ind[up] + 1 == n ? 0 : ind[up] + 1];
		// v[right] = rotate_cw_90(pts[ind[right] + 1 == n ? 0 : ind[right] + 1] - pts[ind[right]]);
		// v[left] = rotate_ccw_90(pts[ind[left] +1 == n ? 0 : ind[left] + 1] - pts[ind[left]]);
		// int min_angle_id = 0;
		// for (int j = 1; j < 4; j++) {
			// if (cross(v[j], v[min_angle_id]) > 0) {
				// min_angle_id = j;
			// }
		// }
		// point V = pts[ind[min_angle_id] + 1 == n ? 0 : ind[min_angle_id] + 1] - pts[ind[min_angle_id]];
		// double w = (double)abs(cross(v, pts[(min_angle_id + 2) % 4] - pts[min_angle_id])) / (double)norm(V);
		// double h = (double)abs(dot(v, pts[(min_angle_id + 1) % 4] - pts[(min_angle_id + 3) % 4])) / (double)norm(V);
		// min_area = min(min_area, w * h);
	// }
	// return min_area;
// }

/// Minkowsky Sum of two convex polygons O(na + nb)
/// Returns a convex polygon
/// ex c = minkowsky_sum(a, b) => every point in c is the sum of two points from a and b respectively
/// and every pair of points from a and b respectively is a point of c (not only the boundary)
/// !!! revisar si esta bien escrito
/// TODO: test
vector <point> minkowsky_sum(vector <point> a, vector <point> b) {
	if (!a.size() || !b.size()) {
		return vector <point> ();
	}
	normalize_convex(a);
	normalize_convex(b);
	vector <point> s;
	s.push_back(a[0] + b[0]);
	int pa = 0;
	int pb = 0;
	while (pa != a.size() && pb != b.size()) {
		point va = a[pa + 1 == a.size() ? 0 : pa + 1] - a[pa];
		point vb = b[pb + 1 == b.size() ? 0 : pb + 1] - b[pb];
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
	while (pa != a.size()) {
		point va = a[pa + 1 == a.size() ? 0 : pa + 1] - a[pa];
		point p = s.back() + va;
		s.push_back(p);
		pa++;
	}
	while (pb != b.size()) {
		point vb = b[pb + 1 == b.size() ? 0 : pb + 1] - b[pb];
		point p = s.back() + vb;
		s.push_back(p);
		pb++;
	}
	assert(s.back() == s[0]);
	normalize_convex(s);
	return s;
}


/// Rotating calipers to find the further pair of points O(n)
/// returns the d ^ 2
/// apply normalize(polygon) before use
/// TODO: test ok
LL convex_diameter_2(vector <point> &polygon) {
	int n = polygon.size();
	int p0 = 0;
	int p1 = 0;
	for (int i = 1; i < n; i++) {
		// < compare y first then x
		if (polygon[i] < polygon[p0]) {
			p0 = i;
		}
		if (polygon[p1] < polygon[i]) {
			p1 = i;
		}
	}
	LL res = dist2(polygon[p0], polygon[p1]);
	int c0 = p0;
	int c1 = p1;
	do {
		point v1 = polygon[p0 + 1 == n ? 0 : p0 + 1] - polygon[p0];
		point v2 = polygon[p1] - polygon[p1 + 1 == n ? 0 : p1 + 1];
		int s = sign(cross(v1, v2));
		if (s == 1) {
			p0 = p0 + 1 == n ? 0 : p0 + 1;
		} else if (s == -1) {
			p1 = p1 + 1 == n ? 0 : p1 + 1;
		} else {
			p0 = p0 + 1 == n ? 0 : p0 + 1;
			p1 = p1 + 1 == n ? 0 : p1 + 1;
		}
		res = max(res, dist2(polygon[p0], polygon[p1]));
	} while (c0 != p0 || c1 != p1);
	return res;
}


double convex_width(vector <point> &polygon) {
	int n = polygon.size();
	int p0 = 0;
	int p1 = 0;
	for (int i = 1; i < n; i++) {
		// < compare y first then x
		if (polygon[i] < polygon[p0]) {
			p0 = i;
		}
		if (polygon[p1] < polygon[i]) {
			p1 = i;
		}
	}
	LL res = dist2(polygon[p0], polygon[p1]);
	int c0 = p0;
	int c1 = p1;
	do {
		point v1 = polygon[p0 + 1 == n ? 0 : p0 + 1] - polygon[p0];
		point v2 = polygon[p1] - polygon[p1 + 1 == n ? 0 : p1 + 1];
		int s = sign(cross(v1, v2));
		if (s == 1) {
			p0 = p0 + 1 == n ? 0 : p0 + 1;
		} else if (s == -1) {
			p1 = p1 + 1 == n ? 0 : p1 + 1;
		} else {
			p0 = p0 + 1 == n ? 0 : p0 + 1;
			p1 = p1 + 1 == n ? 0 : p1 + 1;
		}
		res = max(res, dist2(polygon[p0], polygon[p1]));
	} while (c0 != p0 && c1 != p1);
	return res;
}


//http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_2_B
// test OK
void test_segment_segment_intersection() {
	int q;
	cin >> q;
	while (q--) {
		point A, B, C, D;
		cin >> A.x >> A.y >> B.x >> B.y >> C.x >> C.y >> D.x >> D.y;
		cout << (segment_segment_intersect(A, B, C, D) ? "1" : "0") << "\n";
	}
}

// http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_3_B
// test OK
void test_is_convex() {
	int n;
	cin >> n;
	vector <point> polygon(n);
	for (int i = 0; i < n; i++) {
		cin >> polygon[i].x >> polygon[i].y;
	}
	cout << (is_convex(polygon) ? "1" : "0") << "\n";
}

// http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_3_C
// test OK
void test_point_inside_polygon() {
	int n;
	cin >> n;
	vector <point> polygon(n);
	for (int i = 0; i < n; i++) {
		cin >> polygon[i].x >> polygon[i].y;
	}
	int q;
	cin >> q;
	while (q--) {
		point P;
		cin >> P.x >> P.y;
		cout << point_inside_polygon(P, polygon) << "\n";
	}
}
// https://open.kattis.com/problems/pointinpolygon
// test ???
void test_point_inside_polygon_1() {
	int n;
	while (cin >> n) {
		if (n == 0) {
			break;
		}
		vector <point> polygon(n);
		for (int i = 0; i < n; i++) {
			cin >> polygon[i].x >> polygon[i].y;
		}
		int m;
		cin >> m;
		while (m--) {
			point P;
			cin >> P.x >> P.y;
			int v = point_inside_polygon(P, polygon);
			if (v == 0) {
				cout << "out\n";
			} else if (v == 1) {
				cout << "on\n";
			} else {
				cout << "in\n";
			}
		}
	}
}

// http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_2_A
// test ok
void test_parallel_orthogonal() {
	int q;
	cin >> q;
	while (q--) {
		point A, B, C, D;
		cin >> A.x >> A.y >> B.x >> B.y >> C.x >> C.y >> D.x >> D.y;
		int answer = 0;
		if (cross(B - A, D - C) == 0) {
			answer = 2;
		} else if (dot(B - A, D - C) == 0) {
			answer = 1;
		}
		cout << answer << "\n";
	}
}

// http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_4_A
// test ok
void test_convex_hull() {
	int n;
	cin >> n;
	vector <point> pts(n);
	for (int i = 0; i < n; i++) {
		cin >> pts[i].x >> pts[i].y;
	}
	vector <point> answer = convex_hull_monotone_chain(pts);
	cout << answer.size() << "\n";
	for (int i = 0; i < answer.size(); i++) {
		cout << answer[i].x << " " << answer[i].y << "\n";
	}
}

// https://open.kattis.com/problems/convexhull
// test OK
void test_convex_hull_1() {
	int n;
	while (cin >> n) {
		if (n == 0) {
			break;
		}
		vector <point> pts(n);
		for (int i = 0; i < pts.size(); i++) {
			cin >> pts[i].x >> pts[i].y;
		}
		sort(pts.begin(), pts.end());
		pts.erase(unique(pts.begin(), pts.end()), pts.end());
		pts = convex_hull_monotone_chain(pts);
		cout << pts.size() << "\n";
		for (int i = 0; i < pts.size(); i++) {
			cout << pts[i].x << " " << pts[i].y << "\n";
		}
	}
}

// https://open.kattis.com/problems/convexhull
// test ???
void test_convex_hull_2() {
	int n;
	while (cin >> n) {
		if (n == 0) {
			break;
		}
		vector <point> pts(n);
		for (int i = 0; i < pts.size(); i++) {
			cin >> pts[i].x >> pts[i].y;
		}
		sort(pts.begin(), pts.end());
		pts.erase(unique(pts.begin(), pts.end()), pts.end());
		pts = convex_hull_graham_scan(pts);
		cout << pts.size() << "\n";
		for (int i = 0; i < pts.size(); i++) {
			cout << pts[i].x << " " << pts[i].y << "\n";
		}
	}
}

// http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_1_C
// test OK
void test_counter_clockwise() {
	point A, B;
	cin >> A.x >> A.y >> B.x >> B.y;
	int q;
	cin >> q;
	while (q--) {
		point C;
		cin >> C.x >> C.y;
		int s = sign(cross(B - A, C - A));
		if (s != 0) {
			cout << ((s > 0) ? "COUNTER_CLOCKWISE" : "CLOCKWISE") << "\n";
		} else {
			if (dot(B - A, C - A) < 0) {
				cout << "ONLINE_BACK\n";
			} else {
				cout << (dist2(A, B) < dist2(A, C) ? "ONLINE_FRONT" : "ON_SEGMENT") << "\n";
			}
		}
	}
}

// http://judge.u-aizu.ac.jp/onlinejudge/status.jsp#submit/CGL_7/A
// test OK
void test_count_common_tangents() {
	point C1, C2;
	LL r1, r2;
	cin >> C1.x >> C1.y >> r1 >> C2.x >> C2.y >> r2;
	int answer = count_commont_tangents(C1, r1, C2, r2);
	cout << answer << "\n";
}

// http://codeforces.com/contest/166/problem/B
// test OK
void test_inside_convex() {
	int na;
	cin >> na;
	vector <point> a(na);
	for (int i = 0; i < na; i++) {
		cin >> a[i].x >> a[i].y;
	}
	reverse(a.begin(), a.end()); // points are in cw
	normalize_convex(a);
	int nb;
	cin >> nb;
	vector <point> b(nb);
	for (int i = 0; i < nb; i++) {
		cin >> b[i].x >> b[i].y;
	}
	bool ok = true;
	for (int i = 0; i < nb; i++) {
		if (inside_convex(b[i], a) != IN) {
			ok = false;
			break;
		}
	}
	cout << (ok ? "YES" : "NO") << "\n";
}
// https://open.kattis.com/problems/polygonarea
// test ok
void test_signed_area_2() {
	int n;
	while (cin >> n) {
		if (n == 0) {
			break;
		}
		vector <point> pts(n);
		for (int i = 0; i < n; i++) {
			cin >> pts[i].x >> pts[i].y;
		}
		rotate(pts.begin(), min_element(pts.begin(), pts.end()), pts.end());
		for (int i = 1; i < pts.size(); i++) {
			pts[i] = pts[i] - pts[0];
		}
		pts[0] = point(0, 0);
		LL answer = signed_area_2(pts);
		cout << (answer < 0 ? "CW" : "CCW") << " ";
		answer = abs(answer);
		cout << answer / 2LL << "." << 5LL * (answer % 2LL) << "\n";
	}
}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0);

	// test_segment_segment_intersection();
	// test_is_convex();
	// test_point_inside_polygon();
	test_point_inside_polygon_1();
	// test_parallel_orthogonal();
	// test_convex_hull();
	// test_convex_hull_1();
	// test_convex_hull_2();
	// test_counter_clockwise();
	// test_count_common_tangents();
	// test_inside_convex();
	// test_signed_area_2();
}
