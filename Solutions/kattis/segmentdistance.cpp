// https://open.kattis.com/contests/ffe67x/problems/segmentdistance

#include <bits/stdc++.h>

using namespace std;

typedef long long LL;
typedef long double LD;

inline int sign(const LL x) {
	return (x < 0) ? -1 : (x > 0);
}

struct point_LL {
	LL x, y;
	point_LL() {}
	point_LL(LL _x, LL _y) : x(_x), y(_y) {}
};

point_LL operator - (const point_LL &P, const point_LL &Q) {
	return point_LL(P.x - Q.x, P.y - Q.y);
}

bool operator < (const point_LL &P, const point_LL &Q) {
	if (P.x != Q.x) {
		return P.x < Q.x;
	}
	return P.y < Q.y;
}

bool operator == (const point_LL &P, const point_LL &Q) {
	return !(P < Q) && !(Q < P);
}

inline LL dot_LL(point_LL P, point_LL Q) {
	return P.x * Q.x + P.y * Q.y;
}

inline LL cross_LL(point_LL P, point_LL Q) {
	return P.x * Q.y - P.y * Q.x;
}

inline LD length(point_LL P) {
	return sqrtl((LD)dot_LL(P, P));
}

inline bool inside(point_LL P, point_LL A, point_LL B) {
	return (sign(cross_LL(A - P, B - P)) == 0 && dot_LL(A - P, B - P) <= 0);
}

inline LD point_line_distance_LL(point_LL P, point_LL A, point_LL B) {
	return ((LD)abs(cross_LL(P - A, B - A))) / length(B - A);
}

inline LD point_segment_distance_LL(point_LL P, point_LL A, point_LL B) {
	if (dot_LL(P - A, B - A) <= 0) {
		return length(P - A);
	}
	if (dot_LL(P - A, B - A) >= dot_LL(B - A, B - A)) {
		return length(P - B);
	}
	return point_line_distance_LL(P, A, B);
}

struct point_LD {
	LD x, y;
	point_LD() {}
	point_LD(LD _x, LD _y) : x(_x), y(_y) {}
};

point_LD operator + (point_LD P, point_LD Q) {
	return point_LD(P.x + Q.x, P.y + Q.y);
}

point_LD operator - (point_LD P, point_LD Q) {
	return point_LD(P.x - Q.x, P.y - Q.y);
}

point_LD operator * (point_LD P, double k) {
	return point_LD(P.x * k, P.y * k);
}

inline point_LD get(point_LL P) {
	return point_LD(P.x, P.y);
}

inline LD dot_LD(point_LD P, point_LD Q) {
	return P.x * Q.x + P.y * Q.y;
}

inline LD cross_LD(point_LD P, point_LD Q) {
	return P.x * Q.y - P.y * Q.x;
}

inline point_LD intersect_LD(point_LD A, point_LD B, point_LD C, point_LD D) {
	return A + (B - A) * (cross_LD(C - A, C - D) / cross_LD(B - A, C - D));
}

vector <point_LD> intersect_LL(point_LL A, point_LL B, point_LL C, point_LL D) {
	if (A == B && C == D) {
		if (A == C) {
			return vector <point_LD> (1, get(A));
		}
		return vector <point_LD> ();
	}
	if (A == B) {
		if (inside(A, C, D)) {
			return vector <point_LD> (1, get(A));
		}
		return vector <point_LD> ();
	}
	if (C == D) {
		if (inside(C, A, B)) {
			return vector <point_LD> (1, get(C));
		}
		return vector <point_LD> ();
	}
	if (sign(cross_LL(B - A, D - C)) == 0) {
		if (sign(cross_LL(B - A, C - A)) == 0) {
			if (B < A) {
				swap(A, B);
			}
			if (D < C) {
				swap(C, D);
			}
			if (C < A) {
				swap(A, C);
				swap(B, D);
			}
			if (B < C) {
				return vector <point_LD> ();
			}
			if (B == C) {
				return vector <point_LD> (1, get(B));
			}
			vector <point_LD> res(1, get(C));
			if (D < B) {
				res.push_back(get(D));
			} else {
				res.push_back(get(B));
			}
			return res;
		}
		return vector <point_LD> ();
	}
	
	if (sign(cross_LL(B - A, C - A)) * sign(cross_LL(B - A, D - A)) > 0) {
		return vector <point_LD> ();
	}
	if (sign(cross_LL(D - C, A - C)) * sign(cross_LL(D - C, B - C)) > 0) {
		return vector <point_LD> ();
	}
	return vector <point_LD> (1, intersect_LD(get(A), get(B), get(C), get(D)));
}

LD segment_segment_distance(point_LL A, point_LL B, point_LL C, point_LL D) {
	vector <point_LD> v = intersect_LL(A, B, C, D);
	if (v.size()) {
		return 0.0;
	}
	return min(min(point_segment_distance_LL(A, C, D), point_segment_distance_LL(B, C, D)), min(point_segment_distance_LL(C, A, B), point_segment_distance_LL(D, A, B)));
}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0);
	cout.precision(2);
	
	int n;
	cin >> n;
	
	while (n--) {
		point_LL A, B, C, D;
		cin >> A.x >> A.y >> B.x >> B.y >> C.x >> C.y >> D.x >> D.y;
		cout << fixed << segment_segment_distance(A, B, C, D) << "\n";
	}
}