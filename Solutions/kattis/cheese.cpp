#include <bits/stdc++.h>

using namespace std;

const double EPS = 1e-10;

const int N = 10005;

int n, s;
double r[N], x[N], y[N], z[N];

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0);
	cout.precision(12);
	cin >> n >> s;
	double tvol = 100.0 * 100.0 * 100.0;
	for (int i = 0; i < n; i++) {
		cin >> r[i] >> x[i] >> y[i] >> z[i];
		r[i] /= 1000.0;
		x[i] /= 1000.0;
		y[i] /= 1000.0;
		z[i] /= 1000.0;
		tvol -= 4.0 * M_PI * r[i] * r[i] * r[i] / 3.0;
	}
	double need = tvol / (double)s;
	double curz = 0.0;
	for (int i = 0; i < s - 1; i++) {
		double lo = curz;
		double hi = 100.0;
		for (int it = 0; it < 200; it++) {
			double mid = (lo + hi) * 0.5;
			double vol = 100.0 * 100.0 * (mid - curz);
			for (int j = 0; j < n; j++) {
				if (z[j] - r[j] > mid + EPS) {
					continue;
				}
				if (z[j] + r[j] + EPS < curz) {
					continue;
				}
				double z0 = max(curz, z[j] - r[j]) - z[j];
				double z1 = min(mid, z[j] + r[j]) - z[j];
				vol -= M_PI * (r[j] * r[j] * (z1 - z0) - (z1 * z1 * z1 - z0 * z0 * z0) / 3.0);
			}
			if (vol > need) {
				hi = mid;
			} else {
				lo = mid;
			}
		}
		double cut = (lo + hi) / 2.0 - curz;
		cout << fixed << cut << "\n";
		curz = (lo + hi) / 2.0;
	}
	cout << fixed << 100.0 - curz << "\n";
}
