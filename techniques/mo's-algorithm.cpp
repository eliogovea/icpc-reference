#include <bits/stdc++.h>

using namespace std;

struct query {
	int l, r, id;
};

struct cmp {
	int x;

	cmp(int _x) : x(_x) {}

	bool operator () (const query & a, const query & b) {
		if (a.l / x != b.l / x) {
			return a.l < b.l;
		}
		return a.r < b.r;
	}
};

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0);

	int n, m;
	cin >> n >> m;
	vector <int> a(n);
	for (int i = 0; i < n; i++) {
		cin >> a[i];
		a[i]--;
	}
	vector <query> q(m);
	for (int i = 0; i < m; i++) {
		cin >> q[i].l >> q[i].r;
		q[i].id = i;
	}

	int b = 1 + sqrt(n);
	sort(q.begin(), q.end(), cmp(b));

	int l = 0;
	int r = 0;
	vector <int> f(n);
	vector <set <int> > v(n);
	f[a[0]]++;
	v[f[a[0]]].insert(a[0]);
	int top = 1;

	vector <int> answer(m);

	for (int i = 0; i < m; i++) {
		while (r < q[i].r) {
			r++;
			v[f[a[r]]].erase(a[r]);
			f[a[r]]++;

			if (f[a[r]] > top) {
				top++;
			}

			v[f[a[r]]].insert(a[r]);
		}
		while (r > q[i].r) {
			v[f[a[r]]].erase(a[r]);
			f[a[r]]--;

			if (v[top].size() == 0) {
				top--;
			}

			v[f[a[r]]].insert(a[r]);
			r--;
		}
		while (l < q[i].l) {
			v[f[a[l]]].erase(a[l]);
			f[a[l]]--;

			if (v[top].size() == 0) {
				top--;
			}

			v[f[a[l]]].insert(a[l]);
			l++;
		}

		while (l > q[i].l) {
			l--;
			v[f[a[l]]].erase(a[l]);
			f[a[l]]++;

			if (f[a[l]] > top) {
				top++;
			}

			v[f[a[l]]].insert(a[l]);
		}

		if (f[*v[top].begin()] > q[i].r - q[i].l + 1) {
			answer[q[i].id] = *v[top].begin();
		}
	}

	for (int i = 0; i < m; i++) {
		cout << answer[i] << "\n";
	}
}
