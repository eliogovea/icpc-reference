#include <bits/stdc++.h>

using namespace std;

const double epsilon = 1e-9;

vector <int> gauss(vector <vector <double>> A, vector <double> B, vector <double> & X) {
	int n = A.size();
	int m = A[0].size();
	assert(B.size() == n);
	for (auto & r : A) {
		assert(r.size() == m);
	}
	vector <int> where(m, -1);
	for (int r = 0, c = 0; r < n && c < m; c++) {
		int p = r;
		for (int i = r + 1; i < n; i++) {
			if (abs(A[i][c]) > abs(A[p][c])) {
				p = i;
			}
		}
		if (abs(A[p][c]) < epsilon) {
			continue;
		}
		swap(A[r], A[p]);
		swap(B[r], B[p]);
		where[c] = r;
		double x = 1.0 / A[r][c];
		for (int i = 0; i < n; i++) {
			if (i != r) {
				double y = A[i][c] * x;
				for (int j = c; j < m; j++) {
					A[i][j] = A[i][j] - A[r][j] * y;
				}
				B[i] = B[i] - B[r] * y;
			}
		}
		r++;
	}
	
	X.resize(m, 0);
	for (int i = 0; i < m; i++) {
		if (where[i] != -1) {
			X[i] = B[where[i]] / A[where[i]][i];
		}
	}
	
	for (int i = 0; i < n; i++) {
		double s = 0.0;
		for (int j = 0; j < m; j++) {
			s = s + X[j] * A[i][j];
		}
		if (abs(s - B[i]) >= epsilon) {
			return {};
		}
	}
	for (int i = 0; i < m; i++) {
		if (where[i] != -1) {
			for (int j = 0; j < n; j++) {
				if (j != i && abs(A[where[i]][j]) > epsilon) {
					where[i] = -1;
					break;
				}
			}
		}
	} // !!!!!!!!!!!
	return where;
}


int equationsolver() {
	int n;
	while (cin >> n && n) {
		vector <vector <double> > A(n, vector <double> (n));
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				cin >> A[i][j];
			}
		}
		vector <double> B(n);
		for (int i = 0; i < n; i++) {
			cin >> B[i];
		}
		vector <double> X;
		auto where = gauss(A, B, X);
		if (where.size() == 0) {
			cout << "inconsistent\n";
		} else {
			bool multiple = false;
			for (int i = 0; i < n; i++) {
				if (where[i] == -1) {
					multiple = true;
					break;
				}
			}
			if (multiple) {
				cout << "multiple\n";
			} else {
				for (int i = 0; i < n; i++) {
					cout << fixed << X[i];
					if (i + 1 < n) {
						cout << " ";
					}
				}
				cout << "\n";
			}
		}
	}

}


int equationsolverplus() {
	int n;
	while (cin >> n && n) {
		vector <vector <double> > A(n, vector <double> (n));
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				cin >> A[i][j];
			}
		}
		vector <double> B(n);
		for (int i = 0; i < n; i++) {
			cin >> B[i];
		}
		vector <double> X;
		auto where = gauss(A, B, X);
		if (where.size() == 0) {
			cout << "inconsistent\n";
		} else {
			for (int i = 0; i < n; i++) {
				if (where[i] == -1) {
					cout << "?";
				} else {
					cout << fixed << X[i];
				}
				if (i + 1 < n) {
					cout << " ";
				}
			}
			cout << "\n";
		}
	}

}


int main() {
	ios::sync_with_stdio(false);
	cin.tie(0);
	cout.precision(14);

	equationsolver();
	// equationsolverplus();
}