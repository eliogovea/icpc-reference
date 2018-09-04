// https://open.kattis.com/problems/suffixsorting

#include <bits/stdc++.h>

using namespace std;

vector <int> get_suffix_array(const string &s, const int alpha) {
	int n = s.size() + 1;
	vector <int> suffix_array(n);
	int classes = alpha;
	vector <int> class_id(n);
	for (int i = 0; i < n; i++) {
		suffix_array[i] = i;
		class_id[i] = s[i];
	}
	for (int length = 1; length < 2 * n; length <<= 1) {
		int half_length = length >> 1;
		vector <int> cur_order(n);
		for (int i = 0; i < n; i++) {
			cur_order[i] = suffix_array[i] - half_length;
			if (cur_order[i] < 0) {
				cur_order[i] += n;
			}
		}
		vector <int> cnt(classes);
		for (int i = 0; i < n; i++) {
			cnt[class_id[cur_order[i]]]++;
		}
		for (int i = 1; i < classes; i++) {
			cnt[i] += cnt[i - 1];
		}
		for (int i = n - 1; i >= 0; i--) {
			int p = cur_order[i];
			suffix_array[--cnt[class_id[p]]] = p;
		}
		
		vector <int> new_class_id(n);
		classes = 1;
		new_class_id[suffix_array[0]] = 0;
		for (int i = 1; i < n; i++) {
			if (class_id[suffix_array[i]] != class_id[suffix_array[i - 1]]) {
				classes++;
			} else {
				int pos1 = suffix_array[i] + half_length;
				if (pos1 > n) {
					pos1 -= n;
				}
				int pos2 = suffix_array[i - 1] + half_length;
				if (pos2 > n) {
					pos2 -= n;
				}
				if (class_id[pos1] != class_id[pos2]) {
					classes++;
				}
			}
			new_class_id[suffix_array[i]] = classes - 1;
		}
		class_id = new_class_id;
	}
	for (int i = 0; i + 1 < n; i++) {
		suffix_array[i] = suffix_array[i + 1];
	}
	suffix_array.pop_back();
	return suffix_array;
}

void test(string s) {
	vector <int> sa = get_suffix_array(s, 256);
	for (int i = 0; i < s.size(); i++) {
		cout << sa[i] << ": " << s.substr(sa[i], s.size() - sa[i]) << "\n";
	}
	cout << "\n";
}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0);
	
	const int alpha = 256;
	
	string s, q;
	while (getline(cin, s)) {
		vector <int> suffix_array = get_suffix_array(s, alpha);
		getline(cin, q);
		stringstream in(q);
		int t;
		in >> t;
		while (t--) {
			int p;
			in >> p;
			cout << suffix_array[p];
			if (t != 0) {
				cout << " ";
			}
		}
		cout << "\n";
	}
}