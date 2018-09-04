#include <bits/stdc++.h>

using namespace std;

template <class T>
struct maxStack {
	vector <T> vals;
	vector <T> maxVals;
	int size;

	maxStack(int maxSize) {
		vals.resize(maxSize);
		maxVals.resize(maxSize);
		size = 0;
	}

	bool empty() {
		return size == 0;
	}

	void push(T x) {
		vals[size] = x;
		if (size != 0) {
			x = max(x, vals[size - 1]);
		}
		maxVals[size] = x;
		size++;
	}

	T top() {
		assert(size > 0);
		return vals[--size];
	}

	void pop() {
		if (size > 0) {
			size--;
		}
	}

	T getMax() {
		if (size == 0) {
			return numeric_limits <T> :: min();
		}
		return maxVals[size - 1];
	}
};

template <class T>
struct maxQueue {
	maxStack <T> front;
	maxStack <T> back;

	maxQueue(int maxSize) : front(maxSize), back(maxSize) {}

	void push(T x) {
		back.push(x);
	}

	void pop() {
		if (!front.empty()) {
			front.pop();
		} else {
			while (!back.empty()) {
				front.push(back.top());
				back.pop();
			}
		}
	}

	T getMax() {
		assert(!front.empty() || !back.empty());
		return max(front.getMax(), back.getMax());
	}
};


int main() {
	ios::sync_with_stdio(false);
	cin.tie(0);

	int n, d;
	long long p;

	cin >> n >> p >> d;

	vector <int> w(n);
	for (auto & x : w) {
		cin >> x;
	}
	maxQueue <long long> Q(n);	
	long long sum = 0;
	for (int i = 0; i < d; i++) {
		sum += w[i];
		Q.push(w[i]);
	}
	long long dsum = sum;
	int answer = d;
	for (int i = d, j = 0; i < n; i++) {
		sum += w[i];
		dsum += w[i];
		dsum -= w[i - d];
		Q.push(dsum);
		while (sum - Q.getMax() > p) {
			Q.pop();
			sum -= w[j];
			j++;
		}
		answer = max(answer, i - j + 1);
	}

	cout << answer << "\n";
}
