#include <bits/stdc++.h>

using namespace std;

const int N = 100 * 1000 + 10;

int n;
int x[N], y[N];
int id[N];
pair <int, int> p[N];
vector <int> G[N];

bool cmp_x(const int &a, const int &b) {
	return x[a] < x[b];
}

bool cmp_y(const int &a, const int &b) {
	return y[a] < y[b]; 
}

int scc_count;
int scc_id[N];
int scc_x[N];
int scc_y[N];
int scc_size[N];
int st[N];
int top;
int low[N];
int time_in[N];
int timer;

bool cmp_scc_x(const int &a, const int &b) {
	return scc_x[a] < scc_x[b];
}

void dfs(int u) {
	low[u] = time_in[u] = ++timer;
	st[++top] = u;
	for (auto & v : G[u]) {
		if (time_in[v] == 0) {
			dfs(v);
			low[u] = min(low[u], low[v]);
		} else if (scc_id[v] == -1) {
			low[u] = min(low[u], low[v]);
		}
	}
	if (time_in[u] == low[u]) {
		// cerr << "new scc:\n";
		do {
			// cerr << st[top] << " ";
			scc_id[st[top]] = scc_count;
			scc_size[scc_count]++;
		} while (st[top--] != u);
		// cerr << "\n";
		scc_x[scc_count] = x[u];
		scc_y[scc_count] = y[u];
		scc_count++;
	}
}


int main() {
	ios::sync_with_stdio(false);
	cin.tie(0);
	
	cin >> n;
	for (int i = 0; i < n; i++) {
		cin >> x[i] >> y[i];
	}
	
	for (int i = 0; i < n; i++) {
		p[i] = make_pair(x[i], i);
	}
	
	sort(p, p + n);
	
	for (int i = 1; i < n; i++) {
		G[p[i].second].push_back(p[i - 1].second);
	}
	
	for (int i = 0; i < n; i++) {
		p[i] = make_pair(y[i], i);
	}
	
	sort(p, p + n);

	for (int i = 1; i < n; i++) {
		G[p[i].second].push_back(p[i - 1].second);
	}
	
	for (int i = 0; i < n; i++) {
		time_in[i] = 0;
		scc_id[i] = -1;
	}
	
	timer = 0;
	
	for (int i = 0; i < n; i++) {
		if (!time_in[i]) {
			dfs(i);
		}
	}
	
	for (int i = 0; i < scc_count; i++) {
		id[i] = i;
	}
	
	sort(id, id + scc_count, cmp_scc_x);
	
	for (int i = 1; i < scc_count; i++) {
		scc_size[id[i]] += scc_size[id[i - 1]];
	}
	
	for (int i = 0; i < n; i++) {
		cout << scc_size[scc_id[i]] - 1 << "\n";
	}
}
