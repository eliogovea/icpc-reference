#include <bits/stdc++.h>

using namespace st;

const int N = 100 * 1000 + 10;
const int Q = 50 * 1000 + 10;
const int LN = 20;

int n;
vector <int> G[N];

int q;
int qu[Q], qk[Q];

bool used[N];

int sub_tree_size[N];
int bfs_parent[N];
int bfs_queue[N];
int bfs_head;
int bfs_tail;

inline int get_centroid(int now) {
	bfs_head = 0; bfs_tail = 0;
	bfs_parent[now] = -1;
	bfs_queue[bfs_tail++] = now;
	while (bfs_head < bfs_tail) {
		int u = bfs_queue[bfs_head++];
		sub_tree_size[u] = 1;
		for (auto v : G[u]) {
			if (!used[v] && v != bfs_parent[u]) {
				bfs_parent[v] = u;
				bfs_queue[bfs_tail++] = v;
			}
		}
	}
	int component_size = bfs_head;
	for (int i = bfs_head - 1; i > 0; i--)
		sub_tree_size[bfs_parent[bfs_queue[i]]] += sub_tree_size[bfs_queue[i]];
	int centroid = -1; int best_max;
	for (int i = 0; i < bfs_head; i++) {
		int u = bfs_queue[i];
		int max_sub_tree_size = component_size - sub_tree_size[u];
		for (auto v : G[u])
			if (!used[v] && v != bfs_parent[u])
				max_sub_tree_size = max(max_sub_tree_size, sub_tree_size[v]);
		if (centroid == -1 || max_sub_tree_size < best_max)
			centroid = u, best_max = max_sub_tree_size;
	}
	return centroid;
}

int dist[LN][N];

inline void calc_dist(int now, int _depth) {
	bfs_head = 0; bfs_tail = 0;
	bfs_parent[now] = -1;
	dist[_depth][now] = 0;
	bfs_queue[bfs_tail++] = now;
	while (bfs_head < bfs_tail) {
		int u = bfs_queue[bfs_head++];
		for (auto v : G[u]) {
			if (!used[v] && v != bfs_parent[u]) {
				bfs_parent[v] = u;
				dist[_depth][v] = dist[_depth][u] + 1;
				bfs_queue[bfs_tail++] = v;
			}
		}
	}
}

int depth[N];
int parent[N];
int timer;
int time_in[N];
int time_out[N];

void _centroid_decomposition(int now, int _parent = -1, int _depth = 0) {
	now = get_centroid(now);
	calc_dist(now, _depth);
	parent[now] = _parent;
	depth[now] = _depth;
	time_in[now] = timer++;
	used[now] = true;	
	for (auto to : G[now]) if (!used[to])
		_centroid_decomposition(to, now, _depth + 1);	
	time_out[now] = timer - 1;
}

inline void centroid_decomposition() {
	for (int u = 0; u < n; u++)
		used[u] = false;
	_centroid_decomposition(0);
}

int main() {
    
}
