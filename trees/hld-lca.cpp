#include <bits/stdc++.h>
 
using namespace std;
 
// HLD + Euler Tour Trick (overkill, HLD only, is good enough)
const int N = 100 * 1000 + 10;
int n; vector<int> g[N];
int size[N], depth[N], par[N];
int head[N], chain[N], heavy[N];
int tchains;
vector <int> nodes;
int tin[N], tout[N];
 
void dfs(int u, int p, int d) {
    par[u] = p; depth[u] = d; 
    size[u] = 1; heavy[u] = -1;
    for (auto v : g[u]) if (v != p) {
        dfs(v, u, d + 1); size[u] += size[v];
        if (heavy[u] == -1 || size[v] > size[heavy[u]])
            heavy[u] = v;
    }
}
 
void build_hld(int root, int u, int p) {
    head[u] = root; chain[u] = tchains;
    tin[u] = nodes.size(); nodes.push_back(u);
    if (heavy[u] != -1) build_hld(root, heavy[u], u);
    for (auto v : g[u]) if (v != p && v != heavy[u])
        tchains++, build_hld(v, v, u);
    tout[u] = nodes.size() - 1;
}
 
void hld() { 
    tchains = 0; nodes.clear(); 
    dfs(0, 0, 0); 
    build_hld(0, 0, -1);
}
 
int lca(int u, int v) {
    while (chain[u] != chain[v]) {
        if (depth[head[v]] > depth[head[u]]) swap(u, v);
        u = par[head[u]];
    }
    if (depth[u] > depth[v]) swap(u, v); return u;
}
 
void test_hld() {
    cin >> n;
    for (int i = 0; i < n; i++) {
        g[i].clear();
    }
    for (int i = 0; i < n; i++) {
        int c;
        cin >> c;
        while (c--) {
            int x;
            cin >> x;
            g[i].push_back(x);
            g[x].push_back(i);
        }
    }
    hld();
    int m;
    cin >> m;
    while (m--) {
        int x, y;
        cin >> x >> y;
        cout << lca(x, y) << "\n";
    }
}
 
int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);
     
    test_hld();
}
