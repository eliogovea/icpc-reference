// useful for answer queries on subtrees and paths

const int N = 100 * 1000 + 10;
int n;
vector<int> g[N];
int size[N], depth[N], par[N];
int head[N], chain[N], heavy[N];
int tchains;
vector<int> nodes;
int tin[N], tout[N];

void dfs(int u, int p, int d) {
    par[u] = p;
    depth[u] = d;
    size[u] = 1;
    heavy[u] = -1;
    for (auto v : g[u])
        if (v != p) {
            dfs(v, u, d + 1);
            size[u] += size[v];
            if (heavy[u] == -1 || size[v] > size[heavy[u]]) heavy[u] = v;
        }
}

void build_hld(int root, int u, int p) {
    head[u] = root;
    chain[u] = tchains;
    tin[u] = nodes.size();
    nodes.push_back(u);
    if (heavy[u] != -1) build_hld(root, heavy[u], u);
    for (auto v : g[u])
        if (v != p && v != heavy[u]) tchains++, build_hld(v, v, u);
    tout[u] = nodes.size() - 1;
}

void hld() {
    tchains = 0;
    nodes.clear();
    dfs(0, 0, 0);
    build_hld(0, 0, -1);
}

int queryOnPath(int u, int v) {  // ex. max
    int ans = -INF;
    while (chain[u] != chain[v]) {
        if (depth[head[v]] > depth[head[u]]) swap(u, v);
        int lo = tin[head[u]];
        int hi = lo + depth[u] - depth[head[u]];
        ans = max(ans, query(1, 0, n - 1, lo, hi));
        u = par[head[u]];
    }
    if (depth[u] > depth[v]) swap(u, v);
    // u -> lca!!!!!!!
    int lo = tin[head[u]] + depth[u] - depth[head[u]];
    int hi = tin[head[v]] + depth[v] - depth[head[v]];
    ans = max(ans, query(1, 0, n - 1, lo, hi));
    return ans;
}
