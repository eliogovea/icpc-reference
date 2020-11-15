// Suffix Automaton
struct sa_node {
    int max_length;
    map<char, sa_node*> go;
    sa_node* suffix_link;
    vector<pair<char, sa_node*>> inv_suffix_link;
    int id, pos;
    bool ok, visited;
    sa_node() {}
};

const int N = 500 * 1000 + 10;
const int SIZE = 2 * N + 10;

sa_node all[SIZE];
sa_node* first_free;

inline sa_node* get_new(int max_length, int id, int pos) {
    first_free->max_length = max_length;
    first_free->go = map<char, sa_node*>();
    first_free->suffix_link = nullptr;
    first_free->inv_suffix_link = vector<pair<char, sa_node*>>();
    first_free->id = id;
    first_free->pos = pos;
    first_free->visited = false;
    return first_free++;
}

inline sa_node* sa_init() {
    first_free = all;
    return get_new(0, 0, 0);
}

inline sa_node* get_clone(sa_node* u, int max_length) {
    first_free->max_length = max_length;
    first_free->go = u->go;
    first_free->suffix_link = u->suffix_link;
    first_free->inv_suffix_link = vector<pair<char, sa_node*>>();
    first_free->id = u->id;
    first_free->pos = u->pos;
    first_free->ok = u->ok;
    first_free->visited = false;
    return first_free++;
}

inline sa_node* add(sa_node* root, sa_node* p, char c, int id, int pos) {
    if (p->go.find(c) != p->go.end()) {
        sa_node* q = p->go[c];
        if (p->max_length + 1 == q->max_length) return q;
        sa_node* clone_q = get_clone(q, p->max_length + 1);
        q->suffix_link = clone_q;
        while (p != nullptr && p->go[c] == q) {
            p->go[c] = clone_q, p = p->suffix_link;
        }
        return clone_q;
    }
    sa_node* l = get_new(p->max_length + 1, id, pos);
    while (p != nullptr && p->go.find(c) == p->go.end())
        p->go[c] = l, p = p->suffix_link;
    if (p == nullptr) {
        l->suffix_link = root;
    } else {
        sa_node* q = p->go[c];
        if (p->max_length + 1 == q->max_length) {
            l->suffix_link = q;
        } else {
            auto clone_q = get_clone(q, p->max_length + 1);
            l->suffix_link = q->suffix_link = clone_q;
            while (p != nullptr && p->go[c] == q) {
                p->go[c] = clone_q;
                p = p->suffix_link;
            }
        }
    }
    return l;
}

sa_node* build(vector<string>& s) {
    int n = s.size();
    sa_node* root = sa_init();
    sa_node* last = root;
    for (int i = 0; i < n; i++) {
        reverse(s[i].begin(), s[i].end());
        for (int p = 0; p < (int)s[i].size(); p++)
            last = add(root, last, s[i][p], i, p);
        last = root;
    }
    for (sa_node* now = all; now != first_free; now++) {
        if (now->suffix_link != nullptr) {
            now->suffix_link->inv_suffix_link.push_back(make_pair(
              s[now->id][now->pos - now->suffix_link->max_length], now));
        }
    }
    for (sa_node* now = all; now != first_free; now++)
        sort(now->inv_suffix_link.begin(), now->inv_suffix_link.end());

    // solution

    return root;
}
