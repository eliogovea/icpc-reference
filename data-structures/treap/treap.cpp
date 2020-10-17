// treap
// solution for problem https://www.hackerrank.com/challenges/median/problem

#include <bits/stdc++.h>

using namespace std;

template <class T>
struct treap_node {
    T key;
    int priority;
    int sub_tree_size;
    treap_node<T> *left;
    treap_node<T> *right;
};

template <class T>
int get_sub_tree_size(treap_node<T> *now) {
    return now == nullptr ? 0 : now->sub_tree_size;
}

template <class T>
void update_sub_tree_size(treap_node<T> *now) {
    if (now != nullptr) {
        now->sub_tree_size =
          1 + get_sub_tree_size(now->left) + get_sub_tree_size(now->right);
    }
}

template <class T>
treap_node<T> *merge(treap_node<T> *left, treap_node<T> *right) {
    if (left == nullptr) { return right; }
    if (right == nullptr) { return left; }
    if (left->priority > right->priority) {
        left->right = merge(left->right, right);
        update_sub_tree_size(left);
        return left;
    }
    right->left = merge(left, right->left);
    update_sub_tree_size(right);
    return right;
}

template <class T>
pair<treap_node<T> *, treap_node<T> *> split(treap_node<T> *tree, T key) {
    if (tree == nullptr) { return make_pair(nullptr, nullptr); }
    if (tree->key <= key) {
        pair<treap_node<T> *, treap_node<T> *> right_split =
          split(tree->right, key);
        tree->right = right_split.first;
        update_sub_tree_size(tree);
        return make_pair(tree, right_split.second);
    }
    pair<treap_node<T> *, treap_node<T> *> left_split = split(tree->left, key);
    tree->left = left_split.second;
    update_sub_tree_size(tree);
    return make_pair(left_split.first, tree);
}

template <class T>
struct treap {
    treap_node<T> *root;

    treap() { root = nullptr; }

    treap_node<T> *get_new_node(T key) {
        treap_node<T> *res = new treap_node<T>();
        res->key = key;
        res->priority = rand() * rand();
        res->sub_tree_size = 1;
        res->left = nullptr;
        res->right = nullptr;
        return res;
    }

    treap_node<T> *find(T key) {
        treap_node<T> *now = root;
        while (true) {
            if (now == nullptr) { return nullptr; }
            if (now->key == key) { return now; }
            if (key <= now->key) {
                now = now->left;
            } else {
                now = now->right;
            }
        }
    }

    treap_node<T> *insert(treap_node<T> *now, treap_node<T> *new_node) {
        if (now == nullptr) { return new_node; }
        if (new_node->priority > now->priority) {
            pair<treap_node<T> *, treap_node<T> *> s =
              split(now, new_node->key);
            new_node->left = s.first;
            new_node->right = s.second;
            update_sub_tree_size(new_node);
            return new_node;
        }
        if (new_node->key <= now->key) {
            now->left = insert(now->left, new_node);
        } else {
            now->right = insert(now->right, new_node);
        }
        update_sub_tree_size(now);
        return now;
    }

    void insert(T key) { root = insert(root, get_new_node(key)); }

    treap_node<T> *erase(treap_node<T> *now, T key) {
        if (now->key == key) {
            treap_node<T> *res = merge(now->left, now->right);
            delete now;
            return res;
        }
        if (key <= now->key) {
            now->left = erase(now->left, key);
        } else {
            now->right = erase(now->right, key);
        }
        update_sub_tree_size(now);
        return now;
    }

    bool erase(T key) {
        if (find(key) == nullptr) { return false; }
        root = erase(root, key);
        return true;
    }

    T find_kth(int k) {
        assert(0 < k && k <= get_sub_tree_size(root));
        treap_node<T> *now = root;
        while (true) {
            int left_size = get_sub_tree_size(now->left);
            if (left_size >= k) {
                now = now->left;
            } else {
                k -= left_size;
                if (k == 1) { return now->key; }
                k--;
                now = now->right;
            }
        }
    }

    inline int size() { return get_sub_tree_size(root); }

    void print_treap(treap_node<T> *now) {
        if (now == nullptr) { return; }
        print_treap(now->left);
        cerr << now->key << " ";
        print_treap(now->right);
    }

    void print() {
        print_treap(root);
        cerr << "\n";
    }
};

void test() {
    treap<int> t = treap<int>();
    t.insert(1);
    exit(0);
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);

    // test();

    treap<int> t = treap<int>();

    int n;
    cin >> n;

    while (n--) {
        char c;
        int x;
        cin >> c >> x;

        if (c == 'r') {
            if (!t.erase(x)) {
                cout << "Wrong!\n";
                continue;
            }
        } else {
            t.insert(x);
        }
        if (t.size() == 0) {
            cout << "Wrong!\n";
        } else if (t.size() & 1) {
            cout << t.find_kth((t.size() + 1) >> 1) << "\n";
        } else {
            int a = t.size() >> 1;
            int b = a + 1;
            long long m1 = t.find_kth(a);
            long long m2 = t.find_kth(b);
            if (m1 + m2 == -1LL) {  // !!!
                cout << "-0.5\n";
            } else {
                cout << (m1 + m2) / 2LL;
                if ((m1 + m2) % 2LL) { cout << ".5"; }
                cout << "\n";
            }
        }

        // t.print();
    }
}
