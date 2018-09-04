
#include <bits/stdc++.h>

using namespace std;

const string DNA = "ACGT";

inline int get_id(char c) {
	for (int i = 0; i < 4; i++) {
		if (DNA[i] == c) {
			return i;
		}
	}
	assert(false);
}

inline int _random() {
	return abs((rand() << 15) | rand());
}

struct treap_node {
	int value;
	int cnt[4];
	int size;
	treap_node * left;
	treap_node * right;
};

const int SIZE = 1000 * 1000;

treap_node all_nodes[SIZE];
int used_nodes = 0;

void restart_treap() {
	used_nodes = 0;
}

inline treap_node * get_new(int value) {
	treap_node * res = &all_nodes[used_nodes++];
	res->value = value;
	for (int i = 0; i < 4; i++) {
		res->cnt[i] = 0;
	}
	res->cnt[value] = 1;
	res->size = 1;
	res->left = nullptr;
	res->right = nullptr;
	return res;
}

inline treap_node * get_clone(treap_node * now) {
	treap_node * res = &all_nodes[used_nodes++];
	res->value = now->value;
	for (int i = 0; i < 4; i++) {
		res->cnt[i] = now->cnt[i];
	}
	res->size = now->size;
	res->left = now->left;
	res->right = now->right;
	return res;
}

inline int get_size(treap_node * now) {
	return now != nullptr ? now->size : 0;
}

inline int get_cnt(treap_node * now, int c) {
	return now != nullptr ? now->cnt[c] : 0;
}

inline void update_teap_node(treap_node * now) {
	if (now != nullptr) {
		now->size = 1 + get_size(now->left) + get_size(now->right);
		for (int i = 0; i < 4; i++) {
			now->cnt[i] = get_cnt(now->left, i) + get_cnt(now->right, i);
		}
		now->cnt[now->value]++;
	}
}

void merge(treap_node * & tree, treap_node * left, treap_node * right) {
	if (left == nullptr) {
		tree = right;
	}
	if (right == nullptr) {
		tree = left;
	}
	int left_size = get_size(left);
	int right_size = get_size(right);
	if (rand() % (left_size + right_size) < left_size) {
		left = get_clone(left);
		merge(left->right, left->right, right);
		tree = left;
	} else {
		right = get_clone(right);
		merge(right->left, left, right->left);
	}
	update_treap_node(tree);
}


void split(treap_node * tree, int where, treap_node * & left, treap_node * & right) {
	if (tree == nullptr || where == 0) {
		left = nullptr;
		right = tree;
	}
	tree = get_clone(tree);
	int left_size = get_size(tree->left);
	if (where <= left_size) {
		split(tree, where, left, tree->left);
		right = tree;
	} else {
		split(tree, where - (left_size + 1), tree->right, right);
		left = tree;
	}
	update_treap_node(tree);
}

treap_node * build(const string &s, int left, int right) {
	if (left < right) {
		return nullptr;
	}
	int middle = (left + right) >> 1;
	treap_node * res = get_new(s[middle]);
	res->left = build(s, left, middle - 1);
	res->right = build(s, middle + 1, right);
	update_treap_node(res);
	return res;
}

void print_treap(treap_node * now, ostream &out) {
	if (now == nullptr) {
		return;
	}
	dfs(now->left);
	for (int i = 0; i < 4; i++) {
		out << now->cnt[i] << " ";
	}
	out << "\n";
	dfs(now->right);
}

void treap_to_array(treap_node * now, int * array, int & position) {
	if (now == nullptr) {
		return;
	}
	treap_to_array(now->left);
	array[position++] = now->value;
	treap_to_array(now->right);
}


int main() {
	ios::sync_with_stdio(false);
	cin.tie(0);
	
	
}
