// coj-2418

#include <bits/stdc++.h>

using namespace std;

const int ALPHABET = 26;

struct treap_node {
	// node information
	int value;
	int priority;
	treap_node * left_child;
	treap_node * right_child;

	// sub-tree information
	int size;
	int cnt[ALPHABET];

	// pending operations
	int lazy_reverse;
	int lazy_set;
	
	treap_node(int _value) {
		value = _value;
		priority = (rand() << 15) | rand();
		left_child = nullptr;
		right_child = nullptr;

		size = 1;
		for (int c = 0; c < ALPHABET; c++) {
			cnt[c] = 0;
		}
		cnt[value] = 1;

		lazy_reverse = 0;
		lazy_set = -1;
	}
};

inline int get_size(treap_node * now) {
	return now != nullptr ? now->size : 0;
}

inline int get_cnt(treap_node * now, int c) {
	return now != nullptr ? now->cnt[c] : 0;
}

inline void add_lazy_set(treap_node * now, int c) {
	if (now != nullptr) {
		now->lazy_reverse = 0;
		now->lazy_set = c;
	}
}

inline void add_lazy_reverse(treap_node * now) {
	if (now != nullptr && now->lazy_set == -1) {
		now->lazy_reverse ^= 1;
	}
}

inline void push_lazy(treap_node * now) {
	if (now != nullptr) {
		if (now->lazy_reverse) {
			swap(now->left_child, now->right_child);
			add_lazy_reverse(now->left_child);
			add_lazy_reverse(now->right_child);
			now->lazy_reverse = 0;
		}
		if (now->lazy_set != -1) {
			now->value = now->lazy_set;
			for (int c = 0; c < ALPHABET; c++) {
				now->cnt[c] = 0;
			}
			now->cnt[now->value] = now->size;
			add_lazy_set(now->left_child, now->lazy_set);
			add_lazy_set(now->right_child, now->lazy_set);
			now->lazy_set = -1;
		}
	}
}

void update_treap_node(treap_node * now) {
	if (now != nullptr) {
		// assert(now->lazy_set == -1 && now->lazy_reverse == 0);
		now->size = 1 + get_size(now->left_child) + get_size(now->right_child);
		for (int c = 0; c < ALPHABET; c++) {
			now->cnt[c] = get_cnt(now->left_child, c) + get_cnt(now->right_child, c);
		}
		now->cnt[now->value]++;
	}
}

void merge(treap_node * & root, treap_node * left, treap_node * right) {
	if (left == nullptr) {
		push_lazy(right);
		root = right;
	} else if (right == nullptr) {
		push_lazy(left);
		root = left;
	} else {
		if (left->priority > right->priority) {
			push_lazy(left);
			push_lazy(left->left_child);
			merge(left->right_child, left->right_child, right);
			root = left;
		} else {
			push_lazy(right);
			push_lazy(right->right_child);
			merge(right->left_child, left, right->left_child);
			root = right;
		}
		update_treap_node(root);
	}
}

void split(treap_node * root, int where, treap_node * & left, treap_node * & right) {
	push_lazy(root);
	if (where == 0) {
		left = nullptr;
		right = root;
	} else {
		int left_size = get_size(root->left_child);
		if (where <= left_size) {
			push_lazy(root->right_child);
			split(root->left_child, where, left, root->left_child);
			right = root;
		} else {
			push_lazy(root->left_child);
			split(root->right_child, where - (left_size + 1), root->right_child, right);
			left = root;
		}
		update_treap_node(root);
	}
}

void print_treap(treap_node * now, ostream & out) {
	if (now == nullptr) {
		return;
	}
	push_lazy(now);
	print_treap(now->left_child, out);
	out << char('a' + now->value);
	print_treap(now->right_child, out);
}

void ivan_operations() {
	string s;
	cin >> s;
	
	// cerr << s << "\n";
	
	treap_node * root = nullptr;
	
	for (int i = 0; i < (int)s.size(); i++) {
		merge(root, root, new treap_node(s[i] - 'a'));
	}
	
	// print_treap(root, cout); cerr << "\n";
	
	int q;
	cin >> q;
	
	while (q--) {
		char op;
		int l, r;
		cin >> op >> l >> r;
		
		// cerr << "new query\n";
		// cerr << "array: "; print_treap(root, cerr); cerr << "\n";
		// cerr << op << " " << l << " " << r << "\n";
		
		treap_node * A;
		treap_node * B;
		treap_node * C;
		
		split(root, r, B, C);
		split(B, l - 1, A, B);
		
		// cerr << "( "; print_treap(A, cerr); cerr << ")\n";
		// cerr << "( "; print_treap(B, cerr); cerr << ")\n";
		// cerr << "( "; print_treap(C, cerr); cerr << ")\n";
		
		if (op == 'R') {
			add_lazy_reverse(B);
		} else {
			char c;
			cin >> c;
			
			// cerr << c << "\n";
			
			if (op == 'S') {
				add_lazy_set(B, c - 'a');
			} else {
				push_lazy(B); // !!!
				
				// cout << "-->> "; print_treap(B, cout); cout << "\n";
				// cerr << "answer: " << get_cnt(B, c - 'a') << "\n";
				cout << get_cnt(B, c - 'a') << "\n";
			}
		}
		
		merge(B, A, B);
		merge(root, B, C);
		
		// print_treap(root, cout); cout << "\n";
		
		// cerr << "array: "; print_treap(root, cerr); cerr << "\n\n";
	}
}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0);
	
	ivan_operations();
}
