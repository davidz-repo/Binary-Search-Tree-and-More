#include<iostream>
#include<iomanip>

#include<vector>
#include<deque>
#include<unordered_set>

#include<algorithm>
#include<random>

#include<climits>
#include<cmath>
#include<chrono>

#define SEED 33

using namespace std;

class TP0 {
private:
	struct Node {
		int value;
        int priority;
		Node* left;
		Node* right;
        Node* parent;

		Node(int v, int bf, int prio){
			value = v;
            left = nullptr;
            right = nullptr;
            parent = nullptr;
			priority = prio;
		}
	};

    Node* rotateLeft( Node* A ){
        // get A's right
        Node* B = A->right;
        // break A's right link
        A->right = nullptr;
        // update A's parent's pointer if any
        if(A->parent){
            if(A->parent->left == A)
                A->parent->left = B;
            else
                A->parent->right = B;
            // update B's parent
            B->parent = A->parent;
        }
        else{
            B->parent = nullptr;
        }
        // update A's parent
        A->parent = B;
        // move B's left child as A's right child if any
        if(B->left){
            Node* C = B->left;
            C->parent = A;
            A->right = C;
        }
        B->left = A;
        if( A == root ) root = B;
        return B;
    }
    Node* rotateRight ( Node* A ){
        Node* B = A->left;
        A->left = nullptr;
        // update A's parent's pointer if any
        if(A->parent){
            if(A->parent->left == A)
                A->parent->left = B;
            else
                A->parent->right = B;
            // update B's parent
            B->parent = A->parent;
        }
        else{
            B->parent = nullptr;
        }
        // update A's parent
        A->parent = B;
        // move B's left child as A's right child if any
        if(B->right){
            Node* C = B->right;
            C->parent = A;
            A->left = C;
        }
        B->right = A;
        if( A == root ) root = B;
        return B;
    }
    void siftUp( Node* np ){
        while( np && np!=root && np->priority > np->parent->priority ){
            if( np->parent->left == np )
                np = rotateRight(np->parent);
            else
                np = rotateLeft(np->parent);
        }
    }
	void insert_private(int v, int prio, Node* np){
		Node* node_insert = new Node(v, 0, prio);
		Node* parent = nullptr;
		while(np){
			if(v > np->value){
                parent = np;
				np = np->right;
            }
			else if(v < np->value){
                parent = np;
				np = np->left;
            }
            else{
                cout << "!! node " << v << "already exist!!!" << endl;
                return;
            }
		}
		if(parent->value < v)
			parent->right = node_insert;
		else if(parent->value > v)
			parent->left = node_insert;
        else{
            cout << "!! node " << v << " already exist!!!" << endl;
            return;
        }
        // update the parent
        node_insert->parent = parent;
        // rebalance
        siftUp( node_insert );
        // confirmation message
        cout << "+ node " << v << endl;
	}
	void search_private(int v, Node* np){
		if(np){
			if(np->value > v)
				search_private(v, np->left);
			else if(np->value < v)
				search_private(v, np->right);
		    else if(np->value == v)
				cout << "= node " << v << endl;
		    else
		        cout << "!! NOT FOUND " << v << endl;
		}
		else
			cout << "Search not possible - Tree is empty!!" << endl;
	}
    void siftDown( Node*& np ){
        while( np->left && np->right ){
            if( np->left->priority > np->right->priority )
                np = rotateRight(np)->right;
            else
                np = rotateLeft(np)->left;
        }
        remove_private(np->value, np);
    }
	void remove_private(int v, Node* np){
        if( !root ){
            cout << "!! Empty tree - cannot remove " << v << endl;
            return;
        }
        Node* delete_node = nullptr;
        Node* start = np;
        Node* parent = nullptr;
        while( start ){
            if( v > start->value ){
                parent = start;
                start = start->right;
            }
            else if( v < start->value){
                parent = start;
                start = start->left;
            }
            else{
                delete_node = start;//found the node
                break;
            }
        }
        if( !delete_node ) {
            cout << "!! node with value " << v << "is not found!!!" << endl;
            return;
        }
        else {// found the node

            // if the parent is nullptr, try the delete_node's parent
            if( !parent ) parent = delete_node->parent;

            // get the children
            Node* left_child = delete_node->left;
            Node* right_child = delete_node->right;

            // case - no children
            if( !left_child && !right_child ){
                // no children just delete it
                if( delete_node != root ) {
                    if(parent->left == delete_node) parent->left = nullptr;
                    else if(parent->right == delete_node) parent->right = nullptr;
                    cout << "- node " << delete_node->value << endl;
                    delete delete_node;
                }
                else{
                    cout << "- LAST node " << delete_node->value << endl;
                    delete_node->value = -1;
                    return;
                }
            }

            // only have right child
            else if( !left_child && right_child ){
                // change the pointer to the right child
                if( root == delete_node ){
                    root = right_child;
                    right_child->parent == nullptr;
                    delete delete_node;
                }
                else {
                    if( parent->left == delete_node ) parent->left = right_child;
                    else if( parent->right == delete_node ) parent->right = right_child;
                    right_child->parent = parent;
                    cout << "- node " << delete_node->value << endl;
                    delete delete_node;
                }
            }

            // only have left child
            else if( left_child && !right_child ){// modify modify modify modify modify modify modify
                // change the pointer to the left child
                if( root == delete_node ){
                    root = left_child;
                    left_child->parent == nullptr;
                    delete delete_node;
                }
                else{
                    if( parent->left == delete_node ) parent->left = left_child;
                    else if( parent->right == delete_node ) parent->right = left_child;
                    left_child->parent = parent;
                    cout << "- node " << delete_node->value << endl;
                    delete delete_node;
                }
            }

            // both present - remove the inorder successor
            else {
                delete_node->priority = INT_MIN;
                siftDown( delete_node );
            }
        }
	}

	void postorder(Node* np){
		if (!np) return;
		postorder(np->left);
		postorder(np->right);
		cout << np->value << " ";
	}
	void preorder(Node* np){
		if (!np) return;
		cout << np->value << " ";
		preorder(np->left);
		preorder(np->right);
	}
	void inorder(Node* np){
		if (!np) return;
		inorder(np->left);
		cout << np->value << " ";
		inorder(np->right);
	}
	void padding( char ch, int n ){
        for( ; n>0; --n ){
            cout << ch;
        }
    }
    void print_tree(Node* np){
        if(np && root){
            deque<Node*> d;
            int height = max_height(np);
            cout << "tree height : " << height << endl;
            int max_count = (int)(pow(2, height) - 1);
            vector<int> qv(max_count, -1);
            vector<int> qp(max_count, -1);
            qv[0] = np->value;
            qp[0] = np->priority;
            d.push_back(np);
            int pad = (height+1) * 13;
            while(!d.empty()){
                Node* current = d.front();
                d.pop_front();
                int i = 0;
                for(int pos = 0; pos < qv.size(); pos++)
                    if(qv[pos] == current->value)
                        i = pos;
                if(current->left){
                    qv[2 * i + 1] = current->left->value;
                    qp[2 * i + 1] = current->left->priority;
                    d.push_back(current->left);
                }
                if(current->right){
                    qv[2 * i + 2] = current->right->value;
                    qp[2 * i + 2] = current->right->priority;
                    d.push_back(current->right);
                }
            }
            for(int n = 0; n < max_count; n++){
                if(qv[n] != -1){
                    cout  << "[" << qv[n] << setw(8) << qp[n] << "] ";
                }
                else{
                    cout << "[NULL  NULL] ";
                }
                padding( ' ', pad );
                if(!((n + 2) & (n + 1))){
                    height -= 2;
                    pad = (height+1) * 13;
                    cout << endl << endl;
                }
            }
        }
    }
	void free_tree(Node* np){
		if(np){
            if(np->left){
                free_tree(np->left);
            }
            if(np->right){
                free_tree(np->right);
            }
            cout << "destroyed node " << np->value << endl;
            delete np;
        }
        else{
            cout << "~~~~~Destructor ERROR!!!" << endl;
        }
	}


public:
	~TP0() {
		free_tree(root);
	}
	Node* root;
    int population = 0;
    deque<int> priority_deque;
	int max_height(Node* np){
        if(!np)
            return -1;
        else {
            int hl=-1, hr=-1;
            if( np->left )
                hl = max_height(np->left);
            if( np->right )
                hr = max_height(np->right);
            return 1 + max( hl, hr );
        }
    }
	void create(int count){
        root = nullptr;
        unordered_set<int> priority_gen(count);
        priority_deque.resize(count);
        while( priority_gen.size() < count ){
        	priority_gen.insert( rand() % (count*100) );
        }
        copy( priority_gen.begin(), priority_gen.end(), priority_deque.begin() );
        priority_gen.clear();
        shuffle( priority_deque.begin(), priority_deque.end(), default_random_engine(SEED) );
    }

	void insert(int v, int p){
		if(root)
			insert_private(v, p, root);
		else {
            root = new Node(v, 0, p);
            cout << "+ ROOT " << v << endl;
        }
        ++population;
	}

	void search(int v) {search_private(v, root);}

	void remove(int v){
		remove_private(v, root);
        --population;
	}

	void get_population(){cout << "current population: " << this->population << endl;}

	void display(const string& flag){
        cout << "\n::::::::::::::::: Tree Status ::::::::::::::::::" << endl;
        get_population();
		if(flag == "preorder"){
            cout << "Preorder array: ";
			preorder(root);
        }
		else if(flag == "inorder"){
            cout << "Inorder array: ";
			inorder(root);
        }
		else if(flag == "postorder"){
            cout << "Postorder array: ";
			postorder(root);
        }
        else if(flag == "2d"){
            cout << endl;
            print_tree(root);
        }
        cout << ":::::::::::::::: End Of Status ::::::::::::::::::\n" << endl;
	}

};

class TPf {
/* treap using regular bst insertion (pre-sorted based on priority)
 * treap using regular bst deletion (without rotation) */
private:
	struct Node {
		int value;
        int priority;
        int balance_factor; // delete delete delete delete
		Node* left;
		Node* right;
        Node* parent;

		Node(int v, int bf, int prio){
			value = v;
			left = nullptr;
			right = nullptr;
			parent = nullptr;
			balance_factor = bf;
			priority = prio;
		}
	};
	void insert_private(int v, int prio, Node* np){
		Node* node_insert = new Node(v, 0, prio);
		Node* start = np;
		Node* parent = nullptr;
		while(start){
			if(v > start->value){
                parent = start;
				start = start->right;
            }
			else if(v < start->value){
                parent = start;
				start = start->left;
            }
            else{
                cout << "!! node " << v << " already exist!!!" << endl;
                return;
            }
		}
		if(parent->value < v)
			parent->right = node_insert;
		else if(parent->value > v)
			parent->left = node_insert;
        else{
            cout << "!! node " << v << " already exist!!!" << endl;
            return;
        }
        // update the parent
        node_insert->parent = parent;
        // confirmation message
        cout << "+ node " << v << endl;
	}

	void search_private(int v, Node* np){
		if(np){
			if(np->value > v)
				search_private(v, np->left);
			else if(np->value < v)
				search_private(v, np->right);
		    else if(np->value == v)
				cout << "= node " << v << endl;
		    else
		        cout << "!! NOT FOUND " << v << endl;
		}
		else
			cout << "Search not possible - Tree is empty!!" << endl;
	}
	Node* helper_findmin(Node* np){
        // helper finding the min node starting from np
		Node* min_node = np;
		while( min_node->left )
			min_node = min_node->left;
		return min_node;
	}
    void helper_remove_inorder_successor( Node* np ){
        // helper to remove the node that has both children
        // locate the inorder successor of np's right children
        Node* delete_node = helper_findmin(np->right);
        // swap the values of np and the ssr
        int v = np->value;
        np->value = delete_node->value;
        delete_node->value = v;
        // call remove to remove the ssr starting from the right children
        remove_private( v, np->right );
    }
	void remove_private(int v, Node* np){
        if( !root ){
            cout << "!! Empty tree - cannot remove " << v << endl;
            return;
        }
        Node* delete_node = nullptr;
        Node* start = np;
        Node* parent = nullptr;
        while( start ){
            if( v > start->value ){
                parent = start;
                start = start->right;
            }
            else if( v < start->value){
                parent = start;
                start = start->left;
            }
            else{
                delete_node = start;//found the node
                break;
            }
        }
        if( !delete_node ) {
            cout << "!! node with value " << v << "is not found!!!" << endl;
            return;
        }
        else {// found the node

            // if the parent is nullptr, try the delete_node's parent
            if( !parent ) parent = delete_node->parent;

            // get the children
            Node* left_child = delete_node->left;
            Node* right_child = delete_node->right;

            // case - no children
            if( !left_child && !right_child ){
                // no children just delete it
                if( delete_node != root ) {
                    if(parent->left == delete_node) parent->left = nullptr;
                    else if(parent->right == delete_node) parent->right = nullptr;
                    cout << "- node " << delete_node->value << endl;
                    delete delete_node;
                }
                else{
                    cout << "- LAST node " << delete_node->value << endl;
                    delete_node->value = -1;
                    return;
                }
            }

            // only have right child
            else if( !left_child && right_child ){
                // change the pointer to the right child
                if( root == delete_node ){
                    root = right_child;
                    right_child->parent == nullptr;
                    delete delete_node;
                }
                else {
                    if( parent->left == delete_node ) parent->left = right_child;
                    else if( parent->right == delete_node ) parent->right = right_child;
                    right_child->parent = parent;
                    cout << "- node " << delete_node->value << endl;
                    delete delete_node;
                }
            }

            // only have left child
            else if( left_child && !right_child ){// modify modify modify modify modify modify modify
                // change the pointer to the left child
                if( root == delete_node ){
                    root = left_child;
                    left_child->parent == nullptr;
                    delete delete_node;
                }
                else{
                    if( parent->left == delete_node ) parent->left = left_child;
                    else if( parent->right == delete_node ) parent->right = left_child;
                    left_child->parent = parent;
                    cout << "- node " << delete_node->value << endl;
                    delete delete_node;
                }
            }

            // both present - remove the inorder successor
            else {
                helper_remove_inorder_successor( delete_node );
            }
        }
	}

	void postorder(Node* np){
		if (!np) return;
		postorder(np->left);
		postorder(np->right);
		cout << np->value << " ";
	}
	void preorder(Node* np){
		if (!np) return;
		cout << np->value << " ";
		preorder(np->left);
		preorder(np->right);
	}
	void inorder(Node* np){
		if (!np) return;
		inorder(np->left);
		cout << np->value << " ";
		inorder(np->right);
	}
	void padding( char ch, int n ){
        for( ; n>0; --n ){
            cout << ch;
        }
    }
    void print_tree(Node* np){
        if(np && root){
            deque<Node*> d;
            int height = max_height(np);
            cout << "tree height : " << height << endl;
            int max_count = (int)(pow(2, height) - 1);
            vector<int> qv(max_count, -1);
            vector<int> qp(max_count, -1);
            qv[0] = np->value;
            qp[0] = np->priority;
            d.push_back(np);
            int pad = (height+1) * 13;
            while(!d.empty()){
                Node* current = d.front();
                d.pop_front();
                int i = 0;
                for(int pos = 0; pos < qv.size(); pos++)
                    if(qv[pos] == current->value)
                        i = pos;
                if(current->left){
                    qv[2 * i + 1] = current->left->value;
                    qp[2 * i + 1] = current->left->priority;
                    d.push_back(current->left);
                }
                if(current->right){
                    qv[2 * i + 2] = current->right->value;
                    qp[2 * i + 2] = current->right->priority;
                    d.push_back(current->right);
                }
            }
            for(int n = 0; n < max_count; n++){
                if(qv[n] != -1){
                    cout  << "[" << qv[n] << setw(8) << qp[n] << "] ";
                }
                else{
                    cout << "[NULL  NULL] ";
                }
                padding( ' ', pad );
                if(!((n + 2) & (n + 1))){
                    height -= 2;
                    pad = (height+1) * 13;
                    cout << endl << endl;
                }
            }
        }
    }
	void free_tree(Node* np){
		if(np){
            if(np->left){
                free_tree(np->left);
            }
            if(np->right){
                free_tree(np->right);
            }
            cout << "destroyed node " << np->value << endl;
            delete np;
        }
        else{
            cout << "~~~~~Destructor ERROR!!!" << endl;
        }
	}


public:
	~TPf() {
		free_tree(root);
	}
	Node* root;
    int population = 0;

	int max_height(Node* np){
        if(!np)
            return -1;
        else {
            int hl=-1, hr=-1;
            if( np->left )
                hl = max_height(np->left);
            if( np->right )
                hr = max_height(np->right);
            return 1 + max( hl, hr );
        }
    }
    deque<int> priority_deque;

	void create(int count){
        root = nullptr;
        unordered_set<int> priority_gen;
        priority_deque.resize(count);
        while( priority_gen.size() < count ){
        	priority_gen.insert( rand() % (count*100) );
        }
        copy( priority_gen.begin(), priority_gen.end(), priority_deque.begin() );
        priority_gen.clear();
        shuffle( priority_deque.begin(), priority_deque.end(), default_random_engine(SEED) );
    }

	void insert(int v, int p){
		if(root)
			insert_private(v, p, root);
		else {
            root = new Node(v, 0, p);
            cout << "+ ROOT " << v << endl;
        }
        ++population;
	}

	void search(int v) {search_private(v, root);}

	void remove(int v){
		remove_private(v, root);
        --population;
	}

	void get_population(){cout << "current population: " << this->population << endl;}

	void display(const string& flag){
        cout << "\n::::::::::::::::: Tree Status ::::::::::::::::::" << endl;
        get_population();
		if(flag == "preorder"){
            cout << "Preorder array: ";
			preorder(root);
        }
		else if(flag == "inorder"){
            cout << "Inorder array: ";
			inorder(root);
        }
		else if(flag == "postorder"){
            cout << "Postorder array: ";
			postorder(root);
        }
        else if(flag == "2d"){
            cout << endl;
            print_tree(root);
        }
        cout << ":::::::::::::::: End Of Status ::::::::::::::::::\n" << endl;
	}

};

void test( int size ){
	TP0 tp0 = TP0();
	TPf tpf = TPf();

    int Tsize = size * 1.3, count = 0, c = 0;
    int data[size] = {0}, pool[Tsize] = {0};
    unordered_set<int> record(Tsize);

    pool[0] = 1;
    for( int i=1; i<Tsize; ++i) pool[i] = rand()%20 + pool[i-1];
    shuffle(&pool[0], &pool[Tsize-1], default_random_engine(SEED));
    while( count < size ){
        if( record.find(pool[c]) == record.end() ){
            record.insert(pool[c]);
            data[count++] = pool[c++];
        }
        else c++;
    }
    cout << "\nRandom number collisions: " << c-count << endl;
    cout << "Data population for tree: " << count << endl;
    cout << "Input: "; for(int i=0; i<size; ++i) cout << data[i] << " ";
    cout << endl << endl;

	cout << "\n=================== Create Tree ====================" << endl;
	tp0.create(size);
	tpf.create(size);
	cout << "==================================================\n" << endl;

    cout << "\n================== Insertion start ====================" << endl;
    vector<pair<int,int>> elements(size);
	for(int i=0; i<size; ++i){
        elements[i] = make_pair(data[i],tp0.priority_deque[i]);
    }
    struct mycompare{
        bool operator() (const pair<int,int>& a, const pair<int,int>& b){
            return a.second > b.second;
        }
    };
    for( int i=0; i<size; ++i ) tp0.insert(elements[i].first, elements[i].second);
    sort( elements.begin(), elements.end(), mycompare() );
    for( int i=0; i<size; ++i ) tpf.insert(elements[i].first, elements[i].second);
    cout << "================= Insertion complete ===================\n" << endl;


    cout << "\n=================== Search start ====================" << endl;
	for(int i=0; i<size; ++i){tp0.search(data[i]);}
	for(int i=0; i<size; ++i){tpf.search(data[i]);}
    cout << "==================== Search end =====================\n" << endl;

	cout << "Treap Height: " << tp0.max_height(tp0.root) << endl;
	cout << "Fast Treap H: " << tpf.max_height(tpf.root) << endl << endl;

	tpf.display("postorder");
	tp0.display("postorder");

    cout << "\n=================== Deletion start ====================" << endl;
	for(int i=0; i<size; ++i) tp0.remove(data[i]);
	for(int i=0; i<size; ++i) tpf.remove(data[i]);
    cout << "==================== Deletion end ====================\n" << endl;


    cout << "~~~~ recycle tree (destructor) ~~~ " << endl;
}

int main(){
    int size = 10000;
    test(size);
}
