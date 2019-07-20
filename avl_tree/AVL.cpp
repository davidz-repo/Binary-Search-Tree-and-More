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

class AVL {
private:
	struct Node {
		int value;
        int balance_factor;
		Node* left;
		Node* right;
        Node* parent;

		Node(int v, int bf){
			value = v;
			left = nullptr;
			right = nullptr;
			parent = nullptr;
            balance_factor = bf;
		}
	};
	Node* root;
    int population = 0;

    int max_height(Node* np){
        if(!np)
            return -1;
        else {
            int hl=0, hr=0;
            if( np->left )
                hl = max_height(np->left);
            if( np->right )
                hr = max_height(np->right);
            return 1 + max( hl, hr );
        }
    }
    void update_bf( Node* np){
        int hl=0, hr=0;
        if( np->left )
            hl = max_height(np->left);
        if( np->right )
            hr = max_height(np->right);
        np->balance_factor = hl - hr;
    }
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
			root = B;
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
        // update_bf( B );
        // update_bf( A );
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
			root = B;
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
        // update_bf( B );
        // update_bf( A );
        return B;
    }
    Node* rotateLeftThenRight ( Node* A ){
        A->left = rotateLeft( A->left );
        return rotateRight(A);
    }
    Node* rotateRightThenLeft ( Node* A ){
        A->right = rotateRight( A->right );
        return rotateLeft(A);
    }
    void rebalance_tree( Node* np ){
        update_bf( np );
        if( np->balance_factor == 2 ){
            if( np->left && max_height(np->left->right) > max_height(np->left->left) ){
                np = rotateLeftThenRight( np );
            }
            else {
                np = rotateRight( np );
            }
        }
        else if ( np->balance_factor == -2 ){
            if( (np->right) && max_height(np->right->left) > max_height(np->right->right) ){
                np = rotateRightThenLeft( np );
            }
            else {
                np = rotateLeft( np );
            }
        }

        if (np->parent)
            rebalance_tree( np->parent );
        else if (np!=root)
		// if ( np!=root && !np->parent )
            root = np;
        else; // already a root;
    }
	void insert_private(int v, Node* np){
		Node* node_insert = new Node(v, 0);
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
        // rebalance
        rebalance_tree( parent );
        // confirmation message
        // cout << "+ node " << v << endl;
	}

	void search_private(int v, Node* np){
		if(np){
			if(np->value > v)
				search_private(v, np->left);
			else if(np->value < v)
				search_private(v, np->right);
		    else if(np->value == v){
				// cout << "= node " << v << endl;
			}
		    else{
		        cout << "!! NOT FOUND " << v << endl;
			}
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
                    // cout << "- node " << delete_node->value << endl;
                    delete delete_node;
                    rebalance_tree(parent);
                }
                else{
                    // cout << "- LAST node " << delete_node->value << endl;
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
                    // cout << "- node " << delete_node->value << endl;
                    delete delete_node;
                    rebalance_tree(parent);
                }
            }

            // only have left child
            else if( left_child && !right_child ){
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
                    // cout << "- node " << delete_node->value << endl;
                    delete delete_node;
                    rebalance_tree(parent);
                }
            }

            // both present - remove the inorder successor
            else helper_remove_inorder_successor( delete_node );
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
            int H = max_height(np);
            cout << "tree height : " << H << endl;
            int dimension = (int)(pow((float)2, H+1)-1);
            vector<int> nodes(dimension, -1);
            nodes[0] = np->value;
            d.push_back(np);
            int pad = (H+1) * 7;
            while(!d.empty()){
                Node* current = d.front();
                d.pop_front();
                int position;
                for(int pos = 0; pos < nodes.size(); pos++)
                    if(nodes[pos] == current->value)
                        position = pos;
                if(current->left){
                    nodes[2 * position + 1] = current->left->value;
                    d.push_back(current->left);
                }
                if(current->right){
                    nodes[2 * position + 2] = current->right->value;
                    d.push_back(current->right);
                }
            }
            for(int n = 0; n < dimension; n++){
                if(nodes[n] != -1){
                    cout  << "[" << setw(4) << nodes[n] << "] ";
                }
                else{
                    cout << "[NULL] ";
                }
                padding( ' ', pad );
                if(!((n + 2) & (n + 1))){
                    H -= 2;
                    pad = (H+1) * 7;
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
            // cout << "destroyed node " << np->value << endl;
            delete np;
        }
        else{
            cout << "~~~~~Destructor ERROR!!!" << endl;
        }
	}


public:
	~AVL() {
		free_tree(root);
	}
	int get_height(){
		return max_height(root);
	}
	void create(){
        root = nullptr;
    }

	void insert(int v){
		if(root) insert_private(v, root);
		else {
            root = new Node(v, 0);
            // cout << "+ ROOT " << v << endl;
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

void test_1_2( int size, bool shuf=true ){
	cout << "\n====================================\nInserting " << size << " keys" << endl;
	AVL tree = AVL();

    int Tsize = size * 1.3, count = 0, c = 0;
    int data[size] = {0}, pool[Tsize] = {0};
    unordered_set<int> record(Tsize);

    pool[0] = 1;
    for( int i=1; i<Tsize; ++i) pool[i] = rand()%20 + pool[i-1];
	if( shuf ) shuffle(&pool[0], &pool[Tsize-1], default_random_engine(SEED));
    while( count < size ){
        if( record.find(pool[c]) == record.end() ){
            record.insert(pool[c]);
            data[count++] = pool[c++];
        }
        else c++;
    }
    // cout << "\nRandom number collisions: " << c-count << endl;
    // cout << "Data population for tree: " << count << endl;
    // cout << "Input: "; for(int i=0; i<size; ++i) cout << data[i] << " ";
    // cout << endl << endl;

	cout << "\n  Create Tree: ";
	auto start = std::chrono::steady_clock::now();
	tree.create();
	auto end = std::chrono::steady_clock::now();
	double t = double(std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count()) / 1e3;
	cout << (int)t << endl;


    cout << "  Insertion: ";
	start = std::chrono::steady_clock::now();
	for(int i=0; i<size; ++i){tree.insert(data[i]);}
	end = std::chrono::steady_clock::now();
	double t1 = double(std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count()) / 1e3;
	cout << (int)(t + t1) << endl;
	// tree.get_population();


	cout << "  Height: " << tree.get_height();


    cout << "\n  Search: ";
	start = chrono::steady_clock::now();
	for(int i=0; i<size; ++i){tree.search(data[i]);}
	end = chrono::steady_clock::now();
	t = double(chrono::duration_cast<chrono::nanoseconds>(end-start).count()) / 1e3;
	cout << (int)t << endl;



    cout << "  Deletion: ";
	start = chrono::steady_clock::now();
	for(int i=0; i<size; ++i){tree.remove(data[i]);}
	end = chrono::steady_clock::now();
	t = double(chrono::duration_cast<chrono::nanoseconds>(end-start).count()) / 1e3;
	cout << (int)t << endl << endl;

    cout << "~~~~ recycle tree (destructor) ~~~~" << endl;
}

int main(){
    vector<int> sizes = {10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000};
    for( int i=0; i<10; ++i)
		test_1_2(sizes[i]);
	test_1_2(10000, false);
	test_1_2(30000, false);
	test_1_2(50000, false);
}
