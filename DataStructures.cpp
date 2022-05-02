//
// Created by m7md_nor on 9/4/2021.
//

#include "../../header/DataStructures.h"



namespace __BUFFER {
// pass by refrence....
    void _push_head(BUFFER &buf, void *item,size_t size, int array_length ) {
        if (buf.len == array_length) {
            throw ArrayIsFilled();
        } else {
//            copy_item(buf, item, buf.data + buf.len * buf.size);
            memcpy(buf.data + buf.len * size,item,size);
            buf.len++;
        }
    }

    byte *_pop_head(BUFFER &buf,size_t size ) {
        if (buf.len > 0) {
            buf.len--;
            return buf.data + buf.len * size;
        } else {
            throw ArrayIsEmpty();
        }
    }

    byte *_read_head(BUFFER &buf,size_t size) {

        if (buf.len > 0) {
            return buf.data + (buf.len-1) * size;
        } else {
            throw ArrayIsEmpty();
        }
    }


    void _push_tail(BUFFER &buf, byte *item,size_t size, int array_length ) {
        if (buf.len == array_length) {
            throw ArrayIsFilled();
        } else {
            memcpy(buf.data+size,buf.data, size*buf.len);
            memcpy(buf.data, item, size);
            buf.len++;
        }

    }

    byte *_pop_tail(BUFFER &buf,size_t size, int array_length ) {
        if (buf.len == 0) {
            throw ArrayIsEmpty();
        } else {
            buf.len--;

            byte item[size];
            memcpy(item, buf.data, size);
            memcpy(buf.data, buf.data+size, size*buf.len);
            memcpy(buf.data + (array_length - 1) * size, item,size);

            return buf.data + (array_length - 1) * size;

        }

    }

// pass by pointer....
    void _push_head(BUFFER *buf, void *item,size_t size, int array_length ) {
        if (buf->len > array_length-1) {
            throw ArrayIsFilled();
        } else {
//            copy_item(buf, item, buf.data + buf.len * buf.size);
            memcpy(buf->data + buf->len * size,item,size);
            buf->len++;
        }
    }

    byte *_pop_head(BUFFER *buf,size_t size) {
        if (buf->len > 0) {
            buf->len--;
            return buf->data + buf->len * size;
        } else {
            throw ArrayIsEmpty();
        }
    }

    byte *_read_head(BUFFER *buf,size_t size) {
        if(buf == 0){throw ListIsEmpty();}
        if (buf->len > 0) {
            return buf->data + (buf->len-1) * size;
        } else {
            throw ArrayIsEmpty();
        }
    }


    void _push_tail(BUFFER *buf, byte *item,size_t size, int array_length ) {
        if (buf->len == array_length) {
            throw ArrayIsFilled();
        } else {
            memcpy(buf->data+size,buf->data, size*buf->len);
            memcpy(buf->data, item, size);
            buf->len++;
        }

    }

    byte *_pop_tail(BUFFER *buf,size_t size, int array_length ) {
        if (buf->len == 0) {
            throw ArrayIsEmpty();
        } else {
            buf->len--;

            byte item[size];
            memcpy(item, buf->data, size);
            memcpy(buf->data, buf->data+size, size*buf->len);
            memcpy(buf->data + (array_length - 1) * size, item,size);

            return buf->data + (array_length - 1) * size;

        }
    }


    void *find(BUFFER *buf, void *data, int start, int end, size_t size, hash_value(*hash_function)(void *, void *)) {
        if (end == start) {
            return buf->data + ((start + end) / 2) * size;
        }
        short v = hash_function(data, buf->data + ((start + end) / 2) * size);
//        std::cout<<*(int *)(buf.data +start*buf.size)<<'\t__'<<*(int *)(buf.data +((start + end)/ 2)*buf.size)<<'\t__'<<*(int *)(buf.data +end*buf.size)<<' '<<v<<'\span';
        if (v == gr) {
            return find(buf, data, (start + end) / 2 + 1, end, size, hash_function);
        } else {
            return find(buf, data, start, (start + end) / 2, size, hash_function);
        }
    }

    void ordered_insert(BUFFER *buf, void *data, size_t size, short array_length,
                        hash_value(*hash_function)(void *, void *)) {
        void *place = find(buf, data, 0, buf->len, size, hash_function);
//        printf("%x__  %x__", place, data);
        if (hash_function(data, place) == eq and buf->len != 0) {
            throw ItemDoseExist(place);
        } else {
            if (buf->len == array_length) {
                throw ArrayIsFilled();
            }
            buf->len++;
            memcpy((byte *) place + size, place,
                   (buf->data + (buf->len - 1) * size) - (byte *) place);
            memcpy(place, data, size);
        }
    }

    void ordered_delete(BUFFER *buf, void *data, size_t size, hash_value(*hash_function)(void *, void *)) {
        if (buf->len == 0) {
            throw ArrayIsEmpty();
        }
        void *place = find(buf, data, 0, buf->len, size, hash_function);
//        printf("%x__  %x__", place, data);
        if (place == buf->data + buf->len * size or hash_function(data, place) != eq) {
            throw NotFound();
        } else {
            buf->len--;
            memcpy(place, (byte *) place + size, (buf->data + buf->len * size) - (byte *) place);
        }
    }

}
//--------------------------------------------------------------------------

namespace SingleLinkedList {

    _stack* merge_in_place(_stack* x, _stack* y, size_t size, bool (*hash_function)(void *p1, void *p2))
    {
        _stack* result = new _stack;
        while (true){
            try{
               std::cout<<*(int*)x->_read_head__(size)<<'\t'<<*(int*)y->_read_head__(size)<<'\n';

                if(!hash_function(x->_read_head__(size), y->_read_head__(size))){
                    result->_push(x->_pop(size),size);
                }else{
                    result->_push(y->_pop(size),size);
                }
            } catch (ListIsEmpty) {break;}
        }
        while (true){
            try{result->_push(x->_pop(size),size);} catch (ListIsEmpty) {break;}
        }
        while (true){
            try{result->_push(y->_pop(size),size);} catch (ListIsEmpty) {break;}
        }

        return result;
    }

    _queue* merge_in_place(_queue* x, _queue* y, size_t size, short array_length, bool (*hash_function)(void *p1, void *p2))
    {
        _queue* result = new _queue();
        while (true){
            try{
                std::cout<<*(int*)x->_readTop(size)<<'\t'<<*(int*)y->_readTop(size)<<'\n';

                if(!hash_function(x->_readTop(size), y->_readTop(size))){
                    result->_push(x->_pop(size,array_length),size,array_length);
                }else{
                    result->_push(y->_pop(size,array_length),size,array_length);
                }
            } catch (ListIsEmpty) {break;}
        }
        while (true){
            try{result->_push(x->_pop(size,array_length),size,array_length);} catch (ListIsEmpty) {break;}
        }
        while (true){
            try{result->_push(y->_pop(size,array_length),size,array_length);} catch (ListIsEmpty) {break;}
        }

        return result;
    }


    _stack::_stack() : head(nullptr), array_length(1),len(0) {}

    void _stack::_push(byte *data, size_t size) {
        len++;
        if (head == nullptr) {
            head = new Node;
            head->data = (byte*) malloc ( array_length * size);
            head->len=0;
            head->next=0;
            __BUFFER::_push_head(head, data, size, array_length);
        } else {
            try {
                __BUFFER::_push_head(head, data, size, array_length);
            } catch (ArrayIsFilled&) {

                array_length = short(array_length +1);
                Node * nh= new Node;
                nh->data = (byte*) malloc ( array_length * size);
                nh->next=head;
                nh->len = 0;
                head = nh;
                __BUFFER::_push_head(head, data, size, array_length);
            }
        }
    }

    byte *_stack:: _pop( size_t size) {
        if (head == nullptr) {
            throw ListIsEmpty();
        } else {
            try {
                --len;
                return __BUFFER::_pop_head(head, size);
            } catch (ArrayIsEmpty&) {
//                std::cout<<"head :"<<head<<'\span';
                Node *_next = head;
                head = head->next;
                _next->next = nullptr;
                delete[] _next->data;
                delete _next;
                ++len;
                return _pop(size);
            }
        }
    }

    void _stack::_start_(Node*& node, byte*& reader, size_t size){
        node=head;
        reader = head->data + (head->len-1)*size;
    }

    byte* _stack::_next_(Node*& node, byte* reader, size_t size)
    {
        if(node==0){
            return 0;
        }
        if (reader > node->data) {
            return reader-size;
        } else {
            node = node->next;
            if (node == nullptr) {
                reader= nullptr;
                return reader;
            }
            return node->data + (node->len-1) * size;
        }
    }


    byte * _stack::_read_head__( size_t size) {
        if (head == nullptr) {return nullptr;}
        try {
            return _read_head(head, size);
        } catch (ArrayIsEmpty&) {
            return _read_head(head->next, size);
        }
    }

    int delete_node(Node *head, Node *tail = nullptr) {

        if (head == tail) {
            return 0;
        } else {
            delete[] head->data;
            delete_node(head->next, tail);
            delete head;
            return 0;
        }
    }


    void swap(void * first, void * second, size_t size){
        if(first==second){return;}
        byte item[size];

        memcpy(item, first, size);
        memcpy(first, second, size);
        memcpy(second, item, size);

    }

    void b_sort( Node* node, byte* start, int n, size_t size, bool (*hash_function)(void *p1, void *p2),    byte *  next( Node*& node, byte * reader, size_t size))
    {
        short i;
        bool not_sorted = true;
        byte *reader;
        void *back;
        Node * node1;

        while (not_sorted) {
            not_sorted = false;

            i = 0;
            reader = start;
            back = start;
            node1 = node;

            reader = next(node1, reader, size);

            while (i < n and reader!=0)
            {
                if (hash_function(reader, back)) {
                    swap(back, reader, size);
                    not_sorted = true;
                }
                back = reader;
                reader = next(node1, reader, size);
                i++;
            }
        }
    }

    void print(Node* node, byte* reader, int length, byte *  next( Node*& node, byte * reader, size_t size)){
        int i = 0;
        while (i<length and reader != nullptr){
            std::cout<<*(int*)reader<<'\n';
            reader = next(node, reader, 4);
            i++;
        }
        break_line
    }

    void q_sort(Node * node_i, byte * start , int length, size_t size, bool (*hash_function)(void *p1, void *p2),    byte *  next( Node*& node, byte * reader, size_t size)){
        if(length<4){
            b_sort(node_i, start, length, size, hash_function,next);
            return;
        }

        Node * start_node = node_i;
        Node * node_j = node_i;

        byte * reader_i = start;
        byte * reader_j = start;
        byte * reader_j1 = start;

        reader_j = next(node_j, reader_j, size);
        reader_i = next(node_i, reader_i, size);
        int i =0;
        int j = 0;

        while (i<length and reader_i != nullptr){
            if (hash_function(reader_i , start)) {
                swap(reader_i, reader_j, size);
                reader_j1 = reader_j;
                reader_j = next(node_j, reader_j, size);
                j++;
            }
            reader_i = next(node_i, reader_i, size);
            i++;
        }

        swap(reader_j1, start, size);

        q_sort(start_node, start, j-1, size, hash_function,next);
        q_sort(node_j, reader_j, i - j, size, hash_function,next);
    }

    int delete_node(RNode *head, RNode *tail = nullptr)
    {

        if (head == tail) {
            return 0;
        } else {
            free(head->data);
            delete_node(head->next, tail);
            delete head;
            return 0;
        }
    }

    _stack::~_stack() {
        delete_node(head);
    }

    void _stack::_sort(size_t size, bool (*hash_function)(void *, void *)){

        q_sort(head, _read_head__(size), len, size, hash_function, _stack::_next_);
    }
    //------------------------------------------------------------------------------------------
    // ** queue **

    void _queue::_start_(Node*& node, byte*& reader, size_t size)
    {
        node=tail;
        reader = tail->data ;
    }

    byte* _queue::_next_(Node*& node, byte* reader, size_t size)
    {
        if(node== nullptr){return nullptr;}
        if (node->data + (node->len-1) * size> reader) {
            return reader + size;
        } else {
            node = node->next;
            if (node == nullptr) {
                return nullptr;
            }
            return node->data ;
        }
    }

    byte * _queue::_readTop(size_t size){
        return tail->data;
    }

    _queue::_queue(): head(0), tail(0), len(0)
    {
    }

    _queue::~_queue(){
        delete_node(tail, nullptr);
    }

    void _queue::_push(byte *data, size_t size,short array_length) {
        ++len;
        if (head == nullptr) {
            head = new Node;
            head->len = 0;
            head->data = new byte[array_length * size];
            head->next=0;
            tail = head;
            __BUFFER::_push_head(head, data, size, array_length);
        } else {
            try {
                __BUFFER::_push_head(head, data, size, array_length);
            } catch (ArrayIsFilled &) {
//                std::cout << "new buffer is instantiated of length = " << array_length << std::endl;
                Node *nh = new Node;
                nh->len = 0;
                nh->data = new byte[array_length * size];
                nh->next = nullptr;
                head->next = nh;
                head = nh;
                __BUFFER::_push_head(head, data, size, array_length);
            }
        }
    }

    byte * _queue::_pop( size_t size,short array_length) {
        --len;
        if (head == nullptr) {
            throw ListIsEmpty();
        } else if (head != tail) {
            if (tail->len != 0) {
                tail->len--;
                tail->data += size;
                return tail->data - size;

            } else {

                auto  a = tail->next;
                delete[] (tail->data - array_length * size);
                delete tail;

                tail = a;
                ++len;
                return _pop(size, array_length);
            }
        } else {
            try {
                return __BUFFER::_pop_tail(tail, size, array_length);
            } catch (ArrayIsEmpty &) {
                throw ListIsEmpty();
            }
        }
    }

    int _queue::_length(){return len;}

    void _queue::_sort(size_t size, bool (*hash_function)(void *p1, void *p2)){
        q_sort(tail,tail->data,len,size,hash_function,_queue::_next_);
    }

    //------------------------------------------------------------------------------------------
    // ** ring **
    _ring::_ring(size_t size) : size(size), tail(nullptr) {
    }

    _ring::~_ring() {
        if(tail==0){
            delete tail;
            return;
        }
        delete_node(tail->next, tail);
        free(tail->data);
        delete tail;
    }

    void _ring::_push(byte *data) {
        if (tail == nullptr) {
            tail = new RNode{nullptr, malloc(size)};
            tail->next = tail;
            memcpy(tail->data, data, size);
        } else {
            RNode *node = new RNode{tail->next, malloc(size)};
            tail->next = node;
            memcpy(tail->next->data, data, size);
        }
    }


    void *_ring::_read() {
        if (tail == nullptr) {
            throw ListIsEmpty();
        } else {
            tail = tail->next;
            return tail->data;
        }
    }

    void *_ring::_pop_front() {
        if (tail == nullptr) {
            throw ListIsEmpty();
        } else {
            void *item;
            RNode *index;
            index = tail->next;
            item = tail->next->data;
            tail->next = tail->next->next;
            delete index;

            return item;
        }
    }
}
//--------------------------------------------------------------------------

namespace DoubleLinkedList {
//double linked list used to implement queue


    _list::_list() : head(nullptr), tail(nullptr) {
    }

    _list::_list(size_t size, short array_length) : head(nullptr), tail(nullptr){
    }
    Node * _list::get_head(){
        return head;
    }

    Node * _list::get_tail(){
        return tail;
    }
    void _list::_push(byte *data, size_t size, short array_length) {
        if (head == nullptr) {
            head = new DoubleLinkedList::Node ;
            head->len =0;
            head->data = new byte[array_length * size];
            head->next =0;
            head->back = 0;
            tail = head;

            __BUFFER::_push_head(head, data, size, array_length);
        } else {
            try {
                __BUFFER::_push_head(head, data, size, array_length);
            } catch (ArrayIsFilled&) {
//                std::cout << "new buffer is instantiated of length = " << array_length << std::endl;
//                nh = new DoubleLinkedList::Node{0, new byte[array_length * size],head, nullptr};

                Node *nh = new DoubleLinkedList::Node ;
                nh->len =0;
                nh->data = new byte[array_length * size];
                nh->next =head;
                nh->back = 0;

                head->back = nh;
                head = nh;
                __BUFFER::_push_head(head, data, size, array_length);
            }
        }
    }

    byte *_list::_pop(size_t size, short array_length) {
        if (head == nullptr) {
            throw ListIsEmpty();
        } else if (head != tail) {
            if (tail->len != 0) {
                tail->len--;
                tail->data += size;
                return tail->data-size;

            } else {

                tail = tail->back;
                delete[] (tail->next->data-array_length*size);
                delete tail->next;

                tail->next= nullptr;

                return _pop(size, array_length);
            }
        } else {
            try {
                return __BUFFER::_pop_tail(tail, size, array_length);
            } catch (ArrayIsEmpty&) {
                throw ListIsEmpty();
            }
        }
    }

    int delete_node(Node *head, Node *tail = nullptr) {

        if (head == tail) {
            return 0;
        } else {
            delete[] (head->data);
            delete_node(head->next, tail);
            delete head;
            return 0;
        }
    }

    _list::~_list() {
        delete_node(head);
    }
//------------------------------------------------------------------------------------
// double linked ring

    _ring::_ring(size_t size) : size(size), head(nullptr) {
    }

    _ring::~_ring() {
        if (head == nullptr) {
            return;
        }
        delete_node(head->next, head);
        delete[] head->data;
        delete head;
    }

    void _ring::_push_front(byte *data) {
        if (head == nullptr) {
            head = new Node;
            head->len=0;
            head->data = new byte[size];

            head->next = head;
            head->back = head;
            __BUFFER::_push_head(head, data, size, 1);
        } else {
            Node *node = new Node;
            node->len=0;
            node->data = new byte[size];
            node->next=head->next;
            node->back=head;

            head->next = node;
            head->next->next->back = head->next;
            __BUFFER::_push_head(head->next, data, size, 1);
        }
    }


    byte *_ring::_go_front() {
        if (head == nullptr) {
            throw ListIsEmpty();
        } else {
            head = head->next;
            return head->data;
        }
    }

    byte *_ring::_go_back() {
        if (head == nullptr) {
            throw ListIsEmpty();
        } else {
            head = head->back;
            return head->data;
        }
    }


    byte *_ring::_pop_front() {
        if (head == nullptr) {
            throw ListIsEmpty();
        } else if (head == head->next) {
            byte *h = head->data;
            head = nullptr;
            delete head;

            return h;
        } else {
            byte *item;
            item = head->next->data;
            head->next = head->next->next;
            delete head->next->back;
            head->next->back = head;

            return item;
        }
    }

    byte *_ring::_pop_back() {
        if (head == nullptr) {
            throw ListIsEmpty();
        } else if (head == head->back) {
            byte *h = head->data;
            delete head;
            head = nullptr;
            return h;
        } else {
            byte *item;
            item = head->back->data;
            head->back = head->back->back;
            delete head->back->next;
            head->back->next = head;

            return item;
        }
    }

}
//--------------------------------------------------------------------------


std::string space(int n){
    std::string s = " ";
    for(int i=0; i<n; i++){
        s += ' ';
    }
    return s;
}


namespace ABST {
    using namespace __BUFFER;

    tree::tree(hash_value(*hash_function)(void *, void *), size_t size, short array_length) :
            size(size), tree_root(nullptr), len(0), hash_function(hash_function), array_length(array_length),
        condition(nullptr),
            iterator(new stack) {
    }


    node_ptr min_item(node_ptr parent) {
        // left part
        if (parent->left == nullptr) {
            return parent;
        } else {
            return min_item(parent->left);
        }
    }

    node_ptr max_item(node_ptr parent) {
        // right part
        if (parent->right == nullptr) {
            return parent;
        } else {
            return max_item(parent->right);
        }
    }

// assending reading.....................
    void read_buf(BUFFER *buffer, void **p, void (*func)(void **)) {
        byte *reader = buffer->data;
        size_t size = ((tree *) p[1])->size;
        while (reader < buffer->data + buffer->len * size) {
            p[0] = reader;
            func(p);
            reader += size;
        }
    }

    void asc_read(tree_node *parent, void **p, void (*func)(void **)) {
        if (parent == nullptr) { return; }
        asc_read(parent->left, p, func);
        read_buf(parent, p, func);
        asc_read(parent->right, p, func);

    }

    void tree::ascending_read(void **p, void (*func)(void **)) {
        asc_read(tree_root, p, func);
    }

// descending reading.......
    void des_read_buf(BUFFER *buffer, void **p, node_ptr(*func)(void **)) {
        size_t size = ((tree *) p[1])->size;
        byte *reader = buffer->data + (buffer->len - 1) * size;
        while (reader >= buffer->data) {
            p[1] = reader;
            func(p);
            reader -= size;
        }
    }

    void dsc_read(tree_node *parent, void **p, node_ptr(*func)(void **)) {
        if (parent == nullptr) { return; }

        dsc_read(parent->right, p, func);
        des_read_buf(parent, p, func);
        dsc_read(parent->left, p, func);


    }

    void tree::descending_read(void **p, node_ptr(*func)(void **)) {
        dsc_read(tree_root, p, func);
    }

    void *tree::read_item(tree_node *parent, void *data, void **p, void *(*func_found)(void **),void *(*func_not_found)(void **))
    {
        if (parent == nullptr) {
            p[1] = __BUFFER::find(parent, data, 0, len, size, hash_function);
            return func_not_found(p);
        }

        if(parent->leaf){
            p[1] = __BUFFER::find(parent, data, 0, len, size, hash_function);
            if(hash_function(p[1], data)==eq){
                return func_found(p);
            }
            return func_not_found(p);
        }

        if (hash_function(data, parent->data + (parent->len - 1) * size) == gr) {
            return read_item(parent->right, data, p, func_found, func_not_found);
        } else if (hash_function(data, parent->data) == sma) {
            return read_item(parent->left, data, p, func_found, func_not_found);
        } else {
            p[1] = __BUFFER::find(parent, data, 0, len, size, hash_function);
            return func_found(p);
        }
    }
    // p = {data outside, data in the node, other}
    int height(tree_node *parent) {
        if (parent == nullptr) { return -1; }
        if(parent->leaf){return -1;}
        return parent->height;
    }

    void set_height(tree_node *parent) {
        if (parent == nullptr) { return; }
        if(parent->leaf){return;}
        parent->height = 1 + std::max(height(parent->right), height(parent->left));
    }

    node_ptr right_r(node_ptr x) {
        node_ptr y = x->right;
        node_ptr b = y->left;
        y->left = x;
        x->right = b;
        set_height(b);
        set_height(x);

        return y;
    }

    node_ptr left_r(node_ptr y) {
        node_ptr x = y->left;
        node_ptr b = x->right;
        y->left = b;
        x->right = y;
        set_height(b);
        set_height(y);
        return x;
    }

    int balance(node_ptr parent) {
        if(parent== nullptr){return 0;}
        if(parent->leaf){return 0;}

        return height(parent->left) - height(parent->right);
    }

    node_ptr avl_function(node_ptr parent) {

        if (balance(parent) > 1) {
            if (balance(parent->left) < 0) {
                parent->left = right_r(parent->left);
                parent = left_r(parent);
                set_height(parent);
                return parent;
            }
            parent = left_r(parent);
            set_height(parent);
            return parent;
        } else if (balance(parent) < -1) {
            if (balance(parent->right) > 0) {
                parent->right = (left_r(parent->right));
                parent = right_r(parent);
                set_height(parent);
                return parent;
            }
            parent = right_r(parent);
            set_height(parent);
            return parent;
        } else {
            set_height(parent);
            return parent;
        }
    }


    tree_node *tree::read_node(tree_node *parent, void *data, void **p, node_ptr(*func_found)(node_ptr, void **),
                               node_ptr(*func_not_found)(node_ptr, void **)) {
        if (parent == nullptr) {
            return func_not_found(parent, p);
        }
        if(parent->leaf){
            p[0] = __BUFFER::find(parent, data, 0, parent->len, size, hash_function);
            if(hash_function(p[0], data)==eq){
                return func_found(parent,p);
            }
            return func_not_found(parent,p);
        }
        // nodes
        if (hash_function(data, parent->data + (parent->len-1) * size) == gr and parent->len == array_length) {
            parent->right = avl_function(read_node(parent->right, data, p, func_found, func_not_found));
            return parent;
        } else if (hash_function(data, parent->data) == sma and parent->len == array_length) {
            parent->left = avl_function(read_node(parent->left, data, p, func_found, func_not_found));
            return parent;
        } else {
            return func_found(parent, p);
        }

    }


    tree_node *tree::insert_func(tree_node *parent, void *data) {
        if (parent == nullptr) {
            auto node = new tree_leaf;

            node->len=0;
            node->data = (byte *) malloc(size * array_length);
            node->leaf= true;

            _push_head(node, data, size, array_length);
            return (tree_node*)node;
        }
        if(parent->leaf){
            try{
                ordered_insert(parent,data,size,array_length,hash_function);
                return parent;
            } catch (ArrayIsFilled&) {
                auto node = new tree_node;
//                {array_length, parent->data, false, nullptr, nullptr};

                node->len=0;
                node->data = (byte *) malloc(size * array_length);
                node->leaf= false;
                node->right=0;
                node->left=0;

                insert_func(node, data);
                delete parent;
                return node;
            }
        }
        // nodes
        if (hash_function(data, parent->data + (parent->len-1) * size) == gr ) {
            parent->right = avl_function(insert_func(parent->right, data));
            return parent;
        } else if (hash_function(data, parent->data) == sma ) {
            parent->left = avl_function(insert_func(parent->left, data));
            return parent;
        } else {
            try {
                ordered_insert(parent, data, size, array_length, hash_function);
                return parent;
            } catch (ArrayIsFilled &) {
                parent ->right = avl_function(insert_func(parent->right, _pop_head(parent,size)));
                ordered_insert(parent, data, size, array_length, hash_function);
                return parent;
            }
        }
    }

    void *tree::_find(void *data) {
        void *p[3] = {data,0, this};

        return read_item(tree_root, data, p, [](void **p)->void *{
            tree* tr = (tree*) p[2];
            if(tr->hash_function(p[0], p[1])==eq){return p[1];}
            throw NotFound();
            }, [](void **p)->void *{ throw NotFound();});
    }


    void tree::_insert(void *data){
        tree_root = avl_function(insert_func(tree_root,data));
    }

    tree_node * remove_rec(tree_node *parent, void **p) {

        if (parent->leaf) {
            tree *tr = (tree *) p[1];
            memcpy(p[0], parent->data, tr->size);
            _pop_tail(parent, tr->size, tr->array_length);

            if (parent->len == 0) {
                delete parent->data;
                delete parent;
                return nullptr;
            }
            return parent;
        }
        if (parent->left == 0) {
            if (parent->right == 0) {
                parent->leaf = true;
                return remove_rec(parent, p);
            }
            tree *tr = (tree *) p[1];
            memcpy(p[0], parent->data, tr->size);
            _pop_tail(parent, tr->size, tr->array_length);
            p[0] = _read_head(parent, tr->size);
            parent->right = remove_rec(parent->right, p);
            return parent;
        }
        parent->left = remove_rec(parent->left, p);
        return parent;
    }

    tree_node * removeFound(node_ptr parent, void **p)
    {
        tree* tr=((tree *) p[1]);
        if (parent->leaf)
        {
            ordered_delete(parent, p[0], tr->size, tr->hash_function);
            if(parent->len==0){
                delete[] parent->data;
                delete parent;
                return nullptr;
            }
            return parent;
        } else {
            ordered_delete(parent, p[0],tr->size, tr->hash_function);
            parent->len++;
            p[0] = parent->data + (parent->len-1)*tr->size;
            parent->right=remove_rec(parent->right, p);
            return parent;
        }
    }



    void tree::_remove(void * data)
    {
        void * p[2] = {data, this};
        len--;
        tree_root = avl_function(read_node(tree_root, data, p,removeFound, [](node_ptr, void ** )->node_ptr{throw NotFound();}));
    }



    void delete_rec(node_ptr node){
        if(node== nullptr){return ;}
        if(node->leaf)
        {
            delete[] node->data;
            delete node;
            return ;
        }

        delete_rec(node->right);
        delete_rec(node->left);
        delete[] node->data;
        delete node;
    }


    tree::~tree() {
        delete_rec(tree_root);
    }

    void iterator_start(tree_node* parent, stack* iterator)
    {
        if(parent== nullptr){return;}
        if(parent->leaf){
            iterator->push(parent);
            return;
        }
        iterator->push(parent);
        iterator_start(parent->left,iterator);
    }

    void tree::start()
    {
        end = true;
        iterator_start(tree_root, iterator);
        reader = iterator->top()->data;
        last = iterator->top()->data;
        last += (iterator->top()->len-1)*size;
    }

    void tree::operator++() {
        if (last > reader) {
            reader += size;
            return;
        }
        if (iterator->length() > 0) {
            if (iterator->top()->leaf) {
                iterator->pop();
            } else {
                iterator_start((iterator->pop())->right, iterator);
            }
            reader = iterator->top()->data;
            last = reader;
            last += (iterator->top()->len - 1) * size;

            return;
        }
        end = false;
    }

    void _print_(tree_node *parent, void * tr, int n=0) {
        if (parent == nullptr) { return; }
        void *p[2] = {nullptr, tr};
        n++;

        if(parent->leaf){
            std::cout << space(10*(n-1));
            std::cout <<"h:-1 ";
            read_buf(parent, p, [](void **p) { std::cout << *(int *) p[0] << ' '; });
            std::cout << '\n';
            return;
        }

        _print_(parent->left, p[1], n);
        std::cout << space(10*(n-1));
        std::cout << "height:" << parent->height << ' '<<' ';
        read_buf(parent, p, [](void **p) { std::cout << *(int *) p[0] << ' '; });
        std::cout << '\n';
        _print_(parent->right, p[1], n);

    }

    void tree::_print()
    {
        _print_(tree_root, this);
    }

    node_ptr tree::copy(node_ptr parent)
    {
        if(parent == nullptr){return nullptr;}
        if(parent->leaf){
            auto node = new tree_leaf;
//            {parent->len,( byte*) malloc(array_length*size), true};

            node->len = parent->len;
            node->data = ( byte*) malloc(array_length*size);
            node->leaf = true;

            memcpy(node->data, parent->data, parent->len*size);

            return (tree_node*)node;
        }

        auto node = (tree_node *) malloc(sizeof(tree_node));

        node->height = parent->height;
        node->len = parent->len;

        node->data = new byte[array_length*size];
        memcpy(node->data, parent->data, parent->len*size);

        node->left = copy(parent->left);
        node->right = copy(parent->right);

        return node;
    }

    void c_read_buf(BUFFER *buffer, void **p, void (*func)(void **))
    {
        byte *reader = buffer->data;
        tree* tr = ((tree *) p[1]);
        while (reader < buffer->data + buffer->len * tr->size) {
            p[0] = reader;
            if(tr->condition(p[0])) {
                func(p);
            }
            reader += tr->size;
        }
    }

    void c_read(tree_node *parent, void **p, void (*func)(void **)) {
        if (parent == nullptr) { return; }
        if(parent->leaf){
            c_read_buf(parent, p, func);
            return;
        }
        asc_read(parent->left, p, func);
        c_read_buf(parent, p, func);
        asc_read(parent->right, p, func);

    }

    void tree::conditional_read(void ** p,void (*func)( void **), bool (*_condition)(void *)){
        tree::condition = _condition;
        c_read(tree_root, p,func);
    }

}
