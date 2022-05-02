
#ifndef DATA_STRUCURE_DATA_STRUCTURES_H
#define DATA_STRUCURE_DATA_STRUCTURES_H

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <valarray>

#define break_line     std::cout<<"------------------------------------------------------------------------------------------------------------"<<'\n';

#define node_ptr tree_node*
class NotFound : std::exception {
};

class ItemDoseExist : std::exception {
public:
    void * ptr;
    explicit ItemDoseExist(void * ptr): ptr(ptr){}
};

typedef unsigned char byte;

class EmptyTree : std::exception {
};
class ArrayIsFilled:std::exception{};

class ArrayIsEmpty:std::exception{};

class ListIsEmpty:std::exception{};


namespace __BUFFER {

    typedef struct BUFFER {
        int len;
        byte *data;
    } BUFFER;

    typedef struct BUFFER_ARRAY : BUFFER {

        size_t size;
        int array_length;
        // data are not charcters but they will be treated as chars as a char = 1 byte and
        // copying functionality is simple and copy data bit by bit

    } BUFFER_ARRAY;
// used in short readerStack as buffer based readerStack:
    void _push_head(BUFFER &buf, void *item, size_t size, int array_length);

    byte *_pop_head(BUFFER &buf, size_t size);

    void _push_tail(BUFFER &buf, byte *item, size_t size, int array_length);
// used in short queue and double queue as buffer based queue:
    byte *_pop_tail(BUFFER &buf, size_t size, int array_length);

    byte *_read_head(BUFFER &buf, size_t size);

//-----------------------------------
    void _push_head(BUFFER *buf, void *item, size_t size, int array_length);

    byte *_pop_head(BUFFER *buf, size_t size);

    void _push_tail(BUFFER *buf, void *item, size_t size, int array_length);

    byte *_pop_tail(BUFFER *buf, size_t size, int array_length);

    byte *_read_head(BUFFER *buf, size_t size);

}

template<typename TYPE>
class Array {
    __BUFFER::BUFFER_ARRAY buf;
public:

    explicit Array(short array_length)
    {
        buf.size = sizeof(TYPE);
        buf.array_length = array_length;
        buf.len=0;
        buf.data = (byte *) (new TYPE[array_length]);
    }

    ~Array() {
        delete[] buf.data;
    }

    void push_head(TYPE item) {
        __BUFFER::_push_head(buf, (byte *) (&item), buf.size, buf.array_length);
//        std::cout << *(T*)(char *)(&item)<<'\t__'<< (item) << std::endl;
    }

    TYPE pop_head() {
        return *(TYPE *) (__BUFFER::_pop_head(buf, buf.size));
    }

    void push_tail(TYPE item) {
        __BUFFER::_push_tail(buf, (byte *) (&item), buf.size, buf.array_length);
    }

    TYPE pop_tail() {
        return *(TYPE *) (__BUFFER::_pop_tail(buf, buf.size, buf.array_length));
    }

};


namespace SingleLinkedList {

    typedef struct Node:public __BUFFER::BUFFER {
        Node *next;
    } Node;


    void b_sort( Node* node, byte* start, short n, size_t size, bool (*hash_function)(void *p1, void *p2) , byte *  next( Node*& node, byte * reader, size_t size));

    void q_sort(Node * node_i, byte * start , int length, size_t size, bool (*hash_function)(void *p1, void *p2), byte *  next( Node*& node, byte * reader, size_t size));

    template<typename T>
    class stackIterator;

    class _stack {
    public:
        short array_length;
        int len;
        Node *head;

    public:

        int _length() const{
            return len;
        }

        void _start_(Node*& node, byte*& reader, size_t size);

        static byte * _next_(Node*& node, byte* reader, size_t size);

         _stack();

        ~_stack();

        void _push(byte *data, size_t size);

        byte * _pop( size_t size);

        byte * _read_head__(size_t size);

        void _sort(size_t size, bool (*hash_function)(void *, void *));

        friend _stack* merge_in_place(_stack* x, _stack* y, size_t size, bool (*hash_function)(void *p1, void *p2));
    };

    _stack* merge_in_place(_stack* x, _stack* y, size_t size, bool (*hash_function)(void *p1, void *p2));

    class _queue {
        //  tail (pop) .. -> .. data .. -> .. head (set_item)
    protected:
        Node *head;
        Node * tail;

        int len;
    public:

        void _start_(Node*& node, byte*& reader, size_t size);

        static byte* _next_(Node*& node, byte* reader, size_t size);

        _queue();

        ~_queue();

        void _push(byte *data, size_t size,short array_length);

        byte * _pop( size_t size,short array_length) ;

        byte * _readTop(size_t size);

        int _length();

        void _sort(size_t size, bool (*hash_function)(void *p1, void *p2));

    };

    _queue* merge_in_place(_queue* x, _queue* y, size_t size, short array_length, bool (*hash_function)(void *p1, void *p2));


    template<typename T>
    class stack: protected _stack {
        size_t size;
    public:
        stack() : _stack(), size(sizeof(T))
        {
        }

        void push(T item) {
            _push((byte *) (&item),size);

        }

        T &pop() {
            return *(T *) _pop(size);
        }

        T top() {
            return *(T *) _read_head__(size);
        }

        void sort(){
            _sort(size,[](void *p1, void *p2) {
                if ((*(T *) p1) < (*(T *) p2)) {
                    return true;
                }
                return false;
            });
        }
        int length(){
            return len;
        }
        friend stackIterator<T>;
    };

    template<typename T>
    class stackIterator{
        byte* reader;
        Node * node;
        size_t size;
        _stack* s;
    public:
        stackIterator(){
            s = nullptr;
            reader = nullptr;
            node= nullptr;
            size = sizeof(T);
        }

        explicit stackIterator(stack<T>* stack){
            s = stack;
            reader = s->head->data;
            node=s->head;
            size = sizeof(T);
        }
        void start(stack<T>* stack){
            s = stack;
            reader = s->head->data;
            node=s->head;
            size = sizeof(T);
            s->_start_(node,reader,size);

        }
        void start(stack<T>& stack){
            s = &stack;
            if(s->head== nullptr){
                reader= nullptr;
                node= nullptr;
                return;
            }
            reader = s->head->data;
            node=s->head;
            size = sizeof(T);
            s->_start_(node,reader,size);

        }
        void operator++(){
            reader = stack<T>::_next_(node,reader,size);
        }
        T& read(){
            return *(T *) reader;
        }
        bool end(){
            if(reader== nullptr or node == nullptr){
                return false;
            }
            return true;
        }

    };

    template<typename T>
    class queueIterator;

    template<typename TYPE>
    class queue: protected _queue {

        size_t size;
        short array_length;
    public:
        queue() : _queue(), size(sizeof(TYPE)), array_length(4)
        {
        }

        void push(TYPE item) {
            _push((byte *) (&item),size,array_length);

        }

        TYPE &pop() {
            return *(TYPE *) _pop(size,array_length);
        }

        TYPE top() {
            return *(TYPE *) _readTop(size);
        }
/*
        void operator++(){ reader = _next_(node, reader,size);}

        void start(){_start_(node,reader,size);}

        bool end(){
            if(reader== nullptr){
                return false;
            }
            return true;
        }

        TYPE& read(){
            return *(TYPE *) reader;
        }
        */

        void sort(){
            _sort(size,[](void *p1, void *p2) {
                if (*(TYPE *) p1 > *(TYPE *) p2) {
                    return true;
                }
                return false;
            });
        }
        friend class queueIterator<TYPE>;
    };


    template<typename T>
    class queueIterator{
        byte* reader;
        Node * node;
        size_t size;
        queue<T>* s;
    public:
        queueIterator(){
            s = nullptr;
            reader = nullptr;
            node= nullptr;
            size = sizeof(T);
        }

        explicit queueIterator(queue<T>* stack){
            s = stack;
            reader = s->head->data;
            node=s->head;
            size = sizeof(T);
        }
        void start(queue<T>* stack){
            s = stack;
            reader = s->head->data;
            node=s->head;
            size = sizeof(T);
            s->_start_(node,reader,size);

        }
        void start(queue<T>& stack){
            s = &stack;
            if(s->head== nullptr){
                reader= nullptr;
                node= nullptr;
                return;
            }
            reader = s->head->data;
            node=s->head;
            size = sizeof(T);
            s->_start_(node,reader,size);

        }
        void operator++(){
            reader = _queue::_next_(node,reader,size);
        }
        T& read(){
            return *(T *) reader;
        }
        bool end(){
            if(reader== nullptr or node == nullptr){
                return false;
            }
            return true;
        }
    };

    typedef struct RNode {
        RNode *next;
        void* data;
    } RNode;

    class _ring{

        size_t size;
    protected:
        RNode *tail;

    public:
        explicit _ring(size_t size);

        ~_ring();

        void _push(byte *data);

        void * _pop_front();

        void * _read();

    };

    template<typename TYPE>
    class ring: protected _ring{
    public:
        ring(): _ring(sizeof(TYPE))
        {
        }

        void * get_head(){
            return tail;
        }

        void push(TYPE item){       //will set_item at the beginning
            _push((byte *) (&item));
        }

        void push(TYPE* item){       //will set_item at the beginning
            _push((byte *) (item));
        }

        TYPE& pop(){                 // will pop the last element pushed
            byte* data_ptr = _pop_front();
            TYPE data = *(TYPE*) data_ptr;
            delete[] data_ptr;
            return data;
        }

        TYPE* read(){
            return (TYPE *)_read();
        }
    };

    int delete_node(Node *head, Node *tail );

}

namespace DoubleLinkedList{

    typedef struct Node:public __BUFFER::BUFFER {
        Node *next;
        Node *back;

    } Node;

    class _list {
    protected:
        DoubleLinkedList::Node *head;
        DoubleLinkedList::Node *tail;

    public:

        _list();

        _list(size_t size, short array_length);

        ~_list();

        Node * get_head();

        Node * get_tail();

        void _push(byte *data,size_t size, short array_length);

        byte * _pop(size_t size, short array_length);


    };

    template<typename TYPE>
    class list: protected _list{
        size_t size;
        short array_length;

    public:

        explicit list( short array_length):_list(), array_length(array_length)
        {
        }

        void push(TYPE item){
            _push((byte *) (&item),size, array_length);
        }

        void push(TYPE *item){
            _push((byte *) (item), size, array_length);
        }

        void push(void *item){
            _push((byte *) (item),size, array_length);
        }

        TYPE pop(){
            return *(TYPE *)_pop(size,array_length);
        }


    };


    class _ring{
    private:
        size_t size;

    protected:
        Node *head;

    public:
        explicit _ring(size_t size);

        ~_ring();

        void _push_front(byte *data);

        byte * _pop_front();

        byte * _go_front();

        byte * _go_back();

        byte * _pop_back();

    };
    template<typename TYPE>
    class ring: protected _ring{
    public:

        ring(): _ring(sizeof(TYPE))
        {
        }

        void push_front(TYPE item){       //will set_item at the beginning
            _push_front((byte *) (&item));
        }

//        void push_last(T item){       //will set_item at the beginning
//            _push_((byte *) (&item));
//        }

        TYPE pop_front(){                 // will pop the last element pushed
            byte* data_ptr = _pop_front();
            TYPE data = *(TYPE*) data_ptr;
            delete[] data_ptr;
            return data;
        }

        TYPE pop_back(){                 // will pop the last element pushed
            byte* data_ptr = _pop_back();
            TYPE data = *(TYPE*) data_ptr;
            delete[] data_ptr;
            return data;
        }

        TYPE forward(){
            return *(TYPE *)_go_front();
        }

        TYPE backward(){
            return *(TYPE *)_go_back();
        }
    };


}


enum hash_value {
    sma=0, gr, eq
};


// augmanted data structure
// void **p 0: data     1: this     2: other data

namespace ABST{
    using namespace __BUFFER;

    typedef struct tree_leaf:BUFFER {
        bool leaf;
    }tree_leaf;

    typedef struct tree_node:tree_leaf {
        tree_node* left;
        tree_node* right;
        int height;
    }tree_node;

    typedef SingleLinkedList::stack<node_ptr> stack;

    class tree{
        size_t size;
        short array_length;

        bool (*condition)(void *);

        tree_node * insert_func(tree_node *parent, void *data);

    protected:
        tree_node *tree_root;
        int len;

        hash_value(*hash_function)(void *, void *);  // return true if the first item is greeter than the second argument

        tree_node * read_node(tree_node* parent, void *data, void**p, node_ptr(*func_found)(node_ptr, void **),node_ptr(*func_not_found)(node_ptr, void **));

        void * read_item(tree_node* parent, void *data, void**p, void *(*func_found)( void **),void *(*func_not_found)( void **));

        void descending_read(void ** p,node_ptr(*func)( void **));

        byte *reader;
        byte * last;
        stack* iterator;

    public:

        tree(hash_value(*hash_function)(void *, void *), size_t size, short array_length);

        ~tree();

        byte * read(){return reader;}

        void ascending_read(void ** p,void (*func)( void **));

        void conditional_read(void ** p,void (*func)( void **), bool (*condition)(void *));

        void start();

        void operator++();

        void *_find(void *data);

        void _print();

        void _insert(void *data);

        void _remove(void * data);

        node_ptr copy(node_ptr parent);

        bool end;

// friend functions ............
        friend      void read_buf(BUFFER *buffer, void **p, void (*func)(void **));
        friend      void des_read_buf(BUFFER *buffer, void **p, node_ptr(*func)(void **));
        friend      tree_node *removeFound(node_ptr parent, void **p);
        friend      void c_read_buf(BUFFER *buffer, void **p, void (*func)(void **));
        friend      tree_node *remove_rec(tree_node *parent, void **p);
        };


    template<typename T>
    class Tree: protected tree{

    public:
        Tree(hash_value(*hash_function)(void *, void *), short array_length = 5):
                tree(hash_function, sizeof(T),array_length)
        {

        }

        void insert(T& data){
            _insert(&data);
        }

        void insert(T&& data){
            T d = data;
            _insert(&d);
        }


        void remove(T& data){
            _remove(&data);
        }

        T& find(T& data){
            return _find(&data);
        }

        void print(){
            _print();
        }

    };


}



#endif //DATA_STRUCURE_DATA_STRUCTURES_H
