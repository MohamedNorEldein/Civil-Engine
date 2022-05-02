//
// Created by m7md_nor on 9/21/2021.
//
/* missing  1- clearing memory function based on readerStack
            2- inverse matrix sorting algorithm before inverting*/


#pragma once
#ifndef DATA_STRUCTURES_H_LINEARALGEBRA_H
#define DATA_STRUCTURES_H_LINEARALGEBRA_H

#include "DataStructures.h"
#include <math.h>
#undef node_ptr
#define node_ptr ABST::tree_node *

typedef double TYPE ;

class DimensionError : std::exception {
};

class RedundantSystem : public std::exception {
};

class UnsolvableSystem : public std::exception {
};


namespace short_vector {

    class Vector;

    class Matrix;

    static SingleLinkedList::stack<Vector *> delete_list;

    static SingleLinkedList::stack<Matrix *> mat_delete_list;

    void clear_memory();

    class Vector {
    private:
        TYPE *vector;
        int dimension;

    public:

        explicit Vector(int dimension);

        Vector(int dimension, TYPE vector_[]);

        Vector(const Vector &other);

        ~Vector();

        void print();

        int get_dimension() const;

        TYPE get_item(int i);

        TYPE *get();

        void set_item(int i, TYPE value);

        TYPE& operator[](int i);

        TYPE Norm();

        TYPE modulus();

        // copying involved functions-----------------------------

        void set_vector(Vector &v);

        Vector &operator=(Vector &v2);

        Vector &operator=(double *);

        Vector &operator=(Vector &&v2);

        Vector &operator+(Vector &v2);

        Vector &operator-(Vector &v2);

        Vector &operator*(TYPE scalar);

        TYPE operator*(Vector &v);

        bool operator==(Vector &v);
    };

    Matrix *inverse(Matrix *matrix,Matrix* result = nullptr);

    Matrix identity(int m);

    Matrix *mull(Matrix *mat1, Matrix *mat, Matrix* result= nullptr);

    Matrix *Transpose_mull(Matrix *mat1, Matrix *mat, Matrix* result= nullptr);

    Vector *mull(Matrix *mat1, Vector *v,Vector* result= nullptr);

    class Matrix {
    private:

        int degree, span;
        TYPE *_matrix;
        TYPE det_;

        Matrix *pow_(int p, Matrix &mat);


    public:

        Matrix(int degree, int span);

        Matrix(int degree, int span, TYPE matrix_list[]); // it does not copy data
        // the given location will be treated as the data buffer

        ~Matrix();

        int get_n();

        int get_m();

        TYPE& get_item(int i, int j);

        void set_item(int row, int column, TYPE value);

        void print();

        bool operator==(Matrix &mat);

        Matrix &operator=(const Matrix &mat);

        Matrix &operator+(Matrix &mat);

        Matrix &operator+(Matrix &&mat);

        Matrix &operator-(Matrix &mat);

        Matrix &operator*(TYPE scalar);

        Vector &operator*(Vector &v);

        Matrix &operator*(Matrix &mat);

//        Matrix &operator*=(Matrix &mat);
        Matrix &inverse();

        TYPE det();

        Matrix &pow(int p);

        Matrix &transpose();

        friend Matrix *mull(Matrix *mat1, Matrix *mat, Matrix* result);

        friend Vector *mull(Matrix *mat1, Vector *v, Vector* result);

        friend Matrix *Transpose_mull(Matrix *mat1, Matrix *mat, Matrix* result);

        friend Matrix *inverse(Matrix *matrix, Matrix* result);

        friend class supMatrix;
    };


    class supMatrix{
        Matrix* baseMatrix;
        int rs,re,cs,ce;
        int R,C;
    public:
        supMatrix(Matrix*matrix,int rs, int re,int cs,int ce);
        TYPE &item(int x, int y);
        supMatrix* mull(supMatrix* other, supMatrix* result);
        supMatrix* add(supMatrix* other, supMatrix* result);

    };

    void Recursive_Matrix_Mull(Matrix* matrix1, Matrix* matrix2, Matrix* result);

}

namespace sparse_system {

    class Matrix;
    class AgMatrix;

    static SingleLinkedList::stack<Matrix *> mat_delete_list;
    static SingleLinkedList::stack<short_vector::Vector *> delete_list;

    void clear_memory();

    typedef struct _Item_ {
        int key;
        TYPE value;
        _Item_ *next;
    }item;

    void delete_row(item *row);

    item *q_sort(item *head, item *end);

    item *add_in_place(item *x, item *y);

    item *_add_(item *x, item *y, TYPE , TYPE);

    item *row_matrix_mul(item *row, item *arr[]);

    class Matrix {

    private:
        int span;
        int degree;
        item **array;
        double determinant = 1;
    protected:
        Matrix &copy();

    public:
        Matrix(int span, int degree);
        Matrix(Matrix& matrix);
        void set_item(int row, int column, TYPE value);

        void sort();

        void print();


        Matrix &operator+(Matrix &other);

        Matrix &operator-(Matrix &other);

        Matrix &operator*(TYPE scalar);

        Matrix &operator*(Matrix &other);

        short_vector::Vector & operator*(short_vector::Vector &vector);

        Matrix &inverse();

        ~Matrix();

        Matrix& pow(int p) ;

        Matrix& transpose();

        Matrix& diagonalize_trianglur_matrix(TYPE*& eigenvalues);

        friend Matrix &operator*(TYPE scalar, Matrix &matrix);
        friend Matrix* pow_(int p, Matrix &mat, Matrix* ) ;
        friend sparse_system::AgMatrix;
    };

    Matrix &operator*(TYPE scalar, Matrix &matrix);

    class AgMatrix:public Matrix{
    public:
        inline AgMatrix(int span , int degree):Matrix(span,degree)
        {
        }
        void solveSystem(short_vector::Vector& resultVector);
        inline void setDegree(int i) {
            degree = i;
        }
    };
}


namespace TreeBasedMatrix{

    typedef float type;

    typedef struct node{
        int key;
        type data[5];
        node* left;
        node* right;
        short height;
    }node;

    class Vector{
        node * root;
        int span,counter;
    public:SingleLinkedList::stack<node*> readerStack;
        type* reader;

    public:
        Vector();
        type& operator[](int i);
        void pushItem(type value);
        void print();
        type read();
        void start();
        void operator++();
        int getCounter();
        int getSpan();

        Vector& operator+=(Vector&);

    };



}

#endif //DATA_STRUCTURES_H_LINEARALGEBRA_H
