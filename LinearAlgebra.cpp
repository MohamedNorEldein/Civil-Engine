//
// Created by m7md_nor on 9/21/2021.
//

#include "../..\header\LinearAlgebra.h"

namespace short_vector {

    Vector::Vector(int dimension) : dimension(dimension) {
        vector = (TYPE *) calloc(dimension, sizeof(TYPE));

    }

    Vector::Vector(int dimension, TYPE vector_[]) : dimension(dimension) {
        vector = (TYPE *) calloc(dimension, sizeof(TYPE));
        for (int i = 0; i < dimension; ++i) {
            vector[i] = vector_[i];
        }

    }

    Vector::Vector(const Vector &other) : dimension(other.dimension) {
        vector = (TYPE *) calloc(dimension, sizeof(TYPE));
        for (int i = 0; i < dimension; ++i) {
            vector[i] = other.vector[i];
        }
        clear_memory();
    }

    Vector::~Vector() {
//        std::cout << "delete vector at: " << this << '\t';
//        _print();
        free(vector);
    }

    void Vector::print() {
        std::cout << "< " << vector[0] << " ,";
        for (int i = 1; i < dimension - 1; ++i) {
            std::cout << ' ' << vector[i] << " ,";
        }
        std::cout << ' ' << vector[dimension - 1] << " >" << std::endl;
    }

    int Vector::get_dimension() const { return dimension; }

    TYPE Vector::get_item(int i) { return vector[i]; }

    TYPE *Vector::get() { return vector; }

    void Vector::set_item(int i, TYPE value) { vector[i] = value; }

    TYPE& Vector::operator[](int i) { return vector[i]; }

    TYPE Vector::Norm() {
        TYPE sum = 0;
        for (int i = 0; i < dimension; i++) {
            sum += vector[i] * vector[i];
        }
        return sum;
    }

    TYPE Vector::modulus() { return sqrt(Norm()); }

    void Vector::set_vector(Vector &v) {
        printf("%d, %d\n",dimension,v.dimension);
        if (dimension == v.get_dimension()) {
            for (int i = 0; i < dimension; ++i) {
                vector[i] = v[i];
            }
        } else throw DimensionError();
    }

    Vector &Vector::operator=(Vector &v2) {
        set_vector(v2);
        clear_memory();
        return *this;
    }

    Vector &Vector::operator=(Vector &&v2) {
        set_vector(v2);
        clear_memory();
        return *this;
    }

    Vector &Vector::operator=(double * v){

        for (int i = 0; i < dimension; ++i) {
            vector[i] = v[i];
        }
        return *this;
    }

    Vector &Vector::operator+(Vector &v2) {
        if (v2.get_dimension() == dimension) {
            auto v_sum = new Vector(dimension);
            for (int i = 0; i < dimension; i++) {
                v_sum->set_item(i, vector[i] + v2.get_item(i));
            }
                delete_list.push(v_sum);
            return *v_sum;
        } else {
            throw DimensionError();
        }
    }

    Vector &Vector::operator-(Vector &v2) {
        return *this + v2 * (TYPE(-1));
    }

    Vector &Vector::operator*(TYPE scalar) {
        auto v_sum = new Vector(dimension);
        for (int i = 0; i < dimension; i++) {
            v_sum->vector[i] = scalar * vector[i];
        }
            delete_list.push(v_sum);
        return *v_sum;
    }

    TYPE Vector::operator*(Vector &v) {
        if (v.dimension == dimension) {
            TYPE sum = 0;
            for (int i = 0; i < dimension; i++) {
                sum += v[i] * vector[i];
            }
            return sum;
        } else {
            throw DimensionError();
        }
    }

    bool Vector::operator==(Vector &v) {
        for (int i = 0; i < dimension; i++) {
            if (v[i] != vector[i]) { return false; }

        }
        return true;
    }

    Matrix::Matrix(int degree, int span) : span(span), degree(degree) {
        _matrix = (TYPE *) calloc(span * degree, sizeof(TYPE));

        det_ = 1;
    }

    Matrix::Matrix(int degree, int span, TYPE matrix_list[]) : span(span), degree(degree) {
        _matrix = (TYPE *) calloc(span * degree, sizeof(TYPE));
        _matrix = matrix_list;
        for (int i = 0; i < span * degree; i++) {
            _matrix[i] = matrix_list[i];
        }
//
        det_ = 1;
//        std::cout<< this<<'\t'<< this->_matrix<<'\span';

    }

    Matrix::~Matrix() {
        std::cout<<"delete matrix at"<< this<<'\t'<< this->_matrix <<'\n';
        delete[] _matrix;
    }

    int Matrix::get_n() { return span; }

    int Matrix::get_m() { return degree; }

    TYPE& Matrix::get_item(int i, int j) { return _matrix[i * span + j]; }

    void Matrix::set_item(int i, int j, TYPE value) { _matrix[i * span + j] = value; }

    void Matrix::print() {
        for (int i = 0; i < degree; ++i) {
            std::cout << "[";
            for (int j = 0; j < span; ++j) {
                if (-0.0001 < _matrix[i * span + j] and _matrix[i * span + j] < 0.0001) {
                    std::cout << ' ' << 0 << ' ';
                } else {
                    std::cout << ' ' << _matrix[i * span + j] << ' ';
                }
            }
            std::cout << "]\n";
        }
    }

    bool Matrix::operator==(Matrix &mat) {
        for (int i = 0; i < span * degree; i++) {
            if (_matrix[i] == mat._matrix[i]) {
                return false;
            }
        }
        return true;
    }

    void clear_memory() {
        printf("called\n");
        while (delete_list.length() != 0) {
            delete (delete_list.pop());
        }

        while (mat_delete_list.length() != 0) {
            delete (mat_delete_list.pop());
        }
    }

    Matrix &Matrix::operator=(const Matrix &mat) {
        if(mat.degree!=degree||mat.span!=span){throw DimensionError();}
        for (int i = 0; i < span * degree; i++) {
            _matrix[i] = mat._matrix[i];
        }
        clear_memory();
        return *this;
    }

    Matrix &Matrix::operator+(Matrix &mat) {
        if (degree == mat.degree and span == mat.span) {
            auto mat2 = new Matrix(span, degree);
            for (int i = 0; i < span * span; i++) {
                mat2->_matrix[i] = _matrix[i] + mat._matrix[i];
            }
                mat_delete_list.push(mat2);
            return *mat2;
        } else {
            throw DimensionError();
        }
    }

    Matrix &Matrix::operator+(Matrix &&mat) {
        if (degree == mat.degree and span == mat.span) {
            auto mat2 = new Matrix(span, degree);
            for (int i = 0; i < span * span; i++) {
                mat2->_matrix[i] = _matrix[i] + mat._matrix[i];
            }
            mat_delete_list.push(mat2);
            return *mat2;
        } else {
            throw DimensionError();
        }
    }

    Matrix &Matrix::operator-(Matrix &mat) {
        return *this + mat * (-1);
    }

    Matrix &Matrix::operator*(TYPE scalar) {
        auto mat2 = new Matrix(span, degree);
        for (int i = 0; i < span * span; i++) {
            mat2->_matrix[i] = this->_matrix[i] * scalar;
        }
        mat_delete_list.push(mat2);
        return *mat2;
    }

    Vector &Matrix::operator*(Vector &v) {

        short_vector::Vector* v2 = new Vector(degree);

        for (int i = 0; i < degree; i++) {
            TYPE a = 0;
            for (int j = 0; j < span; j++) {
                a += _matrix[i * span + j] * v[j];
            }
            v2->set_item(i, a);
        }
            delete_list.push(v2);
//        break_line
        return *v2;
    }

    Matrix &Matrix::operator*(Matrix &mat) {
        if (span == mat.degree) {
            auto mat2 = new Matrix(degree, mat.span);
            for (int k = 0; k < mat2->span; k++) {
                for (int i = 0; i < mat2->degree; i++) {
                    TYPE a = 0;
                    for (int j = 0; j < span; j++) {
                        a += _matrix[k * span + j] * mat._matrix[j * mat.span + i];
                    }
                    mat2->_matrix[k * mat2->degree + i] = a;
                }
            }
                mat_delete_list.push(mat2);
            return *mat2;
        } else {
            throw DimensionError();
        }
    }

    Matrix * mull(Matrix* mat1, Matrix* mat, Matrix* mat2){
        int span = mat1->span,degree= mat1->degree;
        TYPE *_matrix=mat1->_matrix;
        if (span == mat->degree) {
            if(!mat2)
                auto mat2 = new Matrix(degree, mat->span);
            for (int k = 0; k < mat2->span; k++) {
                for (int i = 0; i < mat2->degree; i++) {
                    TYPE a = 0;
                    for (int j = 0; j < span; j++) {
                        a += _matrix[k * span + j] * mat->_matrix[j * mat->span + i];
                    }
                    mat2->_matrix[k * mat2->degree + i] = a;
                }
            }
            return mat2;
        } else {
            throw DimensionError();
        }
    }

    Vector * mull(Matrix* mat1, Vector* v, Vector* v2) {
        if(!v2) {
            short_vector::Vector *v2 = new Vector(mat1->degree);
        }
        for (int i = 0; i < mat1->degree; i++) {
            TYPE a = 0;
            for (int j = 0; j < mat1->span; j++) {
                a += mat1->_matrix[i * mat1->span + j] * (*v)[j];
            }
            v2->set_item(i, a);
        }
        return v2;
    }


        Matrix &Matrix::inverse() {
        if (span == degree) {
            auto mat = new Matrix(span, span);
            for (int i = 0; i < degree; i++) {
                mat->set_item(i, i, 1);
            }
            auto mat1 = new Matrix(span, span);
            *mat1 = *this;
            for (int i = 0; i < degree; ++i) {
                if (mat1->_matrix[i * degree + i] == 0) {
                    int k = 0;
                    while (k < degree and  mat1->_matrix[k * degree + i] == 0) {
                        for (int j = 0; j < degree; j++) {
                            mat1->_matrix[k * degree + j] +=  mat1->_matrix[i * degree + j];
                            mat->_matrix[k * degree + j] +=  mat->_matrix[i * degree + j];
                        }
                        k++;
                    }
                }
            }

            for (int i = 0; i < degree; i++) {
                for (int k = 0; k < degree; k++) {
                    if (k != i) {
                        double constant = mat1->_matrix[k * degree + i] / mat1->_matrix[i * degree + i];
                        for (int j = 0; j < degree; j++) {
                            mat1->_matrix[k * degree + j] -= constant * mat1->_matrix[i * degree + j];
                            mat->_matrix[k * degree + j] -= constant * mat->_matrix[i * degree + j];

                        }
                    }
                }
            }

            for (int i = 0; i < degree; i++) {
                det_ *= mat1->_matrix[i * degree + i];
                if (mat1->_matrix[i * degree + i] == 0) { throw UnsolvableSystem(); }
                for (int j = 0; j < degree; j++) {
                    mat->_matrix[i * degree + j] /= mat1->_matrix[i * degree + i];
                }
            }
            delete mat1;
//            mat_delete_list.push(mat);
            mat->degree=degree;
            mat->span=span;
            return *mat;
        } else {
            throw RedundantSystem();
        }

    }

void printf(int i){
    std::printf("%d\n",i);
    }

    Matrix* inverse(Matrix* matrix, Matrix* mat) {
        int span = matrix->span, degree= matrix->degree;
        if (matrix->span == matrix->degree) {
            if(mat== nullptr) {
                auto mat = new Matrix(span, span);
            }
            for (int i = 0; i < degree; i++) {
                mat->set_item(i, i, 1);
            }
            auto mat1 = new Matrix(span, span);
            *mat1 = *matrix;

            for (int i = 0; i < degree; i++) {
                for (int k = 0; k < degree; k++) {
                    if (mat1->_matrix[i * degree + i] == 0) {
                        int s = i+1;
                        while (s < degree and mat1->_matrix[s * degree + i] == 0) {
                            s++;
                        }
                        if (s < degree) {
//                        std::printf("%d %d\n",i,s);
                            for (int j = 0; j < degree; j++) {
                                mat1->_matrix[i * degree + j] += mat1->_matrix[s * degree + j];
                                mat->_matrix[i * degree + j] += mat->_matrix[s * degree + j];
                            }
                        }
                    }
                    if (k != i) {
                        double constant = mat1->_matrix[k * degree + i] / mat1->_matrix[i * degree + i];
                        for (int j = 0; j < degree; j++) {
                            mat1->_matrix[k * degree + j] -= constant * mat1->_matrix[i * degree + j];
                            mat->_matrix[k * degree + j] -= constant * mat->_matrix[i * degree + j];
                        }
                    }
                }
            }

            for (int i = 0; i < degree; i++) {
                if (mat1->_matrix[i * degree + i] == 0) { throw UnsolvableSystem(); }
                for (int j = 0; j < degree; j++) {
                    mat->_matrix[i * degree + j] /= mat1->_matrix[i * degree + i];
                }
            }
            delete mat1;
            mat->degree=degree;
            mat->span=span;
            return mat;
        } else {
            throw RedundantSystem();
        }

    }

    TYPE Matrix::det() {
        this->inverse();
        return det_;
    }

    Matrix &Matrix::pow(int p) {
        if (p == 0) {
            auto mat = new Matrix(span, span);
            for (int i = 0; i < degree; i++) {
                mat->set_item(i, i, 1);
            }
                mat_delete_list.push(mat);
            return *mat;
        } else if (p > 0) {
            clear_memory();
            return *pow_(p, *this);
        } else if (p < 0) { return (pow_(p, *this))->inverse(); }
    }

    Matrix* Matrix::pow_(int p, Matrix &mat) {
        if (p == 1) {
            return &mat;
        }
        Matrix *product = (Matrix *) malloc(sizeof(Matrix));
        product->degree = degree;
        product->span = span;

        if (p == 2) {
            product->_matrix = (mat * mat)._matrix;
            return product;
        }
        if (p % 2 == 0) {
            product->_matrix = pow_(p / 2, mat)->_matrix;
            product->_matrix = ((*product) * (*product))._matrix;
            product = mat_delete_list.pop();
            clear_memory();
            mat_delete_list.push(product);
            return product;
        }
        product->_matrix = pow_(p / 2, mat)->_matrix;
        product->_matrix = ((*product) * (*product) * mat)._matrix;
        product = mat_delete_list.pop();
        clear_memory();
        mat_delete_list.push(product);
        return product;
    }

    Matrix &Matrix::transpose() {
        auto mat_transpose = new Matrix(span,degree);
        for (int i = 0; i < degree; ++i) {
            for (int j = 0; j < span; ++j) {
                mat_transpose->_matrix[j * degree + i] = _matrix[i * span + j];
            }
        }
            mat_delete_list.push(mat_transpose);
        return *mat_transpose;
    }

    Matrix *Transpose_mull(Matrix *mat1, Matrix *mat, Matrix* result) {
        int span = mat1->span,degree= mat1->degree;
        TYPE *_matrix=mat1->_matrix;
        if (span == mat->degree) {

            if(!result) {
                result = new Matrix(degree, mat->span);
            }

            for (int k = 0; k < result->span; k++) {
                for (int i = 0; i < result->degree; i++) {
                    TYPE a = 0;
                    for (int j = 0; j < span; j++) {
                        a += _matrix[j * span + k] * mat->_matrix[j * mat->span + i];
                    }
                    result->_matrix[k * result->degree + i] = a;
                }
            }
            return result;
        } else {
            throw DimensionError();
        }

    }


    Matrix identity(int m) {
        Matrix mat(m, m);
        for (int i = 0; i < m; i++) {
            mat.set_item(i, i, 1);
        }
        return mat;
    }

    supMatrix::supMatrix(Matrix *matrix, int rs, int re, int cs, int ce) :
    baseMatrix(matrix),rs(rs),re(re),cs(cs),ce(ce),R(re-rs),C(ce-cs)
    {
    }

    TYPE &supMatrix::item(int x, int y) {
        return baseMatrix->get_item(x,y);
    }

    Matrix* innerMatrix;

    supMatrix *supMatrix::mull(supMatrix *other, supMatrix *result) {


        return nullptr;
    }


}

namespace sparse_system {


    void clear_memory() {
        while (delete_list.length() != 0) {
            delete (delete_list.pop());
        }

        while (mat_delete_list.length() != 0) {
            delete (mat_delete_list.pop());
        }
    }


    void delete_row(item *row) {
        if (row == 0) { return; }

        delete_row(row->next);
        delete row;
    }

    item *q_sort(item *head, item *end) {
        if (head == end) { return head; }

        item *start = head;
        item *reader = head;

        while (reader->next != end) {
            if (reader->next->key < head->key) {
                item *a = reader->next;
                reader->next = reader->next->next;
                a->next = start;
                start = a;
            }
            else {
                reader = reader->next;
            }
        }
        head->next = q_sort(head->next, end);
        start = q_sort(start, head);
        return start;
    }

    item *add_in_place(item *x, item *y) {
        if (x == nullptr) { return y; }
        if (y == nullptr) { return x; }
        item *z;
        if ((x->key) < (y->key)) {
            x->next = add_in_place(x->next, y);
            return x;
        }
        if ((x->key) == (y->key)) {
            x->value = x->value + y->value;
            x->next = add_in_place(x->next, y->next);
            delete y;
            return x;
        }
        z = y->next;
        y->next = x;
        y->next = add_in_place(x, z);
        delete x;
        return y;
    }

    item *_add_(item *x, item *y, TYPE scalar_x = 1, TYPE scalar_y = 1) {

        if (x == 0 and y == 0) {
            return 0;
        }
        item *result;
        if (x == 0) { return new item{y->key, y->value * scalar_y, _add_(0, y->next, 0, scalar_y)}; }
        if (y == 0) { return new item{x->key, scalar_x * x->value, _add_(x->next, 0, scalar_x, 0)}; }

        if ((x->key) < (y->key)) {
            result = new item{x->key, scalar_x * x->value, 0};
            result->next = _add_(x->next, y, scalar_x, scalar_y);
            return result;
        }
        if ((x->key) > (y->key)) {
            result = new item{y->key, y->value * scalar_y, 0};
            result->next = _add_(x, y->next, scalar_x, scalar_y);
            return result;
        }

        TYPE val = x->value * scalar_x + scalar_y * y->value;
        if (val == 0) {
            result = _add_(x->next, y->next, scalar_x, scalar_y);
            return result;

        }
        result = new item{x->key, x->value * scalar_x + scalar_y * y->value, 0};
        result->next = _add_(x->next, y->next, scalar_x, scalar_y);
        return result;

    }

    item *row_matrix_mul(item *row, item *arr[]) {

        item *result = nullptr;
        item *other = nullptr;
        while (row != nullptr) {
            other = _add_(result, arr[row->key], 1, row->value);
            delete_row(result);
            row = row->next;
            result = other;
        }
        return result;
    }

    int compare(const void* first,const void* second){
//        printf("%d, %x\n",(*(item**)first)->key,(*(item**)second)->key);
        return ((*(item**)first)->key) - ((*(item**)second)->key) ;
    }

    Matrix &Matrix::copy() {
        auto other = new Matrix(degree, span);

        for (int i = 0; i < degree; i++) {
            other->array[i] = _add_(array[i], nullptr);
        }
        mat_delete_list.push(other);
        return *other;
    }

    Matrix::Matrix(int span, int degree) : span(span), degree(degree) {
        array = (item **) calloc(degree, sizeof(double));
    }

    void Matrix::set_item(int row, int column, TYPE value) {
        if (row < degree && column < span)
        {
            if(value !=0)
            {
                item *a = new item{column, value, array[row]};
                array[row] = a;
            }
        } else {
            throw DimensionError();
        }
    }

    void Matrix::sort() {
        for (int i = 0; i < degree; i++) {
            array[i] = q_sort(array[i], nullptr);
        }
//        _print();
//        break_line
        qsort(array,degree,sizeof(void *),compare );
    }

    void Matrix::print() {
        break_line;
        for (int i = 0; i < degree; i++) {
            item *st = array[i];
            std::cout << "row[" << i << ']' << ':' ;
            while (st != nullptr) {
//                if(std::abs(st->value)>0.0001){
                    std::cout << st->key << ':' << st->value << '\t';
//                }
                st = st->next;
            }
            std::cout << '\n';
        }
    }


    Matrix &Matrix::operator+(Matrix &other) {
        if ((degree == other.degree) && (span == other.span)) {
            auto sum = new Matrix(degree, span);
            for (int i = 0; i < degree; i++) {
                sum->array[i] = _add_(array[i], other.array[i]);
            }
            mat_delete_list.push(sum);
            return *sum;
        }
        throw DimensionError();
    }

    Matrix &Matrix::operator-(Matrix &other) {
        if ((degree == other.degree) && (span == other.span)) {
            auto sum = new Matrix(degree, span);
            for (int i = 0; i < degree; i++) {
                sum->array[i] = _add_(array[i], other.array[i], 1, -1);
            }
            mat_delete_list.push(sum);
            return *sum;
        }
        throw DimensionError();
    }

    Matrix &Matrix::operator*(TYPE scalar) {

        auto sum = new Matrix(degree, span);
        for (int i = 0; i < degree; i++) {
            sum->array[i] = _add_(array[i], nullptr, scalar, 0);
        }
        mat_delete_list.push(sum);
        return *sum;

    }


    Matrix &Matrix::operator*(Matrix &other) {
        if (span == other.degree) {
            auto sum = new Matrix(degree, other.span);
            for (int i = 0; i < degree; i++) {
                sum->array[i] = row_matrix_mul(array[i], other.array);
            }
            mat_delete_list.push(sum);

            return *sum;
        }
        throw DimensionError();
    }

    Matrix &Matrix::inverse() {
        Matrix mat1 = copy();

        Matrix *mat2 = new Matrix(degree, span);

        for (int i = 0; i < degree; ++i) {
            mat2->set_item(i, i, 1);
        }
        mat1.sort();
        for (int i = 0; i < degree; ++i) {

            item *x, *y;


            while (mat1.array[i]->key != i) {

                x = _add_(mat1.array[mat1.array[i]->key], mat1.array[i], -1,
                          double(mat1.array[mat1.array[i]->key]->value) / mat1.array[i]->value);
                y = _add_(mat2->array[mat1.array[i]->key], mat2->array[i], -1,
                          double(mat1.array[mat1.array[i]->key]->value) / mat1.array[i]->value);

                delete_row(mat2->array[i]);
                delete_row(mat1.array[i]);
                mat1.array[i] = x;
                mat2->array[i] = y;
            }
            mat1.print();
            break_line
        }
        for (int i = degree - 1; i > -1; --i) {
            item *x, *y, *z;
            z = mat1.array[i]->next;
            while (z != 0) {
                x = _add_(mat1.array[z->key], mat1.array[i], -1, TYPE (mat1.array[z->key]->value) / z->value);
                y = _add_(mat2->array[z->key], mat2->array[i], -1, TYPE (mat1.array[z->key]->value) / z->value);

                delete_row(mat2->array[i]);
                delete_row(mat1.array[i]);
                mat1.array[i] = x;
                mat2->array[i] = y;
                z = mat1.array[i]->next;
            }
        }

        for (int i = 0; i < degree; i++) {
            determinant *= mat1.array[i]->value;
            mat2->array[i] = _add_(mat2->array[i], nullptr, TYPE (1) / mat1.array[i]->value, 0);
        }

        mat_delete_list.push(mat2);

        return *mat2;
    }


    Matrix::~Matrix() {
        std::cout << "delete matrix at : " << this << '\n';
        for (int i = 0; i < degree; i++) {
            delete_row(array[i]);
        }
    }

    Matrix &Matrix::pow(int p) {
        if (span != degree) { throw DimensionError(); }
        if (p == 0) {
            auto mat = new Matrix(span, span);
            for (int i = 0; i < degree; i++) {
                mat->set_item(i, i, 1);
            }
            mat_delete_list.push(mat);
            return *mat;
        } else if (p > 0) {
            clear_memory();
            Matrix *product = (Matrix *) malloc(sizeof(Matrix));
            product->degree = degree;
            product->span = span;
            mat_delete_list.push(product);

            return *pow_(p, *this, product);
        } else if (p < 0) {
            Matrix *product = (Matrix *) malloc(sizeof(Matrix));
            product->degree = degree;
            product->span = span;
            mat_delete_list.push(product);
            return (pow_(p, *this, product))->inverse();
        }
    }

    Matrix *pow_(int p, Matrix &mat, Matrix *product) {

        if (p == 1) {
            return &mat;
        }

        if (p == 2) {
            product->array = (mat * mat).array;
            clear_memory();

            return product;
        }
        if (p % 2 == 0) {
            product->array = pow_(p / 2, mat, product)->array;
            product->array = ((*product) * (*product)).array;
            clear_memory();

            return product;
        }
        product->array = pow_(p / 2, mat, product)->array;
        product->array = ((*product) * (*product) * mat).array;
        clear_memory();

        return product;
    }

    Matrix &Matrix::transpose() {
        auto sum = new Matrix(span, degree);
        for (int i = 0; i < degree; i++) {
            item *st = array[i];
            while (st != 0) {
                sum->set_item(st->key, i, st->value);
                st = st->next;
            }
        }
        sum->sort();
        mat_delete_list.push(sum);
        return *sum;
    }

    Matrix &operator*(TYPE scalar, Matrix &matrix) {

        auto sum = new Matrix(matrix.degree, matrix.span);
        for (int i = 0; i < matrix.degree; i++) {
            sum->array[i] = _add_(matrix.array[i], nullptr, scalar, 0);
        }
        mat_delete_list.push(sum);

        return *sum;

    }

    short_vector::Vector &Matrix::operator*(short_vector::Vector &vector) {
        auto sum = new short_vector::Vector(degree);
        for (int i = 0; i < degree; i++) {
            item *st = array[i];
            TYPE value = 0;
            while (st != 0) {
                value += st->value * vector[st->key];
                st = st->next;
            }
            sum->set_item(i, value);
        }
        delete_list.push(sum);

        return *sum;

    }

    Matrix &Matrix::diagonalize_trianglur_matrix(TYPE *& eigenvalues) {

        auto mat2 = new Matrix(degree, span);
        eigenvalues = (TYPE*) malloc(sizeof(TYPE)*span);
        for (int j = 0; j < degree; ++j)
        {
            eigenvalues[j] = array[j]->value;
            TYPE eigenvector[span];

            for (int i = degree - 1; i > -1; --i)
            {
                item *z;
                z = array[i]->next;
                eigenvector[i] = 0;

                while (z != nullptr) {
                    eigenvector[i] -= z->value * eigenvector[z->key];
                    z = z->next;
                }

                if(array[i]->value==eigenvalues[j])
                {
                    eigenvector[i] = 1;
                }else{
                    eigenvector[i] = eigenvector[i] / (array[i]->value - eigenvalues[j]);
                }

                mat2->set_item(i,j,eigenvector[i]);
            }
        }
        mat2->sort();
        return *mat2;
    }

    Matrix::Matrix(Matrix &matrix1):array(matrix1.array),span(matrix1.span),degree(matrix1.degree)
    {

    }

    void AgMatrix::solveSystem(short_vector::Vector &resultVector) {

        for (int i = 0; i < degree; ++i) {
            item *x, *y;
            while (array[i]->key != i) {
//                printf("%d %f %d %f\n",i,array[i]->key,array[array[i]->key]->value,array[i]->value);
                x = _add_(array[array[i]->key], array[i], -1,
                          TYPE (array[array[i]->key]->value) / array[i]->value);

                delete_row(array[i]);
                array[i] = x;
            }
        }

        for (int i = degree-1; i >-1 ; --i) {
            TYPE sum = 0;
            item * st=array[i]->next;
            while (st!= nullptr) {
                if (st->key == span-1) {
                    sum-=st->value;
                } else {
                    sum += st->value * resultVector[st->key];
                }
                st = st->next;
            }
            resultVector.set_item(i, -sum /array[i]->value);
        }
    }

}

namespace TreeBasedMatrix{
    Vector::Vector() {
        root= nullptr;
        span=0;
    }

    int height(node *parent) {
        if (parent == nullptr) { return -1; }
        return parent->height;
    }

    void set_height(node *parent) {
        if (parent == nullptr) { return; }
        parent->height = 1 + std::max(height(parent->right), height(parent->left));
    }

    node* right_r(node* x) {
        node *y = x->right;
        node *b = y->left;
        y->left = x;
        x->right = b;
        set_height(b);
        set_height(x);

        return y;
    }

    node*  left_r(node*  y) {
        node*  x = y->left;
        node*  b = x->right;
        y->left = b;
        x->right = y;
        set_height(b);
        set_height(y);
        return x;
    }

    int balance(node*  parent) {
        if(parent== nullptr){return 0;}
        return height(parent->left) - height(parent->right);
    }

    node*  avl_function(node*  parent) {

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


    node * push(int key, type value,node* parent){
//        printf("root = %x", parent);

        if(parent== nullptr){
            node * p2= (node *) malloc(sizeof(node));
            p2->key=key;
            p2->left= nullptr;
            p2->right= nullptr;
            p2->data[0] = value;
            p2->height=0;
            return p2;
        }
        if(key - parent->key < 5){
            parent->data[key%5] = value;
            return parent;
        }

        parent->right = push(key,value,parent->right);
        return avl_function(parent);
    }

    void Vector::pushItem(type value) {
        root = push(span,value,root);
        span++;
    }

    void printData(type* data){

        for (int i = 0; i < 1; ++i) {
            printf("%f ",data[i]);
        }
        printf("\t");
    }

    void _print(node* parent){
        if(parent== nullptr){
            return;
        }
        printData(parent->data);

        _print(parent->left);
        _print(parent->right);
    }
    void Vector::print() {
        _print(root);
    }

    type & find( node* parent,int i){
        if(i<parent->key){
            return find(parent->left,i);
        }
        if(i>5+parent->key){
            return find(parent->right,i);
        }
        return parent->data[i%5];
    }

    type & Vector::operator[](int i){
        return find(root,i);
    }

    void start( SingleLinkedList::stack<TreeBasedMatrix::node*>& readerStack, node* root, type *reader) {

        node *p = root;
        while (p != nullptr) {
            readerStack.push(p);
            p = p->left;
        }
        reader = readerStack.top()->data;
    }

    void Vector::start() {
        ::TreeBasedMatrix::start(readerStack,root,reader);
        reader = readerStack.top()->data;
    }

    type Vector::read(){
        return *reader;
    }

    void Vector::operator++() {
        counter ++;
        if(counter%5>0){
            reader++;
            return;
        }
        if (readerStack.top()->right != 0) {
            ::TreeBasedMatrix::start(readerStack, readerStack.pop()->right, reader);
            reader = readerStack.top()->data;
            return;
        }
        readerStack.pop();
        reader = readerStack.top()->data;
    }

    int Vector::getCounter(){
        return counter;
    }
    int Vector::getSpan(){
        return span;
    }
   /*
    Vector& Vector::operator+=(Vector &) {

        while (){

        }
    }
    */
}