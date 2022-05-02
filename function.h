//
// Created by m7mdn on 2/6/2022.
//

#ifndef MAIN_PY_FUNCTION_H
#define MAIN_PY_FUNCTION_H


#include <math.h>
#include "DataStructures.h"

class functionPtr{
    long double (*func)(long double );

public:
    explicit functionPtr(long double (*func)(long double ));
    ~functionPtr()= default;
    long double evaluate(long double input);
};

enum operation{
    add, sub,mull,divide,power ,None
};
enum type{
    constant, node, func,null
};

// left + right
typedef struct funcNode {
    type leftLeaf, rightLeaf;
    operation opr;
    funcNode *left;
    funcNode *right;

    ~funcNode();
}funcNode;


typedef SingleLinkedList::stack<funcNode*> deleteListType;

class function {
private:
    funcNode* root;
    function* in_func;

private:
    function();

    function * binaryOperation(function &function1, operation opr);
    funcNode * selfAssing(function &function1, operation opr) ;

public:

    explicit function(long double (*func1)(long double ));

    explicit function(long double data);

    long double operator() (long double input) ;

    function & operator() (function& input) ;
//------------------------------------------------------------
    function & operator+(function function1);

    function & operator-(function function1);

    function & operator*(function function1);

    function & operator/(function function1);

    function & pow(function function1);

//------------------------------------------------------------
    function & operator += (function function1);

    function & operator -= (function function1);

    function & operator *= (function function1);

    function & operator /= (function function1);

};

void setList(deleteListType * dl);

void deleteSystem();

typedef struct Pitem{
    int power;
    long double value;
}Pitem;

bool operator>(Pitem& p1,Pitem p2);


class Polynomial{
private:
    long double *coff;
    int degree;

    void addRoot(long double root);
    void addScaled(Polynomial& p2, long double c);

public:
    explicit Polynomial(int degree);
    Polynomial(Polynomial&);
    Polynomial(int degree, long double * coff);
    Polynomial(int degree,function& func, long double rangeS, long double rangeE);

    ~Polynomial();
public:
    Polynomial& integrate();
    Polynomial& derivative();

public:
    long double operator()(long double x);
    Polynomial& operator+(Polynomial& p2);
    Polynomial& operator-(Polynomial& p2);
    Polynomial& operator*(Polynomial& p2);
//
    void operator+=(Polynomial& p2);
    void operator-=(Polynomial& p2);
    void operator*=(long double  p2);

    void mullWith(int degree, long double Coff);
    void addWith(int degree, long double Coff);

    friend void printf(Polynomial&);
};

long double integrate(function f);
long double doubleIntegrate(function f);
long double doubleIntegrate(function f, long double& slope);

#endif //MAIN_PY_FUNCTION_H
