//
// Created by m7mdn on 2/6/2022.
//

#include "../../header/function.h"
#include <iostream>

deleteListType *deleteList;


long double evaluate(funcNode* parent, long double x);

long double evaluateLeft(funcNode* parent,long double x) ;

long double evaluateRight(funcNode*parent,long double x) ;


functionPtr::functionPtr(long double (*func)(long double)) : func(func) {
}

long double functionPtr::evaluate(long double input) {
    return func(input);
}

long double evaluateLeft(funcNode *parent, long double x) {
    if(parent->leftLeaf==null){
        return 0;
    }
    if (parent->leftLeaf == func) {
        return ((functionPtr *) parent->left)->evaluate(x);
    }
    if (parent->leftLeaf == constant) {
//        printf("costant %lf\n",*(long double *) parent->left);
        return *(long double *) parent->left;
    }
    if (parent->rightLeaf == node) {
        return evaluate(parent->left, x);
    }
}

long double evaluateRight(funcNode *parent, long double x) {
    if(parent->rightLeaf==null){
        return 0;
    }
    if (parent->rightLeaf == func) {
        return ((functionPtr *) parent->right)->evaluate(x);
    }
    if (parent->rightLeaf == constant) {
        return *(long double *) parent->right;
    }
    if (parent->rightLeaf == node) {
        return evaluate(parent->right, x);
    }
}

long double evaluate(funcNode *parent, long double x) {
    if (parent->opr == None) {
        return evaluateLeft(parent, x);
    }
    if (parent->opr == add) {
        return evaluateLeft(parent, x) + evaluateRight(parent, x);
    }
    if (parent->opr == sub) {
        return evaluateLeft(parent, x) - evaluateRight(parent, x);
    }
    if (parent->opr == mull) {
        return evaluateLeft(parent, x) * evaluateRight(parent, x);
    }
    if (parent->opr == divide) {
        return evaluateLeft(parent, x) / evaluateRight(parent, x);
    }
    if (parent->opr == power) {
        return pow(evaluateLeft(parent, x) , evaluateRight(parent, x));
    }
    throw 0;
}


function::function() {
    root = nullptr;
    in_func= nullptr;
}

function *function::binaryOperation(function &function1, operation opr) {
    function *result = new function;
    if (function1.root->right == nullptr and this->root->right == nullptr) {
        result->root = new funcNode{this->root->leftLeaf, function1.root->leftLeaf, opr, this->root->left,
                                    function1.root->left};
        deleteList->push(result->root);

        return result;
    }
    if (this->root->right == nullptr) {
        result->root = new funcNode{this->root->leftLeaf, node, opr, this->root->left, function1.root};
        deleteList->push(result->root);

        return result;
    }
    if (function1.root->right == nullptr) {
        result->root = new funcNode{node, function1.root->leftLeaf, opr, this->root, function1.root->left};
        deleteList->push(result->root);

        return result;
    }
    result->root = new funcNode{node, node, opr, this->root, function1.root};
    deleteList->push(result->root);

    return result;
}

funcNode *function::selfAssing(function &function1, operation opr) {
    funcNode *result = new funcNode ;
    if (function1.root->rightLeaf == null and this->root->rightLeaf == null) {
        result = new funcNode{this->root->leftLeaf, function1.root->leftLeaf, opr, this->root->left,
                                    function1.root->left};
        deleteList->push(result);
        return result;
    }
    if (this->root->rightLeaf == null) {
        result = new funcNode{this->root->leftLeaf, node, opr, this->root->left, function1.root};
        deleteList->push(result);

        return result;
    }
    if (function1.root->rightLeaf == null) {
        result = new funcNode{node, function1.root->leftLeaf, opr, this->root, function1.root->left};
        deleteList->push(result);

        return result;
    }
    result = new funcNode{node, node, opr, this->root, function1.root};
    deleteList->push(result);

    return result;
}

function::function(long double (*func1)(long double)) {
    root = new funcNode;

    root->leftLeaf = func;

    root->left = (funcNode *) (new functionPtr(func1));
//    printf("%lf\t%lf\n", func1(5), ((functionPtr *) (root->left))->evaluate(5));

    root->opr = None;

    root->rightLeaf = null;
    root->right = nullptr;
    in_func= nullptr;
//    deleteList->push(root);

}

funcNode::~funcNode() {
//    printf("delete node :%x \n", this);

}

function::function(long double data) {
    root = new funcNode;
    in_func=0;

    root->leftLeaf = constant;
    root->left = (funcNode *) (new long double(data));
    root->opr = None;
//    printf("constant node :%x \n", this->root->left);

    root->rightLeaf = constant;
    root->right = nullptr;
    deleteList->push(root);
    deleteList->push(root->left);

//    delete in_func;
}

long double function::operator()(long double input) {
//    printf("%x",in_func);
    if(in_func == nullptr) {
        return ::evaluate(this->root, input);
    }
    long double d = (*in_func)(input);
    return ::evaluate(this->root,d);
}

function& function::operator()(function& input) {

    function *f =  new function;
    f->root = this->root;
    f->in_func = &input;
    return *f;
}

function &function::operator+(function function1) {
    return *binaryOperation(function1, add);
}

function &function::operator-(function function1) {
    return *binaryOperation(function1, sub);
}

function &function::operator*(function function1) {
    return *binaryOperation(function1, mull);
}

function &function::operator/(function function1) {
    return *binaryOperation(function1, divide);
}


function &function::pow(function function1) {
    return *binaryOperation(function1, power);
}

//--------------------------------------------
function& function::operator += (function function1) {
    this->root = selfAssing(function1, add);
    return *this;
}

function& function::operator -= (function function1) {
    this->root = selfAssing(function1, sub);
    return *this;

}

function& function::operator *= (function function1) {
    this->root = selfAssing(function1, mull);
    return *this;
}

function& function::operator /= (function function1) {
    this->root = selfAssing(function1, divide);
    return *this;
}


void setList(deleteListType *dl) {
    deleteList = dl;
}

void deleteSystem() {
    while (true) {
        try {
            delete deleteList->pop();
        } catch (ListIsEmpty) {
            break;
        }
    }
}


Polynomial::Polynomial(int degree) :degree(degree)
{
    coff = (long double *) calloc(degree+1,sizeof(long double ) );
}

Polynomial::Polynomial(int degree, long double *coff) :degree(degree), coff(coff)
{
}

Polynomial &Polynomial::operator+(Polynomial &p2) {
    Polynomial *result = new Polynomial(std::max(degree,p2.degree));
    for (int i = 0; i < degree+1; ++i) {
        result->coff[i] = this->coff[i] + p2.coff[i];
    }
    for (int i = degree+1; i < result->degree+1; ++i) {
        result->coff[i] = - p2.coff[i];
    }
    return *result;
}

Polynomial &Polynomial::operator-(Polynomial &p2) {
    Polynomial *result = new Polynomial(std::max(degree,p2.degree));
    for (int i = 0; i < degree+1; ++i) {
        result->coff[i] = this->coff[i] - p2.coff[i];
    }
    for (int i = degree+1; i < result->degree+1; ++i) {
        result->coff[i] = - p2.coff[i];
    }

    return *result;
}

void Polynomial::operator+=(Polynomial &p2) {
    if(p2.degree>degree){
        long double * newcoff = (long double *) calloc(p2.degree+1,sizeof(long double ) );

        for (int i = 0; i < degree+1; ++i) {
            newcoff[i] = this->coff[i] + p2.coff[i];
        }
        for (int i = degree+1; i < p2.degree+1; ++i) {
            newcoff[i] = + p2.coff[i];
        }
        degree = p2.degree;
        delete[] coff;
        coff = newcoff;
        return;
    }
    int deg=std::min(degree, p2.degree);
    for (int i = 0; i < deg+1; ++i) {
        coff[i] += p2.coff[i] ;
    }

}

void Polynomial::operator-=(Polynomial &p2) {
    if(p2.degree>degree){
        long double * newcoff = (long double *) calloc(p2.degree+1,sizeof(long double ) );

        for (int i = 0; i < degree+1; ++i) {
            newcoff[i] = this->coff[i] - p2.coff[i];
        }
        for (int i = degree+1; i < p2.degree+1; ++i) {
            newcoff[i] = - p2.coff[i];
        }
        degree = p2.degree;
        delete[] coff;
        coff = newcoff;
        return;
    }
    int deg=std::min(degree, p2.degree);
    for (int i = 0; i < deg+1; ++i) {
        coff[i] -= p2.coff[i];
    }
}

void Polynomial::addWith(int degree2, long double Coff) {
    if(degree2>degree){
        throw "out of space";
    }

       coff[degree2] += Coff;
}



long double Polynomial::operator()(long double x) {

    long double result=0;

    for (int i = 0; i < degree +1; ++i) {
        result += coff[i] * pow(x,i);
    }

    return result;
}

Polynomial::Polynomial(Polynomial &p2) :degree(p2.degree), coff(p2.coff)
{

}

Polynomial &Polynomial::operator*(Polynomial &p2) {
    Polynomial *result = new Polynomial(degree+p2.degree);

        for (int i = 0; i < p2.degree+1; ++i) {

            for (int j = 0; j <degree+1; ++j) {
//                printf("%d, %f\n", i-1, p2.coff[i-2]);
                result->coff[j+i] += this->coff[j] * p2.coff[i];
        }

    }

    return *result;

}

void printf(Polynomial& polynomial){
    printf("p = ");
    for (int i = 0; i < polynomial.degree+1; ++i) {
        printf(" +%lf*x^%d ", polynomial.coff[i],i);
    }
    printf("\n");
}

void Polynomial::operator*=(long double p2) {

}

void Polynomial::mullWith(int deg, long double Coff) {

    for (int i =  degree-deg+1; i>-1; --i) {
            coff[i+deg] = Coff*coff[i] ;
    }
    for (int i =  deg-1; i>-1; --i) {
        coff[i] = 0 ;
    }
    degree+= deg;

}

void Polynomial::addRoot(long double root) {

    for (int i =  degree; i>-1; --i) {
        coff[i] = coff[i-1]-root*coff[i] ;
    }
}

Polynomial &Polynomial::integrate() {
    Polynomial* result = new Polynomial(degree+2);
    for (int i =  degree; i>-1; --i) {
        result->coff[i+1] = coff[i]/(i+1);
    }
    coff[0]=0;
    return *result;
}

Polynomial &Polynomial::derivative() {
    Polynomial* result = new Polynomial(degree);
    for (int i =  1; i<degree; ++i) {
        result->coff[i-1] = coff[i]*(i);
    }
    return *result;
}

Polynomial::Polynomial(int deg, function &func, long double rangeS, long double rangeE):degree(int(deg*1.25)) {
    coff = (long double *) calloc(deg+1,sizeof(long double ) );
    long double x = rangeS, step = (rangeE-rangeS)/(deg+1),B;
    long double** a = new long double* [deg];
    a[0] = new long double [deg];
    long double fmax=func(rangeS);
    for (int j = 0; j < deg ; ++j) {
        a[0][j] = func(x);
//        if(std::abs(a[0][j])>fmax){
//            fmax=std::abs(a[0][j]);
//        }
        x+=step;
    }
//    for (int j = 0; j < deg ; ++j) {
//        a[0][j] *= fmax;
//    }
    for (int i = 1; i < deg; ++i) {
        a[i] = new long double [deg-i];
        for (int j = 0; j < deg -i; ++j) {
            a[i][j] = a[i-1][j+1]-a[i-1][j];
            a[i][j] /=(i * step) ;

        }
    }

    Polynomial p(degree+1);
    p.coff[0]=1;
    x = rangeS;
    coff[0]=func(x);
    for (int i = 0; i < deg-1; ++i) {
//
        p.addRoot(x);
//        printf(p);
//        printf("%lf >> %lf\n",x,p(x));

        this->addScaled(p,a[i+1][0]);
        x+=step;
    }

    for (int i = 0; i < deg; ++i) {
        delete[] a[i];
    }
}

void Polynomial::addScaled(Polynomial &p2, long double c) {
    int deg=std::min(degree, p2.degree);
    for (int i = 0; i < deg+1; ++i) {
        coff[i] += p2.coff[i] * c;
    }
}

Polynomial::~Polynomial() {
    delete[] coff;
    printf("deleting\n");
}

long double Gn[90]= {
        0.999646971286638437463248, 0.00090593237121483309373,
        0.998140379938568153561306, 0.00210771877452632989148,
        0.995431812058344663926755, 0.00330886724333601819543,
        0.991523928811062786129147, 0.00450612361367497786414,
        0.986421365057832848734254, 0.005697981560747315260085,
        0.980130251345148385458953, 0.00688298320846328431473,
        0.972658162090193139997465, 0.00805969494462001565867,
        0.964014098171505483393667, 0.00922669695714199094032,
        0.954208473881500336160720, 0.01038258230989321461381,
        0.943253103645357768153575, 0.01152595788914805885059,
        0.931161187500432007005847, 0.01265544583716812886888,
        0.917947295066586383372356, 0.01376968511233709343075,
        0.903627347931302693869986, 0.01486733308804332405038,
        0.888218600434745981298376, 0.01594706715100663901321,
        0.871739618862903434474028, 0.01700758628522267570940,
        0.854210259067071882286021, 0.0180476863446023616405,
        0.835651642533377045564199, 0.01906589303913731842532,
        0.816086130929481056443754, 0.02006120054463959596453,
        0.795537299158248134863601, 0.02103233587872256311706,
        0.774029906950334246806958, 0.02191812889593413383869,
        0.751589869029638468178415, 0.02289743998716318463499,
        0.728244223887390363625800, 0.0237891645252872321010,
        0.704021101202391143555469, 0.02465221883590485293597,
        0.678949687946597146181269, 0.025485571944322848447,
        0.653060193216842191961262, 0.02628821747651458736160,
        0.626383811835045126762609, 0.02705918748154795852161,
        0.598952686760742185888769, 0.02779755327530227515804,
        0.570799870361220978705362, 0.02850242518416141631876,
        0.541959284585913426189305, 0.0291729589210074248656,
        0.512465680093027970988463, 0.0298083346403127548715,
        0.482354594377665692520121, 0.03040179231928695269039,
        0.451662308951869367576379, 0.03097061415408092094594,
        0.420425805628197756101396, 0.03149611781181863607696,
        0.388682721959498206776509, 0.03198367310021857603946,
        0.356471305888567846131189, 0.03243268955425561691179,
        0.323830369662345966511395, 0.03284262714400750457863,
        0.290799243066166651549602, 0.0332129992655131651404,
        0.257417726034420129920888, 0.03354333764112427668293,
        0.223726040694722859268958, 0.03383326624683168725793,
        0.189764782903379019020874, 0.03408242840225399546361,
        0.155574873330529119511405, 0.03429052388637504193170,
        0.121197508153924082968749, 0.03445730196032425617460,
        0.086674109420734770085237, 0.03458256166949689141805,
        0.052046275137206949059279, 0.03466615208568824018827,
        0.017355729146299652461298, 0.03470797248895005792046
};

long double integrate(function f){
    long double sum =0.0;
    for (int i = 0; i <45 ; ++i) {
        sum += Gn[2*i+1] * f((Gn[2*i]));
//        printf("%d , %llf\n",i, (Gn[2*i]));
    }
    return sum;

}


long double doubleIntegrate(function f){
    long double sum =0.0, a=0.0;

    for (int i = 0; i < 45; ++i) {
        sum += Gn[2 * i + 1] * f((Gn[2 * i]));
        a += Gn[2 * i + 1] * f((Gn[2 * i])) * Gn[2 * i];
    }

    return sum - a;

}


long double doubleIntegrate(function f, long double& slope){
    long double sum =0.0, a=0.0;

    for (int i = 0; i < 45; ++i) {
        sum += Gn[2 * i + 1] * f((Gn[2 * i]));
        a += Gn[2 * i + 1] * f((Gn[2 * i])) * Gn[2 * i];
    }
    slope=a;
    return sum - a;

}


/*

using namespace SingleLinkedList;
using namespace std;

bool operator>(Pitem& p1,Pitem p2){
    if(p1.power>p2.power){
        return true;
    }
    return false;
}

Polynomial::Polynomial() : roots(nullptr),degree(0){

}

Polynomial::Polynomial(int degree, const long double *coff):degree(degree) {
    Pitem pi{0,0};
    for (int i = 0; i < degree+1; ++i) {
        pi.value=coff[i];
        pi.power=i;
        this->push(pi);
    }

    this->sort();
}

long double Polynomial::operator()(long double x) {

    queueIterator<Pitem> iterator;
    iterator.start(this);
    long double result=0;

    while (iterator.end()){
        result+= iterator.read().value * pow(x,iterator.read().power);
        ++iterator;
    }

    return result;
}

Polynomial &Polynomial::operator+(Polynomial &p2) {
//    return ;
}

Polynomial &Polynomial::operator-(Polynomial &p2) {
//    return ;
}
*/