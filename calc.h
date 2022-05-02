
/* this program is a machine which takes in a text of an equation or a function and evaluate it
 * first the text is cut down to terms at each '+' or '-' character to represent terms   */

#include <cmath>
#include "function.h"

#ifndef MAIN_PY_CALC_H
#define MAIN_PY_CALC_H

bool strCmpr(char *start, char *end, const char *text);

enum Sign {
    pos , neg
};
/*
typedef struct Quantity {
    Sign state;
    function base;
    function power;

} quantity;

typedef SingleLinkedList::readerStack<quantity*> qstack;

typedef struct Term {

    Sign sign;
    qstack quantities;

} Term;

typedef SingleLinkedList::readerStack<Term*> tstack;
*/
function BuildEqn(std::string eqn);

void delete_eqn();

#endif //MAIN_PY_CALC_H
