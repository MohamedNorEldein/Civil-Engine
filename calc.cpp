//
// Created by m7mdn on 2/5/2022.
//
#include "../../header/calc.h"

SingleLinkedList::queue<function> parallel_bracket;

void print(char * start, char* end){
    while (start<end){
        printf("%c", *start);
        start++;
    }
    printf("\n");
}

bool strCmpr(char* start,char *end,const char * text){
    while (start!=end and *text !='\0'){
        if(*start !=*text){
            return false;
        }
        start++;
        text++;
    }
    if(*text=='\0'){
        return true;
    }
    return true;
}


void strBracket(char* start,char *end){
    *start ='@';
    start++;
    while (start!=end){
        *start='_';
        start++;
    }
}

double func_x(double x){
    return x;
}


function stringObjectConvert(char *start, char *end) {

    if (strCmpr(start, end, "x")) {
        return function(func_x);
    }

    if (*start == '@') {
        return parallel_bracket.pop();
    }
    if (start < end) {
        double v = 0;
//        _print(start,end);
        sscanf(start, "%lf", &v);
        return function(v);
    }
    throw 0;
}


function quanBasePower(char *start, char *end) {
    char *reader = start;
//    _print(start,end);

    while (reader < end and *reader != '^') {
        reader++;
    }
//    printf(reader);
    if (reader == end) {
//        printf("jkxkvdkvjlkvdnall\n");
//        _print(start,reader);
        return stringObjectConvert(start, reader);
    }
    reader++;
//    printf("null\n");

    return (stringObjectConvert(start, reader)).pow(stringObjectConvert(reader,end));
}


function getQuantity(char *termStr, char *end) {
//    break_line
//    _print(termStr,end);

    function f(1);
    char *reader = termStr, *last_read = termStr;
    Sign state = pos;
// divide

    while (reader < end) {
        if (*reader == '*') {
//            _print(last_read,reader);

            if (state == pos) {
                f *= quanBasePower(last_read, reader);
                last_read = reader+1 ;
            }else{
                f /= quanBasePower(last_read, reader);
                last_read = reader+1 ;
            }
            state = pos;

        }
        if (*reader == '/') {
            if (state == pos) {
                f *= quanBasePower(last_read, reader);
                last_read = reader+1 ;
            }else{
                f /= quanBasePower(last_read, reader);
                last_read = reader+1 ;
            }
            state = neg;
        }

        reader++;
    }

    if (state == pos) {
//        _print(last_read,end);
//        printf("%lf\n",quanBasePower(last_read, end)(1));

        f *= quanBasePower(last_read, end);
//        printf("%lf\n",f(1));

    }else{
        f /= quanBasePower(last_read, end);
    }

//    printf("%lf\n",f(0));
//break_line
    return f;
}


function getTerms(char *eqn, char * end) {

    function f(0.0);
    char *reader = eqn, *last_read = eqn;
    Sign sign = pos;

    if (*reader == '-') {
        sign = neg;
        reader++;
        last_read++;
    }

    if (*reader == '+') {
        reader++;
        last_read++;
    }

    while (reader != end) {
        if (*reader == '+') {
//            printf("%lf\n",f(0));

            if (sign == pos) {
                f += getQuantity(last_read, reader);
                last_read = reader + 1;
            }else{
                f -= getQuantity(last_read, reader);
                last_read = reader + 1;
            }
            sign = pos;
        }
        if (*reader == '-') {
//            printf("%lf\n",f(0));

            if (sign == pos) {
                f += getQuantity(last_read, reader);
                last_read = reader + 1;
            }else{
                f -= getQuantity(last_read, reader);
                last_read = reader + 1;
            }
            sign = pos;
        }
        reader++;
    }
    if (sign == pos) {
        f += getQuantity(last_read, reader);

    }else{
        f -= getQuantity(last_read, reader);
    }

    return f;
}

function BuildStrct(char * start, char *end) {
    function sds = getTerms(start+1, end-1);
    strBracket(start,end);
    return sds;
}


function bracket(char * eqn){

    SingleLinkedList::stack<char*> char_bracket ;
    size_t ptrSize = sizeof(char* );
    char * reader = eqn;

    while(*reader!='\0')
    {
        if(*reader=='(')
        {
            char_bracket.push(reader);
        }
        if(*reader==')'){
            function bra = BuildStrct(char_bracket.pop(), reader+1);
            parallel_bracket.push(bra);

        }
        reader++;
    }
    return getTerms(eqn,reader);
}

function BuildEqn(std::string eqn){
    return bracket((char *)eqn.c_str());
}
