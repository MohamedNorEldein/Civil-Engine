//
// Created by m7md_nor on 10/27/2021.
//

#include "../../header/civil_work.h"

namespace Sections {
    Type calculateTriangle(Type &Ix, Type &Iy, Type &Ixy, Type &Tx, Type &Ty, Type &area, Type data[6],
                           int area_state) {

        Tx = (data[0] + data[2] + data[4]) / 3;
        Ty = (data[1] + data[3] + data[5]) / 3;

        area = (data[2] - data[0]) * (data[5] - data[1]);
        area -= (data[4] - data[0]) * (data[3] - data[1]);

        Type a = data[2] - data[0];
        Type c = data[3] - data[1];
        Type b = data[4] - Tx;
        Type d = data[5] - Ty;
        Type r = a * d - c * b;

        Iy = (a * a / 32 + b * b * 3 / 32) * r;
        Ix = (c * c / 32 + d * d * 3 / 32) * r;
        Ixy = (a * c / 32 + d * b * 3 / 32) * r;

        return area * area_state / 2;
    }

    Type arc_calc(Shape &arc_shape) {
        /*
        double r_sqr = (arc_shape.data[2] - arc_shape.data[0])*(arc_shape.data[2] - arc_shape.data[0]);
        r_sqr += (arc_shape.data[3] - arc_shape.data[1])*(arc_shape.data[3] - arc_shape.data[1]);

        double r = sqrt(r_sqr);
    */
        Type r = arc_shape.data[2];
        Type theta = (arc_shape.data[3] + arc_shape.data[4]) / 2;
        Type alpha = (arc_shape.data[4] - arc_shape.data[3]);
        Type h = (4 * r * sin(alpha / 2) / (3 * M_PI));
//    printf("%lf, %f, %f",theta, arc_shape.data[0], arc_shape.data[1]);
        Type rx = h * cos(theta), ry = h * sin(theta);
        arc_shape.cg[0] = arc_shape.data[0] + rx;
        arc_shape.cg[1] = arc_shape.data[1] - ry;

        arc_shape.Ixy = 0;
// problem does exist
        Type r4 = r * r;
        r4 = r4 * r4;
        if (alpha > 0) {
            arc_shape.area = r * r * alpha / 2;
            arc_shape.Ix = (alpha + sin(alpha)) * r4 / 8 - ry * ry * arc_shape.area;
            arc_shape.Iy = (alpha - sin(alpha)) * r4 / 8 - rx * rx * arc_shape.area;

        } else {
            alpha = 2 * M_PI - alpha;
            arc_shape.area = r * r * alpha / 2;
            arc_shape.Ix = (alpha + sin(alpha)) * r4 / 8 - ry * ry * arc_shape.area;
            arc_shape.Iy = (alpha - sin(alpha)) * r4 / 8 - rx * rx * arc_shape.area;

        }
        return 0;
    }

    Type traverseProperties(Shape &traverse1) {
        Type cx = 0, cy = 0;

        Type TAs[traverse1.numOfPoints - 2];
        Type Txs[traverse1.numOfPoints - 2];
        Type Tys[traverse1.numOfPoints - 2];

        traverse1.Ix = 0;
        traverse1.Iy = 0;
        traverse1.Ixy = 0;
        traverse1.area = 0;

        for (int i = 2; i < traverse1.numOfPoints * 2 - 2; i += 2) {
            Type Ix, Iy, Ixy;

            TAs[i / 2 - 1] = calculateTriangle(Ix, Iy, Ixy, Txs[i / 2 - 1], Tys[i / 2 - 1], TAs[i / 2 - 1],
                                               new Type[6]{traverse1.data[0], traverse1.data[1],
                                                           traverse1.data[i], traverse1.data[i + 1],
                                                           traverse1.data[i + 2], traverse1.data[i + 3]},
                                               traverse1.area_state);
//                            std::cout<<"area = "<<TAs[i/2 -1]<<'\n';
            traverse1.area += TAs[i / 2 - 1];
            cx += Txs[i / 2 - 1] * TAs[i / 2 - 1];
            cy += Tys[i / 2 - 1] * TAs[i / 2 - 1];

            traverse1.Ix += Ix;
            traverse1.Iy += Iy;
            traverse1.Ixy += Ixy;

        }
        traverse1.area = abs(traverse1.area);

        if (traverse1.area != 0) {
            cx = cx / traverse1.area;
            cy = cy / traverse1.area;
        }
        traverse1.cg[0] = cx;
        traverse1.cg[1] = cy;

//    std::cout<<traverse1.Ix<<' '<<traverse1.Iy<<' '<<traverse1.Ixy<<'\n';


        for (int i = 0; i < traverse1.numOfPoints - 2; i++) {
            Type Rx = Txs[i] - cx;
            Type Ry = Tys[i] - cy;

            std::cout << Ry * Ry * TAs[i] << '\t' << Rx * Rx * TAs[i] << '\n';
//        std::cout<<traverse1.Ix<<'\t'<<traverse1.Iy<<'\t'<<traverse1.Ixy<<'\n';

            traverse1.Ix = traverse1.Ix + (Ry * Ry) * TAs[i];
            traverse1.Iy = traverse1.Iy + Rx * Rx * TAs[i];
            traverse1.Ixy = traverse1.Ixy + Rx * Ry * TAs[i];

        }

        return traverse1.area;
    }

    void DCShape(Shape &shape1, int id) {
//    printf("Shape type = %d\n", shape1.type);

        if (shape1.type == rectangle) {
            shape1.area = (shape1.data[3] - shape1.data[1]) * (shape1.data[2] - shape1.data[0]);
            shape1.cg[0] = (shape1.data[0] + shape1.data[2]) / 2;
            shape1.cg[1] = (shape1.data[3] + shape1.data[1]) / 2;
            shape1.area = abs(shape1.area) * shape1.area_state;

            //              y intercept                              x intercept
            shape1.Ix = abs(pow((shape1.data[3] - shape1.data[1]), 3) * (shape1.data[2] - shape1.data[0]) / 12);
            shape1.Iy = abs((shape1.data[3] - shape1.data[1]) * pow((shape1.data[2] - shape1.data[0]), 3) / 12);
            shape1.Ixy = 0;


            return;
        }

        if (shape1.type == circle) {
            Type rsq = pow((shape1.data[2] - shape1.data[0]), 2) + pow((shape1.data[3] - shape1.data[1]), 2);
            shape1.area = 3.14159265359 * rsq;
            shape1.area = abs(shape1.area) * shape1.area_state;

            shape1.cg[0] = shape1.data[0];
            shape1.cg[1] = shape1.data[1];

            shape1.Ix = rsq * rsq * M_PI / 4;
            shape1.Iy = rsq * rsq * M_PI / 4;
            shape1.Ixy = 0;

            return;
        }

        if (shape1.type == triangle) {
            calculateTriangle(shape1.Ix, shape1.Iy, shape1.Ixy, shape1.cg[0], shape1.cg[1], shape1.area, shape1.data,
                              shape1.area_state);
            shape1.area *= shape1.area_state;
            shape1.Ix *= shape1.area_state;
            shape1.Ixy *= shape1.area_state;
            shape1.Iy *= shape1.area_state;

            //        printf("area = %f  Ix = %f  Iy = %f  Ixy = %f ",shape1.area, shape1.Ix,shape1.Iy,shape1.Ixy );
            return;
        }

        if (shape1.type == traverse) {
            traverseProperties(shape1);
            return;
        }

        if (shape1.type == arc) {
            arc_calc(shape1);
            return;
        }

    }

    Type Section::sectionProperties() {

        Shape *shape1;
        SingleLinkedList::stackIterator<Shape*> iterator;

        iterator.start(shapeList);
        while (iterator.end()) {
            shape1 = iterator.read();
            DCShape(*shape1, 0);
            area += shape1->area;
            Gx += shape1->area * shape1->cg[0];
            Gy += shape1->area * shape1->cg[1];

            ++iterator;
        }

        if (area == 0) {
//        throw 0;
            Gx = 0;
            Gy = 0;
        } else {
            Gx /= area;
            Gy /= area;
        }

        Ix = 0;
        Iy = 0;
        Ixy = 0;

        iterator.start(shapeList);
        while (iterator.end()) {
            shape1 = iterator.read();
            shape1->Rx = shape1->cg[0] - Gx;
            shape1->Ry = shape1->cg[1] - Gy;

            Ix += shape1->Ix + shape1->Ry * shape1->Ry * shape1->area;
            Iy += shape1->Iy + shape1->Rx * shape1->Rx * shape1->area;
            Ixy += shape1->Ixy + shape1->Rx * shape1->Ry * shape1->area;

           ++iterator;
        }

        Type Iav = (Ix + Iy) / 2;

        if (Iav == Ix) {
            Imin = Ix;
            Imax = Iy;
            return area;
        }

        theta = abs(atan(Ixy / (Ix - Iav)));

        if (theta == 0) {
            Imax = Ix;
            Imin = Iy;
            return area;
        }


        if (Ixy > 0) {
            if (Ix < Iav) {
                theta += M_PI / 4;
            }
        } else {
            if (Ix < Iav) {
                theta += 3 * M_PI / 4;
            } else {
                theta += M_PI / 2;
            }
        }


        Type R = Ixy / sin(theta);

        Imax = Iav + R;
        Imin = Iav - R;

        return area;
    }

    Type Section::calcStressProfile() {

        if (Ix == 0.0 or Iy == 0.0 or area == 0.0) {
            throw std::exception();
        }

        Type ns1 = moment[1] / force[2], ns2 = moment[0] / force[2], a, b, d;
        a = +ns1 * Ix - ns2 * Ixy;
        b = -ns1 * Ixy + ns2 * Iy;
        d = Ix * Iy - Ixy * Ixy;
        a /= d;
        b /= d;

        Type s = 0;
        smax = s;
        smin = s;
        sav = force[2] / area;

        SingleLinkedList::stackIterator<Shape*> iterator;

        iterator.start(shapeList);
        int i=0;
        while (iterator.end()) {

            Type *points = iterator.read()->data;

            for (int j = 0; j < iterator.read()->numOfPoints; ++j) {
                s = sav + a * points[2 * i] + b * points[2 * i + 1];
                if (s > smax) {
                    smax = s;
                }
                if (s < smin) {
                    smin = s;
                }
            }
            ++iterator;
            ++i;
        }
        return 1;
    }

    Section::Section(Type *force, Type *moment) : force(force), moment(moment) {

        area = 0;
        Ix = 0;
        Iy = 0;
        Ixy = 0;
        Imax = 0;
        Imin = 0;
        theta = 0;
        Gx = 0;
        Gy = 0;

        TensionCriticalStress = 0;
        CompressionCriticalStress = 0;

        smax = 0, smin = 0, sav = 0;

        shapeList = new ShapeList();
    }

    void Section::addShape(Shape &shape) const {
        shapeList->push(&shape);
    }

    void Section::addShape(Shape *shape) const {
        shapeList->push(shape);
    }

    Section::Section() : force((Type *) malloc(3 * sizeof(Type))), moment((Type *) malloc(3 * sizeof(Type))) {

        area = 0;
        Ix = 0;
        Iy = 0;
        Ixy = 0;
        Imax = 0;
        Imin = 0;
        theta = 0;
        Gx = 0;
        Gy = 0;

        TensionCriticalStress = 0;
        CompressionCriticalStress = 0;

        smax = 0, smin = 0, sav = 0;
    }

    Section::~Section() {
        Shape *shape1;
        while (shapeList->length() > 0) {
            shape1 = shapeList->pop();
            delete shape1->data;
            delete shape1;
        }
    }

    Shape *Section::addRectangle(Type startx, Type starty, Type endx, Type endy, int areaState) {
        auto shape = new Shape{rectangle, areaState, 2, new Type[4]{startx, starty, endx, endy}};
        addShape(shape);
        return shape;
    }

    Shape *Section::addTraverse(Type *data, int numofpoints, int areastate) {
        auto shape = new Shape{traverse, areastate, numofpoints, data};
        shapeList->push(shape);
        return shape;
    }


}



namespace StructuralSystem::TrussSystem {
    void direction(Joint *startPoint, Joint *endPoint, float *Normal) {
        float length = 0, l;
        for (int i = 0; i < systemDegree; ++i) {
            l = (startPoint->position[i] - endPoint->position[i]);
            l *= l;
            length += l;
        }
        length = sqrt(length);

        for (int i = 0; i < systemDegree; ++i) {
            Normal[i] = (endPoint->position[i] - startPoint->position[i]) / length;
        }

    }

    Joint::Joint(float position_x, float position_y, float position_z, int index, supportType support) : numberOfMembers(0),
                                                                                                                  support(support), index(index) {
        position[0] = position_x;
        position[1] = position_y;

        if (systemDegree == 3) {
            position[2] = position_z;
        }
        externalForce[0] = 0;
        externalForce[0] = 0;
        externalForce[0] = 0;

        iteration = 0;

    }

    Joint::~Joint() {
        printf("deleting joint at %x ...\n", this);
    }

    void Joint::setExternalForce(float px, float py, float pz) {
        externalForce[0] += px;
        externalForce[1] += py;
        if (systemDegree == 3) {
            externalForce[2] += pz;
        }
    }

    void Joint::print() {
        if (systemDegree == 3) {
            printf("joint (%f, %f, %f)........\n", position[0], position[1], position[2]);
        } else {
            printf("joint (%f, %f) jointType type %d node deformation(%e,%e) .......\n", position[0], position[1], support,deformation[0],deformation[1]);
        }
        for (int i = 0; i < numberOfMembers; ++i) {
            float _direction[3];
            direction(this, membersEnd[i], _direction);
            printf("\tend joint (%f, %f) number: %d ---> force = %f \n", membersEnd[i]->position[0],
                   membersEnd[i]->position[1], members[i]->i, members[i]->normalForce);
        }

    }


    void Joint::addMember(trussMember *member1, Joint *joint) {
        members[numberOfMembers] = member1;
        membersEnd[numberOfMembers] = joint;
        numberOfMembers++;
    }

    void Joint::changePoint(float position_x, float position_y, float position_z) {
        position[0] = position_x;
        position[1] = position_y;
        if (systemDegree == 3) {
            position[2] = position_z;
        }
    }

    float &support::getReaction(int i){
        if (reaction == nullptr) {
            throw int(0);
        }
        return reaction[i];
    }

    trussMember::trussMember(int i,Joint *start, Joint*end,float A, float E, float l) : i(i), Area(A), E(E),length(l),start(start),end(end) {
        normalForce = 0;
    }

    float dirDot(float *Normal, float *other) {
        float re = 0;
        for (int i = 0; i < systemDegree; ++i) {
            re += Normal[i] * other[i];
        }
        return re;
    }

    Truss::Truss() : startJoint(nullptr), numOfMembers(0), numOfJoints(0) {
        htmlfile = fopen("file.html","w");
        fprintf(htmlfile,"<html>\n"
                         "<body><h1>canvas</h1>\n"
                         "\n"
                         "\n"
                         "<canvas id=\"0\" width=\"1000\"height=\"500\" style=\"border:1px solid #d3d3d3;\">\n"
                         " </canvas>");
    }

    void Joint::addEqnToMatrix(short_vector::Matrix &matrix, short_vector::Vector& vector,  int Dimn) {
        float member_j[3];
        Joint* joint= this;
        int i = 2*joint->index;
            for (int j = 0; j < joint->numberOfMembers; ++j) {
                direction(joint, joint->membersEnd[j], member_j);

                matrix.set_item(i, joint->members[j]->i, member_j[0]);
                matrix.set_item(i + 1, joint->members[j]->i, member_j[1]);
            }
    }
    void Joint::newLoadSystem(short_vector::Matrix &matrix, short_vector::Vector& vector, int Dimn) {
        float member_j[3];

        vector.set_item(2*index, -externalForce[0]);
        vector.set_item(2*index + 1, -externalForce[1]);

    }

    void Rocker::addEqnToMatrix(short_vector::Matrix &matrix, short_vector::Vector &vector, int Dimn) {
        float member_j[3];
        int i = 2 * index;

        for (int j = 0; j < numberOfMembers; ++j) {
            direction(this, membersEnd[j], member_j);

            matrix.set_item(i, members[j]->i, member_j[0]);
            matrix.set_item(i + 1, members[j]->i, member_j[1]);
        }
        matrix.set_item(i, Dimn - 3, 1);
        matrix.set_item(i + 1, Dimn - 2, 1);

    }

    Rocker::Rocker(float position_x, float position_y, float position_z, int index) :
    support(position_x,position_y,position_z,index,rocker){

    }
    Roller::Roller(float position_x, float position_y, float position_z,float dx, float dy, float dz, int index) :
    support(position_x,position_y,position_z,index,roller){
        directionV[0] =dx;
        directionV[1] =dy;
//        directionV[2] =dz;

    }

    void Roller::addEqnToMatrix(short_vector::Matrix &matrix, short_vector::Vector &vector, int Dimn) {
        printf("roller %f %f \n",directionV[0],directionV[1]);
        float member_j[3];
        int i = 2*index;
        for (int j = 0; j < numberOfMembers; ++j) {
            direction(this, membersEnd[j], member_j);

            matrix.set_item(i, members[j]->i, member_j[0]);
            matrix.set_item(i + 1, members[j]->i, member_j[1]);
        }
        matrix.set_item(i, Dimn - 1, directionV[0]);
        matrix.set_item(i + 1, Dimn - 1, directionV[1]);
    }

    void Truss::iterationSolve() {

    }

    void Truss::MatrixAnalysis() {
        delete vector;
        delete matrix;
        delete DeflectionMatrix;
        delete deformation;

        matrix = new short_vector::Matrix(numOfJoints * 2 , 2 * numOfJoints);
        EA = new short_vector::Matrix(numOfJoints * 2 , 2 * numOfJoints);
        DeflectionMatrix = new short_vector::Matrix(numOfJoints * 2 , 2 * numOfJoints);
        vector = new short_vector::Vector(numOfJoints*2);

       SingleLinkedList::stackIterator<Joint *> iterator0;
       iterator0.start(allJoints);

       while (iterator0.end()) {
           iterator0.read()->addEqnToMatrix(*matrix, *vector, numOfJoints * 2);
           ++iterator0;
       }
       matrix->print();

       short_vector::Matrix* a=short_vector::inverse(matrix);
        delete matrix;
        matrix=a;
       matrix->print();

        trussMember *member;
        SingleLinkedList::stackIterator<trussMember*> iterator;
        iterator.start(allMembers);
        for (int i = 0; i < numOfMembers; i++) {
            member = iterator.read();
            EA->set_item(iterator.read()->i,iterator.read()->i,member->length  / (member->Area * member->E));
            ++iterator;
        }
        *DeflectionMatrix = ( matrix->transpose() * *EA) * *matrix;
        short_vector::clear_memory();
    }


    void Truss::MatrixSolve() {

        SingleLinkedList::stackIterator<Joint *> iterator0;
        iterator0.start(allJoints);

        while (iterator0.end()) {
            iterator0.read()->newLoadSystem(*matrix, *vector, numOfJoints * 2);
            ++iterator0;
        }
        short_vector::Vector *result= short_vector::mull(matrix,vector);

        deformation = short_vector::mull(DeflectionMatrix,vector);

        trussMember *member;
        SingleLinkedList::stackIterator<trussMember*> iterator;
        iterator.start(allMembers);
        for (int i = 0; i < numOfMembers; i++) {
            member = iterator.read();
            member->normalForce = result->get_item(member->i);
            ++iterator;
        }
        SingleLinkedList::stackIterator<Joint*> iterator2;
        iterator2.start(allJoints);
        Joint* j;
        for (int i = 0; i < numOfJoints; i++) {
            j = iterator2.read();
            j->deformation[0] = deformation->get_item(2*j->index);
            j->deformation[1] = deformation->get_item(2*j->index+1);
            ++iterator2;
        }
        short_vector::clear_memory();
        delete result;
    }

    void Truss::print() {

        SingleLinkedList::stackIterator<Joint *> iterator;
        iterator.start(allJoints);

        while (iterator.end()) {
            iterator.read()->print();
            ++iterator;
        }
    }

    Joint *Truss::addJoint(float position_x, float position_y, float position_z) {
        numOfJoints++;

        allJoints.push(new Joint(position_x, position_y, position_z, numOfJoints - 1));
        return allJoints.top();
    }

    Joint *Truss::addRocker(float position_x, float position_y, float position_z) {
        numOfJoints++;
        allJoints.push(new Rocker(position_x, position_y, position_z, numOfJoints - 1));
        return allJoints.top();
    }

    Joint *Truss::addRoller(float position_x, float position_y, float position_z, float dx, float dy, float dz) {
        numOfJoints++;
        allJoints.push(new Roller(position_x, position_y, position_z,dx,dy,dz, numOfJoints - 1));
        return allJoints.top();
    }


    trussMember *Truss::addMember(Joint *start, Joint *end, float A, float E) {

        float length = 0, l;
        for (int i = 0; i < systemDegree; ++i) {
            l = (start->position[i] - end->position[i]);
            l *= l;
            length += l;
        }
        length = sqrt(length);

//        printf("%f, %f, %f\n",length, A,E);
        allMembers.push(new trussMember(numOfMembers,start,end,A,E,length));
        numOfMembers++;
        start->addMember(allMembers.top(), end);
        end->addMember(allMembers.top(), start);
        return allMembers.top();
    }

    void Truss::setStartPoint(Joint *joint) {
        startJoint = joint;
    }


    Truss::~Truss() {
        printf("deleting System at %x ...\n", this);
        SingleLinkedList::stackIterator<Joint *> iterator;
        iterator.start(allJoints);
        while (iterator.end()){
            delete iterator.read();
            ++iterator;
        }
        SingleLinkedList::stackIterator<trussMember *> Titerator;
        Titerator.start(allMembers);
        while (Titerator.end()){
            delete Titerator.read();
            ++Titerator;
        }

        delete vector;
        delete matrix;
        delete DeflectionMatrix;
        delete deformation;

    }




}
