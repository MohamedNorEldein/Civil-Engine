//
// Created by m7mdn on 3/15/2022.
//

#include "../../header/civil_work.h"

double drawFactor=10, factor=100, startF=50;

namespace StructuralSystem::BeamSystem {



        sparse_system::AgMatrix *agMatrix;
//    short_vector::Vector *vec;
        int counter, length, matDegree;

        Type E = 1;


        void shearCondition(Joint *node) {
            stackIterator<Member *> stackIterator;
            stackIterator.start(node->memberList);
            Type shear, moment, slope, deflection;
            Type shearc = 0;

            while (stackIterator.end()) {

                stackIterator.read()->loadActions(node, shear, moment, slope, deflection);
                if (stackIterator.read()->end == node) {
                    shearc += shear;
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id, 1);     // Shear Shear
                } else {
                    shearc -= shear;
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id, -1);     // Shear Shear
                }
                ++stackIterator;
            }
            agMatrix->set_item(matDegree, length - 1, shearc);

            printf("%d id = %d\n", matDegree, node->memberList.top()->id);
            matDegree++;
        }

        void momentCondition(Joint *node) {
            stackIterator<Member *> stackIterator;
            stackIterator.start(node->memberList);
            Type shear, moment, slope, deflection;
            Type momentc = 0;

            while (stackIterator.end()) {
                stackIterator.read()->loadActions(node, shear, moment, slope, deflection);

                if (stackIterator.read()->end == node) {
                    momentc += moment;

                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id,
                                       stackIterator.read()->getLength(node));// moment shear constant
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 1, 1);     // moment

                } else {
                    momentc -= moment;
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id,
                                       -stackIterator.read()->getLength(node));// moment shear constant
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 1, -1);     // moment

                }

                ++stackIterator;
            }

            agMatrix->set_item(matDegree, length - 1, -momentc);

            matDegree++;
        }

        void singleMomentCondition(Joint *node) {
            stackIterator<Member *> stackIterator;
            stackIterator.start(node->memberList);
            Type shear, moment, slope, deflection;
            Type me=0;
            while (stackIterator.end()) {
                stackIterator.read()->loadActions(node, shear, moment, slope, deflection);

                if (stackIterator.read()->end == node) {
                    me+=moment;

                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id,
                                       stackIterator.read()->getLength(node));// moment shear constant
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 1, 1);     // moment

                } else {
                    me+=moment;

                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id,
                                       -stackIterator.read()->getLength(node));// moment shear constant
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 1, -1);     // moment

                }
                agMatrix->set_item(matDegree, length - 1, -me);

                matDegree++;
                ++stackIterator;
            }

        }

        void singleMemberSlope(Joint *node) {
            stackIterator<Member *> stackIterator;
            stackIterator.start(node->memberList);
            Type shear, moment, slope, deflection;

            while (stackIterator.end()) {
                stackIterator.read()->loadActions(node, shear, moment, slope, deflection);

                if (stackIterator.read()->end == node) {

                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id,
                                       pow(stackIterator.read()->getLength(node), 2) / 2 /
                                       (stackIterator.read()->I * E));     // slope  deflection
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 1,
                                       stackIterator.read()->getLength(node) /
                                       (stackIterator.read()->I * E));     // slope deflection
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 2,
                                       1 / (stackIterator.read()->I * E));     // slope deflection

                } else {

                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id,
                                       -pow(stackIterator.read()->getLength(node), 2) / 2 /
                                       (stackIterator.read()->I * E));     // slope  deflection
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 1,
                                       -stackIterator.read()->getLength(node) /
                                       (stackIterator.read()->I * E));     // slope deflection
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 2,
                                       -1 / (stackIterator.read()->I * E));     // slope deflection

                }
                agMatrix->set_item(matDegree + 0, length - 1, 0);
                matDegree++;
                ++stackIterator;
            }

        }

        void slopeCondition(Joint *node) {
            stackIterator<Member *> stackIterator;
            stackIterator.start(node->memberList);
            if (node->memberList.length() == 1) {
                return;
            }
            Type shear, moment, slope, deflection;
            Type slopec = 0;

            while (stackIterator.end()) {
                stackIterator.read()->loadActions(node, shear, moment, slope, deflection);

                if (stackIterator.read()->end == node) {
                    slopec += slope;

                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id,
                                       pow(stackIterator.read()->getLength(node), 2) / 2 /
                                       (stackIterator.read()->I * E));     // slope  deflection
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 1,
                                       stackIterator.read()->getLength(node) /
                                       (stackIterator.read()->I * E));     // slope deflection
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 2,
                                       1 / (stackIterator.read()->I * E));     // slope deflection

                } else {
                    slopec -= slope;

                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id,
                                       -pow(stackIterator.read()->getLength(node), 2) / 2 /
                                       (stackIterator.read()->I * E));     // slope  deflection
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 1,
                                       -stackIterator.read()->getLength(node) /
                                       (stackIterator.read()->I * E));     // slope deflection
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 2,
                                       -1 / (stackIterator.read()->I * E));     // slope deflection

                }

                ++stackIterator;
            }

            agMatrix->set_item(matDegree + 0, length - 1, -slopec);

            matDegree++;
        }

        void deflectionCondition(Joint *node) {
            stackIterator<Member *> stackIterator;
            stackIterator.start(node->memberList);
            Type shear, moment, slope, deflection;
            Type deflectioc = 0;

            while (stackIterator.end()) {
                stackIterator.read()->loadActions(node, shear, moment, slope, deflection);

                if (stackIterator.read()->end == node) {
                    deflectioc += deflection;

                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id,
                                       pow(stackIterator.read()->getLength(node), 3) / 6 /
                                       (stackIterator.read()->I * E));     // deflection
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 1,
                                       pow(stackIterator.read()->getLength(node), 2) / 2 /
                                       (stackIterator.read()->I * E));     // deflection
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 2,
                                       stackIterator.read()->getLength(node) /
                                       (stackIterator.read()->I * E));     // deflection
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 3,
                                       1 / (stackIterator.read()->I * E));     // deflection


                } else {
                    deflectioc -= deflection;

                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id,
                                       -pow(stackIterator.read()->getLength(node), 3) / 6 /
                                       (stackIterator.read()->I * E));     // deflection
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 1,
                                       -pow(stackIterator.read()->getLength(node), 2) / 2 /
                                       (stackIterator.read()->I * E));     // deflection
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 2,
                                       -stackIterator.read()->getLength(node) /
                                       (stackIterator.read()->I * E));     // deflection
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 3,
                                       -1 / (stackIterator.read()->I * E));     // deflection

                }

                ++stackIterator;
            }

            agMatrix->set_item(matDegree + 0, length - 1, -0);

            matDegree++;
        }

        void singleMemberDeflectionCondition(Joint *node) {
            stackIterator<Member *> stackIterator;
            stackIterator.start(node->memberList);
            Type shear, moment, slope, deflection;
            Type deflectioc = 0;

            while (stackIterator.end()) {

                stackIterator.read()->loadActions(node, shear, moment, slope, deflection);

                if (stackIterator.read()->end == node) {
                    deflectioc += deflection;

                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id,
                                       pow(stackIterator.read()->getLength(node), 3) / 6 /
                                       (stackIterator.read()->I * E));     // deflection
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 1,
                                       pow(stackIterator.read()->getLength(node), 2) / 2 /
                                       (stackIterator.read()->I * E));     // deflection
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 2,
                                       stackIterator.read()->getLength(node) /
                                       (stackIterator.read()->I * E));     // deflection
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 3,
                                       1 / (stackIterator.read()->I * E));     // deflection


                } else {
                    deflectioc -= deflection;

                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id,
                                       -pow(stackIterator.read()->getLength(node), 3) / 6 /
                                       (stackIterator.read()->I * E));     // deflection
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 1,
                                       -pow(stackIterator.read()->getLength(node), 2) / 2 /
                                       (stackIterator.read()->I * E));     // deflection
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 2,
                                       -stackIterator.read()->getLength(node) /
                                       (stackIterator.read()->I * E));     // deflection
                    agMatrix->set_item(matDegree, 4 * stackIterator.read()->id + 3,
                                       -1 / (stackIterator.read()->I * E));     // deflection

                }
                agMatrix->set_item(matDegree + 0, length - 1, -deflectioc);
                matDegree++;
                ++stackIterator;
            }
        }


        Force::Force(Type xdir, Type ydir, Type zdir, Type actionPoint, Type magnitude, forceDegree degree)
                : actionPoint(actionPoint), magnitude(magnitude), degree(degree) {
            direction[0] = xdir;
            direction[1] = ydir;
            direction[2] = zdir;

        }

        Type Force::operator()(Type position, int direction_) const {
            if (degree == -1) {
                if (position == actionPoint) {
                    return magnitude * direction[direction_];
                }
                return 0;
            }

            if (position > actionPoint) {
                return direction[direction_] * magnitude * pow(position - actionPoint, degree);
            }
            return 0;
        }

        Type Force::shearAction(Type position, int direction_) {
            if (position >= actionPoint) {
                if (degree == -1) {
                    return direction[direction_] * magnitude;
                }
                return direction[direction_] * magnitude * pow(position - actionPoint, degree + 1) / (degree + 1);
            }
            return 0.0;
        }

        Type Force::momentAction(Type position, int direction_) {
            if (position >= actionPoint) {
                if (degree == -1) {
                    return direction[direction_] * magnitude * pow(position - actionPoint, degree + 2) / (degree + 2);
                }
                return direction[direction_] * magnitude * pow(position - actionPoint, degree + 2) /
                       ((degree + 1) * (degree + 2));
            }
            return 0;

        }

        Type Force::slopeAction(Type position, int direction_) {
            if (position >= actionPoint) {
                if (degree == -1) {
                    return direction[direction_] * magnitude * pow(position - actionPoint, degree + 3) /
                           ((degree + 2) * (degree + 3));
                }
                return direction[direction_] * magnitude * pow(position - actionPoint, degree + 3) /
                       ((degree + 1) * (degree + 2) * (degree + 3));
            }
            return 0;
        }

        Type Force::deflectionAction(Type position, int direction_) {
            if (position >= actionPoint) {
                if (degree == -1) {
                    return direction[direction_] * magnitude * pow(position - actionPoint, degree + 4) /
                           ((degree + 2) * (degree + 3) * (degree + 4));
                }
                return direction[direction_] * magnitude * pow(position - actionPoint, degree + 4) /
                       ((degree + 1) * (degree + 2) * (degree + 3) * (degree + 4));
            }
            return 0;
        }

        Type Force::getMagnitude() const {
            return magnitude;
        }

        Joint::Joint(Type px, Type py) {
            position[0] = px;
            position[1] = py;

        }

        Roller::Roller(Type px1, Type py1) : Support(px1, py1) {

        }

        void Roller::calculateConstants() {
            if (memberList.top()->start == this) {
                memberList.top()->momentConstant = 0.0;
            } else {
                memberList.top()->momentConstant = -memberList.top()->moment(memberList.top()->length);
            }

        }

        Rocker::Rocker(Type px, Type py) : Support(px, py) {

        }

        void Rocker::calculateConstants() {
//        shearCondition(this);
            momentCondition(this);
            slopeCondition(this);
            singleMemberDeflectionCondition(this);
        }

        Support::Support(Type px1, Type py1) : Joint(px1, py1) {

        }

        FixedSupport::FixedSupport(Type px, Type py) : Support(px, py) {

        }

        void FixedSupport::calculateConstants() {
            singleMemberSlope(this);
            singleMemberDeflectionCondition(this);
        }

        Member::Member(Joint *start, Joint *end)
                : start(start), end(end), I(1) {
            double a = 0;
            length = 0;
            for (int i = 0; i < systemDegree; ++i) {
                a = end->position[i] - start->position[i];
                length += a * a;
                direction[i] = a;
            }
//        printf("length = %f\n",length);

            length = sqrt(length);
            for (int i = 0; i < systemDegree; ++i) {
                direction[i] / length;
            }
            start->memberList.push(this);
            end->memberList.push(this);
        }

        void Member::addForce(Type xdir, Type ydir, Type zdir, Type actionPoint, Type magnitude, forceDegree degree) {
            ForceList.push(new Force(xdir, ydir, zdir, actionPoint, magnitude, degree));
        }

        Type Member::getLength() const {
            return length;
        }

        Type Member::shear(float x) {
            stackIterator<Force *> iterator;
            iterator.start(ForceList);
            Type shearX = 0;
            while (iterator.end()) {
                shearX += iterator.read()->shearAction(x, 1);
                ++iterator;
            }
            return -(shearX + shearConstant);
        }

        Type Member::moment(float x) {
            stackIterator<Force *> iterator;
            iterator.start(ForceList);
            Type X = 0;
            while (iterator.end()) {
                X += iterator.read()->momentAction(x, 1);
                ++iterator;

            }

            stackIterator<Moment *> iterator1;

            iterator1.start(MomentList);
            while (iterator1.end()) {
                X += iterator1.read()->momentAction(x, 2);
                ++iterator1;
            }

            return X + shearConstant * x + momentConstant;
        }

        Type Member::slope(float x) {
            stackIterator<Force *> iterator;
            iterator.start(ForceList);
            Type X = 0;
            while (iterator.end()) {
                X += iterator.read()->slopeAction(x, 1);
                ++iterator;

            }

            stackIterator<Moment *> iterator1;

            iterator1.start(MomentList);
            while (iterator1.end()) {
                X += iterator1.read()->slopeAction(x, 2);
                ++iterator1;
            }
            return (X + shearConstant * x * x / 2 + momentConstant * x + slopeConstant) / (E * I);
        }

        Type Member::deflection(float x) {
            stackIterator<Force *> iterator;
            iterator.start(ForceList);
            Type X = 0;
            while (iterator.end()) {
                X += iterator.read()->deflectionAction(x, 1);
                ++iterator;

            }
            stackIterator<Moment *> iterator1;
            iterator1.start(MomentList);
            while (iterator1.end()) {
                X += iterator1.read()->deflectionAction(x, 2);
                ++iterator1;
            }
            return (X + shearConstant * x * x * x / 6 + momentConstant * x * x / 2 + slopeConstant * x +
                    deflectionConstant) / (E * I);
        }

        void Member::actions(Type x, Type &shear, Type &moment, Type &slope, Type &deflection) {
            shear = shearConstant;
            moment = shearConstant * x + momentConstant;
            slope = shearConstant * x * x / 2 + momentConstant * x + slopeConstant;
            deflection =
                    shearConstant * x * x * x / 6 + momentConstant * x * x / 2 + slopeConstant * x + deflectionConstant;

            stackIterator<Force *> iterator;
            iterator.start(ForceList);

            while (iterator.end()) {
                shear += iterator.read()->shearAction(x, 1);
                moment += iterator.read()->momentAction(x, 1);
                slope += iterator.read()->slopeAction(x, 1);
                deflection += iterator.read()->deflectionAction(x, 1);
                ++iterator;
            }

            stackIterator<Moment *> iterator1;
            iterator1.start(MomentList);

            while (iterator1.end()) {
                moment += iterator1.read()->momentAction(x, 2);
                slope += iterator1.read()->slopeAction(x, 2);
                deflection += iterator1.read()->deflectionAction(x, 2);
                ++iterator1;
            }

            shear *= -1;
            slope /= E * I;
            deflection /= E * I;
        }

        int Member::getId() const {
            return id;
        }

        Joint *Member::getStart() const {
            return start;
        }

        Joint *Member::getAnEnd() const {
            return end;
        }

        void Member::setI(Type value) {
            I = value;
        }

/*
    void Member::loadActions(Type x, Type &shear, Type &moment, Type &slope, Type &deflection) {
        stackIterator<Force*> iterator;
        iterator.start(ForceList);
        shear = 0;
        moment = 0;
        slope = 0;
        deflection =0;
        while (iterator.end()) {
            shear += iterator.read()->shearAction(x, 1);
            moment += iterator.read()->momentAction(x, 1);
            slope += iterator.read()->slopeAction(x, 1);
            deflection += iterator.read()->deflectionAction(x, 1);
            ++iterator;
        }

        stackIterator<Moment*> iterator1;
        iterator1.start(MomentList);

        while (iterator1.end()) {
            moment += iterator1.read()->momentAction(x, 2);
            slope += iterator1.read()->slopeAction(x, 2);
            deflection += iterator1.read()->deflectionAction(x, 2);
            ++iterator1;
        }

        shear*=-1;
        slope/=E * I;
        deflection/=E * I;
    }
*/
        void Member::loadActions(Joint *J1, Type &shear, Type &moment, Type &slope, Type &deflection) {
            Type x;
            if (J1 == start) {
                x = 0;
            } else {
                x = this->length;
            }
//        printf("x= %lf\n", x);
            stackIterator<Force *> iterator;
            iterator.start(ForceList);
            shear = 0;
            moment = 0;
            slope = 0;
            deflection = 0;
            while (iterator.end()) {
//            printf("force %f\n", iterator.read()->shearAction(x,1));
                shear += iterator.read()->shearAction(x, 1);
                moment += iterator.read()->momentAction(x, 1);
                slope += iterator.read()->slopeAction(x, 1);
                deflection += iterator.read()->deflectionAction(x, 1);
                ++iterator;
            }

            stackIterator<Moment *> iterator1;
            iterator1.start(MomentList);

            while (iterator1.end()) {
                moment += iterator1.read()->momentAction(x, 2);
                slope += iterator1.read()->slopeAction(x, 2);
                deflection += iterator1.read()->deflectionAction(x, 2);
                ++iterator1;
            }

            shear *= -1;
            slope /= E * I;
            deflection /= E * I;
        }

        Type Member::getLength(Joint *J) {
            if (J == start) {
                return 0;
            }
            return length;
        }

        void Member::addLoad(Type actionPoint, Type magnitude, forceDegree degree) {
            if(std::abs(actionPoint - length) < 0.0001){
                actionPoint -= 0.00001f;
            }
            if(std::abs(actionPoint) < 0.0001){
                actionPoint += 0.00001f;
            }
            addForce(0, -1, 0, actionPoint, magnitude, degree);
        }

        void Member::addMoment(Type actionPoint, Type magnitude) {
            MomentList.push(new Moment(actionPoint, magnitude));
        }

        void Member::print() {
            printf("start  ");
            start->print();
            printf("end  ");
            end->print();

        }

        void Member::printDeflection() {
            Type shear, moment, slope, deflection;

            for (int i = 0; i < length; ++i) {
                actions(i, shear, moment, slope, deflection);
                printf("%d:->shear = %lf moment = %lf slope = %lf deflection %lf\n", i, shear, moment, slope,
                       deflection);
            }
            actions(length, shear, moment, slope, deflection);
            printf("%f:->shear = %lf moment = %lf slope = %lf deflection %lf\n", length, shear, moment, slope,
                   deflection);
        }


    void Member::Draw(FILE* file) {
        fprintf(file, "<script>"

                      "var c = document.getElementById(\"0\");\n"
                      "var ctx = c.getContext(\"2d\");\n"
                      "var cty = c.getContext(\"2d\");\n"

                      "cty.strokeStyle = \"blue\";\n"
                      "cty.lineWidth = 0.5;"
                      "cty.beginPath();\n"
                      "cty.moveTo(%f,%f);\n"
                      "cty.lineTo(%f, %f);\n"

                      "cty.lineTo(%f, %f);\n"
                      "cty.stroke();"
                      "ctx.strokeStyle = \"red\";"
                      "ctx.lineWidth = 1;"
                      "ctx.beginPath()\n"
                      "    ctx.moveTo(%f, %f);\n",
                start->position[0] * factor + startF, 100.0,
                end->position[0] * factor + startF, 100.0,
                end->position[0] * factor + startF, 100.0,
                start->position[0] * factor + startF,
                100.0 - drawFactor * deflection(0.0));

        int n=50;
        float step = length/n;

        for (int i = 0; i < n; i+=3) {

            fprintf(file,
                    " ctx.bezierCurveTo(%f, %f, %f, %f, %f, %f);\n",
//                         "ctx.lineTo(%f, %f)\n"
//                         "ctx.lineTo(%f, %f)\n"
//                         "ctx.lineTo(%f, %f)\n",
                         start->position[0] * factor + startF + i * step * factor, 100 - drawFactor * deflection(i * step)
            , start->position[0] * factor + startF + (i + 1) * step * factor, 100 - drawFactor * deflection((i + 1) * step),
                    start->position[0] * factor + startF + (i + 2) * step * factor, 100 - drawFactor * deflection((i + 2) * step));
//            printf("%f",(i) * step);
        }
        fprintf(file,"    ctx.lineTo(%f, %f);\n"
                     "    ctx.stroke();\n"
                     "</script>"
                     "\n", start->position[0] * factor + startF + length * factor, 100 - drawFactor * deflection(length));

    }

    void Joint::calculateConstants() {
        }

        void Joint::print() {
            printf("position = (%f,%f) ", position[0], position[1]);
            printf("\n");
        }

        void FreeEnd::calculateConstants() {
            shearCondition(this);
            singleMomentCondition(this);
        }

        FreeEnd::FreeEnd(Type px, Type py) : Joint(px, py) {

        }

        System::System() {
            htmlfile = fopen("C:\\Users\\m7mdn\\Desktop\\structure\\newfile.html","w");
            fprintf(htmlfile,"<html>\n"
                    "<body><h1>canvas</h1>\n"
                    "\n"
                    "<body>\n"
                    "\n"
                    "<canvas id=\"0\" width=\"1200\"height=\"400\" style=\"border:1px solid #d3d3d3;\">\n"
                    " </canvas>\n");
        }

        Joint *System::addJoint(Type px, Type py, supportType support) {
            Joint *j1;
            switch (support) {
                case rocker:
                    j1 = new Rocker(px, py);
                    break;
                case roller:
                    j1 = new Roller(px, py);
                    break;
                case None:
                    j1 = 0;
                    break;

                case pin:
                    j1 = 0;
                    break;
                case fixedSupport:
                    j1 = new FixedSupport(px, py);
                    break;
                case freeEnd:
                    j1 = new FreeEnd(px, py);
                    break;
                case pinJoint:
                    j1 = new PinJoint(px, py);
                    break;
                case fixedJoint:
                    j1 = new FixedConnection(px, py);
                    break;
            }
            allJoints.push(j1);
            return j1;
        }

        Member *System::addMember(Joint *j1, Joint *j2) {
            allMembers.push(new Member(j1, j2));
            allMembers.top()->id = allMembers.length() - 1;
            return allMembers.top();
        }

        void System::analyse() {
            length = allMembers.length() * 4 + 1;
            short_vector::Vector vec(length - 1);
            agMatrix = new sparse_system::AgMatrix(length, allJoints.length() * 4);
            counter = 0;
            matDegree = 0;

            stackIterator<Joint *> iterator;
            iterator.start(allJoints);
            while (iterator.end()) {

                iterator.read()->calculateConstants();
                ++iterator;

            }

            agMatrix->setDegree(matDegree);
            agMatrix->print();

            agMatrix->sort();
//        printf("%d\n", matDegree);
            agMatrix->solveSystem(vec);

            stackIterator<Member *> iterator2;

            iterator2.start(allMembers);
            while (iterator2.end()) {
                iterator2.read()->shearConstant = vec[4 * iterator2.read()->id];
                iterator2.read()->momentConstant = vec[4 * iterator2.read()->id + 1];
                iterator2.read()->slopeConstant = vec[4 * iterator2.read()->id + 2];
                iterator2.read()->deflectionConstant = vec[4 * iterator2.read()->id + 3];

                ++iterator2;
            }

        }

        void System::setE_GP(float EG) {
            E = EG;
        }

    void System::Draw() {

        stackIterator<Member*> st;
        st.start(allMembers);

        while (st.end()){
            st.read()->Draw(htmlfile);
            ++st;
        }
        fprintf(htmlfile,"</body>\n"
                         "</body>\n"
                         "</html>");
        fclose(htmlfile);
        system("start C:\\Users\\m7mdn\\Desktop\\structure\\newfile.html");

    }

    void System::print() {
        stackIterator<Member*> st;
        st.start(allMembers);

        while (st.end()){
            st.read()->printDeflection();
            ++st;
        }
    }

    Moment::Moment(Type actionPoint, Type magnitude) : actionPoint(actionPoint), magnitude(magnitude) {
            direction[0] = (0);
            direction[1] = (0);
            direction[2] = (1.0f);
        }

        Type Moment::momentAction(Type position, int direction_) {
//        printf("position = %f ,action point = %f\n", position, actionPoint);

            if (position >= actionPoint) {
                return direction[direction_] * magnitude;
            }
            return 0.0;
        }

        Type Moment::slopeAction(Type position, int direction_) {
            if (position >= actionPoint) {
                return direction[direction_] * magnitude * (position - actionPoint);
            }
            return 0.0;
        }

        Type Moment::deflectionAction(Type position, int direction_) {
            if (position >= actionPoint) {
                return direction[direction_] * magnitude * (position - actionPoint) * (position - actionPoint) / 2.0;
            }
            return 0.0;
        }

        FixedConnection::FixedConnection(Type px, Type py) : Joint(px, py) {

        }

        void FixedConnection::calculateConstants() {
            shearCondition(this);
            momentCondition(this);
            slopeCondition(this);
            deflectionCondition(this);

        }

        void PinJoint::calculateConstants() {
            shearCondition(this);
            singleMomentCondition(this);
            deflectionCondition(this);

        }

    }
