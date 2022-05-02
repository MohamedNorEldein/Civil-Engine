//
// Created by m7md_nor on 10/27/2021.
//

#ifndef UNTITLED1_CIVIL_WORK_H
#define UNTITLED1_CIVIL_WORK_H


#include <iostream>
#include "DataStructures.h"
//#include "calc.h"
#include "LinearAlgebra.h"

typedef double Type;
static const int systemDegree = 2;

namespace Sections {

    const int numOfDimension = 2;

    enum shape_type {
        circle, rectangle, triangle, traverse, arc
    };

    typedef struct Shape {
        shape_type type;
        int area_state;
        int numOfPoints;
        Type *data;
        Type area;
        Type cg[numOfDimension];
        Type Rx, Ry, Ix, Iy, Ixy;
    } Shape;

    typedef SingleLinkedList::stack<Shape *> ShapeList;

    class Section {
    public:

        Type *force;
        // N    Qx    Qy
        Type *moment;
        // Mx   My    t__

        Type area, Ix, Iy, Ixy;
        Type Imax, Imin, theta, Gx, Gy;

        Type TensionCriticalStress, CompressionCriticalStress;
        ShapeList *shapeList;

        Type smax, smin, sav;

        Section(Type *force, Type *moment);

        Section();

        ~Section();

        void addShape(Shape &) const;

        void addShape(Shape *) const;

        Shape *addRectangle(Type startx, Type starty, Type endx, Type endy, int areaState = 1);

        Shape *addTraverse(Type *data, int numofpoints, int areastate = 1);

        Type sectionProperties();

        Type calcStressProfile();

        void drawSection(FILE *file, int id);

//   void readSection()

    };
}
    typedef SingleLinkedList::stack< Sections::Section*> SectionList;
//

//
//void DCShape(Shape&, int id);
//
//void calculateTotalMomentOfInertia(Sections* shapeList, double &Ix, double &Iy, double &Ixy, double Rx, double Ry , double &Imax, double &Imin, double &theta);
//
//calcStressProfile * calcStress(Sections*);
//
//void calcCore(Sections *section);
//
//    int readfile(const char *name);
//}
//}

//-------------------------------------------------------------------
namespace StructuralSystem {
    enum memberType {
        truss, beam, shell, block
    };

    enum forceDegree {
        concentrated = -1, constant, linear, parabolic
    };

    enum supportType {
        None = 0, roller, rocker, pin, fixedSupport, freeEnd, pinJoint, fixedJoint

    };

    namespace TrussSystem {
        class Joint;
        class support;
        class Rocker;
        class Roller;
        class Truss;
        class trussMember;

        class trussMember {
        protected:
            float normalForce, Area,E,length;
            int i,defined;
            Joint *start, *end;

            explicit trussMember(int i,Joint *start,Joint *end,float ,float ,float );
            void draw(FILE*);
            inline ~trussMember() {
                printf("deleting truss member at %x ...\n", this);
            }
            friend Joint;
            friend Truss;
            friend support;
            friend Roller;
            friend Rocker;
        };

        class Joint {
        protected:
            float position[systemDegree]{};
            int numberOfMembers, iteration, index;
            trussMember *members[5]{};
            Joint *membersEnd[5]{};

            supportType support;
            float externalForce[systemDegree]{};
            double deformation[systemDegree]{};

            void addMember(trussMember *member1, Joint *);
            virtual void addEqnToMatrix(short_vector::Matrix &matrix, short_vector::Vector &vector, int Dimn);
            void newLoadSystem(short_vector::Matrix &matrix, short_vector::Vector& vector, int Dimn);

            public:
            Joint(float position_x, float position_y, float position_z, int index, supportType support = None);

            ~Joint();

            void changePoint(float position_x, float position_y, float position_z);

            void setExternalForce(float px, float py, float pz);

            void print();

            friend void direction(Joint *startPoint, Joint *endPoint, float *Normal);
            friend Truss;
            friend trussMember;
        };

        class support:public Joint{
            float reaction[systemDegree]{};
            float &getReaction(int i);
        protected:
            inline support(float position_x, float position_y,float  position_z, int index, supportType support1):
            Joint(position_x,position_y,position_z,index,support1){
            }
        };

        class Rocker:public support{
            float reaction[systemDegree]{};
            void addEqnToMatrix(short_vector::Matrix &matrix, short_vector::Vector &vector, int Dimn) override;

        public:
            Rocker(float position_x, float position_y,float  position_z, int index);

        };

        class Roller:public support{

            float directionV[systemDegree]{}, reaction;
            void addEqnToMatrix(short_vector::Matrix &matrix, short_vector::Vector &vector, int Dimn) override;
        public:
            Roller(float position_x, float position_y,float  position_z,float dx, float dy, float dz, int index);

        };

        class Truss {
            Joint *startJoint;
            int numOfJoints, numOfMembers;
            short_vector::Matrix *matrix, *DeflectionMatrix, *EA;
            short_vector::Vector *vector, *deformation;
            FILE *htmlfile;
        public:

            SingleLinkedList::stack<Joint *> allJoints;
            SingleLinkedList::stack<trussMember *> allMembers;

            Truss();

            ~Truss();

            Joint *addJoint(float position_x, float position_y, float position_z);

            Joint *addRocker(float position_x, float position_y, float position_z);
            Joint *addRoller(float position_x, float position_y, float position_z,float dx, float dy, float dz);

            trussMember *addMember(Joint *start, Joint *end, float A=1, float E=1);

            void setStartPoint(Joint *);

            void iterationSolve();

            void print();

            void MatrixAnalysis();

            void MatrixSolve();
            void Draw();
        };

    }

    namespace BeamSystem {

        using namespace SingleLinkedList;

        class Beam;

        class Joint;

        class Member;

        class System;

        class Force {
            Type direction[3]{}; // x along the line of beam
            Type actionPoint, magnitude;
            forceDegree degree;     // any force action is a polynomial
            // -1 => dirac delta function

        public:
            Force(Type xdir, Type ydir, Type zdir, Type actionPoint, Type magnitude, forceDegree degree);

            Type operator()(Type position, int direction_) const;

            Type shearAction(Type position, int direction_);

            Type momentAction(Type position, int direction_);

            Type slopeAction(Type position, int direction_);

            Type deflectionAction(Type position, int direction_);

            Type getMagnitude() const;

        };

        class Moment {
            Type direction[3]{}; // x along the line of beam
            Type actionPoint, magnitude;

        public:
//        Moment(Type xdir, Type ydir, Type zdir, Type actionPoint, Type magnitude, forceDegree degree);
            Moment(Type actionPoint, Type magnitude);

//        Type operator()(Type position, int direction_) const;

            Type momentAction(Type position, int direction_);

            Type slopeAction(Type position, int direction_);

            Type deflectionAction(Type position, int direction_);

        };

        class Joint {
        protected:
            Type position[systemDegree]{};
            stack<Member *> memberList;
        public:
            Joint(Type px, Type py);

            void print();

            virtual void calculateConstants();

            friend Member;

            friend void shearCondition(Joint *node);

            friend void momentCondition(Joint *node);

            friend void slopeCondition(Joint *node);

            friend void deflectionCondition(Joint *node);

            friend void singleMemberSlope(Joint *node);

            friend void singleMomentCondition(Joint *node);

            friend void singleMemberDeflectionCondition(Joint *node);

        };

        class Support : public Joint {
        public:
            Support(Type px1, Type py1);
        };

        class Member {
        public:
            Joint *start, *end;
            SingleLinkedList::stack<Force *> ForceList;
            stack<Moment *> MomentList;
            Type shearConstant, momentConstant, slopeConstant, deflectionConstant;
            Type length, I, direction[systemDegree];
            int id;
            void Draw(FILE* file);

        public:

            int getId() const;

            void setI(Type inertiaMoment_Dm4);

        public:
            Joint *getStart() const;

            Joint *getAnEnd() const;

            Type getLength(Joint *);

//        void loadActions(Type, Type &, Type &, Type &, Type &);
            void loadActions(Joint *, Type &, Type &, Type &, Type &);

        public:
            Member(Joint *start, Joint *end);

            void printDeflection();

            void print();

            void addForce(Type xdir, Type ydir, Type zdir, Type actionPoint, Type magnitude, forceDegree degree);

            void addLoad(Type actionPoint, Type magnitude, forceDegree degree);

            void addMoment(Type actionPoint, Type magnitude);

            Type getLength() const;

            Type shear(float x);

            Type moment(float x);

            Type slope(float x);

            Type deflection(float x);

            void actions(Type, Type &, Type &, Type &, Type &);

            friend System;
            friend Joint;
        };

        class Roller : public Support {

        public:
//        Type getReaction() const;
//
//        void setReaction(Type reaction);

            Roller(Type px1, Type py1);

//        Type getDirection(int i);
            void calculateConstants() override;

        };

        class Rocker : public Support {
            Type reaction[systemDegree]{};

        public:
            Rocker(Type px, Type py);

//        Type *getReaction();
            void calculateConstants() override;

        };

        class FixedSupport : public Support {
        public:
            FixedSupport(Type px, Type py);

            void calculateConstants() override;
        };

        class PinJoint : public Joint {
        public:
            inline PinJoint(Type px1, Type py1) : Joint(px1, py1) {
            }

            void calculateConstants() override;
        };

// done
        class FreeEnd : public Joint {
        public:
            FreeEnd(Type px, Type py);

            void calculateConstants() override;
        };

        class FixedConnection : public Joint {
        public:
            FixedConnection(Type px, Type py);

            void calculateConstants() override;
        };

        class System {
            stack<Joint *> allJoints;
            stack<Member *> allMembers;
            FILE *htmlfile;
        public:
            System();

            void Draw();

            void print();
            void setE_GP(float);

            void analyse();

            Joint *addJoint(Type px, Type py, supportType);

            Member *addMember(Joint *, Joint *);
        };

    }
}

#endif //UNTITLED1_CIVIL_WORK_H
