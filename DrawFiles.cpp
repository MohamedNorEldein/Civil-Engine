//
// Created by m7md_nor on 10/27/2021.
//

#include "../../header/civil_work.h"

using namespace Sections;

const int width = 400, height = 400, factor = 10;

FILE * draw() {

    std::string text = "<!DOCTYPE html>\n"
                       "<html lang=\"en\">\n"
                       "<head>\n"
                       "    <meta charset=\"UTF-8\">\n"
                       "    <title>program</title>\n"
                       "</head>\n";



    FILE *file = fopen("C:/Users/m7mdn/Desktop/structure/show.html", "w");

    fprintf(file, text.c_str());
    return file;
}



std::string drawRectangle(Shape shape1, int id) {
    char buffer[1000];
    char color[15] = "'#1b8bc2'";
    char color2[15] = "'#ffffff'";
    const char *c = color;

    if(shape1.area<0){
        c  =  color2;
    }

    sprintf(buffer,
            "<script>\n"
            "    var c = document.getElementById(\"%d\");\n"
            "    var ctx = c.getContext(\"2d\");\n"
            "    ctx.beginPath();\n"

            "    ctx.moveTo(%f,%f);\n"
            "    ctx.lineTo(%f,%f);\n"
            "    ctx.lineTo(%f,%f);\n"
            "    ctx.lineTo(%f,%f);\n"
            "    ctx.lineTo(%f,%f);\n"
            "    ctx.closePath();\n"
            "    ctx.fillStyle = %s;\n"
            "    ctx.fill();\n"
            "    ctx.lineWidth = 0;\n"
            "    ctx.strokeStyle = '#1b8bc2';\n"
            "    ctx.stroke();\n"
            "</script>\n",
            id,
            shape1.data[0] * factor, height - shape1.data[1] * factor,       // first point
            shape1.data[0] * factor, height - shape1.data[3] * factor,
            shape1.data[2] * factor, height - shape1.data[3] * factor,
            shape1.data[2] * factor, height - shape1.data[1] * factor,
            shape1.data[0] * factor, height - shape1.data[1] * factor,
            c
    );      // last point
    return std::string(buffer);
}

std::string drawTriangle(Type * data, int id, const char * color = "'#7B8487FF'" , int factor2= factor) {
//    printf("(%f, %f) (%f, %f) (%f, %f)\n", data[0], data[1], data[2], data[3], data[4], data[5] );
    std::string text;
    char buffer[350];

    sprintf(buffer,"<script>\n"
                   "    var c = document.getElementById(\"%d\");\n"
                   "    var ctx = c.getContext(\"2d\");\n"
                   "    ctx.beginPath();\n"

                   "    ctx.moveTo(%f,%f);\n",id, data[4] * factor2,height - data[5] * factor2);   // first point
    text += buffer;

    for (int i = 0; i < 6; i += 2) {
        sprintf(buffer,
                "    ctx.lineTo(%f,%f);\n"
                , data[i] * factor2,height - data[i + 1] * factor2);   // first point
        text += buffer;
    }
    text += "    ctx.fillStyle =";
    text += color;
    text += ";\n"
            "    ctx.fill();\n"
            "    ctx.lineWidth = 0;\n"
            "</script>\n";

    return text;
}

std::string drawTriangle(Shape shape, int id, const char * color = "'#7B8487FF'") {
//    printf("(%f, %f) (%f, %f) (%f, %f)\n", data[0], data[1], data[2], data[3], data[4], data[5] );
    Type * data = shape.data;
    return drawTriangle(data,id,color);
}


std::string drawCircle(Shape shape1, int id) {

    char buffer[1000];
    char color[15] = "'#1b8bc2'";
    char color2[15] = "'#ffffff'";
    const char *c = color;

    if(shape1.area_state==-1){
        c  =  color2;
    }
    double r = sqrt(pow((shape1.data[2] - shape1.data[0]),2)+pow((shape1.data[3] - shape1.data[1]),2));

    sprintf(buffer,
            "<script>\n"
            "    var c = document.getElementById(\"%d\");\n"
            "    var ctx = c.getContext(\"2d\");\n"
            "    ctx.moveTo(%lf, %lf);\n"
            "    ctx.arc(%lf,%lf,%lf,0,3.14159265359 * 2 );\n"
            "    ctx.fillStyle = '#1b8bc2';\n"
            "    ctx.fill();\n"
            "    ctx.lineWidth = 0;\n"
            "    ctx.strokeStyle = '#1b8bc2';\n"
            "    ctx.stroke();\n"
            "</script>\n",
            id, factor * shape1.data[0],height - factor*shape1.data[1], factor * shape1.data[0],height - factor*shape1.data[1],factor*r);
    return std::string(buffer);

}

std::string draw_arc(Shape shape1, int id){
    char buffer[1000];
    char color[15] = "'#1b8bc2'";
    char color2[15] = "'#ffffff'";
    const char *c = color;

    if(shape1.area_state==-1){
        c  =  color2;
    }

    sprintf(buffer,
            "<script>\n"
            "    var c = document.getElementById(\"%d\");\n"
            "    var ctx = c.getContext(\"2d\");\n"
            "    ctx.moveTo(%lf, %lf);\n"
            "    ctx.arc(%lf,%lf,%lf,%lf,%lf);\n"
            "    ctx.fillStyle = '#1b8bc2';\n"
            "    ctx.fill();\n"
            "    ctx.lineWidth = 0;\n"
//            "    ctx.strokeStyle = '#1b8bc2';\n"
//            "    ctx.stroke();\n"
            "</script>\n",
            id, factor * shape1.data[0],height - factor*shape1.data[1], factor * shape1.data[0],height - factor*shape1.data[1],factor*shape1.data[2], shape1.data[3], shape1.data[4]);
    return std::string(buffer);
}

std::string draw_shape(Shape shape1, int id) {
    Type * tr =new Type[6];

    switch (shape1.type) {


        case rectangle:
            return ::drawRectangle(shape1, id);
            break;

        case triangle:
            if(shape1.area<0){
                return ::drawTriangle(shape1, id,"#fff");
            }
            return drawTriangle(shape1, id);
            break;

        case traverse:
        {
            std::string text;


            for (int i = 2; i < shape1.numOfPoints * 2 - 2; i += 2) {

                tr[0] = shape1.data[0];
                tr[1] = shape1.data[1];
                tr[2] = shape1.data[i];
                tr[3] = shape1.data[i + 1];
                tr[4] = shape1.data[i + 2];
                tr[5] = shape1.data[i + 3];

                float area = (tr[2] - tr[0]) * (tr[5] - tr[1]);
                area -= (tr[4] - tr[0]) * (tr[3] - tr[1]);
                if (area > 0) {
                    text += drawTriangle(tr, id);
                } else {
                    text += drawTriangle(tr, id,"'#fff'");
                }
            }
            return text;
        }
            break;

        case circle:
            return drawCircle(shape1, id);
            break;

        case arc:
            return draw_arc(shape1, id);
    }
    delete[] tr;

    return "";
}
/*
std::string drawCoreTraverse(Shape *shape1, int id)
{
    std::string text;
    double * tr =new double[6];

    for (int i = 2; i < shape1->numOfPoints * 2 - 2; i += 2) {

        tr[0] = shape1->ShapeCore[0] + shape1->cg[0];
        tr[1] = shape1->ShapeCore[1] + shape1->cg[1];
        tr[2] = shape1->ShapeCore[i] + shape1->cg[0];
        tr[3] = shape1->ShapeCore[i + 1] + shape1->cg[1];
        tr[4] = shape1->ShapeCore[i + 2] + shape1->cg[0];
        tr[5] = shape1->ShapeCore[i + 3] + shape1->cg[1];

        text += drawTriangle(tr, id, "'rgba(150,94,212,0.5)'");

    }

    delete[] tr;
    return text;
}

std::string drawSectionCore(Sections* section, int id){
    std::string text;
    section->start();
    int i =0;
    while (i < section->length())
        {
        text += drawCoreTraverse(section->read(),id);
        ++i;
        section->operator++();
    }
    return text;
}
*/
void Section::drawSection(FILE* file, int id)
{
    std::string text;

    char text1[400];
    sprintf(text1, "<h1>canvas</h1>\n"
                          "\n"
                          "<body>\n"
                          "\n"
                          "<canvas id=\"%d\" width=\"400\"height=\"400\" style=\"border:1px solid #d3d3d3;\">\n"
                          "    Your browser does not jointType the HTML canvas tag.</canvas>\n",id);

    text+= text1;

    char text2[506];
    Shape *shape1;


    SingleLinkedList::stackIterator<Shape*> iterator;

    iterator.start(shapeList);
    while (iterator.end()) {
        shape1 = (iterator.read());

        text += draw_shape(*shape1, id);

        char type[10];

        switch (shape1->type) {
            case rectangle:
                sprintf(type,"rectangle") ;
                break;

            case triangle:
                sprintf(type,"triangle") ;
                break;

            case circle:
                sprintf(type,"circle") ;
                break;

            case traverse:
                sprintf(type,"traverse") ;
                break;

            case arc:
                sprintf(type,"arc") ;
                break;
        }


        sprintf(text2, "<p><pre>\n"
                       "type : %s   DCShape = %f   center is (%f , %f)    moment of inertia :Ix = %f    Iy = %f    Ixy=%f\n"
                       "</pre></p>", type, shape1->area, shape1->cg[0], shape1->cg[1], shape1->Ix, shape1->Iy,
                shape1->Ixy);


        text += text2;

        ++iterator;
    }

    if (area == 0) {
        Gx = 0;
        Gy = 0;
    } else {
        Gx /= area;
        Gy /= area;
    }

//    calculateTotalMomentOfInertia(section, section->Ix, section->Iy, section->Ixy, section->Gx, section->Gy, section->Imax, section->Imin, section->theta);

    sprintf(text2, "<p><pre>\n"
                   "total Shape : DCShape = %lf     center is (%f , %f)\tIx= %lf\tIy= %lf\tIxy= %lf\tImax = %f\t Imin = %f\ttheta = %f  \n"
                   "</pre></p>\n", area, Gx, Gy, Ix, Iy, Ixy, Imax, Imin, theta * 360);

    text += text2;

//    calcCore(section);
//    text += drawSectionCore(section, id);

    sprintf(text2,"<p><pre>\n"
                  " shear x = %lf, shear y = %lf ,normal = %lf  \n"
                  " moment x = %lf , moment y = %lf ,torque = %lf  \n"
                  "</pre></p>\n"
                  , force[0], force[1], force[2], moment[0], moment[1], moment[2]);
    text += text2;
    fprintf(file, text.c_str());

    calcStressProfile();
    fprintf(file, "<pre>"
                  "<p>"
                  "max normal stress = %f\n"
                  "min normal stress = %f\n"
                  "</pre>"
                  "</p>"
                  "</body>"
                  "</html>",smax, smin);

}

double drawFactor1=10, factor1=100,df=500, startF1=50, start2=400;

void StructuralSystem::TrussSystem::trussMember::draw(FILE *file) {
    fprintf(file, "<script>"

                  "var c = document.getElementById(\"0\");\n"
                  "var ctx = c.getContext(\"2d\");\n"
                  "var cty = c.getContext(\"2d\");\n"

                  "cty.strokeStyle = \"blue\";\n"
                  "cty.lineWidth = 0.5;"
                  "cty.beginPath();\n"
                  "cty.moveTo(%f,%f);\n"
                  "cty.lineTo(%f, %f);\n"

                  "cty.stroke();\n"

                  "ctx.strokeStyle = \"red\";\n"
                  "ctx.lineWidth = 0.5;"
                  "ctx.beginPath();\n"
                  "ctx.moveTo(%f,%f);\n"
                  "ctx.lineTo(%f, %f);\n"
                  "ctx.stroke();\n"
                  "</script>"

            ,
            start->position[0] * factor1 + startF1, -start->position[1] * factor1 + start2,
            end->position[0] * factor1 + startF1, -end->position[1] * factor1 + start2,
            (start->position[0]+df*start->deformation[0]) * factor1 + startF1,
            -(start->position[1]+df*start->deformation[1]) * factor1 + start2,
            (end->position[0]+df*end->deformation[0]) * factor1 + startF1,
            -(end->position[1]+df*end->deformation[1]) * factor1 + start2);


}

void StructuralSystem::TrussSystem::Truss::Draw() {


    SingleLinkedList::stackIterator<trussMember*> iterator;
    iterator.start(allMembers);
    while (iterator.end()){
        iterator.read()->draw(htmlfile);
        ++iterator;
    }
    fprintf(htmlfile,"<pre> the red => the deflected shape\n </pre>"
                     "</body>\n"
                     "</body>\n"
                     "</html>");
    fclose(htmlfile);
    system("start file.html");

}