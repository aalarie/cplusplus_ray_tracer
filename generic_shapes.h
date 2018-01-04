//
//  generic_shapes.h
//
//  Author: Antoine Alarie
//
//  Library providing generic shape objects
//

#ifndef _generic_shapes_h_
#define _generic_shapes_h_

#include "matrix.h"
#include <string>
using namespace std;

// Colour type stores infor about r,g,b channels
//
typedef struct {
    int r, g, b;
} Colour;

//  Class representing any generic shape, with color and shading attributes as well
//  as an associated transformation matrix
//
class GenericShape
{
public:
    Matrix Mi, M;          //Inverse of associated transformation matrix
    Colour colour, specular_colour, diffuse_colour, ambient_colour; //Colour attributes
    double specular_coeff, diffuse_coeff, ambient_coeff; //Shading attributes
    
    GenericShape(){}
    virtual ~GenericShape(){};
    GenericShape(Matrix transformation, Colour col, double ac, double dc, double sc);
    GenericShape(Matrix transformation, Colour a_col, Colour d_col, Colour s_col);
    virtual double intersect(Matrix e, Matrix d){return -1;};
    virtual Matrix normal(Matrix){return Matrix();};
    virtual string type(){return "generic";};
};

//  Class representing a generic sphere
//
class Sphere: public GenericShape
{
public:
    Sphere(){}
    Sphere(Matrix transformation, Colour col, double ac, double dc, double sc) :
        GenericShape(transformation, col, ac, dc, sc){}
    Sphere(Matrix transformation, Colour a_col, Colour d_col, Colour s_col) :
        GenericShape(transformation, a_col, d_col, s_col){}
    double intersect(Matrix e, Matrix d);
    Matrix normal(Matrix);
    string type(){return "sphere";};
};

//  Class representing a generic plane
//
class Plane: public GenericShape
{
public:
    Plane(){}
    Plane(Matrix transformation, Colour col, double ac, double dc, double sc) :
        GenericShape(transformation, col, ac, dc, sc){}
    Plane(Matrix transformation, Colour a_col, Colour d_col, Colour s_col) :
        GenericShape(transformation, a_col, d_col, s_col){}
    double intersect(Matrix e, Matrix d);
    Matrix normal(Matrix);
    string type(){return "plane";};
};

//  Class representing a generic cylinder
//
class Cylinder: public GenericShape
{
public:
    Cylinder(){}
    Cylinder(Matrix transformation, Colour col, double ac, double dc, double sc) :
        GenericShape(transformation, col, ac, dc, sc){}
    Cylinder(Matrix transformation, Colour a_col, Colour d_col, Colour s_col) :
        GenericShape(transformation, a_col, d_col, s_col){}
    double intersect(Matrix e, Matrix d);
    Matrix normal(Matrix);
    string type(){return "cylinder";};
};

//  Class representing a generic cone
//
class Cone: public GenericShape
{
public:
    Cone(){}
    Cone(Matrix transformation, Colour col, double ac, double dc, double sc) :
        GenericShape(transformation, col, ac, dc, sc){}
    Cone(Matrix transformation, Colour a_col, Colour d_col, Colour s_col) :
        GenericShape(transformation, a_col, d_col, s_col){}
    double intersect(Matrix e, Matrix d);
    Matrix normal(Matrix);
    string type(){return "cone";};
};

#endif
