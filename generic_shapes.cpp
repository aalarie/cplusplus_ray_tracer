//
//  generic_shapes.cpp
//
//  Author: Antoine Alarie
//
//  Library providing generic shape objects
//

#include "generic_shapes.h"
#include <math.h>
#include <iostream>

//  Constructs a GenericShape with the transformation Matrix transformation and colour col
//
GenericShape::GenericShape(Matrix transformation, Colour col, double ac, double dc, double sc)
{
    M = transformation;
    Mi = mat_inverse(transformation);
    ambient_colour = col;
    diffuse_colour = col;
    specular_colour = col;
    colour = col;
    ambient_coeff = ac;
    diffuse_coeff = dc;
    specular_coeff = sc;
}

//  Constructs a GenericShape with the transformation Matrix transformation and unique colours for each
//  shading attribute
//
GenericShape::GenericShape(Matrix transformation, Colour a_col, Colour d_col, Colour s_col)
{
    M = transformation;
    Mi = mat_inverse(transformation);
    ambient_colour = a_col;
    diffuse_colour = d_col;
    specular_colour = s_col;
}

///////Intersections

//  Determines the intersect between ray E+D and this sphere
//  Returns intersection value t if they intersect, -1 otherwise
//
double Sphere::intersect(Matrix E, Matrix D)
{
    Matrix Et, Dt;
    
    //Apply inverse transformation to ray
    Et = Mi * E;
    Dt = Mi * D;
    Et.member[3][0] = 0;
    
    //Compute intersection value t
    double a, b, c, dis, t1, t2;
    a = pow(Dt.norm(), 2);
    b = Et.dot(Dt);
    c = pow(Et.norm(), 2)-1;
    
    Et.free_members();
    Dt.free_members();
    
    dis = pow(b, 2) - a * c;
    
    //If discriminant is less than 0 then no intersection, return -1
    if (dis < 0)
        return -1;
    
    t1 = -(b/a) + (sqrt(dis)/a);
    t2 = -(b/a) - (sqrt(dis)/a);
    
    //Else return smallest t value (first intersection)
    if (t1 < t2 && t1 > 0)
        return t1;
    return t2;
}

//  Determines the intersect between ray E+D and this plane
//  Returns intersection value t if they intersect, -1 otherwise
//
double Plane::intersect(Matrix E, Matrix D)
{
    Matrix Et, Dt;
    
    //Apply inverse transformation to ray
    Et = Mi * E;
    Dt = Mi * D;
    
    double dz = Dt.member[2][0];
    double ez = Et.member[2][0];
    
    Dt.free_members();
    Et.free_members();
    
    //If dz is 0 then ray is parallel to plane, no intersection, return -1
    if (dz == 0)
        return -1;
    
    //Else return intersection value t
    return -(ez/dz);
}

//  Determines the intersect between ray E+D and this cylinder
//  Returns intersection value t if they intersect, -1 otherwise
//
double Cylinder::intersect(Matrix E, Matrix D)
{
    Matrix Et, Dt;
    
    //Apply inverse transformation to ray
    Et = Mi * E;
    Dt = Mi * D;
    Et.member[3][0] = 0;
    
    //Compute intersection value t
    double a, b, c, z, dis, t1, t2, t, temp = -1;
    a = pow(Dt.member[0][0], 2) + pow(Dt.member[1][0], 2);
    b = Et.member[0][0] * Dt.member[0][0] + Et.member[1][0] * Dt.member[1][0];
    c = pow(Et.member[0][0], 2) + pow(Et.member[1][0], 2) - 1;
    
    dis = pow(b, 2) - a * c;
    
    //If discriminant is less than 0 then no intersection, return -1
    if (dis < 0)
        return -1;
    
    //Otherwise obtain the intersection points
    t1 = -(b/a) + (sqrt(dis)/a);
    t2 = -(b/a) - (sqrt(dis)/a);
    
    t = min(t1, t2);
    
    //If the smallest intersection is positive, then it is the closest intersection
    if (t > 0)
    {
        //If z is in [-1,1] then the ray intersects the cylinder wall and this is the first intersect
        z = Et.member[2][0] + Dt.member[2][0] * t;
        if (z >= -1 && z <= 1)
            return t;
    }
    
    ////Otherwise check if it intersects the caps (translated generic plane & inside unit circle)
    ////Keep the lowest intersection point
    //Check top cap
    t = (1 - Et.member[2][0])/Dt.member[2][0];
    if (t > 0)
    {
        //Check if intersection is within unit circle if so it intersects
        if ((pow(Et.member[0][0] + Dt.member[0][0] * t, 2) + pow(Et.member[1][0] + Dt.member[1][0] * t, 2)) < 1)
            temp = t;
    }
    
    //Check bottom cap
    t = (-1 - Et.member[2][0])/Dt.member[2][0];
    if (t > 0 && t < temp)
    {
        //Check if intersection is within unit circle if so it intersects
        if ((pow(Et.member[0][0] + Dt.member[0][0] * t, 2) + pow(Et.member[1][0] + Dt.member[1][0] * t, 2)) < 1)
            temp = t;
    }
    
    //If ray intersected a cap then return intersect
    if (temp != -1)
        return temp;
    
    //Otherwise check other wall intersection
    t = max(t1,t2);
    if (t > 0)
    {
        //If z is in [-1,1] then the ray intersects the cylinder wall
        z = Et.member[2][0] + Dt.member[2][0] * t;
        if (z >= -1 && z <= 1)
            return t;
    }
    
    //Otherwise no intersection
    return -1;
}

//  Determines the intersect between ray E+D and this cone
//  Returns intersection value t if they intersect, -1 otherwise
//
double Cone::intersect(Matrix E, Matrix D)
{
    Matrix Et, Dt;
    
    //Apply inverse transformation to ray
    Et = Mi * E;
    Dt = Mi * D;
    
    //Compute intersection value t
    double a, b, c, z, dis, t1, t2, t;
    
    a = pow(Dt.member[0][0], 2) + pow(Dt.member[1][0], 2) - pow(Dt.member[2][0], 2)/4.0;
    b = Et.member[0][0] * Dt.member[0][0] + Et.member[1][0] * Dt.member[1][0] + Dt.member[2][0] * (1.0 - Et.member[2][0])/4.0;
    c = pow(Et.member[0][0], 2) + pow(Et.member[1][0], 2) - pow(1.0 - Et.member[2][0], 2)/4.0;
    
    dis = pow(b, 2) - a * c;
    
    //If discriminant is less than 0 then no intersection, return -1
    if (dis < 0)
        return -1;
    
    //Otherwise obtain the intersection points
    t1 = -(b/a) + (sqrt(dis)/a);
    t2 = -(b/a) - (sqrt(dis)/a);
    
    t = min(t1, t2);
    
    //If the smallest intersection is positive, then it is the closest intersection
    if (t > 0)
    {
        //If z is in [-1,1] then the ray intersects the cone wall and this is the first intersect
        z = Et.member[2][0] + Dt.member[2][0] * t;
        if (z >= -1 && z <= 1)
            return t;
    }
    
    ////Otherwise check if it intersects the base (translated generic plane & inside unit circle)
    //Check base
    t = (-1 - Et.member[2][0])/Dt.member[2][0];
    if (t > 0)
    {
        //Check if intersection is within unit circle if so it intersects
        if ((pow(Et.member[0][0] + Dt.member[0][0] * t, 2) + pow(Et.member[1][0] + Dt.member[1][0] * t, 2)) < 1)
            return t;
    }
    
    //Otherwise check other wall intersection
    t = max(t1,t2);
    if (t > 0)
    {
        //If z is in [-1,1] then the ray intersects the cone wall
        z = Et.member[2][0] + Dt.member[2][0] * t;
        if (z >= -1 && z <= 1)
            return t;
    }
    
    //Otherwise no intersection
    return -1;
}

//////Normals

//  Returns the normal to this sphere
//
Matrix Sphere::normal(Matrix intersect_point)
{
    //If type is sphere, normal vector is intersect point
    Matrix n = Matrix(4,1);
    n.member[0][0] = intersect_point.member[0][0];
    n.member[1][0] = intersect_point.member[1][0];
    n.member[2][0] = intersect_point.member[2][0];
    n.member[3][0] = 0.0;
    
    return n;
}

//  Returns the normal to this plane
//
Matrix Plane::normal(Matrix intersect_point)
{
    //If type is plane, normal vector is [0,0,1,0]
    Matrix n = Matrix(4,1);
    n.member[0][0] = 0.0;
    n.member[1][0] = 0.0;
    n.member[2][0] = 1.0;
    n.member[3][0] = 0.0;
    
    //Transform
    return M * n;
}

//  Returns the normal to this cylinder
//
Matrix Cylinder::normal(Matrix intersect_point)
{
    Matrix n = Matrix(4,1);
    
    //Inverse intersect point
    Matrix Pi = Mi * intersect_point;
    
    //If the intersect is on the bottom cap
    if (Pi.member[2][0] == -1)
    {
        n.member[0][0] = 0.0;
        n.member[1][0] = 0.0;
        n.member[2][0] = -1.0;
        n.member[3][0] = 0.0;
        Pi.free_members();
        //Transform
        return M * n;
    }
    
    //If the intersect is on the top cap
    else if (Pi.member[2][0] == 1)
    {
        n.member[0][0] = 0.0;
        n.member[1][0] = 0.0;
        n.member[2][0] = 1.0;
        n.member[3][0] = 0.0;
        Pi.free_members();
        //Transform
        return M * n;
    }
    
    //Otherwise the intersect is with the cylinder wall
    else
    {
        double denom = sqrt(pow(Pi.member[0][0], 2) + pow(Pi.member[1][0], 2));
        n.member[0][0] = Pi.member[0][0]/denom;
        n.member[1][0] = Pi.member[1][0]/denom;
        n.member[2][0] = 0.0;
        n.member[3][0] = 0.0;
        Pi.free_members();
        //Transform
        return M * n;
    }
}

//  Returns the normal to this cone
//
Matrix Cone::normal(Matrix intersect_point)
{
    Matrix n = Matrix(4,1);
    
    //Inverse intersect point
    Matrix Pi = Mi * intersect_point;
    
    //If the intersect is on the base
    if (Pi.member[2][0] == -1)
    {
        n.member[0][0] = 0.0;
        n.member[1][0] = 0.0;
        n.member[2][0] = -1.0;
        n.member[3][0] = 0.0;
        Pi.free_members();
        //Transform
        return M * n;
    }
    
    //Otherwise the intersect is with the cone wall
    else
    {
        double oper1 = 2.0*pow(Pi.member[0][0], 2) + 2.0*pow(Pi.member[1][0], 2)
            + pow((1.0 - Pi.member[2][0])/2.0, 2);
        oper1 = 1.0/sqrt(oper1);
        n.member[0][0] = oper1 * 2.0 * Pi.member[0][0];
        n.member[1][0] = oper1 * 2.0 * Pi.member[1][0];
        n.member[2][0] = oper1 * ((1.0 - Pi.member[2][0])/2.0);
        n.member[3][0] = 0.0;
        Pi.free_members();
        return M * n;
    }
}
