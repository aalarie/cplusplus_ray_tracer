//
//  camera.cpp
//
//  Author: Antoine Alarie
//
//  Library providing synthetic camera functionality
//

#include "camera.h"
#include <math.h>

////Screen definition

//  Creates a screen with specified height and aspect ratio
//
Screen_t::Screen_t(int height, double aspect)
{
    h = height;
    w = height * aspect;
    ratio = aspect;
}

//  Returns the aspect ratio of this screen
//
double Screen_t::aspect()
{
    return ratio;
}

//  Returns the height of this screen
//
int Screen_t::height()
{
    return h;
}

//  Returns the width of this screen
//
int Screen_t::width()
{
    return w;
}

////Camera definition

//  Creates a camera at origin E gazing at point G in world coordinates with the specified parameters:
//  UP is an indication of the upward direction, N is the distance from the near plane, F is the
//  distance from the far plane, angle is the viewing angle and screen is the size of the screen
//  on which the image will be projected.
//  The near plane and the screen will have the same proportions and are figuratively equivalent
//
Camera::Camera(Matrix E, Matrix G, Matrix UP, double N, double F, double angle, Screen_t screen)
{
    ////Find viewing coordinate system axes
    //Determine n axis
    n = E - G;
    n.normalize();
    n.set_rows(3);
    
    //Determine u axis
    UP.set_rows(3);
    u = cross_product(UP, n);
    u.normalize();
    UP.set_rows(4);
    
    //Determine v axis
    v = cross_product(n, u);
    v.normalize();
    
    //Make u,v,n homogeneous for future use
    u.set_rows(4);
    u.member[3][0] = 0;
    v.set_rows(4);
    v.member[3][0] = 0;
    n.set_rows(4);
    n.member[3][0] = 0;
    
    //Create the Mv transformation matrix (world to viewing system)
    Mv = Matrix(4,4);
    
    Mv.member[0][0] = u.member[0][0];
    Mv.member[0][1] = u.member[1][0];
    Mv.member[0][2] = u.member[2][0];
    Mv.member[0][3] = -(E.member[0][0]*u.member[0][0] + E.member[1][0]*u.member[1][0] + E.member[2][0]*u.member[2][0]);
    
    Mv.member[1][0] = v.member[0][0];
    Mv.member[1][1] = v.member[1][0];
    Mv.member[1][2] = v.member[2][0];
    Mv.member[1][3] = -(E.member[0][0]*v.member[0][0] + E.member[1][0]*v.member[1][0] + E.member[2][0]*v.member[2][0]);
    
    Mv.member[2][0] = n.member[0][0];
    Mv.member[2][1] = n.member[1][0];
    Mv.member[2][2] = n.member[2][0];
    Mv.member[2][3] = -(E.member[0][0]*n.member[0][0] + E.member[1][0]*n.member[1][0] + E.member[2][0]*n.member[2][0]);
    
    Mv.member[3][0] = 0.0;
    Mv.member[3][1] = 0.0;
    Mv.member[3][2] = 0.0;
    Mv.member[3][3] = 1.0;
    
    //Create the Mp transformation matrix (perspective transform)
    float a = -1.0 * (F + N)/(F - N);
    float b = -2.0 * (F * N)/(F - N);
    
    Mp = mat_identity(4);
    Mp.member[0][0] = N;
    Mp.member[1][1] = N;
    Mp.member[2][2] = a;
    Mp.member[2][3] = b;
    Mp.member[3][2] = -1.0;
    Mp.member[3][3] = 0.0;
    
    //Create scaling and translation matrices to fit into canonical view volume and screen coordinates
    float top = N * tan(M_PI/180.0 * angle/2.0);
    float right = screen.aspect() * top;
    float bottom = -top;
    float left = -right;
    
    //T1
    T1 = mat_identity(4);
    T1.member[0][3] = -(right + left)/2.0;
    T1.member[1][3] = -(top + bottom)/2.0;
    
    //S1
    S1 = mat_identity(4);
    S1.member[0][0] = 2.0/(right - left);
    S1.member[1][1] = 2.0/(top - bottom);
    
    //T2
    T2 = mat_identity(4);
    T2.member[0][3] = 1.0;
    T2.member[1][3] = 1.0;
    
    //S2
    S2 = mat_identity(4);
    S2.member[0][0] = screen.width()/2.0;
    S2.member[1][1] = screen.height()/2.0;
    
    //W
    W = mat_identity(4);
    W.member[1][1] = -1.0;
    W.member[1][3] = screen.height();
    
    M = W * S2 * T2 * S1 * T1 * Mp * Mv; //Save composite transformation from world to screen
}

