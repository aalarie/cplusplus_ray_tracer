//
//  camera.h
//
//  Author: Antoine Alarie
//
//  Library providing synthetic camera functionality
//

#ifndef _camera_h_
#define _camera_h_

#include "matrix.h"

//  Class to create a screen object with specified width, height and aspect ratio
//
class Screen_t
{
private:
    int w, h;
    double ratio;
public:
    Screen_t(){}
    Screen_t(int height, double aspect);
    double aspect();
    int height();
    int width();
};

//  Class to create a synthetic camera
//
class Camera
{
public:
    Matrix u, v, n, Mv, Mp, T1, T2, S1, S2, W, M;
    
    Camera(){}
    Camera(Matrix E, Matrix G, Matrix UP, double N, double F, double angle, Screen_t screen);
};

#endif
