//
//  ray_tracer.cpp
//
//  Author: Antoine Alarie
//
//  Simple non-recursive ray tracer that reads scene setup from a text file and produces output
//  on screen. Allows user to navigate through scene with arrows
//

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <memory>
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <X11/Xlib.h>
#include "matrix.h"
#include "camera.h"
#include "generic_shapes.h"
using namespace std;

#define N_CHANNELS 3
#define RED 0
#define GREEN 1
#define BLUE 2

//  Light type specifies light colour and and light position in scene world coordinates
//
typedef struct {
    Matrix light_position;
    Colour light_colour;
    int sun = 0;
} Light;

//  Initializes an X11 window specified by screen
//
Display* InitX(Display* d, Window* w, int* s, Screen_t screen) {
    
    d = XOpenDisplay(NULL);
    if(d == NULL) {
        printf("Cannot open display\n");
        exit(1);
    }
    *s = DefaultScreen(d) ;
    *w = XCreateSimpleWindow(d, RootWindow(d, *s), 0, 0, screen.width(), screen.height(), 1, BlackPixel(d, *s), WhitePixel(d, *s));
    Atom delWindow = XInternAtom(d, "WM_DELETE_WINDOW", 0);
    XSetWMProtocols(d, *w, &delWindow, 1);
    XSelectInput(d, *w, ExposureMask | KeyPressMask);
    XMapWindow(d, *w);
    return(d);
}

//  Set the foreground color of the window
//
void SetCurrentColorX(Display* d, GC* gc, unsigned int r, unsigned int g, unsigned int b) {
    
    XSetForeground(d, *gc, r << 16 | g << 8 | b);
}

//  Set the color of the specified pixel
//
void SetPixelX(Display* d, Window w, int s, int i, int j) {
    
    XDrawPoint(d, w, DefaultGC(d, s), i, j);
}

//  Destroys X11 window
//
void QuitX(Display* d, Window w) {
    
    XDestroyWindow(d,w);
    XCloseDisplay(d);
}

#define I_p 1
#define I_a 1
//  Determines light intensity from a shape, its intersect, and rays casted from it
//  to the light and camera
//
double light_intensity(GenericShape* shape, Matrix intersect, Matrix s, Matrix direction)
{
    Matrix n, r, v;
    
    //Normalize shadow ray
    s.normalize();
    
    //Compute n (normal to intersection point)
    n = shape->normal(intersect);
    n.normalize();
    
    //Compute r (direction of specular reflection)
    r = (s * -1) + n * (2 * s.dot(n)/pow(n.norm(),2));
    //r = (s * -1) + 2 * s.dot(n) * n;
    r.normalize();
    
    //Compute v (vector from intersect to camera)
    v = direction * -1;
    v.normalize();
    
    //Compute Ia, Id, & Is
    double max;
    
    double Ia = shape->ambient_coeff * I_a;
    
    if ((max = s.dot(n)) < 0)
        max = 0;
    double Id = I_p * shape->diffuse_coeff * max;
    
    if ((max = pow(r.dot(v), 10)) < 0)
        max = 0;
    double Is = I_p * shape->specular_coeff * max;
    
    //Return I
    return Ia + Id + Is;
}

//  Reads an input file and sets up the parameters of the scene passed down to it
//
void setup_scene(Camera* camera, vector<GenericShape*>* shapes, vector<Light>* lights, Matrix* E, Screen_t* screen, double* width, double* height, double* N)
{
    Matrix G, UP;
    double ASPECT, THETA, F;
    int H;
    
    double Ex,Ey,Ez,Gx,Gy,Gz,UPx,UPy,UPz,Lx,Ly,Lz;
    *E = Matrix(4,1), G = Matrix(4,1), UP = Matrix(4,1);
    ifstream input; //input stream
    string filename;
    
    //Ask user for input file
    cout << "Enter path to input file:\n";
    cin >> filename;
    
    input.open(filename);
    
    if (input.is_open())
    {
        string token;
        
        //Get point E
        getline(input, token); token = token.substr(3); Ex = stod(token);
        getline(input, token); token = token.substr(3); Ey = stod(token);
        getline(input, token); token = token.substr(3); Ez = stod(token);
        
        //Get point G
        getline(input, token); token = token.substr(3); Gx = stod(token);
        getline(input, token); token = token.substr(3); Gy = stod(token);
        getline(input, token); token = token.substr(3); Gz = stod(token);
        
        //Get point UP
        getline(input, token); token = token.substr(4); UPx = stod(token);
        getline(input, token); token = token.substr(4); UPy = stod(token);
        getline(input, token); token = token.substr(4); UPz = stod(token);
        
        //Set up members
        E->member[0][0] = Ex;    G.member[0][0] = Gx;    UP.member[0][0] = UPx;
        E->member[1][0] = Ey;    G.member[1][0] = Gy;    UP.member[1][0] = UPy;
        E->member[2][0] = Ez;    G.member[2][0] = Gz;    UP.member[2][0] = UPz;
        E->member[3][0] = 1.0;   G.member[3][0] = 1.0;   UP.member[3][0] = 1.0;
        
        //Get near plane
        getline(input, token); token = token.substr(2); *N = stod(token);
        
        //Get far plane
        getline(input, token); token = token.substr(2); F = stod(token);
        
        //Get view angle
        getline(input, token); token = token.substr(6); THETA = stod(token);
        
        //Get screen height
        getline(input, token); token = token.substr(7); H = stoi(token);
        
        //Get screen aspect ratio
        getline(input, token); token = token.substr(7); ASPECT = stod(token);
        
        //Insert lights into vector
        while (1)
        {
            //When SUN line is reached break, all lights have been added
            getline(input, token);
            if (token.substr(0,3) == "SUN")
                break;
            
            //Get light position
            token = token.substr(3); Lx = stod(token);
            getline(input, token); token = token.substr(3); Ly = stod(token);
            getline(input, token); token = token.substr(3); Lz = stod(token);
            
            Light light;
            light.light_position = Matrix(4,1);
            
            //Set up light position
            light.light_position.member[0][0] = Lx;
            light.light_position.member[1][0] = Ly;
            light.light_position.member[2][0] = Lz;
            light.light_position.member[3][0] = 1.0;
            
            //Insert light into vector
            lights->emplace_back(light);
        }
        
        //Determine if sun is desired
        int sun = stoi(token.substr(4));
        if (sun)
        {
            //If so add light at infinity
            Light light;
            light.sun = 1;
            lights->emplace_back(light);
        }
        
        //Insert scene shapes into vector
        while (!input.eof())
        {
            input >> token;
            char type = token.at(0);
            
            input >> token; int r = stoi(token);
            input >> token; int g = stoi(token);
            input >> token; int b = stoi(token);
            
            input >> token; double ac = stod(token);
            input >> token; double dc = stod(token);
            input >> token; double sc = stod(token);
            
            Matrix transform = Matrix(4,4);
            
            for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < 4; ++k) {
                    input >> token;
                    transform.member[j][k] = stod(token);
                }
            }
            
            GenericShape* shape;
            switch (type)
            {
                //Sphere
                case 's':
                    shape = new Sphere(transform, Colour{r,g,b}, ac, dc, sc);
                    shapes->emplace_back(shape);
                    break;
                //Plane
                case 'p':
                    shape = new Plane(transform, Colour{r,g,b}, ac, dc, sc);
                    shapes->emplace_back(shape);
                    break;
                //Cylinder
                case 'l':
                    shape = new Cylinder(transform, Colour{r,g,b}, ac, dc, sc);
                    shapes->emplace_back(shape);
                    break;
                //Cone
                case 'c':
                    shape = new Cone(transform, Colour{r,g,b}, ac, dc, sc);
                    shapes->emplace_back(shape);
                    break;
            }
        }
        input.close();
    }
    else
    {
        cout << "Unable to open file\n";
        exit(EXIT_FAILURE);
    }
    
    //Set up screen and frame buffer
    *screen = Screen_t(H, ASPECT);
    
    //Set up camera
    *camera = Camera(*E, G, UP, *N, F, THETA, *screen);
    
    G.free_members();
    UP.free_members();
    
    //Find screen size in viewing coordinates
    *height = *N * tan(M_PI/180.0 * THETA/2.0);
    *width = *height * ASPECT;
}

//  Simple non-recursive ray tracer that reads scene setup from a text file and produces output
//  on screen. Allows user to navigate through scene with arrows
//
int main()
{
    //X11 parameters
    Display* d;
    Window w;
    XEvent e;
    
    //Setup scene
    double width, height, N;
    Camera camera;
    vector<Light> lights;
    Screen_t screen;
    vector<GenericShape*> shapes;
    Matrix E;
    
    setup_scene(&camera, &shapes, &lights, &E, &screen, &width, &height, &N);
    
    //Initialize framebuffer
    int framebuffer[screen.width()][screen.height()][N_CHANNELS];
    
    Matrix direction = Matrix(4,1); //Direction of casted ray
    //For every pixel, cast a ray into the scene
    for (int i = 0; i < screen.width(); ++i) {
        for (int j = 0; j < screen.height(); ++j) {
            
            //Compute direction component of ray
            direction = (-N * camera.n) + (width * (2.0 * i/screen.width() - 1.0) * camera.u) + (height * (2.0 * j/screen.height() - 1.0) * camera.v);
            direction.member[3][0] = 0.0;
            
            double best_t = -1;                         //Keeps track of smallest t value found
            GenericShape* intersect_shape = nullptr;   //Keeps track of intersected shape associated with best_t
            //Iterate through each shape in the scene & find intersections between them & current ray
            for (vector<GenericShape*>::iterator it = shapes.begin(); it != shapes.end(); ++it)
            {
                //Determine if ray intersects with current shape
                double t = (*it)->intersect(E, direction);
                
                //If so determine if closest intersection yet. If so save shape
                if (t > -1 && (best_t == -1 || t < best_t))
                {
                    best_t = t;
                    intersect_shape = *it;
                }
            }
            
            //If an intersection was found determine pixel colour
            if (best_t != -1)
            {
                framebuffer[i][j][RED] = intersect_shape->colour.r;
                framebuffer[i][j][GREEN] = intersect_shape->colour.g;
                framebuffer[i][j][BLUE] = intersect_shape->colour.b;
                
                for (vector<Light>::iterator light = lights.begin(); light != lights.end(); ++light)
                {
                    bool shadow = false;
                    //Find intersection point
                    Matrix intersect_point = E + (direction * best_t);
                    //Find vector from intersection point to light source
                    Matrix shadow_direction = light->light_position - intersect_point;
                    
                    //For every shape cast a shadow ray
                    for (vector<GenericShape*>::iterator it = shapes.begin(); it != shapes.end(); ++it)
                    {
                        //If this ray intersects any shape then the current point is in shadow
                        double t = (*it)->intersect(intersect_point, shadow_direction);
                        if (t >= 0.01 && t <= 0.99)
                        {
                            shadow = true;
                            break;
                        }
                    }
                
                    //Add pixel colour to framebuffer
                    //If there is shadow, pixel intensity is computed with ambient light only
                    if (shadow)
                    {
                        framebuffer[i][j][RED] *= intersect_shape->ambient_coeff * I_a;
                        framebuffer[i][j][GREEN] *= intersect_shape->ambient_coeff * I_a;
                        framebuffer[i][j][BLUE] *= intersect_shape->ambient_coeff * I_a;
                    }
                    //Else it is given by ambient, diffuse, and specular components
                    else
                    {
                        framebuffer[i][j][RED] *= light_intensity(intersect_shape, intersect_point, shadow_direction, direction);
                        framebuffer[i][j][GREEN] *= light_intensity(intersect_shape, intersect_point, shadow_direction, direction);
                        framebuffer[i][j][BLUE] *= light_intensity(intersect_shape, intersect_point, shadow_direction, direction);
                    }
                }
            }
        }
    }
    
    //Print framebuffer to screen
    int s, r, g, b;
    d = InitX(d, &w, &s, screen);
    while (1) {
        XNextEvent(d, &e);
        if (e.type == Expose) {
            for (int i = 0; i < screen.width(); ++i) {
                for (int j = 0; j < screen.height(); ++j)
                {
                    r = framebuffer[i][j][RED];
                    g = framebuffer[i][j][GREEN];
                    b = framebuffer[i][j][BLUE];
                    SetCurrentColorX(d, &(DefaultGC(d, s)), r, g, b);
                    SetPixelX(d, w, s, i, j);
                }
            }
        }
        if(e.type == KeyPress)
            break;
        if(e.type == ClientMessage)
            break;
    }
    
    //Destroy window
    QuitX(d,w);
}