//
// Created by guizo on 21/01/2024.
//

#ifndef FLUID_FLUIDCUBE_H
#define FLUID_FLUIDCUBE_H

#define N 128
#define heightScreen 512
#define widthScreen 512
//#define IX(x, y) x + (y * N)

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <cstring>

class FluidCube {
public:
    float dt;
    float diff;
    float visc;

    float s[N*N];
    float density[N*N];

    float Vx[N*N], Vy[N*N];
    float Vx0[N*N],Vy0[N*N];

    FluidCube(int diffusion, int viscosity, float dt);
    void AddDensity(int x, int y, float amount);
    void AddVelocity(int x, int y, float amountx, float amounty);
    void Step();
    void Render();

private:

};


#endif //FLUID_FLUIDCUBE_H
