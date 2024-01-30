//
// Created by guizo on 21/01/2024.
//

#include "FluidCube.h"

extern GLFWwindow* window;

int IX(int x, int y) {
    if (x > N * N) x = N * N;
    if (y > N * N) y = N * N;
    if (x < 0) x = 0;
    if (y < 0) y = 0;

    return x + (y * N);
}

void set_bnd(int b, float *x, int n) {
    for (int j = 1; j < n - 1; j++) {
        for (int i = 1; i < n - 1; i++) {
            x[IX(i, j)] = ((b == 3) ? -x[IX(i, j)] : x[IX(i, j)]);
        }
    }

    x[IX(0, 0)] = 0.33f * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, n - 1)] = 0.33f * (x[IX(1, n - 1)] + x[IX(0, n - 2)]);
    x[IX(n - 1, 0)] = 0.33f * (x[IX(n - 2, 0)] + x[IX(n - 1, 1)]);
    x[IX(n - 1, n - 2)] = 0.33f * (x[IX(n - 2, n - 1)] + x[IX(n - 1, n - 2)]);
}

void lin_solve(int b, float *x, float *x0, float a, float c, int iter, int n) {
    float cRecip = 1.0 / c;

    for (int k = 0; k < iter; k++) {
        for (int j = 1; j < n - 1; j++) {
            for (int i = 1; i < n - 1; i++) {
                x[IX(i, j)] =
                        (x0[IX(i, j)] +
                         a *
                         (x[IX(i + 1, j)] +
                          x[IX(i - 1, j)] +
                          x[IX(i, j + 1)] +
                          x[IX(i, j - 1)])) *
                        cRecip;
            }
        }
        set_bnd(b, x, n);
    }
    //x[IX(i, j)] =
    //                        (x0[IX(i, j)] + a * (x[IX(i + 1, j)] + x[IX(i - 1, j)] + x[IX(i, j + 1)] + x[IX(i, j - 1)])) *
    //                        cRecip;
}

void diffuse(int b, float *x, float *x0, float diff, float dt, int iter, int n) {
    float a = dt * diff * (n - 2) * (n - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a, iter, n);
}

void project(float *velocX, float *velocY, float *p, float *div, int iter, int n) {
    for (int j = 1; j < n - 1; j++) {
        for (int i = 1; i < n - 1; i++) {
            div[IX(i, j)] =
                    (-0.5 *
                     (velocX[IX(i + 1, j)] -
                      velocX[IX(i - 1, j)] +
                      velocY[IX(i, j + 1)] -
                      velocY[IX(i, j - 1)])) /
                    N;
            p[IX(i, j)] = 0;

            //{
            //            div[IX(i, j)] = -0.5f * (velocX[IX(i + 1, j)] - velocX[IX(i - 1, j)] + velocY[IX(i, j + 1)] -
            //                                     velocY[IX(i, j - 1)]) / n;
            //            p[IX(i, j)] = 0;
            //        }
        }
    }
    set_bnd(0, div, n);
    set_bnd(0, p, n);
    lin_solve(0, p, div,1, 6, iter, n);

    for (int j = 1; j < n - 1; j++) {
        for (int i = 1; i < n - 1; i++) {
            velocX[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * N;
            velocY[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * N;
            //velocX[IX(i, j)] -= 0.5f * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * n;
            //            velocY[IX(i, j)] -= 0.5f * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * n;
        }
    }
    set_bnd(1, velocX, n);
    set_bnd(2, velocY, n);
}

void advect(int b, float *d, float *d0, float *velocX, float *velocY, float dt, int n) {
    float i0, i1, j0, j1;

    float dtx = dt * (n - 2);
    float dty = dt * (n - 2);

    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;

    float nfloat = n;
    float ifloat, jfloat;
    int i, j;

    for (j = 1, jfloat = 1; j < n - 1; j++, jfloat++) {
        for (i = 1, ifloat = 1; i < n - 2; i++, ifloat++) {
            tmp1 = dtx * velocX[IX(i, j)];
            tmp2 = dty * velocY[IX(i, j)];
            x = ifloat - tmp1;
            y = jfloat - tmp2;

            if (x < 0.5f)
                x = 0.5f;
            if (x > nfloat + 0.5f)
                x = nfloat + 0.5f;
            i0 = floor(x);
            i1 = i0 + 1.0f;

            if (y < 0.5f)
                y = 0.5f;
            if (y > nfloat + 0.5f)
                y = nfloat + 0.5f;
            j0 = floor(y);
            j1 = j0 + 1.0f;

            s1 = x - i0;
            s0 = 1.0f - s1;
            t1 = y - j0;
            t0 = 1.0f - t1;

            int i0i = i0;
            int i1i = i1;
            int j0i = j0;
            int j1i = j1;

            d[IX(i, j)] =
                    s0 * (t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)]) +
                    s1 * (t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]);

            //d[IX(i, j)] = s0 * (t0 * d0[IX(i0i, j0i)]
            //                                + (t1 * d0[IX(i0i, j1i)]))
            //                          + s1 * (t0 * d0[IX(i1i, j0i)]
            //                                  + (t1 * d0[IX(i1i, j1i)]
            //                                  ));

        }
    }
    set_bnd(b, d, n);
}


FluidCube::FluidCube(int diffusion, int viscosity, float dt) {
    this->size = N;
    this->dt = dt;

    this->diff = diffusion;
    this->visc = viscosity;
    memset(this->density, 0, (N * N) * sizeof(*density));
    memset(this->Vx, 0, (N * N) * sizeof(*Vx));
    memset(this->Vy, 0, (N * N) * sizeof(*Vy));
    memset(this->s, 0, (N * N) * sizeof(*s));
    memset(this->Vx0, 0, (N * N) * sizeof(*Vx0));
    memset(this->Vy0, 0, (N * N) * sizeof(*Vy0));
}

void FluidCube::AddDensity(int x, int y, float amount) {
    this->density[IX(x, y)] += amount;
}

void FluidCube::AddVelocity(int x, int y, float amountx, float amounty) {
    int index = IX(x, y);

    this->Vx[IX(x, y)] += amountx;
    this->Vy[IX(x, y)] += amounty;
}

void FluidCube::Step() {
    int n = this->size;
    float visc = this->visc;
    float diff = this->diff;
    float dt = this->dt;
    float *Vx = this->Vx;
    float *Vy = this->Vy;
    float *Vx0 = this->Vx0;
    float *Vy0 = this->Vy0;
    float *s = this->s;
    float *density = this->density;

    std::cout << "step" << std::endl;

    diffuse(1, Vx0, Vx, visc, dt, 4, n);
    diffuse(2, Vy0, Vy, visc, dt, 4, n);

    project(Vx0, Vy0, Vx, Vy, 4, n);

    advect(1, Vx, Vx0, Vx0, Vy0, dt, n);
    advect(2, Vy, Vy0, Vx0, Vy0, dt, n);

    project(Vx, Vy, Vx0, Vy0, 4, n);

    diffuse(0, s, density, diff, dt, 4, n);
    advect(0, density, s, Vx, Vy, dt, n);

}

void FluidCube::Render(){
    for(int x = 0; x < N; x++){
        for(int y = 0; y < N; y++){
            int cube = heightScreen / N;
            int size = heightScreen / N;
            int xo = x * cube;
            int yo = y * cube;

            float dens = (this->density[IX(x,y)]) * 100; //density not being calculated

            //glColor3f(1 * dens, 1 * dens,  1 * dens);
            glColor3ub(255 * dens, 255 * dens,  255 * dens);
            glBegin(GL_QUADS);
            glVertex2i(xo, yo);
            glVertex2i(xo + cube, yo);
            glVertex2i(xo + cube, yo + cube);
            glVertex2i(xo, yo + cube);
            glEnd();

        }
    }
    std::cout << "#### render" << std::endl;
}


