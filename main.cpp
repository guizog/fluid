#include "FluidCube.h"

GLFWwindow* window;

void error_callback(int error, const char* description)
{
    fprintf(stderr, "Error: %s\n", description);
}

int InitGLFW(){
    /* Initialize the library */
    if (!glfwInit())
        return -1;

    /* Create a windowed mode window and its OpenGL context */
    window = glfwCreateWindow(widthScreen, heightScreen, "Hello World", nullptr, nullptr);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    /* Make the window's context current */
    glfwMakeContextCurrent(window);

    if(glewInit() != GLEW_OK)
        std::cout << "GLEW INIT failed" << std::endl;

    glOrtho(0, 512, 512, 0, -1, 1);
    return 1;
}

int main() {

    InitGLFW();

    FluidCube* fluid = new FluidCube(5, 0, 0.00001);

    while (!glfwWindowShouldClose(window))
    {
        fluid->AddDensity(64, 64, 100);
        fluid->AddVelocity(64, 64, 100000, 1000);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        fluid->Render();
        fluid->Step();
        glfwSwapBuffers(window);
    }

    glfwTerminate();
    exit(EXIT_SUCCESS);
}
