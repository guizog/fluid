#include "FluidCube.h"
#define dt 0.00001

GLFWwindow* window;
FluidCube* fluid = new FluidCube(5, 0, dt);

int densAmount = 100;
int amountX = 100000, amountY = 1000;

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

static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos){
    //std::cout << "x: " << xpos << "   y: " << ypos << std::endl;
    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS)
    {
        fluid->AddDensity(xpos/4, ypos/4, densAmount);
        fluid->AddVelocity(xpos/4, ypos/4, amountX, amountY);
    }
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{

}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)//did not work as expected lol
{
    if(action == GLFW_PRESS){
        if(key == GLFW_KEY_UP){
            densAmount += 100;
            std::cout << "density increased to " << densAmount << std::endl;
        }
        if(key == GLFW_KEY_DOWN){
            densAmount -= 100;
            std::cout << "density decreased to " << densAmount << std::endl;
        }
        if(key == GLFW_KEY_LEFT){
            amountY += 5000;
        }
        if(key == GLFW_KEY_RIGHT){
            amountY -= 5000;
        }

    }
}

int main() {

    InitGLFW();

    while (!glfwWindowShouldClose(window))
    {
        fluid->AddDensity(64, 64, densAmount);
        fluid->AddVelocity(64, 64, amountX, amountY);

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        fluid->Render();
        fluid->Step();
        glfwSwapBuffers(window);

        glfwPollEvents();
        //glfwSetMouseButtonCallback(window, mouse_button_callback);
        glfwSetCursorPosCallback(window, cursor_position_callback);
        glfwSetKeyCallback(window, key_callback);
    }

    glfwTerminate();
    exit(EXIT_SUCCESS);
}
