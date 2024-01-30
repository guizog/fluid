cmake_minimum_required(VERSION 3.25)
project(fluid)

set(CMAKE_CXX_STANDARD 17)

find_package(OpenGL REQUIRED)

set(GL_INCLUDE_DIR C:/libs/GL/include/)
set(GL_LIB_DIR C:/libs/GL/lib/Release/Win32)

include_directories(${GL_INCLUDE_DIR})
link_directories(${GL_LIB_DIR})

add_executable(fluid main.cpp
        FluidCube.cpp
        FluidCube.h)
target_link_libraries(${PROJECT_NAME} glfw3 opengl32 glew32)

