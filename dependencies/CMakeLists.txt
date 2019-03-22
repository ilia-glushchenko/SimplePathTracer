# Copyright (C) 2018 by Godlike
# This code is licensed under the MIT license (MIT)
# (http://opensource.org/licenses/MIT)

cmake_minimum_required (VERSION 3.1)

project(Ray_dependencies)

#GLBINDING
include(GlbindingConfig)
add_subdirectory(${GLBINDING_SOURCE_DIR} ${GLBINDING_BINARY_DIR})
include_directories(${GLBINDING_INCLUDE_DIR})

#GLFW
include(GlfwConfig)
add_subdirectory(${GLFW_SOURCE_DIR} ${GLFW_BINARY_DIR})
include_directories(${GLFW_INCLUDE_DIR})

#Imgui
include(ImguiConfig)
include_directories(
    ${PROJECT_SOURCE_DIR}
    ${IMGUI_INCLUDE_DIR}
    ${IMGUI_INCLUDE_DIR}/examples
)

add_library(${PROJECT_NAME} STATIC
    ${IMGUI_SOURCE_DIR}/examples/imgui_impl_glfw.cpp
    ${IMGUI_SOURCE_DIR}/examples/imgui_impl_opengl3.cpp
    ${IMGUI_SOURCE_DIR}/imgui.cpp
    ${IMGUI_SOURCE_DIR}/imgui_demo.cpp
    ${IMGUI_SOURCE_DIR}/imgui_draw.cpp
    ${IMGUI_SOURCE_DIR}/imgui_widgets.cpp
)

target_link_libraries(${PROJECT_NAME}
    ${OPENGL_LIBRARIES}
    ${GLBINDING_LIB}
    ${GLFW_LIB}
)
