# Copyright (C) 2017-2018 by Godlike
# This code is licensed under the MIT license (MIT)
# (http://opensource.org/licenses/MIT)

cmake_minimum_required (VERSION 3.1)

project(Ray)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

set(RAY_ROOT
    ${PROJECT_SOURCE_DIR}
)

set(RAY_INCLUDE_DIR
    ${RAY_ROOT}/include
)

include_directories (
    ${RAY_INCLUDE_DIR}
)

#dependencies (glfw, glbindings, imgui)
add_subdirectory(dependencies)

#IMGUI
include(ImguiConfig)
include_directories(${IMGUI_INCLUDE_DIR})

#STB
include(StbConfig)
include_directories(${STB_INCLUDE_DIR})

set(DEMO_HEADERS
    include/Collision.hpp
    include/Definitions.hpp
    include/Gl.hpp
    include/Globals.hpp
    include/IOHelpers.hpp
    include/Math.hpp
    include/Random.hpp
    include/Renderer.hpp
    include/SceneGenerators.hpp
    include/SingleThreadPathTracer.hpp
    include/TaskBasedPathTracer.hpp
)
set(DEMO_SOURCES 
    Main.cpp
)

add_executable (${PROJECT_NAME}
    ${DEMO_SOURCES}
    ${DEMO_HEADERS}
)

target_link_libraries(${PROJECT_NAME}
    Ray_dependencies
)

target_compile_features(${PROJECT_NAME}
    PRIVATE cxx_std_17
)
