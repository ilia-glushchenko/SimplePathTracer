#pragma once

#include <Math.hpp>
#include <Definitions.hpp>
#include <Gl.hpp>
#include <string>

eGenerationType constexpr g_generationType = eGenerationType::TASK_MULTI_THREAD;
eSceneGenerator constexpr g_sceneGenerator = eSceneGenerator::REFERENCE;

uint32_t constexpr g_maxThreads = 4;
uint32_t constexpr g_bounces = 10;
uint32_t constexpr g_samples = 100;
uint32_t constexpr g_width = 1440;
uint32_t constexpr g_height = 1440;
float constexpr g_ratio = static_cast<float>(g_width) / g_height;
uint8_t constexpr g_stride = 3;
uint32_t constexpr g_size = g_width * g_height * g_stride;
uint8_t *const g_data = (uint8_t *)std::malloc(sizeof(uint8_t) * g_size);

math::Mat4 viewMatrix = math::CreateIdentityMatrix();
math::Vec4 lookAt = {0, 1, 0};
math::Vec4 eyePos = {0, 1, -3};
math::Vec4 constexpr upDir = {0, 1, 0};
math::Vec4 constexpr lightPos = {5, -10, 1};
math::Vec4 constexpr planeNormal = {0, 1, 0};
math::Vec4 constexpr planePoint = {0, -0.5, 0};
math::Vec4 constexpr planeColor = {246, 219, 219};
math::Vec4 constexpr initColor = {137, 207, 240};

std::vector<math::Vec4> g_colors;
std::vector<math::Vec4> g_spheres;
std::vector<float> g_radii;
std::vector<Material> g_materials;
std::vector<math::Vec4> g_attenuations;
std::vector<float> g_diffuses;
uint32_t g_sphereNumber = 10;

GLuint const g_windowWidth = g_width;
GLuint const g_windowHeight = g_height;
int g_frameBufferWidth;
int g_frameBufferHeight;
GLFWwindow *g_window;

GLuint g_program = 0;

GLint g_modelUniform = {};
GLint g_fragColorUniform = {};
GLint g_fragTextureUniform = {};

GLuint g_vertexBufferObject = {};
GLuint g_elementBufferObject = {};
GLuint g_vertexArrayObject = {};

GLchar const *g_texturePath = "textures/TexturesCom_StoneWall2_512_albedo.png";
GLuint g_textureHandle = {};

std::string g_vertexShaderSource;
std::string g_fragmentShaderSource;

GLfloat g_vertices[6 * 8] = {
    //positions          //texture
    -1.0f,
    1.0f,
    0.0f,
    1.0f,
    0.0f,
    1.0f,
    -1.0f,
    0.0f,
    0.0f,
    1.0f,
    -1.0f,
    -1.0f,
    0.0f,
    1.0f,
    1.0f,
    1.0f,
    1.0f,
    0.0f,
    0.0f,
    0.0f,
};

GLuint g_indices[6] = {0, 1, 2, 0, 1, 3};
