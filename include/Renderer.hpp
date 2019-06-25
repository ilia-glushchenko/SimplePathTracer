#pragma once

#include <Math.hpp>
#include <IOHelpers.hpp>
#include <Globals.hpp>
#include <SingleThreadPathTracer.hpp>
#include <TaskBasedPathTracer.hpp>
#include <SceneGenerators.hpp>
#include <GL.hpp>
#include <iostream>
#include <fstream>
#include <atomic>
#include <thread>
#include <condition_variable>

inline void KeyCallback(GLFWwindow *window, int key, int scancode, int action, int mode)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    {
        glfwSetWindowShouldClose(window, 1);
    }
}

inline void ProcessInput(GLFWwindow *window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

inline GLuint CompileShader(GLenum type, GLchar const *source)
{
    static int succes = 0;
    static char infoLog[512];

    GLuint shader_id = glCreateShader(type);
    glShaderSource(shader_id, 1, &source, nullptr);
    glCompileShader(shader_id);

    glGetShaderiv(shader_id, GL_COMPILE_STATUS, &succes);
    if (!succes)
    {
        glGetShaderInfoLog(shader_id, sizeof(infoLog), nullptr, infoLog);
        std::cout << "ERROR::SHADER::COMPILATION_FAILED\n"
                  << infoLog << std::endl;
        return 0;
    }

    return shader_id;
}

inline GLuint CreateShaderProgram(std::vector<GLuint> const &shader_ids)
{
    GLuint const program = glCreateProgram();
    for (auto shader_id : shader_ids)
    {
        glAttachShader(program, shader_id);
    }
    glLinkProgram(program);

    static int succes;
    static char infoLog[512];

    glGetProgramiv(program, GL_LINK_STATUS, &succes);
    if (!succes)
    {
        glGetProgramInfoLog(program, sizeof(infoLog), nullptr, infoLog);
        std::cout << "ERROR::SHADER_PROGRAM::LINKING_FAILED\n"
                  << infoLog << std::endl;
        return 0;
    }

    return program;
}

inline void InitializeShaderProgram()
{
    std::ifstream t1("shaders/shader.vert");
    std::ifstream t2("shaders/shader.frag");
    g_vertexShaderSource = std::string((std::istreambuf_iterator<char>(t1)),
                                       std::istreambuf_iterator<char>());
    g_fragmentShaderSource = std::string((std::istreambuf_iterator<char>(t2)),
                                         std::istreambuf_iterator<char>());

    GLuint const vs_handle = CompileShader(GL_VERTEX_SHADER, g_vertexShaderSource.c_str());
    GLuint const fs_handle = CompileShader(GL_FRAGMENT_SHADER, g_fragmentShaderSource.c_str());

    g_program = CreateShaderProgram({vs_handle, fs_handle});

    glDeleteShader(vs_handle);
    glDeleteShader(fs_handle);
}

inline void InitializeVertexBuffer()
{
    glGenVertexArrays(1, &g_vertexArrayObject);
    glGenBuffers(1, &g_elementBufferObject);
    glGenBuffers(1, &g_vertexBufferObject);

    glBindVertexArray(g_vertexArrayObject);

    glBindBuffer(GL_ARRAY_BUFFER, g_vertexBufferObject);
    glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertices), g_vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_elementBufferObject);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(g_indices), g_indices, GL_STATIC_DRAW);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 5,
                          reinterpret_cast<void *>(0));

    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 5,
                          reinterpret_cast<void *>(sizeof(GLfloat) * 3));
}

inline void InitializeUniforms()
{
    g_fragColorUniform = glGetUniformLocation(g_program, "fragColorShift");
    g_fragTextureUniform = glGetUniformLocation(g_program, "fragTexture");
    g_modelUniform = glGetUniformLocation(g_program, "model");

    if (g_fragColorUniform == -1)
    {
        std::cout << "ERROR::SHADER::UNIFORM_INIT_FAILED\n"
                  << std::endl;
    }
    if (g_fragTextureUniform == -1)
    {
        std::cout << "ERROR::SHADER::UNIFORM_INIT_FAILED\n"
                  << std::endl;
    }
    if (g_modelUniform == -1)
    {
        std::cout << "ERROR::SHADER::UNIFORM_INIT_FAILED\n"
                  << std::endl;
    }
}

inline void InitializeTexture()
{
    glUseProgram(g_program);
    glUniform1i(g_fragTextureUniform, 0);
    glGenTextures(1, &g_textureHandle);

    glBindTexture(GL_TEXTURE_2D, g_textureHandle);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, g_width, g_height, 0,
                 GL_RGB, GL_UNSIGNED_BYTE, g_data);
    glGenerateMipmap(GL_TEXTURE_2D);
}

inline void UpdateTexture()
{
    glBindTexture(GL_TEXTURE_2D, g_textureHandle);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, g_width, g_height, 0,
                 GL_RGB, GL_UNSIGNED_BYTE, g_data);
    glGenerateMipmap(GL_TEXTURE_2D);
}

inline void Startup()
{
    InitializeShaderProgram();
    InitializeUniforms();
    InitializeVertexBuffer();
    InitializeTexture();
}

inline void Shutdown()
{
    glDeleteVertexArrays(1, &g_vertexArrayObject);
    glDeleteBuffers(1, &g_vertexBufferObject);
    glDeleteBuffers(1, &g_elementBufferObject);
}

inline bool Init()
{
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, 0);

    g_window = glfwCreateWindow(g_windowWidth, g_windowHeight, "SimpleRaytracer", nullptr, nullptr);
    if (!g_window)
    {
        glfwTerminate();
        return false;
    }

    glfwMakeContextCurrent(g_window);
    glbinding::Binding::initialize(glfwGetProcAddress);
    glfwSetKeyCallback(g_window, KeyCallback);

    glfwGetFramebufferSize(g_window, &g_frameBufferWidth, &g_frameBufferHeight);
    glViewport(0, 0, g_frameBufferWidth, g_frameBufferHeight);

    return true;
}

inline void Deinit()
{
    glfwTerminate();
}

inline void Render()
{
    static float model[4 * 4] = {
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1};

    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    glUseProgram(g_program);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, g_textureHandle);
    glBindVertexArray(g_vertexArrayObject);
    glUniform4f(g_fragColorUniform, 0.0f, 0.0f, 0.0f, 0.0f);
    glUniformMatrix4fv(g_modelUniform, 1, GL_FALSE, model);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, nullptr);
}

inline RenderSegmentData MakeRenderSegmentData(uint32_t i, uint32_t j, uint32_t segmentWidth, uint32_t segmentHeigth)
{
    uint32_t const yBegin = segmentHeigth * j;
    uint32_t const yEnd = yBegin + segmentHeigth > g_height ? g_height : yBegin + segmentHeigth;
    uint32_t const xBegin = segmentWidth * i;
    uint32_t const xEnd = xBegin + segmentWidth > g_width ? g_width : xBegin + segmentWidth;

    return {yBegin, yEnd, xBegin, xEnd};
}

inline void RenderJob(RenderSegmentData segment, std::atomic<int> &freeThreadsCount, std::condition_variable &freeThreadWaitCondition)
{
    freeThreadsCount.fetch_sub(1);
    if (g_generationType == eGenerationType::MULTI_THREAD)
    {
        RenderSegmentTask(segment);
    }
    else if (g_generationType == eGenerationType::TASK_MULTI_THREAD)
    {
        RenderSegment(segment);
    }
    freeThreadsCount.fetch_add(1);
    freeThreadWaitCondition.notify_one();
}

inline void RenderImageParallelMain()
{
    auto start_time = std::chrono::high_resolution_clock::now();

    uint32_t const hardwearThreads = std::min(std::thread::hardware_concurrency() * 2, g_maxThreads);
    uint32_t const threadCount = (hardwearThreads % 2 != 0 ? hardwearThreads + 1 : hardwearThreads);

    uint32_t const segmentWidth = g_width / threadCount;
    uint32_t const segmentHeigth = g_height / threadCount;
    std::vector<RenderSegmentData> segments;

    for (uint32_t j = 0; j < threadCount; ++j)
    {
        for (uint32_t i = 0; i < threadCount; ++i)
        {
            segments.push_back({MakeRenderSegmentData(i, j, segmentWidth, segmentHeigth)});
        }
    }

    std::atomic<int> freeThreadsCount(static_cast<int>(threadCount));
    std::condition_variable freeThreadWaitCondition;
    std::mutex mutex;
    std::unique_lock<std::mutex> qlock(mutex);
    std::vector<std::thread> threads;

    for (auto segment : segments)
    {
        static int i = 0;
        std::cout << ++i << "/" << segments.size() << std::endl;

        std::thread thread(RenderJob, segment, std::ref(freeThreadsCount), std::ref(freeThreadWaitCondition));
        thread.detach();
        freeThreadWaitCondition.wait(qlock, [&freeThreadsCount]() -> bool { return freeThreadsCount.load() > 0; });
    }

    freeThreadWaitCondition.wait(qlock, [&freeThreadsCount, &threadCount]() -> bool { return freeThreadsCount.load() == threadCount; });

    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;

    std::cout << "Duration: " << std::chrono::duration_cast<std::chrono::milliseconds>(time).count() << " ms\n";
    std::cout << "Duration: " << std::chrono::duration_cast<std::chrono::seconds>(time).count() << " s\n";
    std::cout << "Duration: " << std::chrono::duration_cast<std::chrono::minutes>(time).count() << " m\n";

    io::SaveImage();
}

inline void RenderImage()
{
    RenderSegmentTask({0, g_height, 0, g_width});
    io::SaveImage();
}

inline void RenderImageParallel()
{
    std::thread thread(RenderImageParallelMain);
    thread.detach();
}

inline void MainLoop()
{
    static_assert(g_width % 2 == 0);
    static_assert(g_height % 2 == 0);

    viewMatrix = Transpose(math::CreateCameraBasisMatrix(eyePos, lookAt, upDir));

    switch (g_sceneGenerator)
    {
    case eSceneGenerator::REFERENCE:
        InitSpheres();
        break;
    case eSceneGenerator::RANDOM:
        GenerateSpheres();

    default:
        break;
    }

    switch (g_generationType)
    {
    case eGenerationType::SINGLE_THREAD:
        RenderImage();
        break;

    default:
        RenderImageParallel();
        break;
    }

    while (!glfwWindowShouldClose(g_window))
    {
        ProcessInput(g_window);

        Render();

        glfwSwapBuffers(g_window);
        glfwPollEvents();

        UpdateTexture();
    }
}

inline bool TracePaths()
{
    if (Init())
    {
        Startup();
        MainLoop();
        Shutdown();
        Deinit();

        return true;
    }

    return false;
}
