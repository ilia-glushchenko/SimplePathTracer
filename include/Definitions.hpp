#pragma once

#include <Math.hpp>
#include <cstdint>
#include <vector>

enum class Material : uint8_t
{
    SKYBOX,
    REFLECTIVE,
    REFRACTIVE,
    DIFFUSE,
};

struct RenderSegmentData
{
    uint32_t yBegin = 0;
    uint32_t yEnd = 0;
    uint32_t xBegin = 0;
    uint32_t xEnd = 0;
};

struct RaySampleTask
{
    math::Vec4 origin;
    math::Vec4 direction;
    uint32_t sphereIndex;
    uint32_t bounceCount;
    uint32_t x;
    uint32_t y;
};

struct Tasks
{
    std::vector<RaySampleTask> diffuseTasks;
    std::vector<RaySampleTask> reflectiveTasks;
    std::vector<RaySampleTask> refractiveTasks;
    std::vector<RaySampleTask> skyboxTasks;
};

enum class eGenerationType : uint8_t
{
    SINGLE_THREAD,
    MULTI_THREAD,
    TASK_MULTI_THREAD
};

enum class eSceneGenerator : uint8_t
{
    REFERENCE,
    RANDOM
};