#pragma once

#include <Definitions.hpp>
#include <Random.hpp>
#include <Collision.hpp>
#include <IOHelpers.hpp>
#include <Globals.hpp>

inline void TraceAndSampleColor(RaySampleTask task, Tasks &tasks)
{
    task.sphereIndex = cd::FindClosestIntersectionSphere(task.direction, task.origin);

    if (task.sphereIndex < g_sphereNumber)
    {
        switch (g_materials[task.sphereIndex])
        {
        case Material::DIFFUSE:
            tasks.diffuseTasks.push_back(task);
            return;
        case Material::REFLECTIVE:
            tasks.reflectiveTasks.push_back(task);
            return;
        case Material::REFRACTIVE:
            tasks.refractiveTasks.push_back(task);
            return;
        }
    }

    tasks.skyboxTasks.push_back(task);
}

inline void SwapTasks(Tasks &a, Tasks &b)
{
    a.diffuseTasks.swap(b.diffuseTasks);
    a.reflectiveTasks.swap(b.reflectiveTasks);
    a.refractiveTasks.swap(b.refractiveTasks);
    a.skyboxTasks.swap(b.skyboxTasks);
}

inline void ClearTasks(Tasks &tasks)
{
    tasks.diffuseTasks.clear();
    tasks.reflectiveTasks.clear();
    tasks.refractiveTasks.clear();
    tasks.skyboxTasks.clear();
}

inline uint32_t TasksSize(Tasks &tasks)
{
    return static_cast<uint32_t>(
        tasks.diffuseTasks.size() + tasks.reflectiveTasks.size() + tasks.refractiveTasks.size() + tasks.skyboxTasks.size());
}

void RenderSegmentTask(RenderSegmentData segment)
{
    uint32_t const segmentWidth = segment.xEnd - segment.xBegin;
    uint32_t const segmentHeight = segment.yEnd - segment.yBegin;
    std::vector<math::Vec4> colors(segmentWidth * segmentHeight, {0, 0, 0, 0});
    std::vector<float> samples(segmentWidth * segmentHeight, 0);

    for (uint32_t s = 0; s < g_samples; ++s)
    {
        thread_local Tasks tasks;
        thread_local Tasks nextFrameTasks;
        ClearTasks(tasks);
        ClearTasks(tasks);

        for (uint32_t y = segment.yBegin; y < segment.yEnd; ++y)
        {
            for (uint32_t x = segment.xBegin; x < segment.xEnd; ++x)
            {
                float const u = static_cast<float>(y + random::GenerateUniformRealDist()) / g_width;
                float const v = static_cast<float>(x + random::GenerateUniformRealDist()) / g_height;

                TraceAndSampleColor(
                    RaySampleTask{eyePos, math::Normalize(viewMatrix * math::Vec4{-1.f + 2.f * v, -1.f + 2.f * u, 1.f}), g_sphereNumber, g_bounces, x, y},
                    tasks);
            }
        }

        for (uint32_t pass = 0; pass < 10 && TasksSize(tasks) > 0; ++pass)
        {
            for (uint32_t i = 0; i < tasks.diffuseTasks.size(); ++i)
            {
                uint32_t &bounceCount = tasks.diffuseTasks[i].bounceCount;
                uint32_t &sphereIndex = tasks.diffuseTasks[i].sphereIndex;
                math::Vec4 &direction = tasks.diffuseTasks[i].direction;
                math::Vec4 &origin = tasks.diffuseTasks[i].origin;

                math::Vec4 sampleColor = g_colors[sphereIndex] * 0.5f;
                origin = cd::CalculateRaySphereClosestContactPoint(g_spheres[sphereIndex], g_radii[sphereIndex], origin, direction);
                direction = math::Normalize(cd::CalculateRaySphereContactNormal(origin, g_spheres[sphereIndex]) + random::GenerateUniformDistInsideSphereVector());
                sphereIndex = cd::FindClosestIntersectionSphere(direction, origin);

                while (--bounceCount && sphereIndex < g_sphereNumber)
                {
                    sampleColor = sampleColor * 0.5f;
                    origin = cd::CalculateRaySphereClosestContactPoint(g_spheres[sphereIndex], g_radii[sphereIndex], origin, direction);
                    direction = math::Normalize(origin + cd::CalculateRaySphereContactNormal(origin, g_spheres[sphereIndex]) + random::GenerateUniformDistInsideSphereVector());
                    sphereIndex = cd::FindClosestIntersectionSphere(direction, origin);
                }

                uint32_t const colorIndex = (tasks.diffuseTasks[i].x - segment.xBegin) + (tasks.diffuseTasks[i].y - segment.yBegin) * segmentHeight;
                colors[colorIndex] += sampleColor;
                samples[colorIndex] += 1.f;
            }

            for (uint32_t i = 0; i < tasks.reflectiveTasks.size(); ++i)
            {
                math::Vec4 &direction = tasks.reflectiveTasks[i].direction;
                math::Vec4 &origin = tasks.reflectiveTasks[i].origin;
                uint32_t sphereIndex = tasks.reflectiveTasks[i].sphereIndex;

                origin = cd::CalculateRaySphereClosestContactPoint(g_spheres[sphereIndex], g_radii[sphereIndex], origin, direction);
                math::Vec4 normal = cd::CalculateRaySphereContactNormal(origin, g_spheres[sphereIndex]);
                direction = math::Normalize(math::Reflect(direction, normal) + random::GenerateNormalDistInsideSphereVector() * g_diffuses[sphereIndex]);
                sphereIndex = cd::FindClosestIntersectionSphere(direction, origin);

                TraceAndSampleColor(
                    RaySampleTask{origin, direction, sphereIndex, g_bounces, tasks.reflectiveTasks[i].x, tasks.reflectiveTasks[i].y},
                    nextFrameTasks);
            }

            for (uint32_t i = 0; i < tasks.refractiveTasks.size(); ++i)
            {
                float constexpr nAir = 1.0f;
                float constexpr nGlass = 1.5f;

                math::Vec4 &direction = tasks.refractiveTasks[i].direction;
                math::Vec4 &origin = tasks.refractiveTasks[i].origin;
                uint32_t sphereIndex = tasks.refractiveTasks[i].sphereIndex;
                origin = cd::CalculateRaySphereClosestContactPoint(g_spheres[sphereIndex], g_radii[sphereIndex], origin, direction);
                math::Vec4 n = cd::CalculateRaySphereContactNormal(origin, g_spheres[sphereIndex]);

                float c = math::Dot(-n, direction);
                float r = nAir / nGlass;
                float rSq = pow((nAir - nGlass) / (nAir + nGlass), 2.f);
                float schlick = rSq + (1.f - rSq) * pow(1.f - c, 5.f);

                //External reflection
                if (random::GenerateUniformRealDist(0, 1) < schlick)
                {
                    direction = math::Reflect(direction, n);
                }
                //Internal reflection
                else if (r * sqrt(1.f - c * c) < 1.f)
                {
                    direction = math::Normalize(direction * r + n * (r * c - sqrt(1.f - r * r * (1.f - c * c))));
                    origin = cd::CalculateRaySphereFarthestContactPoint(g_spheres[sphereIndex], g_radii[sphereIndex], origin, direction);
                    n = -cd::CalculateRaySphereContactNormal(origin, g_spheres[sphereIndex]);

                    c = math::Dot(-n, direction);
                    r = nGlass / nAir;

                    rSq = pow((nGlass - nAir) / (nGlass + nAir), 2.f);
                    schlick = rSq + (1.f - rSq) * pow(1.f - c, 5.f);

                    if (random::GenerateUniformRealDist(0, 1) < schlick)
                    {
                        direction = math::Reflect(direction, n);
                    }
                    else if (r * sqrt(1.f - c * c) < 1.f)
                    {
                        direction = math::Normalize(direction * r + n * (r * c - sqrt(1.f - r * r * (1.f - c * c))));
                    }
                    else
                    {
                        direction = math::Reflect(direction, n);
                    }
                }
                else
                {
                    direction = math::Reflect(direction, n);
                }

                sphereIndex = cd::FindClosestIntersectionSphere(direction, origin);
                TraceAndSampleColor(
                    RaySampleTask{origin, direction, sphereIndex, g_bounces, tasks.refractiveTasks[i].x, tasks.refractiveTasks[i].y},
                    nextFrameTasks);
            }

            for (uint32_t i = 0; i < tasks.skyboxTasks.size(); ++i)
            {
                math::Vec4 const sampleColor = initColor * (tasks.skyboxTasks[i].direction.xyzw[1] + 1.f) * 0.5f;

                uint32_t const colorIndex = (tasks.skyboxTasks[i].x - segment.xBegin) + (tasks.skyboxTasks[i].y - segment.yBegin) * segmentHeight;
                colors[colorIndex] += sampleColor;
                samples[colorIndex] += 1.f;
            }

            SwapTasks(tasks, nextFrameTasks);
            ClearTasks(nextFrameTasks);
        } // pass for loop
    }

    for (uint32_t i = 0; i < samples.size(); ++i)
    {
        colors[i] *= 1.f / samples[i];

        uint32_t const x = i % segmentWidth + segment.xBegin;
        uint32_t const y = static_cast<uint32_t>(floor(i / segmentWidth) + segment.yBegin);
        uint32_t const index = g_size - ((g_width - x) * g_stride + y * g_width * g_stride);

        io::WritePixel(index, colors[i]);
    }
}
