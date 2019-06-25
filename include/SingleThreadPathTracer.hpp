#pragma once

#include <Math.hpp>
#include <Random.hpp>
#include <Collision.hpp>
#include <Globals.hpp>
#include <IOHelpers.hpp>

math::Vec4 TraceAndSampleColor(math::Vec4 direction, math::Vec4 origin, uint32_t bounceCount);

inline math::Vec4 SampleColorSkybox(math::Vec4 direction, math::Vec4 origin, uint32_t bounceCount, uint32_t sphereIndex)
{
    return initColor * (direction.xyzw[1] + 1.f) * 0.5f;
}

inline math::Vec4 SampleColorSkybox(math::Vec4 direction, math::Vec4 origin)
{
    return initColor * (direction.xyzw[1] + 1.f) * 0.5f;
}

inline math::Vec4 SampleColorDiffuse(math::Vec4 direction, math::Vec4 origin, uint32_t bounceCount, uint32_t sphereIndex)
{
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

    return sampleColor;
}

inline math::Vec4 SampleColorReflective(math::Vec4 direction, math::Vec4 origin, uint32_t bounceCount, uint32_t sphereIndex)
{
    origin = cd::CalculateRaySphereClosestContactPoint(g_spheres[sphereIndex], g_radii[sphereIndex], origin, direction);
    math::Vec4 normal = cd::CalculateRaySphereContactNormal(origin, g_spheres[sphereIndex]);
    direction = math::Normalize(math::Reflect(direction, normal) + random::GenerateNormalDistInsideSphereVector() * g_diffuses[sphereIndex]);

    return TraceAndSampleColor(direction, origin, bounceCount);
}

inline math::Vec4 SampleColorRefractive(math::Vec4 direction, math::Vec4 origin, uint32_t bounceCount, uint32_t sphereIndex)
{
    float constexpr nAir = 1.0f;
    float constexpr nGlass = 1.5f;

    origin = cd::CalculateRaySphereClosestContactPoint(g_spheres[sphereIndex], g_radii[sphereIndex], origin, direction);
    math::Vec4 n = cd::CalculateRaySphereContactNormal(origin, g_spheres[sphereIndex]);

    float c = math::Dot(-n, direction);
    float r = nAir / nGlass;
    float rSq = pow((nAir - nGlass) / (nAir + nGlass), 2.f);
    float schlick = rSq + (1.f - rSq) * pow(1.f - c, 5.f);

    if (random::GenerateUniformRealDist(0, 1) < schlick)
    {
        return TraceAndSampleColor(math::Reflect(direction, n), origin, bounceCount);
    }

    if (r * sqrt(1.f - c * c) < 1.f)
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
            return TraceAndSampleColor(math::Reflect(direction, n), origin, bounceCount);
        }

        if (r * sqrt(1.f - c * c) < 1.f)
        {
            direction = math::Normalize(direction * r + n * (r * c - sqrt(1.f - r * r * (1.f - c * c))));
            return TraceAndSampleColor(direction, origin, bounceCount);
        }

        return TraceAndSampleColor(math::Reflect(direction, n), origin, bounceCount);
    }

    return TraceAndSampleColor(math::Reflect(direction, n), origin, bounceCount);
}

inline math::Vec4 TraceAndSampleColor(math::Vec4 direction, math::Vec4 origin, uint32_t bounceCount)
{
    uint32_t const sphereIndex = cd::FindClosestIntersectionSphere(direction, origin);

    if (sphereIndex < g_sphereNumber)
    {
        switch (g_materials[sphereIndex])
        {
        case Material::DIFFUSE:
            return SampleColorDiffuse(direction, origin, bounceCount, sphereIndex);
        case Material::REFLECTIVE:
            return SampleColorReflective(direction, origin, bounceCount, sphereIndex);
        case Material::REFRACTIVE:
            return SampleColorRefractive(direction, origin, bounceCount, sphereIndex);
        }
    }

    return SampleColorSkybox(direction, origin, bounceCount, sphereIndex);
}

void RenderSegment(RenderSegmentData segment)
{
    for (uint32_t y = segment.yBegin; y < segment.yEnd; ++y)
    {
        for (uint32_t x = segment.xBegin; x < segment.xEnd; ++x)
        {
            uint32_t const index = g_size - ((g_width - x) * g_stride + y * g_width * g_stride);
            math::Vec4 pixelColor{0, 0, 0, 0};

            for (uint32_t s = 0; s < g_samples; ++s)
            {
                float const u = static_cast<float>(y + random::GenerateUniformRealDist()) / g_width;
                float const v = static_cast<float>(x + random::GenerateUniformRealDist()) / g_height;
                math::Vec4 const primeRayDirection = math::Normalize(viewMatrix * math::Vec4{-1.f + 2.f * v, -1.f + 2.f * u, 1.f});
                math::Vec4 const primeRayOrigin{eyePos};

                pixelColor += TraceAndSampleColor(primeRayDirection, primeRayOrigin, g_bounces);
            }

            pixelColor *= (1.f / static_cast<float>(g_samples));
            io::WritePixel(index, pixelColor);
        }
    }
}