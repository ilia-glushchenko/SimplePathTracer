#pragma once

#include <Math.hpp>
#include <Globals.hpp>

namespace cd
{

inline bool RaySphereIntersection(
    math::Vec4 sphereCenter, float sphereRadius, math::Vec4 rayDirection, math::Vec4 rayOrigin, float threshold = 1e-3f)
{
    sphereCenter -= rayOrigin;
    float const tCenter = Dot(sphereCenter, rayDirection);
    float const distanceSquare = Dot(sphereCenter, sphereCenter) - tCenter * tCenter;

    return (tCenter > threshold) && (sphereRadius * sphereRadius - distanceSquare > threshold);
}

inline float CalculateRaySphereMinIntersectionFactor(
    math::Vec4 raySphere, float sphereRadius, math::Vec4 rayDirection)
{
    float const tCenter = Dot(raySphere, rayDirection);
    float const distanceSquare = Dot(raySphere, raySphere) - tCenter * tCenter;
    float const tDelta = std::sqrtf(sphereRadius * sphereRadius - distanceSquare);

    return tCenter - tDelta;
};

inline float CalculateRaySphereMaxIntersectionFactor(
    math::Vec4 raySphere, float sphereRadius, math::Vec4 rayDirection)
{
    float const tCenter = Dot(raySphere, rayDirection);
    float const distanceSquare = Dot(raySphere, raySphere) - tCenter * tCenter;
    float const tDelta = std::sqrtf(sphereRadius * sphereRadius - distanceSquare);

    return tCenter + tDelta;
};

inline math::Vec4 CalculateRaySphereIntersectionFactors(
    math::Vec4 raySphere, float sphereRadius, math::Vec4 rayDirection)
{
    float const tCenter = Dot(raySphere, rayDirection);
    float const distanceSquare = Dot(raySphere, raySphere) - tCenter * tCenter;
    float const tDelta = std::sqrtf(sphereRadius * sphereRadius - distanceSquare);

    return {tCenter - tDelta, tCenter + tDelta};
};

inline math::Vec4 CalculateRaySphereClosestContactPoint(
    math::Vec4 sphereCenter, float radius, math::Vec4 rayOrigin, math::Vec4 rayDirection)
{
    float const rayFactor = CalculateRaySphereMinIntersectionFactor(
        sphereCenter - rayOrigin, radius, rayDirection);

    return rayOrigin + rayDirection * rayFactor;
}

inline math::Vec4 CalculateRaySphereFarthestContactPoint(
    math::Vec4 sphereCenter, float radius, math::Vec4 rayOrigin, math::Vec4 rayDirection)
{
    float const rayFactor = CalculateRaySphereMaxIntersectionFactor(
        sphereCenter - rayOrigin, radius, rayDirection);

    return rayOrigin + rayDirection * rayFactor;
}

inline math::Vec4 CalculateRaySphereContactNormal(
    math::Vec4 contactPoint, math::Vec4 sphereCenter)
{
    return Normalize(contactPoint - sphereCenter);
}

inline bool RayPlaneIntersection(
    math::Vec4 planeNormal, math::Vec4 planePoint, math::Vec4 rayOrigin, math::Vec4 rayDirection)
{
    return Dot(planeNormal, rayDirection) < 0.0 && Dot(planeNormal, planePoint) < Dot(planeNormal, rayOrigin);
}

inline math::Vec4 CalculateRayPlaneContactPoint(
    math::Vec4 planeNormal, math::Vec4 planePoint, math::Vec4 rayOrigin, math::Vec4 rayDirection)
{
    float const rayPlaneProjection = Dot(planeNormal, rayDirection);
    float const t = Dot(planeNormal, planePoint - rayOrigin) / rayPlaneProjection;
    return rayOrigin + rayDirection * t;
}

inline uint8_t FindClosestIntersectionSphere(math::Vec4 primeRayDirection, math::Vec4 primeRayOrigin)
{
    uint8_t minIndex = g_sphereNumber;
    float minDistanceSq = FLT_MAX;

    for (uint8_t index = 0; index < g_sphereNumber; ++index)
    {
        if (cd::RaySphereIntersection(g_spheres[index], g_radii[index], primeRayDirection, primeRayOrigin))
        {
            math::Vec4 const intersectionPoint = cd::CalculateRaySphereClosestContactPoint(
                g_spheres[index], g_radii[index], primeRayOrigin, primeRayDirection);

            if (math::Dot(primeRayOrigin, primeRayDirection) < math::Dot(intersectionPoint, primeRayDirection))
            {
                float const distanceSq = math::LengthSquared(primeRayOrigin - intersectionPoint);
                minIndex = minDistanceSq > distanceSq ? index : minIndex;
                minDistanceSq = minDistanceSq > distanceSq ? distanceSq : minDistanceSq;
            }
        }
    }

    return minIndex;
}

} // namespace cd