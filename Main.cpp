#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define NOMINMAX
#include <windows.h>
#include <intrin.h>

#include <iostream>
#include <atomic>
#include <condition_variable>
#include <thread>
#include <random>
#include <string>

namespace math
{
struct Vec4
{
    union {
        __m128 m128;
        float xyzw[4];
    };

    inline Vec4 operator-(Vec4 other) const
    {
        other.m128 = _mm_sub_ps(m128, other.m128);

        return other;
    }

    constexpr Vec4 operator-(Vec4 other)
    {
        other.xyzw[0] = xyzw[0] - other.xyzw[0];
        other.xyzw[1] = xyzw[1] - other.xyzw[1];
        other.xyzw[2] = xyzw[2] - other.xyzw[2];
        other.xyzw[3] = xyzw[3] - other.xyzw[3];

        return other;
    }

    inline Vec4 operator+(Vec4 other) const
    {
        other.m128 = _mm_add_ps(other.m128, m128);

        return other;
    }

    constexpr Vec4 operator+(Vec4 other)
    {
        other.xyzw[0] += xyzw[0];
        other.xyzw[1] += xyzw[1];
        other.xyzw[2] += xyzw[2];
        other.xyzw[3] += xyzw[3];

        return other;
    }

    inline Vec4 operator*(float s) const
    {
        return {
            _mm_mul_ps(m128, _mm_set1_ps(s))
        };
    }

    constexpr Vec4 operator*(float s)
    {
        return {
            xyzw[0] * s,
            xyzw[1] * s,
            xyzw[2] * s,
            xyzw[3] * s,
        };
    }

    inline Vec4 operator-() const
    {
        return {
            _mm_xor_ps(m128, _mm_set_ps1(-0.0))
        };
    }

    constexpr Vec4 operator-()
    {
        return {
            -xyzw[0], -xyzw[1], -xyzw[2], -xyzw[3],
        };
    }

    inline Vec4& operator-=(Vec4 other)
    {
        m128 = _mm_sub_ps(m128, other.m128);
        return *this;
    }

    inline Vec4& operator+=(Vec4 other)
    {
        m128 = _mm_add_ps(m128, other.m128);
        return *this;
    }

    inline Vec4& operator*=(float s)
    {
        m128 = _mm_mul_ps(m128, _mm_set1_ps(s));
        return *this;
    }

    inline Vec4& operator/=(float s)
    {
        m128 = _mm_div_ps(m128, _mm_set1_ps(s));
        return *this;
    }
};

inline float Dot(Vec4 a, Vec4 b)
{
    a.m128 = _mm_dp_ps(a.m128, b.m128, 0b11110001);
    return a.xyzw[0];
}

constexpr Vec4 Cross(Vec4 a, Vec4 b)
{
    return {
        a.xyzw[1] * b.xyzw[2] - a.xyzw[2] * b.xyzw[1],
        a.xyzw[2] * b.xyzw[0] - a.xyzw[0] * b.xyzw[2],
        a.xyzw[0] * b.xyzw[0] - a.xyzw[1] * b.xyzw[0],
    };
}

inline float LengthSquared(Vec4 vec)
{
    __m128 vector128f = _mm_load_ps(vec.xyzw);
    vector128f = _mm_mul_ps(vector128f, vector128f);

    vector128f = _mm_hadd_ps(vector128f, vector128f);
    vector128f = _mm_hadd_ps(vector128f, vector128f);

    _mm_store_ss(vec.xyzw, vector128f);

    return vec.xyzw[0];
}

inline float Length(Vec4 vec)
{
    return std::sqrtf(LengthSquared(vec));
}

inline Vec4 Normalize(Vec4 vec)
{
    __m128 lengthVector128f = _mm_load_ps(vec.xyzw);
    __m128 vector128f = lengthVector128f;
    lengthVector128f = _mm_mul_ps(lengthVector128f, lengthVector128f);

    lengthVector128f = _mm_hadd_ps(lengthVector128f, lengthVector128f);
    lengthVector128f = _mm_hadd_ps(lengthVector128f, lengthVector128f);
    lengthVector128f = _mm_sqrt_ps(lengthVector128f);
    vector128f = _mm_div_ps(vector128f, lengthVector128f);

    _mm_store_ps(vec.xyzw, vector128f);

    return vec;
}

inline Vec4 Reflect(Vec4 vec, Vec4 normal)
{
    return vec - normal * math::Dot(vec, normal) * 2.f;
}

struct Mat4
{
    union {
        struct {
            __m128 r0, r1, r2, r3;
        };
        float array[16];
        struct {
            float _00, _01, _02, _03;
            float _10, _11, _12, _13;
            float _20, _21, _22, _23;
            float _30, _31, _32, _33;
        };
    };

    Vec4 operator*(Vec4 vec) const
    {
        return {
            _mm_dp_ps(r0, vec.m128, 0b11110001).m128_f32[0],
            _mm_dp_ps(r1, vec.m128, 0b11110001).m128_f32[0],
            _mm_dp_ps(r2, vec.m128, 0b11110001).m128_f32[0],
            _mm_dp_ps(r3, vec.m128, 0b11110001).m128_f32[0],
        };
    }
};

constexpr Mat4 CreateIdentityMatrix()
{
    return {
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    };
}

Mat4 CreateCameraBasisMatrix(Vec4 eyePos, Vec4 lookAt, Vec4 upDir)
{
    Vec4 viewDir = Normalize(lookAt - eyePos);
    Vec4 right = Normalize(Cross(upDir, viewDir));
    upDir = Cross(viewDir, right);

    return {
        right.m128,
        upDir.m128,
        viewDir.m128,
    };
}

Mat4 Transpose(Mat4 mat)
{
    return {
        mat._00, mat._10, mat._20, mat._30,
        mat._01, mat._11, mat._21, mat._31,
        mat._02, mat._12, mat._22, mat._32,
        mat._03, mat._13, mat._23, mat._33,
    };
}

} // namespace math

namespace
{
class splitmix
{
public:
    using result_type = uint32_t;
    static constexpr result_type(min)() { return 0; }
    static constexpr result_type(max)() { return UINT32_MAX; }

    splitmix() : m_seed(1) {}
    explicit splitmix(uint64_t seed) : m_seed(seed << 31 | seed) {}
    explicit splitmix(std::random_device &rd)
    {
        seed(rd);
    }

    void seed(std::random_device &rd)
    {
        m_seed = uint64_t(rd()) << 31 | uint64_t(rd());
    }

    result_type operator()()
    {
        uint64_t z = (m_seed += UINT64_C(0x9E3779B97F4A7C15));
        z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
        return result_type((z ^ (z >> 31)) >> 31);
    }

    void discard(unsigned long long n)
    {
        for (unsigned long long i = 0; i < n; ++i)
            operator()();
    }

private:
    uint64_t m_seed;
};

class xorshift
{
public:
    using result_type = uint32_t;
    static constexpr result_type(min)() { return 0; }
    static constexpr result_type(max)() { return UINT32_MAX; }

    xorshift() : m_seed(0xc1f651c67c62c6e0ull) {}
    explicit xorshift(uint64_t seed) : m_seed(seed << 31 | seed) {}
    explicit xorshift(std::random_device &rd)
    {
        seed(rd);
    }

    void seed(std::random_device &rd)
    {
        m_seed = uint64_t(rd()) << 31 | uint64_t(rd());
    }

    result_type operator()()
    {
        uint64_t result = m_seed * 0xd989bcacc137dcd5ull;
        m_seed ^= m_seed >> 11;
        m_seed ^= m_seed << 31;
        m_seed ^= m_seed >> 18;
        return uint32_t(result >> 32ull);
    }

    void discard(unsigned long long n)
    {
        for (unsigned long long i = 0; i < n; ++i)
            operator()();
    }

private:
    uint64_t m_seed;
};

__forceinline float GenerateUniformRealDist(float min = -1.f, float max = 1.f)
{
    static thread_local ::splitmix rand(
        static_cast<unsigned int>(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::uniform_real_distribution<float> u(min, max);

    return u(rand);
}

template < typename uint8_t D >
__forceinline math::Vec4 GenerateUnitVectorValue()
{
    auto vec = GenerateUnitVectorValue<D - 1>();
    vec.xyzw[D] = GenerateUniformRealDist();
    return vec;
}

template <>
__forceinline math::Vec4 GenerateUnitVectorValue<0>()
{
    return { GenerateUniformRealDist(), 0, 0, 0 };
}

template < typename uint8_t D = 3 >
__forceinline math::Vec4 GenerateUnitVector()
{
    return math::Normalize(::GenerateUnitVectorValue<D - 1>());
}

__forceinline math::Vec4 GenerateUniformDistInsideSphereVector(float radius = 0.5f)
{
    math::Vec4 result;

    do {
        result = { GenerateUniformRealDist(-radius, radius),
            GenerateUniformRealDist(-radius, radius),
            GenerateUniformRealDist(-radius, radius) };
    } while (math::Length(result) < radius);

    return result;
}

__forceinline math::Vec4 GenerateNormalDistInsideSphereVector(float radius = 0.5f)
{
    math::Vec4 result;

    do {
        result = { GenerateUniformRealDist(-radius, radius),
            GenerateUniformRealDist(-radius, radius),
            GenerateUniformRealDist(-radius, radius) };
    } while (math::Length(result) < radius);

    return result;
}

__forceinline math::Vec4 Attenuate(math::Vec4 color, math::Vec4 attenuation)
{
    color.m128 = _mm_mul_ps(color.m128, attenuation.m128);
    return color;
}
} // namespace ::

namespace collision
{
inline bool RaySphereIntersection(
    math::Vec4 sphereCenter, float sphereRadius, math::Vec4 rayDirection, math::Vec4 rayOrigin, float threshold = 1e-3f
)
{
    sphereCenter -= rayOrigin;
    float const tCenter = Dot(sphereCenter, rayDirection);
    float const distanceSquare = Dot(sphereCenter, sphereCenter) - tCenter * tCenter;

    return (tCenter > threshold) && (sphereRadius * sphereRadius - distanceSquare > threshold);
}

inline float CalculateRaySphereMinIntersectionFactor(
    math::Vec4 raySphere, float sphereRadius, math::Vec4 rayDirection
)
{
    float const tCenter = Dot(raySphere, rayDirection);
    float const distanceSquare = Dot(raySphere, raySphere) - tCenter * tCenter;
    float const tDelta = std::sqrtf(sphereRadius * sphereRadius - distanceSquare);

    return tCenter - tDelta;
};

inline float CalculateRaySphereMaxIntersectionFactor(
    math::Vec4 raySphere, float sphereRadius, math::Vec4 rayDirection
)
{
    float const tCenter = Dot(raySphere, rayDirection);
    float const distanceSquare = Dot(raySphere, raySphere) - tCenter * tCenter;
    float const tDelta = std::sqrtf(sphereRadius * sphereRadius - distanceSquare);

    return tCenter + tDelta;
};

inline math::Vec4 CalculateRaySphereIntersectionFactors(
    math::Vec4 raySphere, float sphereRadius, math::Vec4 rayDirection
)
{
    float const tCenter = Dot(raySphere, rayDirection);
    float const distanceSquare = Dot(raySphere, raySphere) - tCenter * tCenter;
    float const tDelta = std::sqrtf(sphereRadius * sphereRadius - distanceSquare);

    return { tCenter - tDelta, tCenter + tDelta };
};

inline math::Vec4 CalculateRaySphereClosestContactPoint(
    math::Vec4 sphereCenter, float radius, math::Vec4 rayOrigin, math::Vec4 rayDirection
)
{
    float const rayFactor = CalculateRaySphereMinIntersectionFactor(
        sphereCenter - rayOrigin, radius, rayDirection);

    return rayOrigin + rayDirection * rayFactor;
}

inline math::Vec4 CalculateRaySphereFarthestContactPoint(
    math::Vec4 sphereCenter, float radius, math::Vec4 rayOrigin, math::Vec4 rayDirection
)
{
    float const rayFactor = CalculateRaySphereMaxIntersectionFactor(
        sphereCenter - rayOrigin, radius, rayDirection);

    return rayOrigin + rayDirection * rayFactor;
}

inline math::Vec4 CalculateRaySphereContactNormal(
    math::Vec4 contactPoint, math::Vec4 sphereCenter
)
{
    return Normalize(contactPoint - sphereCenter);
}

inline bool RayPlaneIntersection(
    math::Vec4 planeNormal, math::Vec4 planePoint, math::Vec4 rayOrigin, math::Vec4 rayDirection
)
{
    return Dot(planeNormal, rayDirection) < 0.0
        && Dot(planeNormal, planePoint) < Dot(planeNormal, rayOrigin);
}

inline math::Vec4 CalculateRayPlaneContactPoint(
    math::Vec4 planeNormal, math::Vec4 planePoint, math::Vec4 rayOrigin, math::Vec4 rayDirection
)
{
    float const rayPlaneProjection = Dot(planeNormal, rayDirection);
    float const t = Dot(planeNormal, planePoint - rayOrigin) / rayPlaneProjection;
    return rayOrigin + rayDirection * t;
}
} // namespace collision

uint32_t constexpr bounces = 10;
uint32_t constexpr samples = 100;
uint32_t constexpr width  = 1440;
uint32_t constexpr heigth = 1440;
float constexpr ratio = static_cast<float>(width) / heigth;
uint8_t  constexpr stride = 3;
uint32_t constexpr size = width * heigth * stride;
uint8_t* const data = (uint8_t*)std::malloc(sizeof(uint8_t) * size);

math::Mat4 viewMatrix = math::CreateIdentityMatrix();
math::Vec4 constexpr lookAt = { 0, 1, 0 };
math::Vec4 constexpr eyePos = { 0, 1, -3 };
math::Vec4 constexpr upDir  = { 0, 1, 0 };
math::Vec4 constexpr lightPos = { 5, -10, 1 };
math::Vec4 constexpr planeNormal = { 0, 1, 0 };
math::Vec4 constexpr planePoint = { 0, -0.5, 0 };
math::Vec4 constexpr planeColor = { 246, 219, 219 };
math::Vec4 constexpr initColor = { 137, 207, 240 };

enum class Material : uint8_t
{
    SKYBOX,
    REFLECTIVE,
    REFRACTIVE,
    DIFFUSE,
};

struct RenderSegment
{
    uint32_t yBegin = 0;
    uint32_t yEnd = 0;
    uint32_t xBegin = 0;
    uint32_t xEnd = 0;
};

std::vector<math::Vec4> g_colors;
std::vector<math::Vec4> g_spheres;
std::vector<float> g_radii;
std::vector<Material> g_materials;
std::vector<math::Vec4> g_attenuations;
std::vector<float> g_diffuses;
uint32_t g_sphereNumber = 10;

void GenerateSpheres()
{
    g_colors = {{30, 144, 255}};
    g_spheres = {{0, -1e6f, 0}};
    g_radii = { 1e6f };
    g_materials = { Material::DIFFUSE };

    g_colors.push_back({ 0, 0, 0 });
    g_spheres.push_back({ 0, 3, 10 });
    g_radii.push_back(3);
    g_materials.push_back(Material::REFRACTIVE);

    g_colors.push_back({ 0, 0, 0 });
    g_spheres.push_back({ 5, 3, 5 });
    g_radii.push_back(3);
    g_materials.push_back(Material::REFLECTIVE);

    g_colors.push_back({ 223, 55, 132 });
    g_spheres.push_back({ -7, 3, 14 });
    g_radii.push_back(3);
    g_materials.push_back(Material::DIFFUSE);

    float const minSphereDistance = 0.1f;
    float const minR = 0.3f;
    float const maxR = 0.5f;

    for (float z = 0; z < 20; z += 1.25f)
    {
        float const bound = (abs(z)) * 0.85f;
        for (float x = -5 - bound; x < 6 + bound; x += 1.25f)
        {
            if (::GenerateUniformRealDist(0, 1.f) > 0.5f)
            {
                g_radii.push_back(::GenerateUniformRealDist(minR, maxR));
                g_spheres.push_back({ x + ::GenerateUniformRealDist(0, minR), g_radii.back(), z + ::GenerateUniformRealDist(0, minR) });

                if ((math::Length(g_spheres.back() - g_spheres[1]) - g_radii.back() - g_radii[1] < 0.5f)
                    || (math::Length(g_spheres.back() - g_spheres[2]) - g_radii.back() - g_radii[2] < 0.5f)
                    || (math::Length(g_spheres.back() - g_spheres[3]) - g_radii.back() - g_radii[3] < 0.5f))
                {
                    g_radii.pop_back();
                    g_spheres.pop_back();
                    continue;
                }

                g_colors.push_back({ ::GenerateUniformRealDist(0, 255), ::GenerateUniformRealDist(0, 255), ::GenerateUniformRealDist(0, 255) });
                g_materials.push_back({ static_cast<Material>(static_cast<uint8_t>(std::min(std::round(::GenerateUniformRealDist(0.5f, 6.0f)), 3.0f))) });
            }
        }
    }
    g_sphereNumber = static_cast<uint32_t>(g_radii.size());

    g_attenuations.resize(g_sphereNumber, { 1, 1, 1 });
    for (auto& v : g_attenuations)
        if (::GenerateUniformRealDist(0, 1) > 0.2f)
            v = GenerateUnitVector();

    g_diffuses.resize(g_sphereNumber, 0.0f);
    for (auto& s : g_diffuses)
        if (::GenerateUniformRealDist(0, 1) > 0.2f)
            s = GenerateUniformRealDist(0, 1);
    g_diffuses[2] = 0.01f;
}

void InitSpheres()
{
    g_colors = {
        { 30, 144, 255},
        { 10, 255, 110}, {110,  10, 255}, {255, 100, 230},
        {200, 255, 110}, {210,  10, 255}, {255, 100, 150},
        { 50, 255, 200}, { 10, 210, 255}, {255, 100, 220},
    };

    g_spheres = {
        {0, -1e3f - 0.5f, 0},
        {-1, 0, 0}, {0, 0, 0}, {1, 0, 0},
        {-1, 1, 0}, {0, 1, 0}, {1, 1, 0},
        {-1, 2, 0}, {0, 2, 0}, {1, 2, 0},
    };

    g_radii = {
        1e3f,
        0.5f, 0.5f, 0.5f,
        0.5f, 0.5f, 0.5f,
        0.5f, 0.5f, 0.5f,
    };

    g_materials = {
        Material::DIFFUSE,
        Material::DIFFUSE, Material::REFLECTIVE, Material::DIFFUSE,
        Material::DIFFUSE, Material::REFRACTIVE, Material::DIFFUSE,
        Material::DIFFUSE, Material::REFLECTIVE, Material::DIFFUSE,
    };

    g_attenuations.resize(g_sphereNumber, {1, 1, 1});
    for (auto& v : g_attenuations)
        if (::GenerateUniformRealDist(0, 1) > 0.3f)
            v = GenerateUnitVector();

    g_diffuses.resize(g_sphereNumber, 0.01f);
    for (auto& s : g_diffuses)
        if (::GenerateUniformRealDist(0, 1) > 0.3f)
            s = GenerateUniformRealDist(0, 1);

    g_diffuses[2] = 0;
}

inline void WritePixel(uint32_t index, math::Vec4 color)
{
    data[index + 0] = static_cast<uint8_t>(std::round(std::sqrt(color.xyzw[0] / 255.f) * 255.f));
    data[index + 1] = static_cast<uint8_t>(std::round(std::sqrt(color.xyzw[1] / 255.f) * 255.f));
    data[index + 2] = static_cast<uint8_t>(std::round(std::sqrt(color.xyzw[2] / 255.f) * 255.f));
}

inline uint8_t FindClosestIntersectionSphere(math::Vec4 primeRayDirection, math::Vec4 primeRayOrigin)
{
    uint8_t minIndex = g_sphereNumber;
    float minDistanceSq = FLT_MAX;

    for (uint8_t index = 0; index < g_sphereNumber; ++index)
    {
        if (collision::RaySphereIntersection(g_spheres[index], g_radii[index], primeRayDirection, primeRayOrigin))
        {
            math::Vec4 const intersectionPoint = collision::CalculateRaySphereClosestContactPoint(
                g_spheres[index], g_radii[index], primeRayOrigin, primeRayDirection
            );

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

math::Vec4 SampleColor(math::Vec4 direction, math::Vec4 origin, uint32_t bounceCount);

template < typename Material M >
inline math::Vec4 SampleColor(math::Vec4 direction, math::Vec4 origin, uint32_t bounceCount, uint32_t sphereIndex);

template <>
inline math::Vec4 SampleColor<Material::SKYBOX>(math::Vec4 direction, math::Vec4 origin, uint32_t bounceCount, uint32_t sphereIndex)
{
    return initColor * (direction.xyzw[1] + 1.f) * 0.5f;
}

template <>
inline math::Vec4 SampleColor<Material::DIFFUSE>(math::Vec4 direction, math::Vec4 origin, uint32_t bounceCount, uint32_t sphereIndex)
{
    math::Vec4 sampleColor = g_colors[sphereIndex] * 0.5f;
    origin = collision::CalculateRaySphereClosestContactPoint(g_spheres[sphereIndex], g_radii[sphereIndex], origin, direction);
    direction = math::Normalize(collision::CalculateRaySphereContactNormal(origin, g_spheres[sphereIndex]) + GenerateUniformDistInsideSphereVector());
    sphereIndex = FindClosestIntersectionSphere(direction, origin);

    while (--bounceCount && sphereIndex < g_sphereNumber && g_materials[sphereIndex] != Material::REFRACTIVE) {
        sampleColor = sampleColor * 0.5f;
        origin = collision::CalculateRaySphereClosestContactPoint(g_spheres[sphereIndex], g_radii[sphereIndex], origin, direction);
        direction = math::Normalize(origin + collision::CalculateRaySphereContactNormal(origin, g_spheres[sphereIndex]) + GenerateUniformDistInsideSphereVector());
        sphereIndex = FindClosestIntersectionSphere(direction, origin);
    }

    return sampleColor;
}

template <>
inline math::Vec4 SampleColor<Material::REFLECTIVE>(math::Vec4 direction, math::Vec4 origin, uint32_t bounceCount, uint32_t sphereIndex)
{
    origin = collision::CalculateRaySphereClosestContactPoint(g_spheres[sphereIndex], g_radii[sphereIndex], origin, direction);
    math::Vec4 normal = collision::CalculateRaySphereContactNormal(origin, g_spheres[sphereIndex]);
    direction = math::Normalize(math::Reflect(direction, normal) + GenerateNormalDistInsideSphereVector() * g_diffuses[sphereIndex]);

    return SampleColor(direction, origin, bounceCount);
}

template <>
inline math::Vec4 SampleColor<Material::REFRACTIVE>(math::Vec4 direction, math::Vec4 origin, uint32_t bounceCount, uint32_t sphereIndex)
{
    float constexpr nAir = 1.0f;
    float constexpr nGlass = 1.5f;

    origin = collision::CalculateRaySphereClosestContactPoint(g_spheres[sphereIndex], g_radii[sphereIndex], origin, direction);
    math::Vec4 n = collision::CalculateRaySphereContactNormal(origin, g_spheres[sphereIndex]);

    float c = math::Dot(-n, direction);
    float r = nAir / nGlass;
    float rSq = pow((nAir - nGlass) / (nAir + nGlass), 2.f);
    float schlick = rSq + (1.f - rSq) * pow(1.f - c, 5.f);

    if (::GenerateUniformRealDist(0, 1) < schlick)
    {
        return SampleColor(math::Reflect(direction, n), origin, bounceCount);
    }

    if (r*sqrt(1.f - c*c) < 1.f)
    {
        direction = math::Normalize(direction*r + n * (r*c - sqrt(1.f - r*r * (1.f - c*c))));
        origin = collision::CalculateRaySphereFarthestContactPoint(g_spheres[sphereIndex], g_radii[sphereIndex], origin, direction);
        n = -collision::CalculateRaySphereContactNormal(origin, g_spheres[sphereIndex]);

        c = math::Dot(-n, direction);
        r = nGlass / nAir;

        rSq = pow((nGlass - nAir) / (nGlass + nAir), 2.f);
        schlick = rSq + (1.f - rSq) * pow(1.f - c, 5.f);
        if (::GenerateUniformRealDist(0, 1) < schlick)
        {
            return SampleColor(math::Reflect(direction, n), origin, bounceCount);
        }

        if (r*sqrt(1.f - c*c) < 1.f)
        {
            direction = math::Normalize(direction*r + n * (r*c - sqrt(1.f - r*r * (1.f - c*c))));
            return SampleColor(direction, origin, bounceCount);
        }

        return SampleColor(math::Reflect(direction, n), origin, bounceCount);
    }

    return SampleColor(math::Reflect(direction, n), origin, bounceCount);
}

math::Vec4 SampleColor(math::Vec4 direction, math::Vec4 origin, uint32_t bounceCount)
{
    uint32_t const sphereIndex = FindClosestIntersectionSphere(direction, origin);

    if (sphereIndex < g_sphereNumber)
    {
        switch (g_materials[sphereIndex])
        {
        case Material::DIFFUSE:
            return SampleColor<Material::DIFFUSE>(direction, origin, bounceCount, sphereIndex);
        case Material::REFLECTIVE:
            return SampleColor<Material::REFLECTIVE>(direction, origin, bounceCount, sphereIndex);
        case Material::REFRACTIVE:
            return SampleColor<Material::REFRACTIVE>(direction, origin, bounceCount, sphereIndex);
        }
    }

    return SampleColor<Material::SKYBOX>(direction, origin, bounceCount, sphereIndex);
}

math::Vec4 SampleColor(math::Vec4 direction, math::Vec4 origin, uint32_t sphereIndex, uint32_t bounceCount)
{
    switch (g_materials[sphereIndex])
    {
    case Material::DIFFUSE:
        return SampleColor<Material::DIFFUSE>(direction, origin, bounceCount, sphereIndex);
    case Material::REFLECTIVE:
        return SampleColor<Material::REFLECTIVE>(direction, origin, bounceCount, sphereIndex);
    case Material::REFRACTIVE:
        return SampleColor<Material::REFRACTIVE>(direction, origin, bounceCount, sphereIndex);
    }

    return { 255, 0, 0 };
}

void Render(RenderSegment segment)
{
    for (uint32_t y = segment.yBegin; y < segment.yEnd; ++y)
    {
        for (uint32_t x = segment.xBegin; x < segment.xEnd; ++x)
        {
            uint32_t const index = size - ((width - x) * stride + y * width * stride);
            math::Vec4 pixelColor{ 0, 0, 0, 0 };

            for (uint32_t s = 0; s < samples; ++s)
            {
                float const u = static_cast<float>(y + GenerateUniformRealDist()) / width;
                float const v = static_cast<float>(x + GenerateUniformRealDist()) / heigth;
                math::Vec4 const primeRayDirection = math::Normalize(viewMatrix * math::Vec4{ -1.f + 2.f * v, -1.f + 2.f * u, 1.f } );
                math::Vec4 const primeRayOrigin{ eyePos };

                uint32_t const sphereIndex = FindClosestIntersectionSphere(primeRayDirection, primeRayOrigin);
                if (sphereIndex < g_sphereNumber)
                {
                    pixelColor += SampleColor(primeRayDirection, primeRayOrigin, sphereIndex, bounces);
                }
                else
                {
                    pixelColor += SampleColor<Material::SKYBOX>(primeRayDirection, primeRayOrigin, sphereIndex, bounces);
                }
            }

            pixelColor *= (1.f / static_cast<float>(samples));
            WritePixel(index, pixelColor);
        }
    }
}

void SaveImage()
{
    stbi_write_bmp((std::string("output") + std::to_string(samples) + "s" + std::to_string(bounces) + "b" + ".bmp").c_str(), width, heigth, stride, data);
}

RenderSegment MakeRenderSegment(uint32_t i, uint32_t j, uint32_t segmentWidth, uint32_t segmentHeigth)
{
    uint32_t const yBegin = segmentHeigth * j;
    uint32_t const yEnd = yBegin + segmentHeigth > heigth ? heigth : yBegin + segmentHeigth;
    uint32_t const xBegin = segmentWidth * i;
    uint32_t const xEnd = xBegin + segmentWidth > width ? width : xBegin + segmentWidth;

    return { yBegin, yEnd, xBegin, xEnd };
}

void RenderJob(RenderSegment segment, std::atomic<int>& freeThreadsCount, std::condition_variable& freeThreadWaitCondition)
{
    freeThreadsCount.fetch_sub(1);
    Render(segment);
    freeThreadsCount.fetch_add(1);
    freeThreadWaitCondition.notify_one();
}

void RenderImageParallel()
{
    uint32_t const hardwearThreads = std::thread::hardware_concurrency();
    uint32_t const threadCount = (hardwearThreads % 2 != 0 ? hardwearThreads + 1 : hardwearThreads);

    uint32_t const segmentWidth = width / threadCount;
    uint32_t const segmentHeigth = heigth / threadCount;
    std::vector<RenderSegment> segments;

    for (uint32_t j = 0; j < threadCount; ++j) {
        for (uint32_t i = 0; i < threadCount; ++i) {
            segments.push_back({ MakeRenderSegment(i, j, segmentWidth, segmentHeigth) });
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
        freeThreadWaitCondition.wait(qlock, [&freeThreadsCount]()->bool { return freeThreadsCount.load() > 0; });
    }

    freeThreadWaitCondition.wait(qlock, [&freeThreadsCount, &threadCount]()->bool { return freeThreadsCount.load() == threadCount; });
}

void RenderImage()
{
    Render({ 0, heigth, 0, width });
}

int main()
{
    static_assert(width % 2 == 0);
    static_assert(heigth % 2 == 0);

    viewMatrix = Transpose(math::CreateCameraBasisMatrix(eyePos, lookAt, upDir));
    InitSpheres();
    RenderImageParallel();
    SaveImage();
}
