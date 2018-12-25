#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

#define NOMINMAX
#include <windows.h>
#include <intrin.h>

#include <algorithm>
#include <thread>
#include <random>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
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
	_mm_store_ss(a.xyzw, _mm_dp_ps(_mm_load_ps(a.xyzw), _mm_load_ps(b.xyzw), 0b11110001));
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
} // namespace math

namespace
{
inline float GenerateUniformRealDist(float min = -1.f, float max = 1.f)
{
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<float> dis(min, max);

    return dis(gen);
}

inline float GenerateNormalRealDist(float min = -1.f, float max = 1.f)
{
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::normal_distribution<float> dis(min, max);

    return dis(gen);
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
uint32_t constexpr samples = 1000;
uint32_t constexpr width  = 1440;
uint32_t constexpr heigth = 1440;
float constexpr ratio = static_cast<float>(width) / heigth;
uint8_t  constexpr stride = 3;
uint32_t constexpr size = width * heigth * stride;
uint8_t* const data = (uint8_t*)std::malloc(sizeof(uint8_t) * size);

math::Vec4 constexpr lookAt = { 0, 0, 0 };
math::Vec4 constexpr eyePos = { 0, 3, -5 };
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

std::vector<math::Vec4> g_colors;
std::vector<math::Vec4> g_spheres;
std::vector<float> g_radii;
std::vector<Material> g_materials;
uint32_t g_sphereNumber = 10;

void GenerateSpheres()
{
    g_colors = { { 30, 144, 255 } };
    g_spheres = { {0, -1000, 0} };
    g_radii = { 1000.f };
    g_materials = { Material::DIFFUSE };

    float const minSphereDistance = 0.1f;
    float const minR = 0.3f;
    float const maxR = 0.5f;

    for (float z = -3; z < 20; z += 1.25f)
    {
        for (float x = -5; x < 6; x += 1.25f)
        {
            if (::GenerateUniformRealDist(0, 1.f) > 0.1f)
            {
                g_radii.push_back(::GenerateUniformRealDist(minR, maxR));
                g_spheres.push_back({ x + ::GenerateUniformRealDist(0, minR), g_radii.back(), z + ::GenerateUniformRealDist(0, minR) });
                g_colors.push_back({ ::GenerateUniformRealDist(0, 255), ::GenerateUniformRealDist(0, 255), ::GenerateUniformRealDist(0, 255) });
                g_materials.push_back({ static_cast<Material>(static_cast<uint8_t>(std::min(std::round(::GenerateUniformRealDist(0.5f, 6.0f)), 3.0f))) });
            }
        }
    }

    g_sphereNumber = g_radii.size();
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
        {0, -1000.5f, 0},
        {-1, 0, 0}, {0, 0, 0}, {1, 0, 0},
        {-1, 1, 0}, {0, 1, 0}, {1, 1, 0},
        {-1, 2, 0}, {0, 2, 0}, {1, 2, 0},
    };
    g_radii = {
        1000.f,
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
}

void InitImage()
{
	for (uint32_t i = 0; i < size; i += stride)
	{
		data[i + 0] = static_cast<uint8_t>(initColor.xyzw[0] * (float(size - i) / size));
		data[i + 1] = static_cast<uint8_t>(initColor.xyzw[1] * (float(size - i) / size));
		data[i + 2] = static_cast<uint8_t>(initColor.xyzw[2] * (float(size - i) / size));
	}
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
inline math::Vec4 SampleColor(math::Vec4 direction, math::Vec4 origin, math::Vec4 sampleColor, uint32_t bounceCount, uint32_t sphereIndex);

template <>
inline math::Vec4 SampleColor<Material::SKYBOX>(math::Vec4 direction, math::Vec4 origin, math::Vec4 sampleColor, uint32_t bounceCount, uint32_t sphereIndex)
{
    return initColor * (direction.xyzw[1] + 1.f) * 0.5f;
}

template <>
inline math::Vec4 SampleColor<Material::DIFFUSE>(math::Vec4 direction, math::Vec4 origin, math::Vec4 sampleColor, uint32_t bounceCount, uint32_t sphereIndex)
{  
    sampleColor = g_colors[sphereIndex] * 0.5f;
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
inline math::Vec4 SampleColor<Material::REFLECTIVE>(math::Vec4 direction, math::Vec4 origin, math::Vec4 sampleColor, uint32_t bounceCount, uint32_t sphereIndex)
{
    origin = collision::CalculateRaySphereClosestContactPoint(g_spheres[sphereIndex], g_radii[sphereIndex], origin, direction);
    math::Vec4 normal = collision::CalculateRaySphereContactNormal(origin, g_spheres[sphereIndex]);
    direction = math::Normalize(math::Reflect(direction, normal) + GenerateNormalDistInsideSphereVector() * 0.01f);

    sampleColor = SampleColor(direction, origin, bounceCount);

    sampleColor.xyzw[0] *= 0.8f;
    sampleColor.xyzw[1] *= 0.8f;
    sampleColor.xyzw[2] *= 0.7f;

    return sampleColor;
}

template <>
inline math::Vec4 SampleColor<Material::REFRACTIVE>(math::Vec4 direction, math::Vec4 origin, math::Vec4 sampleColor, uint32_t bounceCount, uint32_t sphereIndex)
{
    float constexpr nAir = 1.0f;
    float constexpr nGlass = 1.5f;

    origin = collision::CalculateRaySphereClosestContactPoint(g_spheres[sphereIndex], g_radii[sphereIndex], origin, direction);
    math::Vec4 n = collision::CalculateRaySphereContactNormal(origin, g_spheres[sphereIndex]);

    float c = math::Dot(-n, direction);
    float r = nAir / nGlass;
    float schlick = pow((nAir - nGlass) / (nAir + nGlass), 2.f) + (1.f - pow((nAir - nGlass) / (nAir + nGlass), 2)) * pow(1.f - c, 5.f);

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

        schlick = pow((nGlass - nAir) / (nGlass + nAir), 2.f) + (1.f - pow((nGlass - nAir) / (nGlass + nAir), 2.f)) * pow(1.f - c, 5.f);
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
            return SampleColor<Material::DIFFUSE>(direction, origin, { 255.f, 255.f, 255.f, 0.f }, bounceCount, sphereIndex);
        case Material::REFLECTIVE:
            return SampleColor<Material::REFLECTIVE>(direction, origin, { 255.f, 255.f, 255.f, 0.f }, bounceCount, sphereIndex);
        case Material::REFRACTIVE:
            return SampleColor<Material::REFRACTIVE>(direction, origin, { 255.f, 255.f, 255.f, 0.f }, bounceCount, sphereIndex);
        }
    }

    return SampleColor<Material::SKYBOX>(direction, origin, { 255.f, 255.f, 255.f, 0.f }, bounceCount, sphereIndex);
}

void Render(uint32_t yBegin, uint32_t yEnd, uint32_t xBegin, uint32_t xEnd)
{
	math::Vec4 const primeRayOrigin{ eyePos };

	for (uint32_t y = yBegin; y < yEnd; ++y)
	{
		for (uint32_t x = xBegin; x < xEnd; ++x)
		{
            uint32_t const index = size - ((width - x) * stride + y * width * stride);
            math::Vec4 pixelColor{ 0, 0, 0, 0 };

            for (uint32_t s = 0; s < samples; ++s)
            {
                math::Vec4 sampleColor{ (float)data[index + 0], (float)data[index + 1], (float)data[index + 2] };
			    float const u = static_cast<float>(y + GenerateUniformRealDist()) / width;
			    float const v = static_cast<float>(x + GenerateUniformRealDist()) / heigth;
			    math::Vec4 const primeRayDirection = math::Normalize({ -1.f + 2.f * v, -1.f + 2.f * u, 1.f } );

                pixelColor += SampleColor(primeRayDirection, primeRayOrigin, bounces);
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

void RenderImageParallel()
{
    uint32_t const hardwearThreads = std::thread::hardware_concurrency();
    uint32_t const threadCount = (hardwearThreads % 2 != 0 ? hardwearThreads + 1 : hardwearThreads) - 1;
    uint32_t const segmentWidth = width / threadCount;
    uint32_t const segmentHeigth = heigth / threadCount;
    std::vector<std::thread> threads(threadCount);

    for (uint32_t j = 0; j < threadCount; ++j)
    {
        for (uint32_t i = 0; i < threadCount; ++i)
        {
            uint32_t yBegin = segmentHeigth * j;
            uint32_t yEnd = yBegin + segmentHeigth > heigth ? heigth : yBegin + segmentHeigth;
            uint32_t xBegin = segmentWidth * i;
            uint32_t xEnd = xBegin + segmentWidth > width ? width : xBegin + segmentWidth;
            threads[i] = std::thread(Render, yBegin, yEnd, xBegin, xEnd);
        }
        for (auto& thread : threads) 
            if (thread.joinable()) 
                thread.join();
        SaveImage();
    }
}

void RenderImage()
{
    Render(0, heigth, 0, width);
}

int main()
{
    static_assert(width % 2 == 0);
    static_assert(heigth % 2 == 0);

    GenerateSpheres();
	InitImage();
    RenderImageParallel();
	SaveImage();
}