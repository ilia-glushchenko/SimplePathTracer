#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

#include <windows.h>
#include <intrin.h>

#include <random>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>

namespace 
{
constexpr float Clamp(float s)
{
	return s < 0 ? 0 : s;
}
}

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
} // namespace math

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

inline float CalculateRaySphereIntersectionFactors(
	math::Vec4 raySphere, float sphereRadius, math::Vec4 rayDirection
)
{
	float const tCenter = Dot(raySphere, rayDirection);
	float const distanceSquare = Dot(raySphere, raySphere) - tCenter * tCenter;
	float const tDelta = std::sqrtf(sphereRadius * sphereRadius - distanceSquare);

	return tCenter - tDelta;
};

inline math::Vec4 CalculateRaySphereContactPoint(
	math::Vec4 sphereCenter, float radius, math::Vec4 rayOrigin, math::Vec4 rayDirection
)
{
	float const rayFactor = CalculateRaySphereIntersectionFactors(
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
uint32_t constexpr samples = 12;
uint32_t constexpr width  = 1080;
uint32_t constexpr heigth = 1080;
float constexpr ratio = static_cast<float>(width) / heigth;
uint8_t  constexpr stride = 3;
uint32_t constexpr size = width * heigth * stride;
uint8_t* const data = (uint8_t*)std::malloc(sizeof(uint8_t) * size);

math::Vec4 constexpr eyePos = { 0, 1, 0 };
math::Vec4 constexpr eyeDir = { 0, 0, 1 };
math::Vec4 constexpr upDir  = { 0, 1, 0 };
float constexpr nearPlane = 1;
float constexpr farPlane  = 100;
math::Vec4 constexpr lightPos = { 5, -10, 1 };
math::Vec4 constexpr planeNormal = { 0, 1, 0 };
math::Vec4 constexpr planePoint = { 0, -0.5, 0 };
math::Vec4 constexpr planeColor = { 246, 219, 219 };
math::Vec4 constexpr initColor = { 137, 207, 240 };

enum class Material : uint8_t
{
    DIFFUSE = 0,
    REFLECTIVE = 1,
};

math::Vec4 constexpr colors[10] = {
    { 30, 144, 255 },
    {10,  255, 110}, {110, 10, 255}, {255, 100, 230}, 
	{200, 255, 110}, {210, 10, 255}, {255, 100, 150},
	{50,  255, 200}, {10, 210, 255}, {255, 100, 220},
};
math::Vec4 constexpr sphere[10] = {
    {0, -1000.5f, 0},
    { 0, 0, 3}, {-1, 0, 3}, {1, 0, 3},
	{-1, 1, 3}, {0, 1, 3}, {1, 1, 3}, 
	{-1, 2, 3}, {0, 2, 3}, {1, 2, 3}, 
};
float constexpr radius[10] = {
    1000.f,
    0.5f, 0.5f, 0.5f,
	0.5f, 0.5f, 0.5f,
	0.5f, 0.5f, 0.5f,
};
Material constexpr materials[10] = {
    Material::DIFFUSE, 
    Material::REFLECTIVE, Material::DIFFUSE, Material::DIFFUSE,
    Material::DIFFUSE, Material::DIFFUSE, Material::DIFFUSE,
    Material::DIFFUSE, Material::DIFFUSE, Material::DIFFUSE,
};
uint32_t constexpr sphereNumber = 10;

inline float GenerateUniformRealDist(float min = -1.f, float max = 1.f)
{
    static thread_local std::random_device rd;
    static thread_local std::mt19937_64 gen(rd());
    std::uniform_real_distribution<float> dis(min, max);

    return dis(gen);
}

namespace
{
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
} // namespace ::

template < typename uint8_t D = 3 >
__forceinline math::Vec4 GenerateUnitVector()
{
    return math::Normalize(::GenerateUnitVectorValue<D - 1>());
}

__forceinline math::Vec4 GenerateInsideSphereVector(float radius = 0.5f)
{
    math::Vec4 result;
    
    do {
        result = { GenerateUniformRealDist(-radius, radius), 
            GenerateUniformRealDist(-radius, radius), 
            GenerateUniformRealDist(-radius, radius) };
    } while (math::Length(result) < radius);

    return result;
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
    uint8_t minIndex = sphereNumber;
    float minDistanceSq = FLT_MAX;

    for (uint8_t index = 0; index < sphereNumber; ++index)
    {
        if (collision::RaySphereIntersection(sphere[index], radius[index], primeRayDirection, primeRayOrigin))
        {
            math::Vec4 const intersectionPoint = collision::CalculateRaySphereContactPoint(
                sphere[index], radius[index], primeRayOrigin, primeRayDirection
            );
            float const distanceSq = math::LengthSquared(primeRayOrigin - intersectionPoint);
            minIndex = minDistanceSq > distanceSq ? index : minIndex;
            minDistanceSq = minDistanceSq > distanceSq ? distanceSq : minDistanceSq;
        }
    }

    return minIndex;
}

template < typename Material M >
inline math::Vec4 SampleColor(math::Vec4 direction, math::Vec4 origin, math::Vec4 sampleColor, uint32_t bounceCount);

template <>
inline math::Vec4 SampleColor<Material::DIFFUSE>(math::Vec4 direction, math::Vec4 origin, math::Vec4 sampleColor, uint32_t bounceCount)
{  
    uint8_t sphereIndex = FindClosestIntersectionSphere(direction, origin);
    if (sphereIndex < sphereNumber)
    {
        sampleColor = colors[sphereIndex] * 0.5f;
        origin = collision::CalculateRaySphereContactPoint(sphere[sphereIndex], radius[sphereIndex], origin, direction);
        direction = math::Normalize(collision::CalculateRaySphereContactNormal(origin, sphere[sphereIndex]) + GenerateInsideSphereVector());
        sphereIndex = FindClosestIntersectionSphere(direction, origin);

        while (--bounceCount && sphereIndex < sphereNumber) {
            sampleColor = sampleColor * 0.5f;
            origin = collision::CalculateRaySphereContactPoint(sphere[sphereIndex], radius[sphereIndex], origin, direction);
            direction = math::Normalize(origin + collision::CalculateRaySphereContactNormal(origin, sphere[sphereIndex]) + GenerateInsideSphereVector());
            sphereIndex = FindClosestIntersectionSphere(direction, origin);
        }
    }
    else
    {
        sampleColor = initColor * (direction.xyzw[1] + 1.f) * 0.5f;
    }

    return sampleColor;
}

template <>
inline math::Vec4 SampleColor<Material::REFLECTIVE>(math::Vec4 direction, math::Vec4 origin, math::Vec4 sampleColor, uint32_t bounceCount)
{
    uint8_t sphereIndex = FindClosestIntersectionSphere(direction, origin);
    if (sphereIndex < sphereNumber)
    {
        origin = collision::CalculateRaySphereContactPoint(sphere[sphereIndex], radius[sphereIndex], origin, direction);
        math::Vec4 normal = collision::CalculateRaySphereContactNormal(origin, sphere[sphereIndex]);
        direction = math::Normalize(direction - normal * math::Dot(direction, normal) * 2.f);// +GenerateInsideSphereVector(0.3f));

        sampleColor = SampleColor<Material::DIFFUSE>(direction, origin, { 255.f, 255.f, 255.f, 0.f }, 10);
        sampleColor.xyzw[0] *= 0.8f;
        sampleColor.xyzw[1] *= 0.8f;
        sampleColor.xyzw[2] *= 0.7f;
    }

    return sampleColor;
}

void RenderImage()
{
	math::Vec4 const primeRayOrigin{ eyePos };

	for (uint32_t y = 0; y < heigth; ++y)
	{
		for (uint32_t x = 0; x < width; ++x)
		{
            uint32_t const index = size - ((width - x) * stride + y * width * stride);
            math::Vec4 pixelColor{ 0, 0, 0, 0 };

            for (uint8_t s = 0; s < samples; ++s)
            {
                math::Vec4 sampleColor{ (float)data[index + 0], (float)data[index + 1], (float)data[index + 2] };
			    float const u = static_cast<float>(y + GenerateUniformRealDist()) / width;
			    float const v = static_cast<float>(x + GenerateUniformRealDist()) / heigth;
			    math::Vec4 const primeRayDirection = math::Normalize({ -1.f + 2.f * v, -1.f + 2.f * u, 1.f } );

                uint8_t sphereIndex = FindClosestIntersectionSphere(primeRayDirection, primeRayOrigin);
                if (sphereIndex < sphereNumber)
                {
                    switch (sphereIndex)
                    {                            
                        case 1:
                            pixelColor += SampleColor<Material::REFLECTIVE>(primeRayDirection, primeRayOrigin, sampleColor, bounces);
                            break;
                        default:
                            pixelColor += SampleColor<Material::DIFFUSE>(primeRayDirection, primeRayOrigin, sampleColor, bounces);
                            break;
                    }
                }
                else
                {
                    pixelColor += SampleColor<Material::DIFFUSE>(primeRayDirection, primeRayOrigin, sampleColor, bounces);
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

int main()
{
	InitImage();
	RenderImage();
	SaveImage();
}