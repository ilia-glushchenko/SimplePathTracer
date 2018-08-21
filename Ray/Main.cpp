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
	math::Vec4 sphereCenter, float sphereRadius, math::Vec4 rayDirection, math::Vec4 rayOrigin
)
{
    sphereCenter -= rayOrigin;
	float const tCenter = Dot(sphereCenter, rayDirection);
	float const distanceSquare = Dot(sphereCenter, sphereCenter) - tCenter * tCenter;

	return (tCenter > 0) && (sphereRadius * sphereRadius - distanceSquare > 0);
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

uint32_t constexpr samples = 4;
uint32_t constexpr width  = 1080;
uint32_t constexpr heigth = 1080;
float constexpr ratio = static_cast<float>(width) / heigth;
uint8_t  constexpr stride = 3;
uint32_t constexpr size = width * heigth * stride;
uint8_t* const data = (uint8_t*)std::malloc(sizeof(uint8_t) * size);

math::Vec4 constexpr eyePos = { 0, 0, 0 };
math::Vec4 constexpr eyeDir = { 0, 0, 1 };
math::Vec4 constexpr upDir  = { 0, 1, 0 };
float constexpr nearPlane = 1;
float constexpr farPlane  = 100;
math::Vec4 constexpr lightPos = { 5, 10, 1 };
math::Vec4 constexpr planeNormal = { 0, 1, 0 };
math::Vec4 constexpr planePoint = { 0, -0.5, 0 };
math::Vec4 constexpr planeColor = { 246, 219, 219 };
math::Vec4 constexpr initColor = { 137, 207, 240 };

math::Vec4 constexpr colors[9] = {
	{10,  255, 110}, {110, 10, 255}, {255, 100, 230}, 
	{200, 255, 110}, {210, 10, 255}, {255, 100, 150},
	{50,  255, 200}, {10, 210, 255}, {255, 100, 220},
};
math::Vec4 constexpr sphere[9] = {
	{-1, 0, 3}, {0, 0, 3}, {1, 0, 3}, 
	{-1, 1, 3}, {0, 1, 3}, {1, 1, 3}, 
	{-1, 2, 3}, {0, 2, 3}, {1, 2, 3}, 
};
float constexpr radius[9] = {
	0.5f, 0.5f, 0.5f, 
	0.5f, 0.5f, 0.5f,
	0.5f, 0.5f, 0.5f,
};
uint32_t constexpr sphereNumber = 9;

float GenerateUniformRealDist()
{
    static std::random_device rd;
    static std::mt19937 gen(rd()); 
    static std::uniform_real_distribution<float> dis(0.f, 1.0f);

    return dis(gen);
}

namespace
{
template < typename uint8_t D >
math::Vec4 GenerateUnitVectorValue()
{
    auto vec = GenerateUnitVectorValue<D - 1>();
    vec.xyzw[D] = GenerateUniformRealDist();
    return vec;
}

template <>
math::Vec4 GenerateUnitVectorValue<0>()
{
    return { GenerateUniformRealDist(), 0, 0, 0 };
}
} // namespace ::

template < typename uint8_t D = 3 >
math::Vec4 GenerateUnitVector()
{
    return math::Normalize(::GenerateUnitVectorValue<D - 1>());
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
    data[index + 0] = static_cast<uint8_t>(std::round(color.xyzw[0]));
    data[index + 1] = static_cast<uint8_t>(std::round(color.xyzw[1]));
    data[index + 2] = static_cast<uint8_t>(std::round(color.xyzw[2]));
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

inline math::Vec4 SampleColor(math::Vec4 direction, math::Vec4 origin, math::Vec4 sampleColor, uint32_t bounces)
{  
    uint8_t const sphereIndex = FindClosestIntersectionSphere(direction, origin);

    if (sphereIndex < sphereNumber)
    {
        math::Vec4 const shadowRayOrigin = collision::CalculateRaySphereContactPoint(sphere[sphereIndex], radius[sphereIndex], origin, direction);
        math::Vec4 const shadowRayDirection = math::Normalize(lightPos - shadowRayOrigin);
        math::Vec4 const uint = GenerateUnitVector();

        uint8_t const shadowSphereIndex = FindClosestIntersectionSphere(shadowRayDirection, shadowRayOrigin);
        if (shadowSphereIndex < sphereNumber)
        {
            sampleColor = colors[sphereIndex] * 0.01f;
        }
        else
        {
            math::Vec4 const sphereContactNormal = collision::CalculateRaySphereContactNormal(shadowRayOrigin, sphere[sphereIndex]);
            float const light = ::Clamp(math::Dot(sphereContactNormal, shadowRayDirection));
            sampleColor = colors[sphereIndex] * light;
        }
    }
    else if (collision::RayPlaneIntersection(planeNormal, planePoint, origin, direction))
    {
        math::Vec4 const shadowRayOrigin = collision::CalculateRayPlaneContactPoint(planeNormal, planePoint, origin, direction);
        math::Vec4 const shadowRayDirection = Normalize(lightPos - shadowRayOrigin);

        uint8_t const shadowSphereIndex = FindClosestIntersectionSphere(shadowRayDirection, shadowRayOrigin);
        float const light = shadowSphereIndex < sphereNumber ? 0.1f : ::Clamp(Dot(planeNormal, shadowRayDirection));
        sampleColor = planeColor * light;
    }

    return sampleColor;
}

void RenderImage()
{
	math::Vec4 const primeRayOrigin{ 0, 1, 0, 0 };

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

                pixelColor += SampleColor(primeRayDirection, primeRayOrigin, sampleColor, 10);
            }

            pixelColor *= (1.f / static_cast<float>(samples));
            WritePixel(index, pixelColor);
		}
	}
}

void SaveImage()
{
	stbi_write_bmp("output.bmp", width, heigth, stride, data);
}

int main()
{
	InitImage();
	RenderImage();
	SaveImage();
}