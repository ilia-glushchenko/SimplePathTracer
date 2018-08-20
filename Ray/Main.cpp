#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

#include <windows.h>
#include <intrin.h>

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
	math::Vec4 raySphere, float sphereRadius, math::Vec4 rayDirection
)
{
	float const tCenter = Dot(raySphere, rayDirection);
	float const distanceSquare = Dot(raySphere, raySphere) - tCenter * tCenter;

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
math::Vec4 constexpr planeColor = { 246,219,219 };

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

void InitImage()
{
	for (uint32_t i = 0; i < size; i += stride)
	{
		data[i + 0] = static_cast<uint8_t>(137 * (float(size - i) / size));
		data[i + 1] = static_cast<uint8_t>(207 * (float(size - i) / size));
		data[i + 2] = static_cast<uint8_t>(240 * (float(size - i) / size));
	}
}

void RenderImage()
{
	math::Vec4 const primeRayOrigin = { 0, 1, 0 };

	for (uint32_t y = 0; y < heigth; ++y)
	{
		for (uint32_t x = 0; x < width; ++x)
		{
			float const u = static_cast<float>(y) / width;
			float const v = static_cast<float>(x) / heigth;
			math::Vec4 const primeRayDirection = math::Normalize({ -1.f + 2.f * v, -1.f + 2.f * u, 1.f } );

			bool isEmpty = true;
			for (uint8_t i = 0; i < sphereNumber; ++i)
			{
				if (isEmpty = collision::RaySphereIntersection(sphere[i] - primeRayOrigin, radius[i], primeRayDirection))
				{
					math::Vec4 const shadowRayOrigin = collision::CalculateRaySphereContactPoint(
						sphere[i], radius[i], primeRayOrigin, primeRayDirection);
					math::Vec4 const shadowRayDirection = math::Normalize(lightPos - shadowRayOrigin);

					bool isInShadow = false;
					for (uint8_t j = 0; j < sphereNumber; ++j)
					{
						isInShadow = (i != j) && collision::RaySphereIntersection(
							sphere[j] - shadowRayOrigin, radius[j], shadowRayDirection);
						if (isInShadow)
							break;
					}

					uint32_t const index = size - ((width - x) * stride + y * width * stride);
					if (!isInShadow)
					{
						math::Vec4 const sphereContactNormal 
							= collision::CalculateRaySphereContactNormal(shadowRayOrigin, sphere[i]);

						float const light = ::Clamp(math::Dot(sphereContactNormal, shadowRayDirection));
						data[index + 0] = static_cast<uint8_t>(colors[i].xyzw[0] * light);
						data[index + 1] = static_cast<uint8_t>(colors[i].xyzw[1] * light);
						data[index + 2] = static_cast<uint8_t>(colors[i].xyzw[2] * light);
					}
					else
					{
						data[index + 0] = static_cast<uint8_t>(colors[i].xyzw[0] * 0.01f);
						data[index + 1] = static_cast<uint8_t>(colors[i].xyzw[1] * 0.01f);
						data[index + 2] = static_cast<uint8_t>(colors[i].xyzw[2] * 0.01f);
					}

					break;
				}
			}

			if (!isEmpty)
			{
				if (collision::RayPlaneIntersection(planeNormal, planePoint, primeRayOrigin, primeRayDirection))
				{
					math::Vec4 const shadowRayOrigin = collision::CalculateRayPlaneContactPoint(
						planeNormal, planePoint, primeRayOrigin, primeRayDirection);
					math::Vec4 const shadowRayDirection = Normalize(lightPos - shadowRayOrigin);

					bool isInShadow = false;
					for (uint8_t i = 0; i < sphereNumber; ++i)
					{
						isInShadow = collision::RaySphereIntersection(sphere[i] - shadowRayOrigin, radius[i], shadowRayDirection);
						if (isInShadow)
							break;
					}

					float const light = isInShadow ? 0.1f : ::Clamp(Dot(planeNormal, shadowRayDirection));
					uint32_t const index = size - ((width - x) * stride + y * width * stride);
					data[index + 0] = static_cast<uint8_t>(planeColor.xyzw[0] * light);
					data[index + 1] = static_cast<uint8_t>(planeColor.xyzw[1] * light);
					data[index + 2] = static_cast<uint8_t>(planeColor.xyzw[2] * light);
				}
			}
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