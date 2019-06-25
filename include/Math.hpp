#pragma once

#include <intrin.h>
#include <cmath>

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
            _mm_mul_ps(m128, _mm_set1_ps(s))};
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
            _mm_xor_ps(m128, _mm_set_ps1(-0.0))};
    }

    constexpr Vec4 operator-()
    {
        return {
            -xyzw[0],
            -xyzw[1],
            -xyzw[2],
            -xyzw[3],
        };
    }

    inline Vec4 &operator-=(Vec4 other)
    {
        m128 = _mm_sub_ps(m128, other.m128);
        return *this;
    }

    inline Vec4 &operator+=(Vec4 other)
    {
        m128 = _mm_add_ps(m128, other.m128);
        return *this;
    }

    inline Vec4 &operator*=(float s)
    {
        m128 = _mm_mul_ps(m128, _mm_set1_ps(s));
        return *this;
    }

    inline Vec4 &operator/=(float s)
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
    return vec - normal * Dot(vec, normal) * 2.f;
}

struct Mat4
{
    union {
        struct
        {
            __m128 r0, r1, r2, r3;
        };
        float array[16];
        struct
        {
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
        0, 0, 0, 1};
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
        mat._00,
        mat._10,
        mat._20,
        mat._30,
        mat._01,
        mat._11,
        mat._21,
        mat._31,
        mat._02,
        mat._12,
        mat._22,
        mat._32,
        mat._03,
        mat._13,
        mat._23,
        mat._33,
    };
}
} // namespace math