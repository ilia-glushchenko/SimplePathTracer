#pragma once

#include <Math.hpp>
#include <cstdint>
#include <random>
#include <chrono>

namespace random
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

inline __forceinline float GenerateUniformRealDist(float min = -1.f, float max = 1.f)
{
    static thread_local splitmix rand(
        static_cast<unsigned int>(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::uniform_real_distribution<float> u(min, max);

    return u(rand);
}

template <typename uint8_t D>
inline __forceinline math::Vec4 GenerateUnitVectorValue()
{
    auto vec = GenerateUnitVectorValue<D - 1>();
    vec.xyzw[D] = GenerateUniformRealDist();
    return vec;
}

template <>
inline __forceinline math::Vec4 GenerateUnitVectorValue<0>()
{
    return {GenerateUniformRealDist(), 0, 0, 0};
}

template <typename uint8_t D = 3>
inline __forceinline math::Vec4 GenerateUnitVector()
{
    return math::Normalize(GenerateUnitVectorValue<D - 1>());
}

inline __forceinline math::Vec4 GenerateUniformDistInsideSphereVector(float radius = 0.5f)
{
    math::Vec4 result;

    do
    {
        result = {GenerateUniformRealDist(-radius, radius),
                  GenerateUniformRealDist(-radius, radius),
                  GenerateUniformRealDist(-radius, radius)};
    } while (math::Length(result) < radius);

    return result;
}

inline __forceinline math::Vec4 GenerateNormalDistInsideSphereVector(float radius = 0.5f)
{
    math::Vec4 result;

    do
    {
        result = {GenerateUniformRealDist(-radius, radius),
                  GenerateUniformRealDist(-radius, radius),
                  GenerateUniformRealDist(-radius, radius)};
    } while (math::Length(result) < radius);

    return result;
}

inline __forceinline math::Vec4 Attenuate(math::Vec4 color, math::Vec4 attenuation)
{
    color.m128 = _mm_mul_ps(color.m128, attenuation.m128);
    return color;
}

} // namespace rand