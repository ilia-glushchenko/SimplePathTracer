#pragma once

#include <Random.hpp>
#include <Globals.hpp>

inline void GenerateSpheres()
{
    g_colors = {{30, 144, 255}};
    g_spheres = {{0, -1e6f, 0}};
    g_radii = {1e6f};
    g_materials = {Material::DIFFUSE};

    g_colors.push_back({0, 0, 0});
    g_spheres.push_back({0, 3, 10});
    g_radii.push_back(3);
    g_materials.push_back(Material::REFRACTIVE);

    g_colors.push_back({0, 0, 0});
    g_spheres.push_back({5, 3, 5});
    g_radii.push_back(3);
    g_materials.push_back(Material::REFLECTIVE);

    g_colors.push_back({223, 55, 132});
    g_spheres.push_back({-7, 3, 14});
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
            if (random::GenerateUniformRealDist(0, 1.f) > 0.5f)
            {
                g_radii.push_back(random::GenerateUniformRealDist(minR, maxR));
                g_spheres.push_back({x + random::GenerateUniformRealDist(0, minR), g_radii.back(), z + random::GenerateUniformRealDist(0, minR)});

                if ((math::Length(g_spheres.back() - g_spheres[1]) - g_radii.back() - g_radii[1] < 0.5f) || (math::Length(g_spheres.back() - g_spheres[2]) - g_radii.back() - g_radii[2] < 0.5f) || (math::Length(g_spheres.back() - g_spheres[3]) - g_radii.back() - g_radii[3] < 0.5f))
                {
                    g_radii.pop_back();
                    g_spheres.pop_back();
                    continue;
                }

                g_colors.push_back({random::GenerateUniformRealDist(0, 255), random::GenerateUniformRealDist(0, 255), random::GenerateUniformRealDist(0, 255)});
                g_materials.push_back({static_cast<Material>(static_cast<uint8_t>(std::min(std::round(random::GenerateUniformRealDist(0.5f, 6.0f)), 3.0f)))});
            }
        }
    }
    g_sphereNumber = static_cast<uint32_t>(g_radii.size());

    g_attenuations.resize(g_sphereNumber, {1, 1, 1});
    for (auto &v : g_attenuations)
        if (random::GenerateUniformRealDist(0, 1) > 0.2f)
            v = random::GenerateUnitVector();

    g_diffuses.resize(g_sphereNumber, 0.0f);
    for (auto &s : g_diffuses)
        if (random::GenerateUniformRealDist(0, 1) > 0.2f)
            s = random::GenerateUniformRealDist(0, 1);
    g_diffuses[2] = 0.01f;
}

inline void InitSpheres()
{
    g_colors = {
        {30, 144, 255},
        {10, 255, 110},
        {110, 10, 255},
        {255, 100, 230},
        {200, 255, 110},
        {210, 10, 255},
        {255, 100, 150},
        {50, 255, 200},
        {10, 210, 255},
        {255, 100, 220},
    };

    g_spheres = {
        {0, -1e3f - 0.5f, 0},
        {-1, 0, 0},
        {0, 0, 0},
        {1, 0, 0},
        {-1, 1, 0},
        {0, 1, 0},
        {1, 1, 0},
        {-1, 2, 0},
        {0, 2, 0},
        {1, 2, 0},
    };

    g_radii = {
        1e3f,
        0.5f,
        0.5f,
        0.5f,
        0.5f,
        0.5f,
        0.5f,
        0.5f,
        0.5f,
        0.5f,
    };

    g_materials = {
        Material::DIFFUSE,
        Material::DIFFUSE,
        Material::REFLECTIVE,
        Material::DIFFUSE,
        Material::DIFFUSE,
        Material::REFRACTIVE,
        Material::DIFFUSE,
        Material::DIFFUSE,
        Material::REFLECTIVE,
        Material::DIFFUSE,
    };

    g_attenuations.resize(g_sphereNumber, {1, 1, 1});
    for (auto &v : g_attenuations)
        if (random::GenerateUniformRealDist(0, 1) > 0.3f)
            v = random::GenerateUnitVector();

    g_diffuses.resize(g_sphereNumber, 0.01f);
    for (auto &s : g_diffuses)
        if (random::GenerateUniformRealDist(0, 1) > 0.3f)
            s = random::GenerateUniformRealDist(0, 1);

    g_diffuses[2] = 0;
}
