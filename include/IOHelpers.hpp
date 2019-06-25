#pragma once

#include <Math.hpp>
#include <Globals.hpp>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <cstdint>
#include <cmath>

namespace io
{

inline void WritePixel(uint32_t index, math::Vec4 color)
{
    g_data[index + 0] = static_cast<uint8_t>(std::round(std::sqrt(color.xyzw[0] / 255.f) * 255.f));
    g_data[index + 1] = static_cast<uint8_t>(std::round(std::sqrt(color.xyzw[1] / 255.f) * 255.f));
    g_data[index + 2] = static_cast<uint8_t>(std::round(std::sqrt(color.xyzw[2] / 255.f) * 255.f));
}

void SaveImage()
{
    stbi_write_bmp((std::string("output") + std::to_string(g_samples) + "s" + std::to_string(g_bounces) + "b" + ".bmp").c_str(), g_width, g_height, g_stride, g_data);
}

} // namespace io
