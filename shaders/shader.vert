#version 460 core

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec2 aTexCoords;

uniform mat4 model;
out vec2 vsTexCoordsOut;

void main()
{
    gl_Position = model * vec4(aPos.xyz, 1.0);
    vsTexCoordsOut = aTexCoords;
}
