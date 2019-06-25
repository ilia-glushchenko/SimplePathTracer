#version 460 core

in  vec2 vsTexCoordsOut;
out vec4 fragColor;

uniform vec4 fragColorShift;
uniform sampler2D fragTexture;

void main()
{
   fragColor = texture(fragTexture, vsTexCoordsOut) + fragColorShift;
}
