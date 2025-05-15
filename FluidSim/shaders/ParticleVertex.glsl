#version 410 core

layout (location = 0) in vec4 pos;
layout (location = 1) in vec3 colour;
layout(location = 2) in float density; // Add this line
out vec3 particle_colour;

uniform mat4 MVP;
out float v_density; // Add this line

void main()
{
    gl_Position = MVP * vec4(pos.xyz, 1.0);
    gl_PointSize= pos.w;
    particle_colour = colour;
    v_density = density;

}