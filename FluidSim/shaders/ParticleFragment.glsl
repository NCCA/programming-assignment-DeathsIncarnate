#version 410 core

// Input from vertex shader
in float v_density;
in vec3 particle_colour;  // Make sure this matches your vertex shader output

// Fragment shader output
layout(location = 0) out vec4 fragColour;

void main()
{
    // Circular point sprite
    vec2 circlecoord = 2.0 * gl_PointCoord - 1.0;
    if(dot(circlecoord, circlecoord) > 1.0)
    {
        discard;
    }

    // Normalize density (0 to 1 range)
    float normalizedDensity = clamp(v_density / 2.0, 0.0, 1.0);

    // Create color gradient from original color to red
    vec3 baseColor = vec3(0.0, 0.0, 1.0);  // Using the input color
    vec3 redColor = vec3(1.0, 0.0, 0.0);

    // Mix between base color and red based on density
    vec3 finalColor = mix(baseColor, redColor, normalizedDensity);

    // Add brightness to dense particles
    finalColor = mix(finalColor, vec3(1.0), pow(normalizedDensity, 3.0) * 0.5);

    fragColour = vec4(finalColor, 1.0);
}