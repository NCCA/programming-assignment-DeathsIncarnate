#include "Emitter.h"
#include <iostream>
#include <ngl/Random.h>
#include <algorithm>
#include <fstream>
#include <ngl/VAOPrimitives.h>
#include <ngl/Mat4.h>
#include <ngl/ShaderLib.h>
#include <ngl/Transformation.h>
#include <ngl/MultiBufferVAO.h>
#include <ngl/Text.h>




Emitter::Emitter(size_t _num, size_t _maxAlive, size_t _numPerFrame, ngl::Vec3 _pos):
m_maxParticles{_num}, m_maxAlive{_maxAlive}, m_numPerFrame{_numPerFrame}, m_ppos{_pos}
{
    m_maxParticles = _num;
    m_ppos.resize(m_maxParticles);
    m_pcolour.resize(m_maxParticles);
    m_pdir.resize(m_maxParticles);
    m_psize.resize(m_maxParticles);
    m_plife.resize(m_maxParticles);
    m_pstate.resize(m_maxParticles);
    m_densities.resize(m_maxParticles);
    m_showDensity = false;
    m_showSmoothing = false;
    m_maxDensity = 0.0f;

    // for (size_t i = 0; i < m_maxParticles; ++i)
    // {
    //     resetParticle(i);
    // }

    // Initialize SPH parameters
    const float m_particleSpacing = 2.0f;
    m_targetDensity = m_restDensity;
    m_smoothingRadius = 2.0f * m_particleSpacing; // Typically 2x particle spacing

    initializeParticles();

    //m_vao = ngl::VAOFactory::createVAO(ngl::simpleVAO, GL_POINTS);
    m_vao = ngl::vaoFactoryCast<ngl::MultiBufferVAO>(ngl::VAOFactory::createVAO(ngl::multiBufferVAO, GL_POINTS));
    m_vao->bind();
    m_vao->setData(ngl::MultiBufferVAO::VertexData(0,0)); // index 0 points
    m_vao->setData(ngl::MultiBufferVAO::VertexData(0,0)); // index 1 colours
    m_vao->setData(ngl::MultiBufferVAO::VertexData(0,0)); // index 3 Densities
    m_vao->unbind();
}
size_t Emitter::size() const
{
    return m_maxParticles;
}

// void Emitter::birthParticles()
// {
//     auto births = 0 + static_cast<int>(ngl::Random::randomPositiveNumber(m_numPerFrame));
//     for(size_t i = 0; i < births; ++i)
//     {
//         for(size_t p = 0; p < m_maxParticles; ++p)
//         {
//             if(m_pstate[p] == ParticleState::Dead)
//             {
//                 resetParticle(p);
//                 m_pstate[p] = ParticleState::Active;
//                 break;
//             }
//         }
//     }
// }

// static float SmoothingKernel(float radius, float dst)
// {
//     float volume = M_PI * pow(radius, 8)/4;
//     float value = std::max(0.0f, std::max(radius - dst * dst, 0.0f));
//
//     return value * value * value/volume;
//
//     // float value = std::max(0.0f, std::max(radius - dst * dst, 0.0f));
//     //
//     // return value * value * value;
// }
//
// static float SmoothingKernelDerivative(float radius, float dst)
// {
//     if (dst >= radius) return 0;
//     float f = radius * radius - dst * dst;
//     float scale = -24 / (M_PI * pow(radius, 8));
//     return scale * dst * f * f;
// }
//
// float Emitter::CalculateDensity(const ngl::Vec3 samplePoint)
// {
//     float density = 0.0f;
//     const float mass = 1.0f;
//     float smoothingRadius = 0.5f;
//
//     for (size_t i = 0; i < m_maxParticles; ++i)
//     {
//         float dst = (m_ppos[i] - samplePoint).length(); //length = magnitude;
//         float influence = SmoothingKernel(smoothingRadius, dst);
//         density += influence * mass;
//     }
//
//     return density;
// }
// Standard cubic spline kernel (Müller et al.)
// Optimized cubic spline kernel for 3D SPH
static float SmoothingKernel(float radius, float dst)
{
    if (dst >= radius) return 0.0f;

    const float h = radius;
    const float h3 = h * h * h;
    const float q = dst / h;

    // Precomputed normalization constant
    const float norm = 8.0f / (M_PI * h3);

    if (q <= 0.5f) {
        const float q2 = q * q;
        const float q3 = q2 * q;
        return norm * (1.0f - 6.0f * q2 + 6.0f * q3);
    }
    else {
        const float term = 2.0f - q;
        return norm * (0.5f * term * term * term);
    }
}

// Optimized derivative of cubic spline kernel
static float SmoothingKernelDerivative(float radius, float dst)
{
    if (dst >= radius || dst == 0.0f) return 0.0f;

    const float h = radius;
    const float h4 = h * h * h * h;
    const float q = dst / h;

    // Precomputed normalization constant
    const float norm = 8.0f / (M_PI * h4);

    if (q <= 0.5f) {
        return norm * (3.0f * q - 4.5f * q * q);
    }
    else {
        const float term = 2.0f - q;
        return norm * (-1.5f * term * term);
    }
}

float Emitter::CalculateDensity(const ngl::Vec3& samplePoint)
{
    float density = 0.0f;
    const float radiusSq = m_smoothingRadius * m_smoothingRadius;

    for (size_t i = 0; i < m_maxParticles; ++i)
    {
        const float dx = m_ppos[i].m_x - samplePoint.m_x;
        const float dy = m_ppos[i].m_y - samplePoint.m_y;
        const float dz = m_ppos[i].m_z - samplePoint.m_z;
        const float dstSq = dx*dx + dy*dy + dz*dz;

        if (dstSq >= radiusSq) continue;

        const float dst = std::sqrt(dstSq);
        density += m_particleMass * SmoothingKernel(m_smoothingRadius, dst);
    }
    std::cout <<"Density" << density << "\n";
    return density;
}

float Emitter::CalculateProperty(const ngl::Vec3& samplePoint,
                               const std::vector<float>& particleProperties)
{
    float property = 0.0f;
    const float radiusSq = m_smoothingRadius * m_smoothingRadius;

    for (int i = 0; i < m_maxParticles; ++i)
    {
        const float dx = m_ppos[i].m_x - samplePoint.m_x;
        const float dy = m_ppos[i].m_y - samplePoint.m_y;
        const float dz = m_ppos[i].m_z - samplePoint.m_z;
        const float dstSq = dx*dx + dy*dy + dz*dz;

        if (dstSq >= radiusSq) continue;

        const float dst = std::sqrt(dstSq);
        const float influence = SmoothingKernel(m_smoothingRadius, dst);
        const ngl::Vec3 position3D(m_ppos[i].m_x, m_ppos[i].m_y, m_ppos[i].m_z);
        const float density = CalculateDensity(position3D);

        // Avoid division by zero
        if (density > 0.0f) {
            property += particleProperties[i] * m_particleMass * influence / density;
        }
    }

    return property;
}

ngl::Vec3 Emitter::CalculatePropertyGradient(const ngl::Vec3& samplePoint,
                                           const std::vector<float>& particleProperties)
{
    ngl::Vec3 gradient(0.0f, 0.0f, 0.0f);
    const float radiusSq = m_smoothingRadius * m_smoothingRadius;

    for (size_t i = 0; i < m_maxParticles; ++i)
    {
        const ngl::Vec3 displacement = ngl::Vec3(m_ppos[i].m_x, m_ppos[i].m_y, m_ppos[i].m_z) - samplePoint;
        const float dstSq = displacement.lengthSquared();

        if (dstSq >= radiusSq || dstSq == 0.0f)
            continue;

        const float dst = std::sqrt(dstSq);
        const ngl::Vec3 dir = displacement / dst;  // Normalized direction vector
        const float slope = SmoothingKernelDerivative(m_smoothingRadius, dst);
        const float density = CalculateDensity(ngl::Vec3(m_ppos[i].m_x, m_ppos[i].m_y, m_ppos[i].m_z));

        // Avoid division by zero and accumulate gradient
        if (density > 0.0f)
        {
            gradient += -particleProperties[i] * dir * slope * m_particleMass / density;
        }
    }

    return gradient;
}
void Emitter::update(float _dt)
{
    // Calculate densities first
    calculateAllDensities();

    // Apply visualizations if enabled
    if(m_showDensity)
        visualizeDensities();
    else if(m_showSmoothing)
        visualizeSmoothing();
    // else
    // {
    //     // Default appearance when not visualizing
    //     for(size_t i = 0; i < m_maxParticles; ++i)
    //     {
    //         m_pcolour[i] = ngl::Random::getRandomColour3();
    //         m_psize[i] = 2.0f;
    //     }
    // }


    const ngl::Vec3 gravity(0.0f, -9.81f, 0.0f);

    auto numAlive = std::count_if(std::begin(m_pstate), std::end(m_pstate), [](auto p){return p == ParticleState::Active;});

    // if (numAlive < m_maxAlive)
    // {
    //     birthParticles();
    // }
//#pragma omp parallel for
    for(size_t i = 0; i < m_maxParticles; ++i)
    {
        // if(m_pstate[i] == ParticleState::Dead)
        //     continue;
        if (!m_simulate)
        {

            m_psize[i] += 0.1f;
            m_psize[i] = std::clamp(m_psize[i], 0.0f, 4.0f);
            m_ppos[i].m_w = m_psize[i];
            m_pcolour[i] = ngl::Random::getRandomColour3();
        }
        else
        {
            ngl::Vec3 emitDir(1.0f, 2.0f, 0.0f);
            m_pdir[i] += gravity * _dt * 0.5f;
            m_ppos[i] += m_pdir[i] * 0.5f;
            m_psize[i] += 0.1f;
            m_psize[i] = std::clamp(m_psize[i], 0.0f, 5.0f);
            m_ppos[i].m_w = m_psize[i];
            m_pcolour[i] = ngl::Random::getRandomColour3();
            resolveCollisions(i);
        }
        // if (--m_plife[i] <= 0 || m_ppos[i].m_y < 0.0f)
        // {
        //     resetParticle(i);
        // }
    }

}



void Emitter::resolveCollisions(size_t _i)
{
    ngl::Vec3 boundSize(50.0f, 50.0f, 50.0f); // Proper bounding box
    float collisionDamping = 0.7f; // Controls bounce effect
    ngl::Vec3 halfBoundSize = boundSize * 0.5f - ngl::Vec3(1.0f, 1.0f, 1.0f) * m_psize[_i];

    // X-axis collision
    if (std::abs(m_ppos[_i].m_x) > halfBoundSize.m_x)
    {
        m_ppos[_i].m_x = halfBoundSize.m_x * (m_ppos[_i].m_x > 0 ? 1.0f : -1.0f);
        m_pdir[_i].m_x *= -collisionDamping;
    }

    // Floor collision at y = 0
     if (m_ppos[_i].m_y <= 0.0f)
     {
         m_ppos[_i].m_y = 0.0f;  // Stick to the floor
         m_pdir[_i].m_y *= -collisionDamping;  // Reverse velocity

    //     // Small upward push to avoid getting stuck
    //     if (std::abs(m_pdir[_i].m_y) < 0.1f)
    //     {
    //         m_pdir[_i].m_y += 0.5f;
    //     }
     }

    // Floor collision at y = 0
    //if (m_ppos[_i].m_y <= 0.0f)
    //{
    //   m_ppos[_i].m_y = 0.0f;

        // Reverse velocity and apply damping
     //   m_pdir[_i].m_y = -std::abs(m_pdir[_i].m_y) * collisionDamping;

        // // Prevent sticking by setting a minimum bounce velocity
        // if (std::abs(m_pdir[_i].m_y) < 1.0f)
        // {
        //     m_pdir[_i].m_y = 1.0f; // minimum bounce impulse
        // }
    //}

    // Ceiling collision at upper bound
    if (m_ppos[_i].m_y > halfBoundSize.m_y)
    {
        m_ppos[_i].m_y = halfBoundSize.m_y;
        m_pdir[_i].m_y *= -collisionDamping;
    }

    // Z-axis collision
    if (std::abs(m_ppos[_i].m_z) > halfBoundSize.m_z)
    {
        m_ppos[_i].m_z = halfBoundSize.m_z * (m_ppos[_i].m_z > 0 ? 1.0f : -1.0f);
        m_pdir[_i].m_z *= -collisionDamping;
    }
}

void Emitter::calculateAllDensities()
{
    m_maxDensity = 0.0f;
    m_densities.resize(m_maxParticles);

    for(size_t i = 0; i < m_maxParticles; ++i)
    {
        m_densities[i] = CalculateDensity(ngl::Vec3(m_ppos[i].m_x, m_ppos[i].m_y, m_ppos[i].m_z));
        if(m_densities[i] > m_maxDensity)
        {
            m_maxDensity = m_densities[i];
        }
    }
}

void Emitter::visualizeDensities()
{
    if(m_densities.empty() || m_maxDensity <= 0.0f)
        return;

    for(size_t i = 0; i < m_maxParticles; ++i)
    {
        // Normalize density for visualization
        float normalizedDensity = m_densities[i] / m_maxDensity;

        // Color based on normalized density (blue = low, red = high)
        m_pcolour[i].m_r = normalizedDensity;   // Red channel (higher density = more red)
        m_pcolour[i].m_g = 0.2f;                // Green channel (constant low)
        m_pcolour[i].m_b = 1.0f - normalizedDensity; // Blue channel (higher density = less blue)

        // Optional: Scale size by density to make the high-density particles bigger
        m_psize[i] = 2.0f + 8.0f * normalizedDensity; // Size increases with density
    }
}

void Emitter::visualizeSmoothing()
{
    const float radius = m_smoothingRadius;

    for(size_t i = 0; i < m_maxParticles; ++i)
    {
        // Color particles based on smoothing influence
        float selfInfluence = SmoothingKernel(radius, 0.0f);
        float normalizedInfluence = selfInfluence / SmoothingKernel(radius, 0.0f);

        m_pcolour[i].m_r = normalizedInfluence;
        m_pcolour[i].m_g = normalizedInfluence;
        m_pcolour[i].m_b = 1.0f;

        // Show smoothing radius
        m_psize[i] = radius * 0.1f; // Scale for visibility
    }
}

float Emitter::convertDensityToPressure(float _density) const
{
    const float densityError = _density - m_targetDensity;
    return densityError * m_pressureMultiplier;
}

float Emitter::calculateSharedPressure(float densityA, float densityB) const
{
    float pressureA = convertDensityToPressure(densityA);
    float pressureB = convertDensityToPressure(densityB);
    return (pressureA + pressureB) / 2.0f;
}

ngl::Vec3 Emitter::calculatePressureForce(size_t _particleIndex) const
{
    ngl::Vec3 pressureForce(0.0f, 0.0f, 0.0f);
    const float radiusSq = m_smoothingRadius * m_smoothingRadius;

    // Convert current particle position to Vec3
    const ngl::Vec3 currentPos(m_ppos[_particleIndex].m_x,
                              m_ppos[_particleIndex].m_y,
                              m_ppos[_particleIndex].m_z);

    for(size_t otherIndex = 0; otherIndex < m_maxParticles; ++otherIndex)
    {
        if(_particleIndex == otherIndex) continue;

        // Convert other particle position to Vec3
        const ngl::Vec3 otherPos(m_ppos[otherIndex].m_x,
                                m_ppos[otherIndex].m_y,
                                m_ppos[otherIndex].m_z);

        const ngl::Vec3 offset = otherPos - currentPos;
        const float dstSq = offset.lengthSquared();

        if(dstSq >= radiusSq || dstSq == 0.0f) continue;

        const float dst = std::sqrt(dstSq);
        const ngl::Vec3 dir = offset / dst;
        const float slope = SmoothingKernelDerivative(m_smoothingRadius, dst);
        const float density = m_densities[otherIndex];
        const float sharedPressure = calculateSharedPressure(density, m_densities[_particleIndex]);

        pressureForce += sharedPressure * dir * slope * m_particleMass / density;
    }

    return pressureForce;
}

void Emitter::simulationStep(float _dt)
{
    // 1. Update densities
    calculateAllDensities();

    // 2. Calculate forces and update velocities
    for(size_t i = 0; i < m_maxParticles; ++i)
    {
        // Gravity
        const ngl::Vec3 gravityForce(0.0f, -9.81f * m_particleMass, 0.0f);

        // Pressure force
        const ngl::Vec3 pressureForce = calculatePressureForce(i);

        // Total force
        const ngl::Vec3 totalForce = pressureForce + gravityForce;

        // Update velocity (F = ma => a = F/m)
        m_pdir[i] += (totalForce / m_particleMass) * _dt;
    }

    // 3. Update positions and handle collisions
    for(size_t i = 0; i < m_maxParticles; ++i)
    {
        m_ppos[i] += m_pdir[i] * _dt;
        resolveCollisions(i);

        // Visual updates
        m_psize[i] = std::clamp(2.0f + m_densities[i]/1000.0f, 1.0f, 5.0f);
        m_ppos[i].m_w = m_psize[i];
    }
}

// void Emitter::resetParticle(size_t _i)
// {
//     ngl::Vec3 emitDir(1.0f, 2.0f, 0.0f);
//     float spread = 5.5f;
//     m_ppos[_i].set(m_pos.m_x, m_pos.m_y ,m_pos.m_z, 0.0f);
//     m_pdir[_i] = emitDir; // * ngl::Random::randomPositiveNumber() + randomVectorOnSphere() * spread;
//     m_pdir[_i].m_y = std::abs(m_pdir[_i].m_y);
//     m_psize[_i] = 2.0f;
//     m_plife[_i] = 20 + static_cast<int>(ngl::Random::randomPositiveNumber(100));
//     m_pcolour[_i] = ngl::Random::getRandomColour3();
//     m_pstate[_i] = ParticleState::Dead;
//
//
//
// }

// ngl::Vec3 Emitter::randomVectorOnSphere(float _radius)
// {
//     auto phi = ngl::Random::randomPositiveNumber(M_PI * 2.0f);
//     auto costheta = ngl::Random::randomNumber();
//     auto u = ngl::Random::randomNumber();
//     auto theta = std::acos(costheta);
//     auto r = _radius * std::cbrt(u);
//     return ngl::Vec3(r * std::sin(theta) * cos(phi), r * std::sin(theta) * sin(phi), r * std::cos(theta));
// }

// void Emitter::initializeParticles()
// {
//     int particleSpacing = 1.0f; // Space between particles
//     int jitterStrength = 0.1f;
//     int particlesPerRow = static_cast<int>(cbrt(m_maxParticles));   // number of particles per row
//     int particlesPerYColumn = (m_maxParticles - 1) / (particlesPerRow + 1);
//     int particlesPerZColumn = (m_maxParticles - 1) / (particlesPerRow + 1);// determine column count
//     float spacing = m_psize[0] * 2 + particleSpacing;   // define spacing between particles
//
//     int i = 0;
//
//     // Iterate over particles per row, column, and axis
//     for (int x = 0; x < particlesPerRow && i < m_maxParticles; ++x)
//     {
//         for (int y = 0; y < particlesPerColumn && i < m_maxParticles; ++y)
//         {
//             for (int z = 0; z < particlesPerColumn && i < m_maxParticles; ++z)
//             {
//                 // Calculate position with some random jitter and spacing
//                 float tx = (x / (float)(particlesPerRow - 1)) - 0.5f;
//                 float ty = (y / (float)(particlesPerColumn - 1)) - 0.5f;
//                 float tz = (z / (float)(particlesPerColumn - 1)) - 0.5f;
//
//                 // Incorporate spacing between particles by adding the spacing to positions
//                 float px = tx * (spacing) + m_ppos[i].m_x;
//                 float py = ty * (spacing) + m_ppos[i].m_y;
//                 float pz = tz * (spacing) + m_ppos[i].m_z;
//
//                 // Adding some random jitter to the particle position
//                 ngl::Vec3 jitter = ngl::Vec3(
//                     ngl::Random::randomNumber() * jitterStrength,  // X-axis jitter
//                     ngl::Random::randomNumber() * jitterStrength,  // Y-axis jitter
//                     ngl::Random::randomNumber() * jitterStrength   // Z-axis jitter
//                 );
//
//                 m_ppos[i] = ngl::Vec3(px, py, pz) + jitter;
//
//                 // Assigning initial velocity (could be randomized or predefined)
//                 m_pdir[i] = ngl::Vec3(
//                     ngl::Random::randomNumber(1.0f), // Small lateral jitter for X velocity
//                     ngl::Random::randomPositiveNumber(2.0f) + 2.0f, // Upward velocity for Y
//                     ngl::Random::randomNumber(1.0f)  // Small lateral jitter for Z velocity
//                 );
//
//                 m_psize[i] = 2.0f;  // Assign a default size for the particle
//                 m_plife[i] = 100;   // Can be randomized if desired
//                 m_pcolour[i] = ngl::Random::getRandomColour3();  // Random color for each particle
//                 m_pstate[i] = ParticleState::Active;  // Mark particle as active
//
//                 ++i;  // Increment the index to move to the next particle
//             }
//         }
//     }
// }
void Emitter::initializeParticles()
{
    const float cubeSideLength = 25.0f; // Fixed 1000x1000x1000 volume
    const float jitterStrength = 0.1f;
    const float initialLife = 100.0f;
    const float baseVelocityY = 2.0f;
    const float particleSpacing = 2.0f; // Adjusts density


    // Calculate particles per axis based on spacing and size
    const float effectiveParticleDiameter = m_psize[0] * 2 + particleSpacing;
    const int particlesPerAxis = static_cast<int>(cubeSideLength / effectiveParticleDiameter);

    // Calculate exact particle count that fits in the volume
    const size_t particlesThatFit = particlesPerAxis * particlesPerAxis * particlesPerAxis;

    // Update max particles to match what actually fits
    m_maxParticles = std::min(m_maxParticles, particlesThatFit);

    // Recalculate actual side length to maintain exact volume
    const float actualSideLength = particlesPerAxis * effectiveParticleDiameter;
    const float halfLength = actualSideLength * 0.5f;

    // Set uniform particle size
    for (size_t i = 0; i < m_maxParticles; ++i) {
        m_psize[i] = 4.0f; // Or your preferred size
    }

    // Place particles in perfect grid
    size_t particleIndex = 0;
    for (int x = 0; x < particlesPerAxis && particleIndex < m_maxParticles; ++x) {
        for (int y = 0; y < particlesPerAxis && particleIndex < m_maxParticles; ++y) {
            for (int z = 0; z < particlesPerAxis && particleIndex < m_maxParticles; ++z) {
                // Calculate grid position (centered around origin)
                float px = x * effectiveParticleDiameter - halfLength;
                float py = y * effectiveParticleDiameter - halfLength;
                float pz = z * effectiveParticleDiameter - halfLength;

                // Add slight randomness
                m_ppos[particleIndex] = ngl::Vec3(
                    px + ngl::Random::randomNumber() * jitterStrength,
                    py + ngl::Random::randomNumber() * jitterStrength,
                    pz + ngl::Random::randomNumber() * jitterStrength
                );

                // Initial velocity (mostly upward)
                m_pdir[particleIndex] = ngl::Vec3(
                    ngl::Random::randomNumber(1.0f),
                    ngl::Random::randomPositiveNumber(2.0f) + baseVelocityY,
                    ngl::Random::randomNumber(1.0f)
                );

                m_plife[particleIndex] = initialLife;
                m_pcolour[particleIndex] = ngl::Random::getRandomColour3();
                m_pstate[particleIndex] = ParticleState::Active;

                particleIndex++;
            }
        }
    }

    CalculateDensity(ngl::Vec3(0.0f, 0.0f, 0.0f));
    CalculateProperty(ngl::Vec3(0.0f,0.0f,0.0f), m_densities);


    // No leftover particles will be placed - m_maxParticles now equals exact grid count
}


void Emitter::render(int _width, int _height) const
{
    //std::cout<<"Render\n";
    // auto view = ngl::lookAt({0,20,20},{0,0,0}, {0,1,0});
    // auto project = ngl::perspective(45.0f, 1.0f, 0.1f, 200.0f);
    // //ngl::ShaderLib::use(ngl::nglColourShader);
    // ngl::ShaderLib::setUniform("MVP", project * view);
    // glPointSize(4.0f);
    m_vao->bind();
    m_vao->setData(0, ngl::MultiBufferVAO::VertexData(m_ppos.size() * sizeof(ngl::Vec3), m_ppos[0].m_x));
    m_vao->setVertexAttributePointer(0, 4, GL_FLOAT, 0, 0);

    m_vao->setData(1, ngl::MultiBufferVAO::VertexData(m_pcolour.size() * sizeof(ngl::Vec3), m_pcolour[0].m_x));
    //m_vao->setVertexAttributePointer(1, 3, GL_FLOAT, sizeof(Particle), 6);
    m_vao->setVertexAttributePointer(1, 3, GL_FLOAT, 0, 0);

    // Density data (location = 2) - Corrected version
    if(!m_densities.empty())
    {
        m_vao->setData(2, ngl::MultiBufferVAO::VertexData(
            m_densities.size() * sizeof(float),
            m_densities[0]  // Pass first element by reference
        ));
        m_vao->setVertexAttributePointer(2, 1, GL_FLOAT, 0, 0);
    }

    m_vao->setNumIndices(m_maxParticles);
    glEnable(GL_PROGRAM_POINT_SIZE);
    m_vao->draw();
    glDisable(GL_PROGRAM_POINT_SIZE);


    m_vao->unbind();
    //ngl::Transformation tx;
    // for(auto &p : m_particles)
    // {
    //     tx.setPosition(p.pos);
    //     //tx.setScale(p.size, p.size, p.size);
    //     ngl::ShaderLib::setUniform("MVP", project * view * tx.getMatrix());
    //     ngl::ShaderLib::setUniform("Colour", p.colour.m_r,  p.colour.m_g, p.colour.m_b, 1.0f);
    //     ngl::VAOPrimitives::draw("sphere", GL_POINTS);
    // }

    // Add debug text

    if(m_showDensity)
    {
        ngl::Text text("Arial", 20); // Font and size
        text.setColour(1,1,1);
        text.renderText(10, 50, "Density: " + std::to_string(m_maxDensity));
    }
    else if(m_showSmoothing)
    {
        ngl::Text text("Arial", 20);
        text.setColour(1,1,1);
        text.renderText(10, 50, "Radius: " + std::to_string(m_smoothingRadius));
    }
}


void Emitter::move(float _dx, float _dy, float _dz)
{
    m_pos.m_x += _dx;
    m_pos.m_y += _dy;
    m_pos.m_z += _dz;
}



