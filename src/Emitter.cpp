#include "Emitter.h"
#include <iostream>
#include <ngl/Random.h>
#include <algorithm>
#include <fstream>
#include <ngl/VAOPrimitives.h>
#include <ngl/VAOFactory.h>
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

    //initBoundingBoxVAO();
}
size_t Emitter::size() const
{
    return m_maxParticles;
}

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
    //std::cout <<"Density" << density << "\n";
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
//
// ngl::Vec3 Emitter::calculateViscosityForce(size_t _i) const
// {
//     ngl::Vec3 viscosityForce(0.0f, 0.0f, 0.0f);
//     const float radiusSq = m_smoothingRadius * m_smoothingRadius;
//
//
//     for (size_t j = 0; j < m_maxParticles; ++j)
//     {
//         if (_i == j) continue;
//
//         ngl::Vec3 offset = ngl::Vec3(m_ppos[j].m_x, m_ppos[j].m_y, m_ppos[j].m_z) -
//                    ngl::Vec3(m_ppos[_i].m_x, m_ppos[_i].m_y, m_ppos[_i].m_z);
//
//         float dstSq = offset.lengthSquared();
//         if (dstSq >= radiusSq || dstSq == 0.0f) continue;
//
//         float dst = std::sqrt(dstSq);
//         float influence = SmoothingKernel(m_smoothingRadius, dst);
//
//
//         viscosityForce += (m_pdir[j] - m_pdir[_i]) * (m_viscosityStrength * influence * m_particleMass / m_densities[j]);
//     }
//
//     return viscosityForce;
// }

ngl::Vec3 Emitter::calculateViscosityForce(size_t _particleIndex) const
{
    ngl::Vec3 viscosityForce(0.0f, 0.0f, 0.0f);
    const float radiusSq = m_smoothingRadius * m_smoothingRadius;

    const ngl::Vec3 currentPos(m_ppos[_particleIndex].m_x,
                          m_ppos[_particleIndex].m_y,
                          m_ppos[_particleIndex].m_z);

    for(size_t otherIndex = 0; otherIndex < m_maxParticles; ++otherIndex)
    {
        if(_particleIndex == otherIndex) continue;

        const ngl::Vec3 otherPos(m_ppos[otherIndex].m_x,
                                m_ppos[otherIndex].m_y,
                                m_ppos[otherIndex].m_z);

        const ngl::Vec3 offset = otherPos - currentPos;
        const float dstSq = offset.lengthSquared();

        if(dstSq >= radiusSq || dstSq == 0.0f) continue;

        const float dst = std::sqrt(dstSq);
        const float influence = SmoothingKernel(m_smoothingRadius, dst);

        // Velocity difference
        const ngl::Vec3 velocityDiff = m_pdir[otherIndex] - m_pdir[_particleIndex];

        // The viscous force that THIS particle exerts on THE OTHER particle
        const ngl::Vec3 force = velocityDiff * (m_viscosityStrength * influence * m_particleMass / m_densities[otherIndex]);

        // Newton's Third Law application:
        viscosityForce -= force;  // Our particle gets the opposite force
        // Note: The other particle will get +force when its turn comes
    }

    return viscosityForce;
}

void Emitter::update(float _dt)
{
    // Calculate densities first
    calculateAllDensities();

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

    #pragma omp parallel for
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
            ngl::Vec3 emitDir(0.0f, 2.0f, 0.0f);
            //m_pdir[i] += gravity * _dt * 0.5f;
            m_ppos[i] += m_pdir[i] * 0.5f;
            m_psize[i] += 0.1f;
            m_psize[i] = std::clamp(m_psize[i], 0.0f, 5.0f);
            m_ppos[i].m_w = m_psize[i];

            //ngl::Vec3 cursorForce = calculateCursorForce(i);

            m_pcolour[i] = ngl::Random::getRandomColour3();
            ngl::Vec3 pressureForce = calculatePressureForce(i); /// 10.0f;
            ngl::Vec3 viscosityForce = calculateViscosityForce(i); /// 10.0f;
            ngl::Vec3 gravityForce(0.0f, -9.81f * m_particleMass, 0.0f);
            ngl::Vec3 totalForce = pressureForce + gravityForce + viscosityForce * _dt; //+ (cursorForce*2000000) * _dt;
            m_pdir[i] = totalForce / m_particleMass / 200.0f; /// 10000.0f;

            resolveCollisions(i);

        }
        // if (--m_plife[i] <= 0 || m_ppos[i].m_y < 0.0f)
        // {
        //     resetParticle(i);
        // }
    }

}
//
//  void Emitter::update(float _dt)
//  {
// //     if (!m_simulate)
// //     {
// //         // Non-physics visualization mode
// //         for(size_t i = 0; i < m_maxParticles; ++i)
// //         {
// //             m_psize[i] += 0.1f;
// //             m_psize[i] = std::clamp(m_psize[i], 0.0f, 4.0f);
// //             m_ppos[i].m_w = m_psize[i];
// //             m_pcolour[i] = ngl::Random::getRandomColour3();
// //         }
// //         return;
// //     }
// //
//     if(m_simulate)
//     {
//         simulationStep(_dt);
//     }
//
//     // Visualization updates
//     if(m_showDensity)
//         visualizeDensities();
//     else if(m_showSmoothing)
//         visualizeSmoothing();
//  }



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
         // // Reverse velocity and apply damping
         // m_pdir[_i].m_y = -std::abs(m_pdir[_i].m_y) * collisionDamping;
         m_pdir[_i].m_y *= -collisionDamping;  // Reverse velocity

         // Small upward push to avoid getting stuck
         if (std::abs(m_pdir[_i].m_y) < 0.1f)
         {
             m_pdir[_i].m_y += 0.5f;
             // m_pdir[_i].m_y = 1.0f; // minimum bounce impulse
         }
     }


    // // Ceiling collision at upper bound
    // if (m_ppos[_i].m_y > halfBoundSize.m_y)
    // {
    //     m_ppos[_i].m_y = halfBoundSize.m_y;
    //     m_pdir[_i].m_y *= -collisionDamping;
    // }

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

    #pragma omp parallel for reduction(max:m_maxDensity)
    for(size_t i = 0; i < m_maxParticles; ++i)
    {
        m_densities[i] = CalculateDensity(ngl::Vec3(m_ppos[i].m_x, m_ppos[i].m_y, m_ppos[i].m_z));
        if(m_densities[i] > m_maxDensity)
        {
            m_maxDensity = m_densities[i];
        }
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
        const float sharedPressure = calculateSharedPressure(m_densities[_particleIndex], m_densities[otherIndex]);

        //The force that THIS particle exerts on THE OTHER particle
        const ngl::Vec3 force = sharedPressure * dir * slope * m_particleMass / m_densities[otherIndex];

        // Newton's Third Law application:
        pressureForce -= force;  // Our particle gets the opposite force
        // Note: The other particle will get +force when its turn comes in the simulation loop
    }

    return pressureForce;
}

// void Emitter::initBoundingBoxVAO()
// {
//     // Match your collision bounds from resolveCollisions()
//     const float width = 100.0f, height = 50.0f, depth = 100.0f;
//     const float hw = width/2, hh = height/2, hd = depth/2;
//
//     // Define the 8 corners of the box
//     std::vector<ngl::Vec3> vertices = {
//         // Bottom face
//         {-hw, 0, -hd}, {hw, 0, -hd},
//         {hw, 0, -hd}, {hw, 0, hd},
//         {hw, 0, hd}, {-hw, 0, hd},
//         {-hw, 0, hd}, {-hw, 0, -hd},
//
//         // Top face
//         {-hw, height, -hd}, {hw, height, -hd},
//         {hw, height, -hd}, {hw, height, hd},
//         {hw, height, hd}, {-hw, height, hd},
//         {-hw, height, hd}, {-hw, height, -hd},
//
//         // Vertical edges
//         {-hw, 0, -hd}, {-hw, height, -hd},
//         {hw, 0, -hd}, {hw, height, -hd},
//         {hw, 0, hd}, {hw, height, hd},
//         {-hw, 0, hd}, {-hw, height, hd}
//     };
//
//     // Correct VAO creation - no need for reset()
//     m_boxVAO = ngl::vaoFactoryCast<ngl::MultiBufferVAO>(ngl::VAOFactory::createVAO(ngl::multiBufferVAO, GL_LINES));
//
//     m_boxVAO->bind();
//     m_boxVAO->setNumIndices(vertices.size());
//     m_boxVAO->setData(2, ngl::SimpleVAO::VertexData(vertices.size()*sizeof(ngl::Vec3), vertices[0].m_x));
//     m_boxVAO->unbind();
// }
//
// void Emitter::renderBoundingBox() const
// {
//
//     glLineWidth(2.0f);
//     m_boxVAO->bind();
//     m_boxVAO->setVertexAttributePointer(2, 3, GL_FLOAT, 0, 0);
//     m_boxVAO->draw();
//     m_boxVAO->unbind();
//     glLineWidth(1.0f);
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
//
// void Emitter::setCursorPos(const ngl::Vec3 &_pos, float _radius, float _strength)
// {
//     m_cursorPos = _pos;
//     m_cursorRadius = _radius;
//     m_cursorStrength = _strength;
//     m_cursorInteraction = true; // ensure it's active
// }
//
// ngl::Vec3 Emitter::calculateCursorForce(size_t particleIdx) const
// {
//     if (!m_cursorInteraction) return ngl::Vec3(0.0f, 0.0f, 0.0f);
//
//     const ngl::Vec3 particlePos(m_ppos[particleIdx].m_x,
//                               m_ppos[particleIdx].m_y,
//                               m_ppos[particleIdx].m_z);
//
//     const ngl::Vec3 toCursor = particlePos - m_cursorPos;
//     const float dist = toCursor.length();
//
//     if (dist > m_cursorRadius || dist < 0.001f)
//         return ngl::Vec3(0.0f, 0.0f, 0.0f);
//
//     // Inverse square law for smooth repulsion
//     const float falloff = 1.0f - (dist / m_cursorRadius);
//     const float forceMagnitude = m_cursorStrength * falloff * falloff;
//
//     ngl::Vec3 forceDir = toCursor;
//     forceDir.normalize();
//
//     ngl::Vec3 force = forceDir * forceMagnitude;
//
//     // Debug output
//     std::cout << "Particle " << particleIdx << " dist=" << dist
//               << " force=(" << force.m_x << ", " << force.m_y << ", " << force.m_z << ")\n";
//
//
//     return force;
// }
//

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
}


void Emitter::move(float _dx, float _dy, float _dz)
{
    m_pos.m_x += _dx;
    m_pos.m_y += _dy;
    m_pos.m_z += _dz;
}



