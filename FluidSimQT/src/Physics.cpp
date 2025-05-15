//
// Created by s5616052 on 13/05/25.
//
#include "Physics.h"

#include <iostream>
#include <ngl/Random.h>
#include <algorithm>
#include <ngl/VAOPrimitives.h>
#include <ngl/ShaderLib.h>
#include <ngl/Transformation.h>

Physics::Physics(std::vector<ngl::Vec4>& m_ppos, std::vector<ngl::Vec3>& m_pdir)
  :m_ppos(m_ppos), m_pdir(m_pdir)
{
}


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

float Physics::CalculateDensity(const ngl::Vec3& samplePoint, size_t m_maxParticles)
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

ngl::Vec3 Physics::calculateViscosityForce(size_t _particleIndex, size_t m_maxParticles) const
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
        viscosityForce += force;  // Our particle gets the opposite force
        // Note: The other particle will get +force when its turn comes
    }

    return viscosityForce;
}

void Physics::calculateAllDensities(size_t m_maxParticles)
{
    m_maxDensity = 0.0f;
    m_densities.resize(m_maxParticles);

    #pragma omp parallel for reduction(max:m_maxDensity)
    for(size_t i = 0; i < m_maxParticles; ++i)
    {
        m_densities[i] = CalculateDensity(ngl::Vec3(m_ppos[i].m_x, m_ppos[i].m_y, m_ppos[i].m_z), m_maxParticles);
        if(m_densities[i] > m_maxDensity)
        {
            m_maxDensity = m_densities[i];
        }
    }
}

float Physics::convertDensityToPressure(float _density) const
{
    const float densityError = _density - m_targetDensity;
    return densityError * m_pressureMultiplier;
}

// float Physics::convertNearDensityToNearPressure(float _density) const
// {
//     const float densityError = _density - m_targetDensity;
//     return nearDensity * m_nearPressureMultiplier;
// }


// float Emitter::convertDensityToPressure(float _density) const {
//     const float ratio = _density / m_restDensity;
//     return m_pressureMultiplier * (pow(ratio, 7) - 1.0f); // Tait equation
// }

float Physics::calculateSharedPressure(float densityA, float densityB) const
{
    float pressureA = convertDensityToPressure(densityA);
    float pressureB = convertDensityToPressure(densityB);
    return (pressureA + pressureB) / 2.0f;
}

ngl::Vec3 Physics::calculatePressureForce(size_t _particleIndex, size_t m_maxParticles) const
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

