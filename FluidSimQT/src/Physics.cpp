//
// Created by s5616052 on 13/05/25.
//
#include "Physics.h"
#include <iostream>
#include <ngl/Random.h>
#include <ngl/ShaderLib.h>

Physics::Physics(std::vector<ngl::Vec4>& m_ppos, std::vector<ngl::Vec3>& m_pdir)
  :m_ppos(m_ppos), m_pdir(m_pdir)
{
}


static float SmoothingKernel(float _radius, float _dst)
{
    if (_dst >= _radius) return 0.0f;

    const float c_h = _radius;
    const float c_h3 = c_h * c_h * c_h;
    const float c_q = _dst / c_h;

    // Precomputed normalization constant
    const float c_norm = 8.0f / (M_PI * c_h3);

    if (c_q <= 0.5f)
    {
        const float c_q2 = c_q * c_q;
        const float c_q3 = c_q2 * c_q;
        return c_norm * (1.0f - 6.0f * c_q2 + 6.0f * c_q3);
    }
    else
    {
        const float c_term = 2.0f - c_q;
        return c_norm * (0.5f * c_term * c_term * c_term);
    }
}

// Optimized derivative of cubic spline kernel
static float SmoothingKernelDerivative(float _radius, float _dst)
{
    if (_dst >= _radius || _dst == 0.0f) return 0.0f;

    const float c_h = _radius;
    const float c_h4 = c_h * c_h * c_h * c_h;
    const float c_q = _dst / c_h;

    // Precomputed normalization constant
    const float c_norm = 8.0f / (M_PI * c_h4);

    if (c_q <= 0.5f)
    {
        return c_norm * (3.0f * c_q - 4.5f * c_q * c_q);
    }
    else
    {
        const float c_term = 2.0f - c_q;
        return c_norm * (-1.5f * c_term * c_term);
    }
}

float Physics::CalculateDensity(const ngl::Vec3& c_samplePoint, size_t m_maxParticles)
{
    float density = 0.0f;
    const float c_radiusSq = m_smoothingRadius * m_smoothingRadius;

    for (size_t i = 0; i < m_maxParticles; ++i)
    {
        const float c_dx = m_ppos[i].m_x - c_samplePoint.m_x;
        const float c_dy = m_ppos[i].m_y - c_samplePoint.m_y;
        const float c_dz = m_ppos[i].m_z - c_samplePoint.m_z;
        const float c_dstSq = c_dx*c_dx + c_dy*c_dy + c_dz*c_dz;

        if (c_dstSq >= c_radiusSq) continue;

        const float c_dst = std::sqrt(c_dstSq);
        density += m_particleMass * SmoothingKernel(m_smoothingRadius, c_dst);
    }
    //std::cout <<"Density" << density << "\n";
    return density;
}

ngl::Vec3 Physics::calculateViscosityForce(size_t _particleIndex, size_t m_maxParticles) const
{
    ngl::Vec3 viscosityForce(0.0f, 0.0f, 0.0f);
    const float c_radiusSq = m_smoothingRadius * m_smoothingRadius;

    const ngl::Vec3 c_currentPos(m_ppos[_particleIndex].m_x,
                          m_ppos[_particleIndex].m_y,
                          m_ppos[_particleIndex].m_z);

    for(size_t otherIndex = 0; otherIndex < m_maxParticles; ++otherIndex)
    {
        if(_particleIndex == otherIndex) continue;

        const ngl::Vec3 c_otherPos(m_ppos[otherIndex].m_x,
                                m_ppos[otherIndex].m_y,
                                m_ppos[otherIndex].m_z);

        const ngl::Vec3 c_offset = c_otherPos - c_currentPos;
        const float c_dstSq = c_offset.lengthSquared();

        if(c_dstSq >= c_radiusSq || c_dstSq == 0.0f) continue;

        const float c_dst = std::sqrt(c_dstSq);
        const float c_influence = SmoothingKernel(m_smoothingRadius, c_dst);

        // Velocity difference
        const ngl::Vec3 c_velocityDiff = m_pdir[otherIndex] - m_pdir[_particleIndex];

        // The viscous force that THIS particle exerts on THE OTHER particle
        const ngl::Vec3 c_force = c_velocityDiff * (m_viscosityStrength * c_influence * m_particleMass / m_densities[otherIndex]);

        // Newton's Third Law application:
        viscosityForce += c_force;  // Our particle gets the opposite force
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
    const float c_densityError = _density - m_targetDensity;
    return c_densityError * m_pressureMultiplier;
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
    const float c_radiusSq = m_smoothingRadius * m_smoothingRadius;

    // Convert current particle position to Vec3
    const ngl::Vec3 c_currentPos(m_ppos[_particleIndex].m_x,
                              m_ppos[_particleIndex].m_y,
                              m_ppos[_particleIndex].m_z);

    for(size_t otherIndex = 0; otherIndex < m_maxParticles; ++otherIndex)
    {
        if(_particleIndex == otherIndex) continue;

        // Convert other particle position to Vec3
        const ngl::Vec3 c_otherPos(m_ppos[otherIndex].m_x,
                                m_ppos[otherIndex].m_y,
                                m_ppos[otherIndex].m_z);

        const ngl::Vec3 c_offset = c_otherPos - c_currentPos;
        const float c_dstSq = c_offset.lengthSquared();

        if(c_dstSq >= c_radiusSq || c_dstSq == 0.0f) continue;

        const float c_dst = std::sqrt(c_dstSq);
        const ngl::Vec3 c_dir = c_offset / c_dst;
        const float c_slope = SmoothingKernelDerivative(m_smoothingRadius, c_dst);
        const float c_sharedPressure = calculateSharedPressure(m_densities[_particleIndex], m_densities[otherIndex]);

        //The force that THIS particle exerts on THE OTHER particle
        const ngl::Vec3 c_force = c_sharedPressure * c_dir * c_slope * m_particleMass / m_densities[otherIndex];

        // Newton's Third Law application:
        pressureForce -= c_force;  // Our particle gets the opposite force
        // Note: The other particle will get +force when its turn comes in the simulation loop
    }

    return pressureForce;
}

