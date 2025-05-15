//
// Created by s5616052 on 13/05/25.
//

#ifndef PHYSICS_H
#define PHYSICS_H


#include <ngl/VAOPrimitives.h>


#pragma once
#include <ngl/Vec3.h>
#include <vector>

class Physics
{
public:
    Physics(std::vector<ngl::Vec4>& m_ppos, std::vector<ngl::Vec3>& m_pdir);
    ~Physics() = default; // Destructor
    float CalculateDensity(const ngl::Vec3& samplePoint, size_t m_maxParticles);
    ngl::Vec3 calculateViscosityForce(size_t _particleIndex, size_t m_maxParticles) const;
    void calculateAllDensities(size_t m_maxParticles);
    float convertDensityToPressure(float _density) const;
    float calculateSharedPressure(float densityA, float densityB) const;
    ngl::Vec3 calculatePressureForce(size_t _particleIndex, size_t m_maxParticles) const;

    float m_particleMass = 1.0f;      // Typically 1.0 for simplicity
    float m_restDensity = 1000.0f;    // Water-like density
    std::vector<float> m_densities;
    float m_maxDensity = 0.0f;
    // Add these to your Emitter.h in the private section:
    float m_targetDensity = 1000.0f; // Typical water density kg/m³
    float m_pressureMultiplier = 0.15f; // Stiffness constant
    float m_viscosityStrength = 0.115f;
    float m_damping = 0.01f;
    float m_particleSpacing = 1.5f;


    float m_smoothingRadius = 2 * m_particleSpacing;   // Adjust based on scale

    size_t m_maxParticles;


private:

    std::vector<ngl::Vec4>& m_ppos; // Reference to particle positions
    std::vector<ngl::Vec3>& m_pdir; // Reference to particle directions

};

#endif //PHYSICS_H
