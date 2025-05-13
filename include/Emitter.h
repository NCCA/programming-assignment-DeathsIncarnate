#ifndef EMITTER_H_
#define EMITTER_H_
#include <vector>
#include "Particle.h"
#include <ngl/Vec3.h>
#include <ngl/AbstractVAO.h>
#include <memory>
#include <ngl/MultiBufferVAO.h>
#include <ngl/VAOFactory.h>
#include <ngl/Vec4.h>
#include "BoundingBox.h"
#include "Physics.h"

class Emitter
{
    public:
    Emitter(size_t _num, size_t _maxALive, size_t _numPerFrame, ngl::Vec3 _pos);
    ~Emitter() = default;

    size_t size() const;
    void update(float _dt);
    void render(int _width, int _height) const;
    void move(float _dx, float _dy, float _dz);
    // float CalculateDensity(const ngl::Vec3 &samplePoint);
    float CalculateDensity(const ngl::Vec3& samplePoint);
    bool m_simulate = false; // default: no simulation

    void calculateAllDensities();
    ngl::Vec3 calculatePressureForce(size_t _particleIndex) const;
    float convertDensityToPressure(float _density) const;
    float calculateSharedPressure(float densityA, float densityB) const;
    ngl::Vec3 calculateViscosityForce(size_t particleIndex) const;
    void simulationStep(float _dt);

    void setCursorPos(const ngl::Vec3 &_pos, float _radius, float _strength);
    void toggleCursorInteraction() { m_cursorInteraction = !m_cursorInteraction; }




    private:
    BoundingBox m_boundingBox;
    Physics m_physics;

    //std::vector<Particle> m_particles;
    std::vector<ngl::Vec4> m_ppos;
    std::vector<ngl::Vec3> m_pdir;
    std::vector<ngl::Vec3> m_pcolour;
    std::vector<float> m_psize;
    std::vector<int> m_plife;
    enum class ParticleState : bool {Active, Dead};
    std::vector<ParticleState> m_pstate;
    size_t m_maxParticles;
    size_t m_maxAlive = 1000;
    size_t m_numPerFrame = 120;

    void initializeParticles();

    std::unique_ptr<ngl::MultiBufferVAO> m_vao;

    // Particle data


    std::vector<float> m_pdensity;     // Densities


    // // SPH parameters
    // float m_smoothingRadius = 3.0f;   // Adjust based on scale
    // float m_particleMass = 1.0f;      // Typically 1.0 for simplicity
    // float m_restDensity = 1000.0f;    // Water-like density
    // std::vector<float> m_densities;
    // float m_maxDensity = 0.0f;
    // // Add these to your Emitter.h in the private section:
    // float m_targetDensity = 1000.0f; // Typical water density kg/m³
    // float m_pressureMultiplier = 0.1f; // Stiffness constant
    // float m_viscosityStrength = 0.01f;
    // float m_damping = 0.01f;
    // const float m_particleSpacing = 1.5f;


    // float m_pressureStiffness = 100.0f;   // Adjusts fluid compressibility
    // float m_viscosity = 0.1f;            // For damping oscillations

    ngl::Vec3 m_cursorPos;
    float m_cursorRadius = 5.0f;
    float m_cursorStrength = 10.0f;
    bool m_cursorInteraction = false;

    ngl::Vec3 calculateCursorForce(size_t particleIdx) const;


    struct SpatialEntry {
        size_t particleIndex;
        size_t cellKey;
    };

    std::vector<SpatialEntry> m_spatialLookup;
    std::vector<size_t> m_startIndices;
    float m_gridCellSize;

    // Spatial hashing methods
    void UpdateSpatialLookup();
    std::array<int, 3> PositionToCellCoord(const ngl::Vec3& point) const;
    size_t HashCell(int cellX, int cellY, int cellZ) const;
    void ForEachNeighbor(size_t particleIdx, std::function<void(size_t)> callback);
    size_t GetKeyFromHash(size_t hash) const;



    //         // Physical properties
    //         float particleMass = 0.02f;
    //         float restDensity = 1000.0f;
    //         float smoothingRadius = 0.1f;
    //         float particleSpacing = 0.05f;
    //
    //         // Force coefficients
    //         float pressureMultiplier = 100.0f;
    //         float viscosityStrength = 0.1f;
    //         float surfaceTension = 0.072f;
    //
    //         // Simulation control
    //         float timeStep = 0.001f;
    //         int maxIterations = 5;
    //         float boundaryDamping = 0.5f;
    //     };


};

#endif