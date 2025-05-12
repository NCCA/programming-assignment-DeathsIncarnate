#ifndef EMITTER_H_
#define EMITTER_H_
#include <vector>
#include <string_view>
#include "Particle.h"
#include <ngl/Vec3.h>
#include <ngl/AbstractVAO.h>
#include <memory>
#include <ngl/MultiBufferVAO.h>
#include <ngl/VAOFactory.h>
#include <ngl/Vec4.h>

class Emitter
{
    public:
    Emitter(size_t _num, size_t _maxALive, size_t _numPerFrame, ngl::Vec3 _pos);
    size_t size() const;
    void update(float _dt);
    void render(int _width, int _height) const;
    void move(float _dx, float _dy, float _dz);
    // float CalculateDensity(const ngl::Vec3 &samplePoint);
    float CalculateDensity(const ngl::Vec3 &samplePoint);
    float CalculateProperty(const ngl::Vec3& samplePoint,
                           const std::vector<float>& particleProperties);
    ngl::Vec3 CalculatePropertyGradient(const ngl::Vec3& samplePoint,
                                               const std::vector<float>& particleProperties);
    bool m_simulate = false; // default: no simulation

    void calculateAllDensities();
    void visualizeDensities();
    void visualizeSmoothing();
    void toggleDensityVisualization() { m_showDensity = !m_showDensity; }
    void toggleSmoothingVisualization() { m_showSmoothing = !m_showSmoothing; }
    ngl::Vec3 calculatePressureForce(size_t _particleIndex) const;
    float convertDensityToPressure(float _density) const;
    float calculateSharedPressure(float densityA, float densityB) const;
    ngl::Vec3 calculateViscosityForce(size_t particleIndex) const;
    void simulationStep(float _dt);
    void renderBoundingBox() const;

    void setCursorPos(const ngl::Vec3 &_pos, float _radius, float _strength);
    void toggleCursorInteraction() { m_cursorInteraction = !m_cursorInteraction; }



    private:
    ngl::Vec3 m_pos;
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

    void resolveCollisions(size_t _i);
    void resetParticle(size_t _i);
    void birthParticles();
    void initializeParticles();
    ngl::Vec3 randomVectorOnSphere(float _radius = 1.0f);
    std::unique_ptr<ngl::MultiBufferVAO> m_vao;

    // Particle data
    std::vector<ngl::Vec3> m_pvel;    // Velocities
    std::vector<float> m_pdensity;     // Densities
    std::vector<float> m_ppressure;    // Pressures

    // SPH parameters
    float m_smoothingRadius = 100.0f;   // Adjust based on scale
    float m_particleMass = 1.0f;      // Typically 1.0 for simplicity
    float m_restDensity = 1000.0f;    // Water-like density
    std::vector<float> m_densities;
    bool m_showDensity = false;
    bool m_showSmoothing = false;
    float m_maxDensity = 0.0f;
    // Add these to your Emitter.h in the private section:
    float m_targetDensity = 1000.0f; // Typical water density kg/m³
    float m_pressureMultiplier = 1.0f; // Stiffness constant
    float m_viscosityStrength = 0.1f;
    float m_damping = 0.01f;

    std::unique_ptr<ngl::MultiBufferVAO> m_boxVAO;
    void initBoundingBoxVAO();

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




};

#endif