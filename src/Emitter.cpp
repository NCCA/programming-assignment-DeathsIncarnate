#include "Emitter.h"
#include <iostream>
#include <ngl/Random.h>
#include <algorithm>
#include <ngl/VAOPrimitives.h>
#include <ngl/VAOFactory.h>
#include <ngl/ShaderLib.h>
#include <ngl/Transformation.h>
#include <ngl/MultiBufferVAO.h>
#include "Physics.h"




Emitter::Emitter(size_t _num, size_t _maxAlive, size_t _numPerFrame, ngl::Vec3 _pos):
m_maxParticles{_num}, m_maxAlive{_maxAlive}, m_numPerFrame{_numPerFrame}, m_ppos{_pos}, m_boundingBox()
{
    m_physics = std::make_unique<Physics>(m_ppos, m_pdir);
    m_maxParticles = _num;
    m_ppos.resize(m_maxParticles);
    m_pcolour.resize(m_maxParticles);
    m_pdir.resize(m_maxParticles);
    m_psize.resize(m_maxParticles);
    m_plife.resize(m_maxParticles);
    m_pstate.resize(m_maxParticles);


    m_physics->m_densities.resize(m_maxParticles);
    m_physics->m_maxDensity = 0.0f;



    // for (size_t i = 0; i < m_maxParticles; ++i)
    // {
    //     resetParticle(i);
    // }

    m_startIndices.resize(m_maxParticles);
    std::fill(m_startIndices.begin(), m_startIndices.end(), std::numeric_limits<size_t>::max());

    // Initialize SPH parameters
    m_physics->m_targetDensity = m_physics->m_restDensity;
    //m_smoothingRadius = 2.0f * m_particleSpacing; // Typically 2x particle spacing

    // Initialize spatial hashing
    m_gridCellSize = m_physics->m_smoothingRadius * 2.0f;
    m_spatialLookup.resize(m_maxParticles);
    m_startIndices.resize(m_maxParticles, std::numeric_limits<size_t>::max());

    initializeParticles();

    //m_vao = ngl::VAOFactory::createVAO(ngl::simpleVAO, GL_POINTS);
    m_vao = ngl::vaoFactoryCast<ngl::MultiBufferVAO>(ngl::VAOFactory::createVAO(ngl::multiBufferVAO, GL_POINTS));
    m_vao->bind();
    m_vao->setData(ngl::MultiBufferVAO::VertexData(0,0)); // index 0 points
    m_vao->setData(ngl::MultiBufferVAO::VertexData(0,0)); // index 1 colours
    m_vao->setData(ngl::MultiBufferVAO::VertexData(0,0)); // index 2 Densities
    m_vao->unbind();

    m_boundingBox.initBoundingBoxVAO();
}

// void Emitter::UpdateSpatialLookup()
// {
// #pragma omp parallel for
//     for(size_t i = 0; i < m_maxParticles; ++i) {
//         if(i >= m_ppos.size()) continue; // Safety check
//
//         auto cellCoord = PositionToCellCoord(ngl::Vec3(m_ppos[i].m_x, m_ppos[i].m_y, m_ppos[i].m_z));
//         size_t cellKey = HashCell(cellCoord[0], cellCoord[1], cellCoord[2]);
//         if(i < m_spatialLookup.size()) {
//             m_spatialLookup[i] = {i, cellKey};
//         }
//     }
//
//     // Sort safely
//     std::sort(m_spatialLookup.begin(), m_spatialLookup.begin() + std::min(m_maxParticles, m_spatialLookup.size()),
//         [](const SpatialEntry& a, const SpatialEntry& b) {
//             return a.cellKey < b.cellKey;
//         });
//
//     // Reset start indices
//     std::fill(m_startIndices.begin(), m_startIndices.end(), std::numeric_limits<size_t>::max());
//
//     // Build start indices with bounds checking
//     for(size_t i = 0; i < m_spatialLookup.size() && i < m_maxParticles; ++i) {
//         size_t key = m_spatialLookup[i].cellKey % m_startIndices.size();
//         if(i == 0 || m_spatialLookup[i].cellKey != m_spatialLookup[i-1].cellKey) {
//             if(key < m_startIndices.size()) {
//                 m_startIndices[key] = i;
//             }
//         }
//     }
// }
//
// std::array<int, 3> Emitter::PositionToCellCoord(const ngl::Vec3& point) const
// {
//     return {
//         static_cast<int>(std::floor(point.m_x / m_gridCellSize)),
//         static_cast<int>(std::floor(point.m_y / m_gridCellSize)),
//         static_cast<int>(std::floor(point.m_z / m_gridCellSize))
//     };
// }
//
// size_t Emitter::HashCell(int cellX, int cellY, int cellZ) const
// {
//     const size_t p1 = 73856093;
//     const size_t p2 = 19349663;
//     const size_t p3 = 83492791;
//     return (cellX * p1) ^ (cellY * p2) ^ (cellZ * p3);
// }
//
// size_t Emitter::GetKeyFromHash(size_t hash) const
// {
//     return hash % m_maxParticles;
// }
//
// void Emitter::ForEachNeighbor(size_t particleIdx, std::function<void(size_t)> callback)
// {
//     if(particleIdx >= m_ppos.size()) return;
//
//     const ngl::Vec3 pos(m_ppos[particleIdx].m_x, m_ppos[particleIdx].m_y, m_ppos[particleIdx].m_z);
//     auto centerCell = PositionToCellCoord(pos);
//     const float radiusSq = m_smoothingRadius * m_smoothingRadius;
//
//     for(int dx = -1; dx <= 1; ++dx) {
//         for(int dy = -1; dy <= 1; ++dy) {
//             for(int dz = -1; dz <= 1; ++dz) {
//                 size_t cellKey = HashCell(centerCell[0]+dx, centerCell[1]+dy, centerCell[2]+dz);
//                 size_t key = cellKey % m_startIndices.size();
//                 size_t startIdx = m_startIndices[key];
//
//                 if(startIdx == std::numeric_limits<size_t>::max()) continue;
//
//                 for(size_t i = startIdx; i < m_spatialLookup.size() &&
//                     i < m_maxParticles &&
//                     m_spatialLookup[i].cellKey == cellKey; ++i)
//                 {
//                     size_t otherIdx = m_spatialLookup[i].particleIndex;
//                     if(otherIdx >= m_ppos.size()) continue;
//                     if(particleIdx == otherIdx) continue;
//
//                     float distSq = (ngl::Vec3(m_ppos[otherIdx].m_x, m_ppos[otherIdx].m_y, m_ppos[otherIdx].m_z) - pos).lengthSquared();
//                     if(distSq <= radiusSq) {
//                         callback(otherIdx);
//                     }
//                 }
//             }
//         }
//     }
// }


size_t Emitter::size() const
{
    return m_maxParticles;
}



// float Emitter::CalculateDensity(size_t particleIdx)
// {
//     // Safety check
//     if(particleIdx >= m_ppos.size()) return 0.0f;
//
//     float density = 0.0f;
//     const ngl::Vec3 pos(m_ppos[particleIdx].m_x,
//                        m_ppos[particleIdx].m_y,
//                        m_ppos[particleIdx].m_z);
//
//     // Always include self-density
//     density += m_particleMass * SmoothingKernel(m_smoothingRadius, 0.0f);
//
//     // Process neighbors using spatial hashing
//     ForEachNeighbor(particleIdx, [&](size_t neighborIdx) {
//         if(neighborIdx >= m_ppos.size()) return; // Safety check
//
//         const ngl::Vec3 neighborPos(m_ppos[neighborIdx].m_x,
//                                   m_ppos[neighborIdx].m_y,
//                                   m_ppos[neighborIdx].m_z);
//         float dist = (neighborPos - pos).length();
//
//         // Only consider particles within smoothing radius
//         if(dist <= m_smoothingRadius) {
//             density += m_particleMass * SmoothingKernel(m_smoothingRadius, dist);
//         }
//     });
//
//     // Ensure non-zero density to prevent division by zero later
//     return std::max(density, 0.0001f);
// }


void Emitter::update(float _dt)
{
    // Calculate densities first
    m_physics->calculateAllDensities(m_maxParticles);

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
            // ngl::Vec3 emitDir(0.0f, 2.0f, 0.0f);
            //m_pdir[i] += gravity * _dt * 0.5f;
            m_psize[i] += 0.1f;
            m_psize[i] = std::clamp(m_psize[i], 0.0f, 5.0f);
            m_ppos[i].m_w = m_psize[i];

            //ngl::Vec3 cursorForce = calculateCursorForce(i);

            m_pcolour[i] = ngl::Random::getRandomColour3();
            ngl::Vec3 pressureForce = m_physics->calculatePressureForce(i, m_maxParticles); /// 10.0f;
            ngl::Vec3 viscosityForce = m_physics->calculateViscosityForce(i, m_maxParticles); /// 10.0f;
            ngl::Vec3 gravityForce(0.0f, -9.81f * m_physics->m_particleMass, 0.0f);
            ngl::Vec3 totalForce = pressureForce + viscosityForce + gravityForce ; //+ (cursorForce*2000000) * _dt;
            ngl::Vec3 acceleration = totalForce / m_physics->m_particleMass;// / 2.0f * _dt; /// 10000.0f;
            m_pdir[i] += acceleration * _dt;
            m_ppos[i] += m_pdir[i] * (1.0f / 50.0f);// * _dt; // * 0.5;

            m_boundingBox.resolveCollisions(i, m_ppos, m_pdir, m_psize);

        }
        // if (--m_plife[i] <= 0 || m_ppos[i].m_y < 0.0f)
        // {
        //     resetParticle(i);
        // }
    }

}



// void Emitter::update(float _dt)
// {
//     // 1. Update spatial lookup
//     UpdateSpatialLookup();
//
//     // 2. Calculate densities
//     #pragma omp parallel for
//     for(size_t i = 0; i < m_maxParticles; ++i) {
//         m_densities[i] = 0.0f;
//         ForEachNeighbor(i, [&](size_t j) {
//             float dist = (ngl::Vec3(m_ppos[j].m_x, m_ppos[j].m_y, m_ppos[j].m_z) -
//                          ngl::Vec3(m_ppos[i].m_x, m_ppos[i].m_y, m_ppos[i].m_z)).length();
//             m_densities[i] += m_particleMass * SmoothingKernel(m_smoothingRadius, dist);
//         });
//     }
//
//     // 3. Calculate forces and update
//     #pragma omp parallel for
//     for(size_t i = 0; i < m_maxParticles; ++i) {
//         if (!m_simulate) {
//             // Visual-only updates
//             m_psize[i] = std::clamp(m_psize[i] + 0.1f, 0.0f, 4.0f);
//             m_ppos[i].m_w = m_psize[i];
//             m_pcolour[i] = ngl::Random::getRandomColour3();
//             continue;
//         }
//
//         ngl::Vec3 pressureForce(0.0f, 0.0f, 0.0f);
//         ngl::Vec3 viscosityForce(0.0f, 0.0f, 0.0f);
//         ngl::Vec3 pos(m_ppos[i].m_x, m_ppos[i].m_y, m_ppos[i].m_z);
//
//         ForEachNeighbor(i, [&](size_t j) {
//             ngl::Vec3 dir = ngl::Vec3(m_ppos[j].m_x, m_ppos[j].m_y, m_ppos[j].m_z) - pos;
//             float dist = dir.length();
//             if(dist > 0.0f) {
//                 // Pressure force
//                 dir /= dist;
//                 float slope = SmoothingKernelDerivative(m_smoothingRadius, dist);
//                 float sharedPressure = calculateSharedPressure(m_densities[i], m_densities[j]);
//                 pressureForce -= dir * sharedPressure * slope * m_particleMass / m_densities[j];
//
//                 // Viscosity force
//                 float influence = SmoothingKernel(m_smoothingRadius, dist);
//                 viscosityForce += (m_pdir[j] - m_pdir[i]) *
//                                 (m_viscosityStrength * influence * m_particleMass / m_densities[j]);
//             }
//         });
//
//         // Update particle
//         ngl::Vec3 gravityForce(0.0f, -9.81f * m_particleMass, 0.0f);
//         ngl::Vec3 totalForce = pressureForce + gravityForce + viscosityForce * _dt;
//         m_pdir[i] += totalForce / m_particleMass / 200.0f;
//         m_ppos[i] += ngl::Vec4(m_pdir[i].m_x, m_pdir[i].m_y, m_pdir[i].m_z, 0.0f) * 0.5f;
//
//         // Visual updates
//         m_psize[i] = std::clamp(m_psize[i] + 0.1f, 0.0f, 5.0f);
//         m_ppos[i].m_w = m_psize[i];
//         m_pcolour[i] = ngl::Random::getRandomColour3();
//
//         resolveCollisions(i);
//     }
// }



void Emitter::initializeParticles()
{
    const float cubeSideLength = 25.0f; // Fixed 1000x1000x1000 volume
    const float jitterStrength = 0.0f;
    const float initialLife = 100.0f;
    const float baseVelocityY = 2.0f;
    const float particleSpacing = m_physics->m_particleSpacing; // Adjusts density
    const float yOffset = 20.0f;


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
                float py = y * effectiveParticleDiameter - halfLength + yOffset;
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
    if(!m_physics->m_densities.empty())
    {
        m_vao->setData(2, ngl::MultiBufferVAO::VertexData(
            m_physics->m_densities.size() * sizeof(float),
            m_physics->m_densities[0]  // Pass first element by reference
        ));
        m_vao->setVertexAttributePointer(2, 1, GL_FLOAT, 0, 0);
    }

    m_vao->setNumIndices(m_maxParticles);
    glEnable(GL_PROGRAM_POINT_SIZE);
    m_vao->draw();
    glDisable(GL_PROGRAM_POINT_SIZE);



    m_vao->unbind();

    m_boundingBox.renderBoundingBox();

}


// void Emitter::move(float _dx, float _dy, float _dz)
// {
//     m_pos.m_x += _dx;
//     m_pos.m_y += _dy;
//     m_pos.m_z += _dz;
// }



