//
// Created by s5616052 on 13/05/25.
//
// BoundingBox.cpp
#include "BoundingBox.h"
#include <ngl/VAOFactory.h>
#include <ngl/ShaderLib.h>
#include <ngl/VAOPrimitives.h>
#include <ngl/VAOFactory.h>
#include <ngl/ShaderLib.h>
#include <ngl/Transformation.h>
#include <ngl/MultiBufferVAO.h>

BoundingBox::BoundingBox()
{
}


void BoundingBox::initBoundingBoxVAO()
{
    // Match your collision bounds from resolveCollisions()
    const float hw = m_width/2, hd = m_depth/2;

    // Define the 8 corners of the box
    std::vector<ngl::Vec3> vertices = {
        // Bottom face
        {-hw, 0, -hd}, {hw, 0, -hd},
        {hw, 0, -hd}, {hw, 0, hd},
        {hw, 0, hd}, {-hw, 0, hd},
        {-hw, 0, hd}, {-hw, 0, -hd},

        // Top face
        {-hw, m_height, -hd}, {hw, m_height, -hd},
        {hw, m_height, -hd}, {hw, m_height, hd},
        {hw, m_height, hd}, {-hw, m_height, hd},
        {-hw, m_height, hd}, {-hw, m_height, -hd},

        // Vertical edges
        {-hw, 0, -hd}, {-hw, m_height, -hd},
        {hw, 0, -hd}, {hw, m_height, -hd},
        {hw, 0, hd}, {hw, m_height, hd},
        {-hw, 0, hd}, {-hw, m_height, hd}
    };

    // Create and bind VAO
    m_boxVAO = ngl::vaoFactoryCast<ngl::MultiBufferVAO>(
        ngl::VAOFactory::createVAO(ngl::multiBufferVAO, GL_LINES)
    );
    m_boxVAO->bind();

    // Set the vertex data (location = 0)
    m_boxVAO->setData(ngl::MultiBufferVAO::VertexData(
        vertices.size()*sizeof(ngl::Vec3),
        vertices[0].m_x
    ));

    // Set the vertex attribute pointer
    m_boxVAO->setVertexAttributePointer(
        0,                  // attribute location
        3,                  // number of components
        GL_FLOAT,           // type
        GL_FALSE,          // normalized?
        0,                 // stride
        0                  // offset
    );

    // Set number of indices (vertices) to draw
    m_boxVAO->setNumIndices(vertices.size());

    m_boxVAO->unbind();
}

void BoundingBox::renderBoundingBox() const
{
    if(!m_boxVAO) return;  // Safety check

    ngl::ShaderLib::use("nglColourShader");
    ngl::ShaderLib::setUniform("Colour", 1.0f, 0.0f, 0.0f, 1.0f); // Red color

    glLineWidth(2.0f);
    m_boxVAO->bind();
    m_boxVAO->draw();
    m_boxVAO->unbind();
    glLineWidth(1.0f);
}


void BoundingBox::resolveCollisions(size_t _i, std::vector<ngl::Vec4>& m_ppos, std::vector<ngl::Vec3>& m_pdir, std::vector<float>& m_psize)
{
    ngl::Vec3 boundSize(m_width, m_height, m_depth); // Proper bounding box
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
