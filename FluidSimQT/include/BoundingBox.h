//
// Created by s5616052 on 13/05/25.
//

#ifndef BOUNDINGBOX_H
#define BOUNDINGBOX_H

#pragma once
#include <vector>
#include <ngl/Vec3.h>
#include <ngl/AbstractVAO.h>
#include <memory>
#include <ngl/MultiBufferVAO.h>
#include <ngl/Vec4.h>

class BoundingBox
{
public:
    BoundingBox();

    void initBoundingBoxVAO();
    void renderBoundingBox() const;
    void resolveCollisions(size_t _i, std::vector<ngl::Vec4>& io_ppos, std::vector<ngl::Vec3>& io_pdir, std::vector<float>& io_psize);

    float m_width = 100.0f;
    float m_height = 50.0f;
    float m_depth = 50.0f;


private:
    std::unique_ptr<ngl::MultiBufferVAO> m_boxVAO;
};

#endif //BOUNDINGBOX_H
