//
//  pmmvps.cpp
//  PMMVPS
//
//  Created by KaiWu on Oct/24/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#include "pmmvps.hpp"
#include <stdio.h>

PmMvps::PmMvps() : m_patchManager(*this), m_dnInit(*this), m_propagate(*this), m_optim(*this), m_filter(*this) {
}

PmMvps::~PmMvps() {
}

void PmMvps::init(const Option& option) {
    m_images.clear();
    m_images = option.m_images;
    m_nimages = option.m_nimages;
    m_nillums = option.m_nillums;
    
    m_prefix = option.m_prefix;
    m_level = option.m_level;
    m_csize = option.m_csize;
    m_nccThreshold = option.m_nccThreshold;
    m_wsize = option.m_wsize;
    m_minImageNumThreshold = option.m_minImageNum;
    m_visdata = option.m_visdata;
    m_visdata2 = option.m_visdata2;
    m_tau = std::min(option.m_minImageNum * 2, m_nimages);
    m_depth = 0;
    
    // We set m_level + 3, to use multi-resolutional texture grabbing
    m_photoSet.init(m_images, m_prefix, m_nimages, m_nillums, m_level + 3, m_wsize, 1);
    
    // set pairwise distance based on distance of camera centers and angular diatance
    m_photoSet.setDistances();
    
    // Initialize each core member
    // m_patchManager should be initialize first
    m_patchManager.init();
    // init depths, normals, and patches
    m_dnInit.init(m_prefix, m_nimages + 1);
    // init propagation
    m_propagate.init();
    // init optimization
    m_optim.init();
    // init filter
    m_filter.init();
    
    // Initialize thresholds;
    m_angleThreshold0 = 60.0f * M_PI / 180.0f;
    m_angleThreshold1 = 60.0f * M_PI / 180.0f;
    
    m_countThreshold1 = 4;
    
    m_neighborThreshold = 0.5f;
    m_neighborThreshold1 = 1.0f;
    m_neighborThreshold2 = 1.0f;
    
    m_nccThresholdBefore = m_nccThreshold - 0.3f; // 0.7 - 0.3 = 0.4;
    
    m_maxAngleThreshold = option.m_maxAngleThreshold;
    
    m_quadThreshold = option.m_quadThreshold;
}

void PmMvps::updateThreshold() {
    m_nccThreshold -= 0.05f;
    m_nccThresholdBefore -= 0.05f;
    m_countThreshold1 = 2;
}

void PmMvps::run() {
    time_t tv;
    time(&tv);
    time_t curtime = tv;
    
    //----------------------------------------------------------------------
    // initialize depth and normal, and create initial patches
    //----------------------------------------------------------------------
    m_dnInit.createPatches(); // target = 0
    ++m_depth; // m_depth = 1 after init, before propag/optim
    
    //----------------------------------------------------------------------
    // patch expansion and optimization
    //----------------------------------------------------------------------
    const int ITER = 3;
    for (int iter = 0; iter < ITER; ++iter) {
        cerr << "\n---------------------" << endl
             << "Iteration: " << iter << endl
             << "---------------------" << endl;
        m_propagate.run(iter);
        // write to ply after each iteration
        cerr << "\nWriting Iter " << iter << " to file..." << endl;
        string name = m_prefix + "ply/refined_patches_before_refine_" + std::to_string(iter);
        m_patchManager.writePatches(name, true, false, false);
        
        m_filter.run();
        
        updateThreshold();
        
        ++m_depth;
        // write to ply after each iteration
        cerr << "\nWriting Iter " << iter << " to file..." << endl;
        name = m_prefix + "ply/refined_patches_" + std::to_string(iter);
        m_patchManager.writePatches(name, true, false, false);
    }
    
    time(&tv);
    cerr << "---- Total: " << (tv - curtime) / CLOCKS_PER_SEC << " secs ----" << endl;
}

// used to determine if two patches are spatial neighbors
int PmMvps::isNeighbor(const Patch& lhs, const Patch& rhs, const float neighborThreshold) const {
    const float hunit = (m_optim.getUnit(lhs.m_images[0], lhs.m_coord) + m_optim.getUnit(rhs.m_images[0], rhs.m_coord)) / 2.0f * m_csize; // (horizontal diff per pixel) * (pixel per cell) = (horizontal diff per cell)
    
    return isNeighbor(lhs, rhs, hunit, neighborThreshold);
}

int PmMvps::isNeighbor(const Patch& lhs, const Patch& rhs, const float hunit, const float neighborThreshold) const {
    if (lhs.m_normal.dot(rhs.m_normal) < cosf(120.0f / M_PI * 180.0f)) {
        return 0;
    }
    
    const Vector4f diff = lhs.m_coord - rhs.m_coord;
    const float vunit = lhs.m_dscale + rhs.m_dscale;
    const float f0 = lhs.m_normal.dot(diff);
    const float f1 = rhs.m_normal.dot(diff);
    float ftmp = (fabsf(f0) + fabsf(f1)) / 2.0f;
    ftmp /= vunit;
    
    const float hsize = (diff - f0 * lhs.m_normal + diff - f1 * rhs.m_normal).norm() / 2.0f / hunit;
    
    if (1.0f < hsize) {
        ftmp /= std::min(2.0f, hsize);
    }
    
    if (ftmp < neighborThreshold) {
        return 1;
    }
    else {
        return 0;
    }
}

int PmMvps::isNeighborRadius(const Patch& lhs, const Patch& rhs, const float hunit, const float neighborThreshold, const float radius) const {
    if (lhs.m_normal.dot(rhs.m_normal) < cos(120.0f * M_PI / 180.0f)) {
        return 0;
    }
    const Vector4f diff = rhs.m_coord - lhs.m_coord;
    const float vunit = lhs.m_dscale + rhs.m_dscale; // that's why we compute the vertical distance diff for 1/2 pixel diff
    // compute the vector's (diff) vertical component
    const float f0 = lhs.m_normal.dot(diff); // shouldn't we use line of sight instead
    const float f1 = rhs.m_normal.dot(diff);
    float ftmp = (fabsf(f0) + fabsf(f1)) / 2.0f;
    ftmp /= vunit;
    
    // compute the vector's (diff) horizontal component;
    const float hsize = (2 * diff - lhs.m_normal * f0 - rhs.m_normal * f1).norm() / 2.0f / hunit; // (horizontal diff)/(horizontal diff per cell) = cells
//    const float hsize = (diff - f0 * lhs.m_normal + diff - f1 * rhs.m_normal).norm() / 2.0f / hunit; // rewrite of the above statement
    
    // radius check
    if (radius / hunit < hsize) {
        return 0;
    }
    
    if (1.0f < hsize) {
        ftmp /= std::min(2.0f, hsize);
    }
    
    if (ftmp < neighborThreshold) {
        return 1;
    }
    else {
        return 0;
    }
}

