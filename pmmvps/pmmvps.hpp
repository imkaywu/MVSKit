//
//  pmmvps.hpp
//  PMMVPS
//
//  Created by KaiWu on Oct/24/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#ifndef pmmvps_hpp
#define pmmvps_hpp

#include "option.hpp"
#include "depth_normal_init.hpp"
#include "../image/photoSet.hpp"
#include "patch_manager.hpp"
#include "propagate.hpp"
#include "optim.hpp"
#include "filter.hpp"
#include <vector>
#include <string>

using std::vector;
using std::string;

class PmMvps {
public:
    PmMvps();
    virtual ~PmMvps();
    
    void init(const Option& option);
    void run();
    int isNeighborRadius(const Patch& lhs, const Patch& rhs, const float hunit, const float neighborThreshold, const float radius) const;
    int isNeighbor(const Patch& lhs, const Patch& rhs, const float hunit, const float neighborThreshold) const;
    int isNeighbor(const Patch& lhs, const Patch& rhs, const float neighborThreshold) const;
    
    // number of images/viewpoints
    int m_nimages;
    // number of illuminations
    int m_nillums;
    // total images
    vector<int> m_images;
    
    // root directory
    string m_prefix;
    // number of levels
    int m_level;
    // cell size
    int m_csize;
    // ncc threshold
    float m_nccThreshold;
    // window size
    int m_wsize;
    // minimum image number threshold
    int m_minImageNumThreshold;
    // visdata from SfM, m_nimages x m_nimages matrix
    vector<vector<int> > m_visdata;
    // an array of relevant images
    vector<vector<int> > m_visdata2;
    // threshold on filterQuad
    float m_quadThreshold;
    // Maximum number of images used in the optimization
    int m_tau;
    // ??? if patches are dense or not, that is, if we use check(patch) after patch optimization
    int m_depth;
    
    //----------------------------------------------------------------------
    // Thresholds
    //----------------------------------------------------------------------
    // for first feature matching. Images within this angle are used in matching
    float m_angleThreshold0;
    // tigher angle
    float m_angleThreshold1;
    
    // number of counts, expansion can be tried
    int m_countThreshold1;
    
    // Parameter for isNeighbor in findemptyblocks
    float m_neighborThreshold;
    // parameter for isNeighbor in filterOutside;
    float m_neighborThreshold1;
    // parameter for filterNeighbor
    float m_neighborThreshold2;
    
    // ncc threshold before optim
    float m_nccThresholdBefore;
    // Maximum angle of images must be at least as large as this
    float m_maxAngleThreshold;
    
    //----------------------------------------------------------------------
    // Core memebers
    //----------------------------------------------------------------------
    // photo set from different viewpoints, under different illuminations
    PhotoSet m_photoSet;
    // depth, normal, and patch initialization
    DepthNormInit m_dnInit;
    // patch manager
    PatchManager m_patchManager;
    // patch propagation
    Propagate m_propagate;
    // patch optimization
    Optim m_optim;
    // patch filter
    Filter m_filter;
    
protected:
    void updateThreshold();
};

#endif /* pmmvps_hpp */
