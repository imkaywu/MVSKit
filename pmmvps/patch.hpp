//
//  patch.hpp
//  PMMVPS
//
//  Created by KaiWu on Oct/19/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#ifndef patch_hpp
#define patch_hpp

#include <iostream>
#include <vector>
#include <memory>
#include "Eigen/Dense"

using std::endl;
using std::vector;
using std::shared_ptr;
using Eigen::Vector2i;
using Eigen::Vector4f;

class Patch
{
public:
    Patch();
    
//    float score (const float threshold) const;
    float score2 (const float threshold) const;
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    // 3D coordinate of the center of the patch
    Vector4f m_coord;
    // normal vector
    Vector4f m_normal;
    
    // potentially visible views, first image is the reference
    vector<int> m_images;
    vector<Vector2i> m_grids;
    // truely visible views
    vector<int> m_vimages;
    vector<Vector2i> m_vgrids;
    
    // average ncc score
    float m_ncc;
    // number of images in m_images;
    int m_nimages;
    // flag of PatchMatch stereo
    int m_iter;
    // flag for collectPatches
    // 0: not collected, 1: collected
    int m_collected;
    // flat for propagation
    // 0: not visited, 1: visited
    int m_flag;
    // for directional flag, not understood
    unsigned char m_dflag;
    // fixed patch or not, not understood, probably not used
    int m_fix;
    // id number in m_ppatches
    int m_id;
    // scaling factor corresponding to one pixel difference -- not understood
    float m_dscale; // distance/depth scaling
    float m_ascale; // angle scaling
    // used in priority_queue
    float m_tmp;
};

typedef shared_ptr<Patch> Ppatch;

std::istream& operator >>(std::istream& istr, Patch& rhs);
std::ostream& operator <<(std::ostream& ostr, const Patch& rhs);
std::istream& operator >>(std::istream& istr, Vector4f& v);
std::ostream& operator <<(std::ostream& ostr, const Vector4f& v);

#endif /* patch_hpp */
