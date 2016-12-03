//
//  depth_normal_init.hpp
//  PMMVPS
//
//  Created by KaiWu on Nov/2/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#ifndef depth_normal_init_hpp
#define depth_normal_init_hpp

#include <iostream>
#include "patch.hpp"
#include "Eigen/Dense"

using Eigen::Vector3f;

using std::vector;
using std::string;

class PmMvps;

class DepthNormInit {
public:
    DepthNormInit(PmMvps& pmmvps);
    ~DepthNormInit();
    
    void init(const string prefix, const int nfiles);
    void createPatches();
    void readDepths(vector<Vector3f>& coords);
    void readNormals(vector<vector<Vector3f> >& normals);
    void readDepthsNorms(vector<Vector3f>& coords, vector<Vector3f>& normals);
    
    // the plys are named as xxxxxxxx.ply, each x represents a decimal number
    // the first ply files stores the depth info, the rest of the file stores the
    // normal info for each viewpoint;
    string m_prefix;
    int m_nplys;
    
protected:
    PmMvps& m_pmmvps;
};



#endif /* depth_normal_init_hpp */
