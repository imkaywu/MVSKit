//
//  propagate.hpp
//  PMMVPS
//
//  Created by KaiWu on Oct/24/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#ifndef propagate_hpp
#define propagate_hpp

#include <iostream>
#include <vector>
#include <queue>
#include "patch_manager.hpp"
#include "Eigen/Dense"

using std::priority_queue;
using Eigen::Vector3f;
using Eigen::Vector4f;

class PmMvps;

class Propagate {
public:
    Propagate(PmMvps& pmmvps);
    ~Propagate();
    
    void init();
    void run(const int iter);
    void propagate(const int iter); // PatchMatch
    void propagateTest(const int iter); // PMVS
    void findViewNeighbors(Ppatch& ppatch, vector<int>& images, vector<int>& indexes);
    int propagatePatch(Ppatch& ppatch, vector<int>& indexes);
    int propagatePatch(const Ppatch& ppatch, const int image, const int index);
    // generate new patch from old patch and randomly generated icood
    Ppatch generatePatch(const Ppatch& patch, const Vector3f& icoord) const;
    /*
    int propagatePatch(const int image, const int x, const int y, const int xinc, const int yinc);
    void propagatePatch(Ppatch& ppatch, vector<int>& images, vector<int>& indexes);
    void findNeighbors(Ppatch& ppatch, const int image, const int x, const int y, const int xinc, const int yinc, vector<int>& images, vector<int>& indexes);
     */
    
    // for testing purpose
    float computeRadius(const Patch& patch);
    void findEmptyBlocks(const Ppatch& ppatch, vector<vector<Vector4f, Eigen::aligned_allocator<Vector4f> > >& canCoords);
    void ortho(const Vector4f& z, Vector4f& x, Vector4f& y);
    int expandSub(const Ppatch& orgppatch, const Vector4f& canCoord);
    int updateCounts(const Patch& patch);
    int checkCounts(Patch& patch);
    
    // maximum times of propagation of a patch to a cell
    int MAX_NUM_OF_PATCHES;
    int MAX_NUM_OF_PROPAG;
    
protected:
    priority_queue<Ppatch, vector<Ppatch>, PatchCmp> m_queue;
    PmMvps& m_pmmvps;
    
    // propagation related
    // number of trials;
    int m_ecount;
    // number of failures in the preProcess
    int m_fcount0;
    // number of failures in the postProcess
    int m_fcount1;
    // number of passes
    int m_pcount;
};

#endif /* propagate_hpp */
