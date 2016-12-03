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

using std::priority_queue;

class PmMvps;

class Propagate {
public:
    Propagate(PmMvps& pmmvps);
    ~Propagate();
    
    void init();
    void run(const int iter);
    void propagate(const int iter);
    void findViewNeighbors(pPatch& ppatch, vector<int>& images, vector<int>& indexes);
    int propagatePatch(pPatch& ppatch, vector<int>& indexes);
    /*
    int propagatePatch(const int image, const int x, const int y, const int xinc, const int yinc);
    void propagatePatch(pPatch& ppatch, vector<int>& images, vector<int>& indexes);
    void findNeighbors(pPatch& ppatch, const int image, const int x, const int y, const int xinc, const int yinc, vector<int>& images, vector<int>& indexes);
     */
    
    // maximum times of propagation of a patch to a cell
    int MAX_NUM_OF_PATCHES;
    
protected:
    priority_queue<pPatch, vector<pPatch>, PatchCmp> m_queue;
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
