//
//  filter.hpp
//  PMMVPS
//
//  Created by KaiWu on Nov/24/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#ifndef filter_hpp
#define filter_hpp

#include "patch.hpp"
#include "Eigen/Dense"
#include <vector>
#include <list>

using std::vector;
using std::list;
using Eigen::Vector2i;

class PmMvps;

class Filter {
    
public:
    Filter(PmMvps& pmmvps);
    
    void init();
    void run();
    float computeGain(const Patch& patch);
    int filterQuad(const Patch& patch, const vector<Ppatch>& neighbors) const;
    
protected:
    void filterOutside();
    void filterOutsideSub();
    
    void filterExact();
    void filterExactSub();
    
    void filterNeighbor(const int times);
    void filterNeighborSub();
    
    void filterSmallGroups();
    void filterSmallGroupsSub(const int pid, const int id, vector<int>& label, list<int>& ltmp) const;
    
    void setDepthMaps();
    void setDepthMapsSub();
    void setDepthMapsVGridsVPGridsAddPatchV(const int additive);
    void setVGridsVPGrids();
    void addPatchV();
    
    vector<float> m_gains;
    vector<vector<int> > m_newimages, m_removeimages;
    vector<vector<Vector2i> > m_newgrids, m_removegrids;
    int m_time;
    vector<int> m_rejects;
    
    PmMvps& m_pmmvps;
    
private:
    void ortho(const Vector4f& z, Vector4f& x, Vector4f& y) const;
    void lls(const vector<vector<float> >& A, const vector<float>& b, vector<float>& x) const;
    
};

#endif /* filter_hpp */
