//
//  patch_manager.hpp
//  PMMVPS
//
//  Created by KaiWu on Oct/19/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#ifndef patch_manager_hpp
#define patch_manager_hpp

#include <vector>
#include <queue>
#include "patch.hpp"
#include "Eigen/Dense"

using std::string;
using std::vector;
using std::priority_queue;
using Eigen::Vector3i;

class PmMvps;

class PatchCmp {
public:
    bool operator()(const Ppatch& lhs, const Ppatch& rhs) const {
        return lhs->m_tmp < rhs->m_tmp;
    }
};

class PatchManager
{
public:
    PatchManager(PmMvps& pmmvps);
    void init();
    // change the contents of m_images from images to indexes;
    void image2index(Patch& patch);
    // change the contents of m_images from indexes to images;
    void index2image(Patch& patch);
    //
    void collectPatches(const int target = 0);
    // collect all patches in m_pgrids to a priority-queue
    void collectPatches(priority_queue<Ppatch, vector<Ppatch>, PatchCmp>& pqPatches);
    // collect all patches visible in one image (m_pgrids[image]) to a priority-queue
    void collectPatches(const int index, priority_queue<Ppatch, vector<Ppatch>, PatchCmp>& pqPatches);
    // reset m_counts
    void clearCounts();
    // reset m_flag
    void clearFlags();
    // set m_pgrids, m_vpgrids and m_dpgrids from Ppatch
    void addPatch(Ppatch& ppatch);
    // update the m_dpgrids, called by addPatch()
    void updateDepthMaps(Ppatch& ppatch);
    // set Patch (m_images and m_grids) from a list of images
    void setGridsImages(Patch& patch, vector<int>& images) const;
    // set Patch (m_grids)
    void setGrids(Patch& patch);
	// set Patch (m_vgrids);
	void setVGrids(Patch& patch);
    // set Patch (m_vimages and m_vgrids)
    void setVImagesVGrids(Ppatch& ppatch);
    void setVImagesVGrids(Patch& patch);
    // remove a patch from m_pgrids and m_vpgrids
    void removePatch(const Ppatch& ppatch);
    // check is a patch is visible in an image
    int isVisible0(const Patch& patch, const int image, int& ix, int& iy, const float strict);
    int isVisible(const Patch& patch, const int image, const int& ix, const int& iy, const float strict);
    
    // methods used in propagation
    void findNeighbors(const Patch& patch, vector<Ppatch>& neighbors, const float scale = 1.0f, const int margin = 1, const int skipvis = 0);
    // there is a similar one in Optim, seem not used
    float computeUnit(const Patch& patch) const;
    // set m_dscale and m_ascale
    void setScales(Patch& patch) const;
    
    // compute ncc score, added by Kai, probably not used
    void computeNcc(Patch& patch) const;
    // sort a vector of patches based on ncc score
    void sortPatches(vector<Ppatch>& ppatches, const int ascend = 1) const;
    
    // read patches
    void readPatches();
	void readPatches(const int iter);
    // write results
    void writePatches(const string prefix, bool bExportPLY, bool bExportPatch, bool bExportPSet);
    void writePly(const vector<Ppatch>& ppatches, const string filename, const vector<Vector3i>& colors);
    void writePly(const vector<Ppatch>& ppatches, const string filename);
    
    // heights of the grids for all images
    vector<int> m_gheights;
    // widths of the grids for all images
    vector<int> m_gwidths;
    
    // image grid: m_nimages x (m_gheight * m_gwidth) x m_npatches
    // grids of patches for m_images
    vector<vector<vector<Ppatch> > > m_pgrids;
    // grids of patches for m_vimages
    vector<vector<vector<Ppatch> > > m_vpgrids;
    // grids of patches with closest depths, m_nimages x (m_gheight * m_gwidth)
    vector<vector<Ppatch> > m_dpgrids;
    // all the patches in the current level of m_pgrids
    vector<Ppatch> m_ppatches;
    // check how many times patch optimization was performed for expansion
    vector<vector<unsigned char> > m_counts;
    
    static Ppatch m_MAXDEPTH;
    static Ppatch m_BACKGROUND;
    
protected:
    PmMvps& m_pmmvps;
    
private:
    int m_nimages;
};

#endif /* patch_manager_hpp */
