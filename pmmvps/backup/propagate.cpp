//
//  propagate.cpp
//  PMMVPS
//
//  Created by KaiWu on Oct/24/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#include "propagate.hpp"
#include "pmmvps.hpp"

using std::cerr;
using std::endl;
using std::flush;
using Eigen::Vector3f;

Propagate::Propagate(PmMvps& pmmvps) : m_pmmvps(pmmvps){
}

Propagate::~Propagate() {
}

void Propagate::init() {
    MAX_NUM_OF_PATCHES = m_pmmvps.m_csize * m_pmmvps.m_csize;
}

void Propagate::run(const int iter) {
//    m_pmmvps.m_patchManager.collectPatches(m_queue);
    
    cerr << "Expanding patches..." << flush;
    
    propagate(iter);
}

void Propagate::propagate(const int iter) {
//    while (1) {
//        pPatch ppatch;
//        int isEmpty = 0;
//        if (m_queue.empty()) {
//            isEmpty = 1;
//        }
//        else {
//            ppatch = m_queue.top();
//            m_queue.pop();
//        }
//        
//        if (isEmpty) {
//            break;
//        }
//        // check if this patches has been visited before
//        if  (ppatch->m_iter <= iter) { // originally !=
//            continue;
//        }
//        
//        int image = ppatch->m_images[0];
        for (int image = 0; image < m_pmmvps.m_nviews; ++image) {
            int gwidth = m_pmmvps.m_patchManager.m_gwidths[image];
            int gheight = m_pmmvps.m_patchManager.m_gheights[image];
            int start = 0, end = gwidth * gheight - 1, inc = 1;
            if (iter % 2 == 1) {
                start = gwidth * gheight - 1;
                end = 0;
                inc = -1;
            }

            for (int index = start; index != end; index += inc) {
                vector<pPatch> ppatches = m_pmmvps.m_patchManager.m_pgrids[image][index];
                m_pmmvps.m_patchManager.sortPatches(ppatches, 0); // descending order
                int npatches = std::min(static_cast<int>(ppatches.size()), MAX_NUM_OF_PATCHES);
                if (npatches == 0) {
                    continue;
                }
                
                for (int n = 0; n < npatches; ++n) {
                    if (ppatches[n]->m_images[0] != image) {
                        continue;
                    }
                    
                    vector<int> images;
                    vector<int> indexes;
                    // spatial propagation
                    int lr = index + inc;
                    int tb = index + gwidth * inc;
                    image = ppatches[n]->m_images[0];
                    images.push_back(image);    indexes.push_back(lr);
                    images.push_back(image);    indexes.push_back(tb);
                    
                    // propagate patch to neighboring cells
                    cerr << "source: cell index: " << index << ", total patches: " << npatches << ", patch index: " << n + 1 << endl;
                    int isContinue = propagatePatch(ppatches[n], indexes);
                    if (isContinue == -1) {
                        break;
                    }
                }
            }
        }
//    }
}

// not used
void Propagate::findViewNeighbors(pPatch& ppatch, vector<int>& images, vector<int>& indexes) {
    for (int i = 1; i < ppatch->m_nimages; ++i) {
        Vector3f icoord = m_pmmvps.m_photoSet.project(ppatch->m_images[i], ppatch->m_coord, m_pmmvps.m_level);
        int cx = ((int)floorf(icoord(0) + 0.5f)) / m_pmmvps.m_csize; // used in setGrids, setGridsImage
        int cy = ((int)floorf(icoord(1) + 0.5f)) / m_pmmvps.m_csize;
        int gwidth = m_pmmvps.m_patchManager.m_gwidths[ppatch->m_images[i]];
        images.push_back(ppatch->m_images[i]);
        indexes.push_back(cy * gwidth + cx);
    }
}

// need to be cautious that the original patch should not be changed
int Propagate::propagatePatch(pPatch& ppatch, vector<int>& indexes) {
    int refImage = ppatch->m_images[0];
    int ncells = static_cast<int>(indexes.size()); // number of neighboring cells to propagate to
    int isSuccess = 0; // at least one successful propagation
    
    for (int i = 0; i < ncells; ++i) {
        vector<pPatch> ppatches = m_pmmvps.m_patchManager.m_pgrids[refImage][indexes[i]]; // patches in the neighboring cell
        m_pmmvps.m_patchManager.sortPatches(ppatches, 0); // descending order
        int npatches = std::min(static_cast<int>(ppatches.size()), MAX_NUM_OF_PATCHES); // number of patches to get refined in the neighboring cell
        int niter = std::max(1, npatches);
        
        vector<pPatch>& ppatchesRef = m_pmmvps.m_patchManager.m_pgrids[refImage][indexes[i]]; // used to remove patches from grid
        
        //----------initialize a new patch for each cell
        Patch patch(*ppatch);
        // estimate depth of the current patch
        const float depth = m_pmmvps.m_photoSet.m_photos[refImage].m_oaxis.dot(ppatch->m_coord);
        // backproject to get the coord for the neighbouring pixel
        const int cx = indexes[i] % m_pmmvps.m_patchManager.m_gwidths[refImage];
        const int cy = indexes[i] / m_pmmvps.m_patchManager.m_gwidths[refImage];
        const float fx = (m_pmmvps.m_csize * (2 * cx + 1) - 1) / 2.0f;
        const float fy = (m_pmmvps.m_csize * (2 * cy + 1) - 1) / 2.0f;
        Vector3f icoord;
        icoord << std::max(0.0f, std::min((float)m_pmmvps.m_photoSet.getWidth(refImage, m_pmmvps.m_level), fx)),
                  std::max(0.0f, std::min((float)m_pmmvps.m_photoSet.getWidth(refImage, m_pmmvps.m_level), fy)),
                  1.0f; // probably won't be outside of the image, but just in case
        icoord = depth * icoord;
        patch.m_coord = m_pmmvps.m_photoSet.m_photos[refImage].unproject(icoord, m_pmmvps.m_level);
        // set grids, the grid of the reference SHOULD be fixed since only the depth is updated
        m_pmmvps.m_patchManager.setGrids(patch);
        // recompute the ncc score
        m_pmmvps.m_patchManager.computeNcc(patch);
        // check if the grid coord to image coord is correct
        if (patch.m_grids[0](0) != cx || patch.m_grids[0](1) != cy) {
            cerr << "!!!pixel coordinate wrong" << endl;
        }
        
        for (int n = 0; n < niter; ++n) {
            cerr << "target: cell index: " << indexes[i] << ", total patches: " << std::max(1, npatches) << ", patch index: " << n + 1 << endl;
            // a NON-empty cell
            if (npatches > 0) {
                // they're the SAME patch, and the ncc score is not higher
                if (ppatch == ppatches[n] || ppatch->m_ncc <= ppatches[n]->m_ncc) {
                    cerr << "patch ncc is too lower to get propagated any more" << endl;
                    return -1;
                }
            }
            
            // generate a new patch for each patch in this cell
            Patch newpatch(patch);
            
            //----------patch optimization
            if (m_pmmvps.m_optim.preProcess(newpatch) == -1) {
                ++m_fcount0;
                cerr << "failed before optimization" << endl;
                break; // this patch won't work in this cell
            }
            
            m_pmmvps.m_optim.refinePatch(newpatch, 100);
            
            // the optimized patch is lower than the original one
            if (npatches > 0 && newpatch.m_ncc <= ppatches[n]->m_ncc) {
                cerr << "the optimized patch has lower ncc score" << endl;
                continue;
            }
            
            if (m_pmmvps.m_optim.postProcess(newpatch) == -1) {
                ++m_fcount1;
                cerr << "failed after optimization" << endl;
                break; // this patch won't work in this cell
            }
            
            ++m_pcount;
            isSuccess = 1;
            
            // patches that is currently being refined need to be removed from the grid
//            cerr << "# patches before removing: " << npatches << endl;
            ppatchesRef.erase(remove(ppatchesRef.begin(), ppatchesRef.end(), ppatches[n]), ppatchesRef.end()); // ppatches is ordered
//            cerr << "# patches after removing: " << static_cast<int>(ppatchesRef.size()) << endl;
//            cerr << "# patches after removing (double checking): " << static_cast<int>(m_pmmvps.m_patchManager.m_pgrids[refImage][indexes[i]].size()) << endl;

            pPatch newppatch(new Patch(newpatch));
            
            m_pmmvps.m_patchManager.addPatch(newppatch);

        }
        // to make sure each cell has most 2 * MAX_NUM_OF_PATCHES patches
        m_pmmvps.m_patchManager.sortPatches(ppatchesRef, 0); // descending order
        for (int n = static_cast<int>(ppatchesRef.size()) - 1; n >= 2 * MAX_NUM_OF_PATCHES; --n) {
            pPatch ppatch = ppatches[n];
            m_pmmvps.m_patchManager.removePatch(ppatch);
        }
    }
    if (isSuccess == 0) {
        return 0;
    }
//    ++ppatch->m_iter;
    return -1; // current patch successfully propagated, the lower ncc patch is less likely to succeed
}

//-------------------------------------------------
// trash code
//-------------------------------------------------
/*
void Propagate::propagate() {
    vector<int> images = m_pmmvps.m_images;
    int wsize = m_pmmvps.m_wsize;
    int hwsize = wsize / 2;
    int level = m_pmmvps.m_level;
    
    for (int index = 0; index < images.size(); ++index) {
        int cwidth = m_pmmvps.m_patchManager.m_gwidths[level];
        int cheight = m_pmmvps.m_patchManager.m_gheights[level];
        
        for (int iter = 0; iter < 3; ++iter) {
            int xstart = hwsize - 1, xend = cwidth - hwsize + 1, xinc = 1;
            int ystart = hwsize - 1, yend = cheight - hwsize + 1, yinc = 1;
            if (iter % 2 == 1) {
                xstart = cwidth - hwsize + 1; xend = hwsize - 1; xinc = -1;
                ystart = cheight - hwsize + 1; yend = hwsize - 1; yinc = -1;
            }
            for (int y = ystart; y < yend; y += yinc) {
                for (int x = xstart; x < xend; x += xinc) {
                    int ind = y * cwidth + x;
                    if (m_pmmvps.m_patchManager.m_vgrids[index][ind].empty()) {
                        continue;
                    }
                    
                    vector<pPatch> ppatches = m_pmmvps.m_patchManager.m_vgrids[index][ind];
                    int nPatches = static_cast<int>(ppatches.size());
                    
                    for (int n = 0; n < nPatches; ++n) {
                        pPatch ppatch = ppatches[n];
                        // ---------spatial propagation
                        // propagation along the x-axis
                        ind = y * cwidth + x + xinc;
                        vector<pPatch>& ppatchesNeighbor = m_pmmvps.m_patchManager.m_vgrids[index][ind];
                        int nPatchesNeighbor = static_cast<int>(ppatchesNeighbor.size());
                        if (ppatchesNeighbor.empty()) {
                            ppatchesNeighbor.push_back(ppatch);
                        }
                        else {
                            for (int m = 0; m < nPatchesNeighbor; ++m) {
                                pPatch ppatchNeighbor = ppatchesNeighbor[m];
                                // calculate the photo-consistency score of ppatch
                                float score0 = ppatch->score(m_pmmvps.m_nccThreshold);
                                // calculate the ncc score of ppatchNeighbor
                                float score1 = ppatchNeighbor->score(m_pmmvps.m_nccThreshold);
                                
                                if (score0 < score1) {
                                    ppatchesNeighbor[m] = ppatch;
                                }
//                                vector<vector<Vector3f> > texs;
//                                m_pmmvps.m_photoSet.m_photos[index].getTex(ppatch->m_vgrids[ind], m_pmmvps.m_level, Vector2f(1.0f, 0.0f), Vector2f(0.0f, 1.0f), wsize, texs[0], 1);
//                                m_pmmvps.m_photoSet.m_photos[index].getTex(ppatchNeighbor->m_vgrids[ind], m_pmmvps.m_level, Vector2f(1.0f, 0.0f), Vector2f(0.0f, 1.0f), wsize, texs[1], 1);
                            }
                        }
                        
                        // propagation along the y-axis
                        ind = (y + yinc) & cwidth + x;
                        ppatchesNeighbor = m_pmmvps.m_patchManager.m_vgrids[index][ind];
                        if (ppatchesNeighbor.empty()) {
                            ppatchesNeighbor.push_back(ppatch);
                        }
                        else {
                            for (int m = 0; m < nPatchesNeighbor; ++m) {
                                pPatch ppatchNeighbor = ppatchesNeighbor[m];
                                // calculate the photo-consistency score of ppatch
                                float score0 = ppatch->score(m_pmmvps.m_nccThreshold);
                                // calculate the photo-consistency score of ppatchNeighbor
                                float score1 = ppatchNeighbor->score(m_pmmvps.m_nccThreshold);
                                
                                if (score0 < score1) {
                                    ppatchesNeighbor[m] = ppatch;
                                }
                            }
                        }
                        
                        // ---------view propagation
                        
                    }
                }
            }
        }
    }
}

void Propagate::propagate() {
    while (1) {
        pPatch ppatch;
        int isEmpty = 0;
        if (m_queue.empty()) {
            isEmpty = 1;
        }
        else {
            ppatch = m_queue.top();
            m_queue.pop();
        }
        
        if (isEmpty) {
            break;
        }
        int wsize = m_pmmvps.m_wsize;
        int hwsize = wsize / 2;
        int image = ppatch->m_images[0];
        int cwidth = m_pmmvps.m_patchManager.m_gwidths[image];
        int cheight = m_pmmvps.m_patchManager.m_gheights[image];
        
        for (int iter = 0; iter < 3; ++iter) {
            int xstart = hwsize - 1, xend = cwidth - hwsize + 1, xinc = 1;
            int ystart = hwsize - 1, yend = cheight - hwsize + 1, yinc = 1;
            if (iter % 2 == 1) {
                xstart = cwidth - hwsize + 1; xend = hwsize - 1; xinc = -1;
                ystart = cheight - hwsize + 1; yend = hwsize - 1; yinc = -1;
            }
            for (int y = ystart; y < yend; y += yinc) {
                for (int x = xstart; x < xend; x += xinc) {
                    propagatePatch(image, x, y, xinc, yinc);
                }
            }
        }
    }
}

int Propagate::propagatePatch(const int image, const int x, const int y, const int xinc, const int yinc) {
    int width = m_pmmvps.m_patchManager.m_gwidths[image];
    int index = y * width + x;
    if (m_pmmvps.m_patchManager.m_pgrids[image][index].empty()) {
        cerr << "No patches in the current grid cell" << endl;
        return -1;
    }
    
    vector<pPatch> ppatches = m_pmmvps.m_patchManager.m_pgrids[image][index];
    int npatches = static_cast<int>(ppatches.size());
    
    for (int n = 0; n < npatches; ++n) {
        vector<int> images;
        vector<int> indexes;
        findNeighbors(ppatches[n], image, x, y, xinc, yinc, images, indexes);
        propagatePatch(ppatches[n], images, indexes);
    }
    
    return 0;
}

void Propagate::findNeighbors(pPatch& ppatch, const int image, const int x, const int y, const int xinc, const int yinc, vector<int>& images, vector<int>& indexes) {
    const int width = m_pmmvps.m_patchManager.m_gwidths[image];
    
    // spatial propagation
    images.clear();
    indexes.clear();
    images.push_back(image);
    indexes.push_back(y * width + x + xinc);
    indexes.push_back((y + yinc) * width + x);
    
    // view propagation
    for (int i = 1; i < static_cast<int>(ppatch->m_images.size()); ++i) {
        Vector3f icoord = m_pmmvps.m_photoSet.project(ppatch->m_images[i], ppatch->m_coord, m_pmmvps.m_level);
        // !!!transform to the cell coordinates, need to test
        int cx = (icoord(0) + m_pmmvps.m_csize - 1) / m_pmmvps.m_csize;
        int cy = (icoord(1) + m_pmmvps.m_csize - 1) / m_pmmvps.m_csize;
        images.push_back(ppatch->m_images[i]);
        indexes.push_back(cy * width + cx);
    }
    
}
*/